#!/usr/bin/env python3
"""Validate external inputs and build local manifests without downloading data."""
from __future__ import annotations

import argparse
import csv
import gzip
from functools import lru_cache
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CATALOG_COLUMNS = {"accession", "trait", "trait_display", "genome_build", "sample_size", "download_url"}
METABOLIC_COLUMNS = {
    "variant_id", "effect_allele", "other_allele", "beta", "standard_error",
    "effect_allele_frequency", "neg_log_10_p_value",
}
PGC_COLUMNS = {"CHR", "BP", "SNP", "A1", "A2", "OR", "SE", "P", "Nca", "Nco"}
FINNGEN_COLUMNS = [
    "#chrom", "pos", "ref", "alt", "rsids", "nearest_genes", "pval", "mlogp",
    "beta", "sebeta", "af_alt", "af_alt_cases", "af_alt_controls",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def resolve(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def load_config(path: Path) -> dict[str, Path]:
    rows = read_tsv(path)
    if not rows or set(rows[0]) != {"key", "path"}:
        raise ValueError("configuration must contain exactly the columns key and path")
    values = {row["key"].strip(): resolve(row["path"].strip()) for row in rows}
    required = {
        "catalog", "metabolic_download_dir", "rsid_map_dir", "pgc_main", "pgc_noukbb",
        "finngen", "ld_bfile", "plink", "bcftools", "dbsnp", "dbsnp_index",
    }
    missing = sorted(required - set(values))
    if missing:
        raise ValueError("configuration is missing keys: " + ", ".join(missing))
    return values


def load_catalog(path: Path) -> list[dict[str, str]]:
    rows = read_tsv(path)
    if len(rows) != 249 or set(rows[0]) != CATALOG_COLUMNS:
        raise ValueError("metabolic catalog must contain 249 rows and the documented six columns")
    accessions = [row["accession"] for row in rows]
    expected = [f"GCST{number}" for number in range(90451106, 90451355)]
    if accessions != expected or len({row["trait"] for row in rows}) != 249:
        raise ValueError("metabolic catalog must contain the ordered contiguous accession range and 249 unique traits")
    if any(row["genome_build"] != "GRCh38" or row["sample_size"] != "599249" for row in rows):
        raise ValueError("metabolic catalog build or sample-size metadata differs from the recorded meta_EUR series")
    return rows


@lru_cache(maxsize=None)
def file_index(directory: Path) -> dict[str, list[Path]]:
    index: dict[str, list[Path]] = {}
    for path in directory.rglob("*"):
        if path.is_file():
            index.setdefault(path.name, []).append(path)
    return index


def find_one(directory: Path, filename: str) -> Path:
    hits = file_index(directory).get(filename, [])
    if len(hits) != 1:
        raise ValueError(f"expected one {filename} below {directory}; found {len(hits)}")
    return hits[0]


def find_map(directory: Path, trait: str) -> Path:
    names = (f"{trait}_rsid_map.tsv", f"{trait}_forward_ivs_rsid_b157.tsv")
    hits = [path for name in names for path in file_index(directory).get(name, [])]
    if len(hits) != 1:
        raise ValueError(f"expected one rsID map for {trait} below {directory}; found {len(hits)}")
    return hits[0]


def header(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="strict") as handle:
        line = handle.readline().strip()
    return line.split("\t") if "\t" in line else line.split()


def require_file(path: Path, label: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"{label} is missing or empty: {path}")


def inventory(config: dict[str, Path], level: str) -> tuple[list[dict[str, str]], dict[str, Path]]:
    catalog = load_catalog(config["catalog"])
    files: dict[str, Path] = {}
    for row in catalog:
        source = find_one(config["metabolic_download_dir"], f"{row['accession']}.tsv.gz")
        missing = METABOLIC_COLUMNS - set(header(source))
        if missing:
            raise ValueError(f"{row['accession']} is missing columns: {', '.join(sorted(missing))}")
        files[row["trait"]] = source

    for key in ("pgc_main", "pgc_noukbb", "finngen", "plink", "bcftools", "dbsnp", "dbsnp_index"):
        require_file(config[key], key)
    for suffix in (".bed", ".bim", ".fam"):
        require_file(Path(str(config["ld_bfile"]) + suffix), f"LD reference {suffix}")

    for key in ("pgc_main", "pgc_noukbb"):
        fields = set(header(config[key]))
        missing = PGC_COLUMNS - fields
        if missing or not any(name.startswith("FRQ_A_") for name in fields) or not any(name.startswith("FRQ_U_") for name in fields):
            raise ValueError(f"{key} does not match the documented PGC schema")
    if header(config["finngen"]) != FINNGEN_COLUMNS:
        raise ValueError("FinnGen file does not match the locked 13-column R13 schema")

    if level == "analysis":
        for row in catalog:
            mapping = find_map(config["rsid_map_dir"], row["trait"])
            fields = set(header(mapping))
            if not {"variant_id", "rsid", "rsid_status"}.issubset(fields):
                raise ValueError(f"rsID map has an invalid schema: {mapping}")
    return catalog, files


def manifest_path(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return path.resolve().as_posix()


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build_manifests(config: dict[str, Path], out_dir: Path) -> None:
    catalog, files = inventory(config, "analysis")
    exposure_rows = [{"trait": row["trait"], "source_file": manifest_path(files[row["trait"]])} for row in catalog]
    source_rows = [
        {
            "trait": row["trait"],
            "source_file": manifest_path(files[row["trait"]]),
            "source_type": "metabolite",
            "rsid_map_file": manifest_path(find_map(config["rsid_map_dir"], row["trait"])),
        }
        for row in catalog
    ]
    source_rows.append({
        "trait": "Major_depressive_disorder",
        "source_file": manifest_path(config["pgc_main"]),
        "source_type": "mdd",
        "rsid_map_file": "",
    })
    write_tsv(out_dir / "exposure_manifest.local.tsv", ["trait", "source_file"], exposure_rows)
    write_tsv(out_dir / "source_manifest.local.tsv", ["trait", "source_file", "source_type", "rsid_map_file"], source_rows)
    print(f"PASS: wrote 249 exposure rows and 250 source-manifest rows to {out_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    check = sub.add_parser("check", help="validate downloaded or analysis-ready inputs")
    check.add_argument("--config", type=Path, required=True)
    check.add_argument("--level", choices=("sources", "analysis"), default="sources")
    build = sub.add_parser("build-manifests", help="write complete local manifests after analysis-level checks")
    build.add_argument("--config", type=Path, required=True)
    build.add_argument("--out-dir", type=Path, default=ROOT / "analysis" / "manifests")
    args = parser.parse_args()
    config = load_config(args.config.resolve())
    if args.command == "check":
        catalog, _ = inventory(config, args.level)
        print(f"PASS: {args.level} inputs validated for {len(catalog)} metabolic traits")
    else:
        build_manifests(config, args.out_dir.resolve())


if __name__ == "__main__":
    main()
