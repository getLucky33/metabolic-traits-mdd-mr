#!/usr/bin/env python3
"""Build allele-validated metabolic-trait variant_id-to-rsID maps from dbSNP 157."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REQUIRED_SOURCE_COLUMNS = {
    "variant_id", "effect_allele", "other_allele", "beta", "standard_error",
    "effect_allele_frequency", "neg_log_10_p_value",
}
EXPECTED_DBSNP_MD5 = "6a6f313e92a39c337571174dad12cfe1"
CONTIGS = {
    "1": "NC_000001.11", "2": "NC_000002.12", "3": "NC_000003.12",
    "4": "NC_000004.12", "5": "NC_000005.10", "6": "NC_000006.12",
    "7": "NC_000007.14", "8": "NC_000008.11", "9": "NC_000009.12",
    "10": "NC_000010.11", "11": "NC_000011.10", "12": "NC_000012.12",
    "13": "NC_000013.11", "14": "NC_000014.9", "15": "NC_000015.10",
    "16": "NC_000016.10", "17": "NC_000017.11", "18": "NC_000018.10",
    "19": "NC_000019.10", "20": "NC_000020.11", "21": "NC_000021.9",
    "22": "NC_000022.11", "X": "NC_000023.11", "23": "NC_000023.11",
    "Y": "NC_000024.10", "24": "NC_000024.10", "MT": "NC_012920.1",
}
RSID = re.compile(r"rs[0-9]+$")


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
    required = {"catalog", "metabolic_download_dir", "rsid_map_dir", "bcftools", "dbsnp", "dbsnp_index"}
    missing = sorted(required - set(values))
    if missing:
        raise ValueError("configuration is missing keys: " + ", ".join(missing))
    return values


def find_one(directory: Path, filename: str) -> Path:
    hits = [path for path in directory.rglob(filename) if path.is_file()]
    if len(hits) != 1:
        raise ValueError(f"expected one {filename} below {directory}; found {len(hits)}")
    return hits[0]


def file_md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_variant_id(value: str) -> tuple[str, int, str, str] | None:
    parts = value.split("_")
    if len(parts) != 4 or not parts[1].isdigit():
        return None
    chrom, position, other, effect = parts
    return chrom.upper(), int(position), other.upper(), effect.upper()


def valid_allele(value: str) -> bool:
    return bool(value) and value not in {".", "*", "-"} and "," not in value and set(value) <= set("ACGT")


def allele_set_match(other: str, effect: str, ref: str, alt: str) -> bool:
    return (other == ref and effect == alt) or (other == alt and effect == ref)


def collect_candidates(
    source: Path,
    p_threshold: float,
) -> tuple[list[str], dict[str, tuple[str, int, str, str] | None], dict[str, str | int]]:
    digest = hashlib.sha256()
    variants: dict[str, tuple[str, int, str, str] | None] = {}
    ordered: list[str] = []
    rows_scanned = 0
    with gzip.open(source, "rb") as handle:
        header_raw = handle.readline()
        digest.update(header_raw)
        fields = header_raw.decode("utf-8").rstrip("\r\n").split("\t")
        missing = REQUIRED_SOURCE_COLUMNS - set(fields)
        if missing:
            raise ValueError(f"{source} is missing columns: {', '.join(sorted(missing))}")
        index = {name: fields.index(name) for name in REQUIRED_SOURCE_COLUMNS}
        max_index = max(index.values())
        for raw in handle:
            digest.update(raw)
            if not raw.strip():
                continue
            rows_scanned += 1
            values = raw.decode("utf-8").rstrip("\r\n").split("\t")
            if len(values) <= max_index:
                raise ValueError(f"truncated source row {rows_scanned} in {source}")
            try:
                beta = float(values[index["beta"]])
                se = float(values[index["standard_error"]])
                eaf = float(values[index["effect_allele_frequency"]])
                neglog = float(values[index["neg_log_10_p_value"]])
            except ValueError:
                continue
            if not (math.isfinite(beta) and math.isfinite(se) and se > 0 and
                    math.isfinite(eaf) and 0 < eaf < 1 and math.isfinite(neglog)):
                continue
            if neglog <= -math.log10(p_threshold):
                continue
            variant_id = values[index["variant_id"]].strip()
            if not variant_id or variant_id in variants:
                raise ValueError(f"empty or duplicate candidate variant_id in {source}: {variant_id or '<empty>'}")
            parsed = parse_variant_id(variant_id)
            if parsed is not None:
                _, _, other, effect = parsed
                source_other = values[index["other_allele"]].strip().upper()
                source_effect = values[index["effect_allele"]].strip().upper()
                if (other, effect) != (source_other, source_effect):
                    raise ValueError(f"variant_id allele order differs from source columns: {variant_id}")
            variants[variant_id] = parsed
            ordered.append(variant_id)
    if not ordered:
        raise ValueError(f"no variants passed the documented mapping-stage significance rule: {source}")
    return ordered, variants, {
        "source_rows_scanned": rows_scanned,
        "candidate_rows": len(ordered),
        "source_decompressed_sha256": digest.hexdigest(),
    }


def query_dbsnp(
    bcftools: Path,
    dbsnp: Path,
    variants: dict[str, tuple[str, int, str, str] | None],
    work_dir: Path,
) -> dict[tuple[str, int], list[tuple[str, str, str]]]:
    positions = sorted({
        (CONTIGS[value[0]], value[1])
        for value in variants.values()
        if value is not None and value[0] in CONTIGS and valid_allele(value[2]) and valid_allele(value[3])
    })
    if not positions:
        raise ValueError("no valid GRCh38 positions remained for the dbSNP query")
    bed = work_dir / "exposure_rsid_positions.bed"
    with bed.open("w", encoding="utf-8", newline="") as handle:
        for contig, position in positions:
            handle.write(f"{contig}\t{position - 1}\t{position}\n")
    hits = work_dir / "dbsnp_hits.tsv"
    with hits.open("w", encoding="utf-8", newline="") as output:
        result = subprocess.run(
            [str(bcftools), "query", "-R", str(bed), "-f", "%CHROM\t%POS\t%REF\t%ALT\t%ID\n", str(dbsnp)],
            stdout=output, stderr=subprocess.PIPE, text=True,
        )
    if result.returncode != 0:
        raise RuntimeError(f"bcftools query failed ({result.returncode}): {result.stderr[:2000]}")
    records: dict[tuple[str, int], list[tuple[str, str, str]]] = {}
    with hits.open(encoding="utf-8") as handle:
        for line in handle:
            values = line.rstrip("\r\n").split("\t")
            if len(values) != 5:
                raise ValueError("bcftools returned a malformed dbSNP row")
            try:
                position = int(values[1])
            except ValueError as error:
                raise ValueError("bcftools returned a non-numeric dbSNP position") from error
            records.setdefault((values[0], position), []).append((values[2].upper(), values[3].upper(), values[4]))
    return records


def assign_rsids(
    variants: dict[str, tuple[str, int, str, str] | None],
    records: dict[tuple[str, int], list[tuple[str, str, str]]],
) -> dict[str, tuple[str, str]]:
    mapped: dict[str, tuple[str, str]] = {}
    for variant_id, value in variants.items():
        if value is None or value[0] not in CONTIGS or not valid_allele(value[2]) or not valid_allele(value[3]):
            mapped[variant_id] = (".", "skipped")
            continue
        chrom, position, other, effect = value
        position_records = records.get((CONTIGS[chrom], position), [])
        if not position_records:
            mapped[variant_id] = (".", "no_dbsnp_entry")
            continue
        matched: set[str] = set()
        for ref, alts, identifiers in position_records:
            for alt in alts.split(","):
                if allele_set_match(other, effect, ref, alt):
                    matched.update(identifier for identifier in re.split(r"[;,]", identifiers) if RSID.fullmatch(identifier))
        if not matched:
            mapped[variant_id] = (".", "no_allele_match")
        elif len(matched) == 1:
            mapped[variant_id] = (next(iter(matched)), "matched")
        else:
            mapped[variant_id] = (",".join(sorted(matched)), "ambiguous")
    return mapped


def write_outputs(
    output_dir: Path,
    trait_variants: dict[str, list[str]],
    mapped: dict[str, tuple[str, str]],
    source_qc: dict[str, dict[str, str | int]],
) -> list[dict[str, str | int | float]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {trait: output_dir / f"{trait}_rsid_map.tsv" for trait in trait_variants}
    qc_rows: list[dict[str, str | int | float]] = []
    for trait, variant_ids in trait_variants.items():
        path = paths[trait]
        temporary = path.with_suffix(path.suffix + ".tmp")
        counts = {status: 0 for status in ("matched", "ambiguous", "no_allele_match", "no_dbsnp_entry", "skipped")}
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(("variant_id", "rsid", "rsid_status"))
            for variant_id in variant_ids:
                rsid, status = mapped[variant_id]
                counts[status] += 1
                writer.writerow((variant_id, rsid, status))
        temporary.replace(path)
        content = path.read_bytes()
        total = len(variant_ids)
        qc_rows.append({
            "trait": trait,
            **source_qc[trait],
            **counts,
            "match_fraction": f"{counts['matched'] / total:.12g}",
            "map_sha256": hashlib.sha256(content).hexdigest(),
        })
    return qc_rows


def write_qc(output_dir: Path, qc_rows: list[dict[str, str | int | float]]) -> None:
    qc_path = output_dir / "exposure_rsid_mapping_qc.tsv"
    with qc_path.open("w", encoding="utf-8", newline="") as handle:
        columns = list(qc_rows[0])
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(qc_rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--traits", help="comma-separated trait identifiers; omit only after a one-trait smoke test")
    parser.add_argument("--p-threshold", type=float, default=5e-8)
    parser.add_argument("--batch-size", type=int, default=8, help="traits per indexed dbSNP query")
    parser.add_argument("--expected-dbsnp-md5", default=EXPECTED_DBSNP_MD5)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if not 0 < args.p_threshold < 1:
        parser.error("require 0 < p-threshold < 1")
    if args.batch_size < 1:
        parser.error("batch-size must be positive")

    config = load_config(args.config.resolve())
    for key in ("catalog", "metabolic_download_dir", "bcftools", "dbsnp", "dbsnp_index"):
        if not config[key].exists():
            raise FileNotFoundError(f"missing {key}: {config[key]}")
    adjacent_indexes = {Path(str(config["dbsnp"]) + ".tbi"), Path(str(config["dbsnp"]) + ".csi")}
    if config["dbsnp_index"].resolve() not in {path.resolve() for path in adjacent_indexes}:
        raise ValueError("dbSNP index must be the adjacent .tbi or .csi used by bcftools")
    if file_md5(config["dbsnp"]) != args.expected_dbsnp_md5.lower():
        raise ValueError("dbSNP file MD5 differs from the frozen build-157 source")

    catalog = read_tsv(config["catalog"])
    if len(catalog) != 249 or len({row.get("trait") for row in catalog}) != 249:
        raise ValueError("catalog must contain 249 unique metabolic traits")
    requested = set(filter(None, (args.traits or "").split(",")))
    if requested:
        unknown = sorted(requested - {row["trait"] for row in catalog})
        if unknown:
            raise ValueError("unknown trait identifiers: " + ", ".join(unknown))
        catalog = [row for row in catalog if row["trait"] in requested]

    output_dir = (args.out_dir or config["rsid_map_dir"]).resolve()
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    final_maps = [output_dir / f"{row['trait']}_rsid_map.tsv" for row in catalog]
    existing = [path for path in final_maps if path.exists()]
    if existing and not args.overwrite:
        raise FileExistsError(f"refusing to overwrite {len(existing)} existing rsID map(s); use --overwrite explicitly")

    qc_rows: list[dict[str, str | int | float]] = []
    total_unique = 0
    with tempfile.TemporaryDirectory(prefix="exposure-rsid-stage-", dir=output_dir.parent) as staging_name:
        staging = Path(staging_name)
        for start in range(0, len(catalog), args.batch_size):
            batch = catalog[start:start + args.batch_size]
            trait_variants: dict[str, list[str]] = {}
            source_qc: dict[str, dict[str, str | int]] = {}
            batch_variants: dict[str, tuple[str, int, str, str] | None] = {}
            for row in batch:
                trait = row["trait"]
                if Path(trait).name != trait:
                    raise ValueError(f"unsafe trait identifier: {trait}")
                source = find_one(config["metabolic_download_dir"], f"{row['accession']}.tsv.gz")
                ordered, variants, qc = collect_candidates(source, args.p_threshold)
                for variant_id, value in variants.items():
                    if variant_id in batch_variants and batch_variants[variant_id] != value:
                        raise ValueError(f"conflicting variant_id definitions across traits: {variant_id}")
                    batch_variants[variant_id] = value
                trait_variants[trait] = ordered
                source_qc[trait] = {"accession": row["accession"], **qc}
                print(f"scanned {trait}: {qc['source_rows_scanned']} rows; {qc['candidate_rows']} candidates", flush=True)
            total_unique += len(batch_variants)
            with tempfile.TemporaryDirectory(prefix="dbsnp-query-", dir=staging) as query_name:
                records = query_dbsnp(config["bcftools"], config["dbsnp"], batch_variants, Path(query_name))
            mapped = assign_rsids(batch_variants, records)
            qc_rows.extend(write_outputs(staging, trait_variants, mapped, source_qc))

        write_qc(staging, qc_rows)
        run = {
            "catalog": str(config["catalog"]),
            "dbsnp_build": "157 GRCh38.p14",
            "dbsnp_md5": args.expected_dbsnp_md5.lower(),
            "traits": len(catalog),
            "candidate_variant_occurrences": sum(int(row["candidate_rows"]) for row in qc_rows),
            "batch_unique_candidate_variants": total_unique,
            "p_threshold": args.p_threshold,
            "batch_size": args.batch_size,
            "maps": len(qc_rows),
        }
        (staging / "exposure_rsid_mapping_run.json").write_text(
            json.dumps(run, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n"
        )
        output_dir.mkdir(parents=True, exist_ok=True)
        for path in staging.iterdir():
            if path.is_file():
                path.replace(output_dir / path.name)
    print(f"PASS: wrote {len(qc_rows)} allele-validated rsID map(s) to {output_dir}")


if __name__ == "__main__":
    main()
