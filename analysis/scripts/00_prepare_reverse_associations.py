#!/usr/bin/env python3
"""Prepare MDD-instrument associations in metabolite GWAS files.

The program scans each provider file once. It does not download data and does
not embed credentials. A dbSNP GRCh38 VCF/BCF and bcftools are used only when a
precomputed rsID-coordinate table is not supplied.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import subprocess
import tempfile
from pathlib import Path


CONTIG = {
    f"NC_{i:06d}.{v}": str(i)
    for i, v in {
        1: 11, 2: 12, 3: 12, 4: 12, 5: 10, 6: 12, 7: 14, 8: 11,
        9: 12, 10: 11, 11: 10, 12: 12, 13: 11, 14: 9, 15: 10,
        16: 10, 17: 11, 18: 10, 19: 10, 20: 11, 21: 9, 22: 11,
    }.items()
}


def canonical_chrom(value: str) -> str:
    value = value.removeprefix("chr")
    return CONTIG.get(value, value)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def require_unique(rows: list[dict[str, str]], key: str, label: str) -> None:
    values = [row.get(key, "").strip() for row in rows]
    if any(not value for value in values):
        raise ValueError(f"{label} contains a blank {key}")
    if len(values) != len(set(values)):
        raise ValueError(f"{label} must contain unique {key} values")


def build_coords(mdd: dict[str, dict[str, str]], bcftools: Path, dbsnp: Path) -> dict[str, tuple[str, int]]:
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", delete=False, suffix=".txt") as handle:
        handle.write("\n".join(sorted(mdd)) + "\n")
        ids_file = Path(handle.name)
    try:
        cmd = [str(bcftools), "query", "-i", f"ID=@{ids_file}", "-f", "%ID\t%CHROM\t%POS\n", str(dbsnp)]
        run = subprocess.run(cmd, check=True, capture_output=True, text=True)
    finally:
        ids_file.unlink(missing_ok=True)
    coords: dict[str, tuple[str, int]] = {}
    for line in run.stdout.splitlines():
        rid, chrom, pos = line.split("\t")[:3]
        chrom = canonical_chrom(chrom)
        if rid in coords and coords[rid] != (chrom, int(pos)):
            raise RuntimeError(f"conflicting dbSNP positions for {rid}")
        coords[rid] = (chrom, int(pos))
    missing = sorted(set(mdd) - set(coords))
    if missing:
        raise RuntimeError(f"dbSNP mapping failed for {len(missing)} instruments; first={missing[0]}")
    return coords


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mdd-iv", type=Path, required=True)
    parser.add_argument("--exposure-manifest", type=Path, required=True,
                        help="TSV columns: trait,source_file")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--coords", type=Path,
                        help="Optional TSV columns: rsid,chr,pos")
    parser.add_argument("--bcftools", type=Path)
    parser.add_argument("--dbsnp", type=Path)
    args = parser.parse_args()

    mdd_rows = read_tsv(args.mdd_iv)
    required_mdd = {"rsid", "A1", "A2", "beta", "se", "p", "eaf"}
    if not mdd_rows or not required_mdd.issubset(mdd_rows[0]):
        raise ValueError(f"MDD IV file needs columns: {sorted(required_mdd)}")
    require_unique(mdd_rows, "rsid", "MDD IV file")
    mdd = {row["rsid"]: row for row in mdd_rows}

    if args.coords:
        coords_rows = read_tsv(args.coords)
        if not coords_rows or not {"rsid", "chr", "pos"}.issubset(coords_rows[0]):
            raise ValueError("coordinate table needs rsid, chr and pos columns")
        require_unique(coords_rows, "rsid", "coordinate table")
        coords = {r["rsid"]: (canonical_chrom(r["chr"]), int(r["pos"])) for r in coords_rows}
        if set(coords) != set(mdd):
            raise ValueError("coordinate table must contain exactly the MDD instrument rsIDs")
    else:
        if not args.bcftools or not args.dbsnp:
            parser.error("provide --coords or both --bcftools and --dbsnp")
        coords = build_coords(mdd, args.bcftools, args.dbsnp)

    by_position: dict[tuple[str, int], list[str]] = {}
    for rid, value in coords.items():
        by_position.setdefault(value, []).append(rid)

    manifest = read_tsv(args.exposure_manifest)
    if not manifest or not {"trait", "source_file"}.issubset(manifest[0]):
        raise ValueError("exposure manifest needs trait and source_file columns")
    require_unique(manifest, "trait", "exposure manifest")

    fields = [
        "rsid", "mdd_A1", "mdd_A2", "mdd_beta", "mdd_se", "mdd_p", "mdd_eaf",
        "trait", "met_effect", "met_other", "met_beta", "met_se", "met_eaf", "met_p",
    ]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    seen: dict[tuple[str, str], tuple[str, ...]] = {}
    with args.out.open("w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for item in manifest:
            trait = item["trait"]
            source = Path(item["source_file"])
            opener = gzip.open if source.suffix == ".gz" else open
            with opener(source, "rt", encoding="utf-8", errors="replace", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                required = {
                    "chromosome", "base_pair_location", "effect_allele", "other_allele",
                    "beta", "standard_error", "effect_allele_frequency", "neg_log_10_p_value",
                }
                if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                    raise ValueError(f"{source} is missing required metabolite-GWAS columns")
                for row in reader:
                    key = (canonical_chrom(row["chromosome"]), int(row["base_pair_location"]))
                    for rid in by_position.get(key, []):
                        exposure = mdd[rid]
                        if {exposure["A1"], exposure["A2"]} != {row["effect_allele"], row["other_allele"]}:
                            continue
                        signature = tuple(row[k] for k in sorted(required))
                        pair = (rid, trait)
                        if pair in seen:
                            if seen[pair] != signature:
                                raise RuntimeError(f"conflicting duplicate association for {rid} x {trait}")
                            continue
                        seen[pair] = signature
                        writer.writerow({
                            "rsid": rid, "mdd_A1": exposure["A1"], "mdd_A2": exposure["A2"],
                            "mdd_beta": exposure["beta"], "mdd_se": exposure["se"],
                            "mdd_p": exposure["p"], "mdd_eaf": exposure["eaf"], "trait": trait,
                            "met_effect": row["effect_allele"], "met_other": row["other_allele"],
                            "met_beta": row["beta"], "met_se": row["standard_error"],
                            "met_eaf": row["effect_allele_frequency"],
                            "met_p": 10.0 ** (-float(row["neg_log_10_p_value"])),
                        })

    expected_pairs = len(mdd) * len(manifest)
    print(f"wrote {len(seen)} associations; possible rsID-trait pairs={expected_pairs}; output={args.out}")


if __name__ == "__main__":
    main()
