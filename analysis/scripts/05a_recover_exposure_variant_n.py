#!/usr/bin/env python3
"""Recover source-GWAS per-variant exposure sample sizes for selected IVs."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def require_manifest(path: Path, columns: set[str], label: str) -> list[dict[str, str]]:
    rows = read_tsv(path)
    if not rows or not columns.issubset(rows[0]):
        raise ValueError(f"{label} needs columns: {sorted(columns)}")
    traits = [row["trait"].strip() for row in rows]
    if any(not trait for trait in traits) or len(traits) != len(set(traits)):
        raise ValueError(f"{label} must contain unique, non-empty trait values")
    return rows


def read_iv_targets(path: Path) -> dict[bytes, str]:
    rows = read_tsv(path)
    required = {"rsid", "variant_id"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"{path}: missing columns {sorted(required)}")
    variants = [row["variant_id"].strip() for row in rows]
    rsids = [row["rsid"].strip() for row in rows]
    if (any(not value for value in variants + rsids) or
            len(variants) != len(set(variants)) or len(rsids) != len(set(rsids))):
        raise ValueError(f"{path}: rsid and variant_id must be complete and unique")
    return {variant.encode(): rsid for variant, rsid in zip(variants, rsids, strict=True)}


def extract_one(job: tuple[str, str, str, int]) -> tuple[list[dict[str, str]], dict[str, str]]:
    trait, source_name, iv_name, max_n = job
    source, iv_path = Path(source_name), Path(iv_name)
    targets = read_iv_targets(iv_path)
    found: dict[bytes, int] = {}
    started = time.perf_counter()
    rows_scanned = 0
    opener = gzip.open if source.suffix.lower() == ".gz" else open
    with opener(source, "rb") as handle:
        header = handle.readline().decode("utf-8-sig").rstrip("\r\n").split("\t")
        try:
            variant_idx, n_idx = header.index("variant_id"), header.index("n")
        except ValueError as exc:
            raise ValueError(f"{source}: variant_id/n column missing") from exc
        last_idx = max(variant_idx, n_idx)
        for line_number, line in enumerate(handle, start=2):
            rows_scanned += 1
            parts = line.rstrip(b"\r\n").split(b"\t", last_idx + 1)
            if len(parts) <= last_idx:
                raise ValueError(f"{source}:{line_number}: truncated row")
            variant = parts[variant_idx]
            if variant not in targets:
                continue
            try:
                value = int(parts[n_idx])
            except ValueError as exc:
                raise ValueError(f"{source}:{line_number}: invalid n for {variant!r}") from exc
            if value <= 0 or value > max_n:
                raise ValueError(f"{source}:{line_number}: out-of-range n={value} for {variant!r}")
            if variant in found:
                raise ValueError(f"{source}:{line_number}: duplicate target variant_id {variant!r}")
            found[variant] = value
    missing = sorted(set(targets) - set(found))
    if missing:
        preview = ", ".join(value.decode(errors="replace") for value in missing[:5])
        raise ValueError(f"{trait}: {len(missing)} IV variant_id values not found; first: {preview}")
    output = [
        {
            "trait": trait,
            "rsid": targets[variant],
            "variant_id": variant.decode(),
            "n_exposure": str(found[variant]),
            "source_gwas_file": source.resolve().as_posix(),
            "source_n_field": "n",
        }
        for variant in sorted(found)
    ]
    values = list(found.values())
    qc = {
        "trait": trait,
        "target_iv_n": str(len(targets)),
        "matched_iv_n": str(len(found)),
        "missing_iv_n": "0",
        "n_min": str(min(values)),
        "n_max": str(max(values)),
        "n_unique": str(len(set(values))),
        "n_values": ";".join(map(str, sorted(set(values)))),
        "source_rows_scanned": str(rows_scanned),
        "elapsed_seconds": f"{time.perf_counter() - started:.3f}",
        "source_gwas_file": source.resolve().as_posix(),
    }
    return output, qc


def write_tsv_atomic(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exposure-manifest", type=Path, required=True,
                        help="TSV columns: trait,source_file")
    parser.add_argument("--iv-manifest", type=Path, required=True,
                        help="TSV columns: trait,iv_file")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--qc-out", type=Path)
    parser.add_argument("--traits", help="Optional comma-separated trait names for a smoke test")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--max-n", type=int, default=599249)
    args = parser.parse_args()
    if not 1 <= args.workers <= 8:
        parser.error("--workers must be between 1 and 8")
    if args.max_n < 1:
        parser.error("--max-n must be positive")

    exposures = require_manifest(
        args.exposure_manifest, {"trait", "source_file"}, "exposure manifest"
    )
    instruments = require_manifest(args.iv_manifest, {"trait", "iv_file"}, "IV manifest")
    source_by_trait = {row["trait"].strip(): row["source_file"].strip() for row in exposures}
    iv_by_trait = {row["trait"].strip(): row["iv_file"].strip() for row in instruments}
    selected = ({value.strip() for value in args.traits.split(",") if value.strip()}
                if args.traits else set(iv_by_trait))
    if not selected:
        parser.error("no traits selected")
    missing = sorted(selected - set(source_by_trait) | selected - set(iv_by_trait))
    if missing:
        raise ValueError(f"selected traits are absent from a manifest: {missing}")

    jobs: list[tuple[str, str, str, int]] = []
    for trait in sorted(selected):
        source, iv_path = Path(source_by_trait[trait]), Path(iv_by_trait[trait])
        if not source.is_file() or not iv_path.is_file():
            raise FileNotFoundError(f"{trait}: missing source or IV file: {source}; {iv_path}")
        jobs.append((trait, str(source), str(iv_path), args.max_n))

    all_rows: list[dict[str, str]] = []
    qc_rows: list[dict[str, str]] = []
    with ProcessPoolExecutor(max_workers=min(args.workers, len(jobs))) as pool:
        futures = {pool.submit(extract_one, job): job[0] for job in jobs}
        for future in as_completed(futures):
            rows, qc = future.result()
            all_rows.extend(rows)
            qc_rows.append(qc)
            print(f"{qc['trait']}: {qc['matched_iv_n']}/{qc['target_iv_n']} IVs; "
                  f"N={qc['n_values']}; {qc['elapsed_seconds']} s", flush=True)
    all_rows.sort(key=lambda row: (row["trait"], row["rsid"], row["variant_id"]))
    qc_rows.sort(key=lambda row: row["trait"])
    qc_out = args.qc_out or args.out.with_name(f"{args.out.stem}_qc.tsv")
    write_tsv_atomic(args.out, all_rows)
    write_tsv_atomic(qc_out, qc_rows)
    print(f"wrote {len(all_rows)} IV rows for {len(qc_rows)} traits to {args.out}")


if __name__ == "__main__":
    main()
