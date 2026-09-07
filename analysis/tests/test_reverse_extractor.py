#!/usr/bin/env python3
"""Synthetic fail-closed checks for reverse-association input keys."""
from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "analysis" / "scripts" / "00_prepare_reverse_associations.py"


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def expect_failure(work: Path, mdd: list[dict[str, object]], coords: list[dict[str, object]],
                   manifest: list[dict[str, object]], message: str) -> None:
    mdd_path = work / "mdd.tsv"
    coords_path = work / "coords.tsv"
    manifest_path = work / "manifest.tsv"
    write_tsv(mdd_path, mdd)
    write_tsv(coords_path, coords)
    write_tsv(manifest_path, manifest)
    run = subprocess.run(
        [sys.executable, str(SCRIPT), "--mdd-iv", str(mdd_path), "--coords", str(coords_path),
         "--exposure-manifest", str(manifest_path), "--out", str(work / "out.tsv")],
        cwd=ROOT, capture_output=True, text=True,
    )
    assert run.returncode != 0 and message in run.stderr, run.stderr


def main() -> None:
    mdd = [{"rsid": "rs1", "A1": "A", "A2": "G", "beta": 0.1, "se": 0.02,
            "p": 1e-8, "eaf": 0.2}]
    coords = [{"rsid": "rs1", "chr": 1, "pos": 100}]
    manifest = [{"trait": "Trait_1", "source_file": "missing.tsv"}]
    with tempfile.TemporaryDirectory(prefix="reverse_key_test_") as tmp:
        work = Path(tmp)
        expect_failure(work, mdd + mdd, coords, manifest, "MDD IV file must contain unique rsid values")
        expect_failure(work, mdd, coords + coords, manifest,
                       "coordinate table must contain unique rsid values")
        expect_failure(work, mdd, coords, manifest + manifest,
                       "exposure manifest must contain unique trait values")
    print("PASS: reverse extractor rejects duplicate MDD, coordinate and trait keys")


if __name__ == "__main__":
    main()
