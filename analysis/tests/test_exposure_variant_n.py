#!/usr/bin/env python3
"""Synthetic gzip test for exposure-side variant sample-size recovery."""
from __future__ import annotations

import csv
import gzip
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "analysis" / "scripts" / "05a_recover_exposure_variant_n.py"


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def write_gzip(path: Path, rows: list[dict[str, object]]) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="variant_n_test_") as tmp:
        work = Path(tmp)
        source_a, source_b = work / "a.tsv.gz", work / "b.tsv.gz"
        iv_a, iv_b = work / "a_iv.tsv", work / "b_iv.tsv"
        write_gzip(source_a, [
            {"variant_id": "1:90:A:G", "n": 599249, "beta": 0.01},
            {"variant_id": "1:100:A:C", "n": 413897, "beta": 0.02},
            {"variant_id": "1:200:G:T", "n": 599249, "beta": -0.03},
        ])
        write_gzip(source_b, [
            {"variant_id": "2:100:C:T", "n": 599249, "beta": 0.04},
            {"variant_id": "2:110:A:G", "n": 413897, "beta": 0.05},
        ])
        write_tsv(iv_a, [
            {"rsid": "rs2", "variant_id": "1:200:G:T"},
            {"rsid": "rs1", "variant_id": "1:100:A:C"},
        ])
        write_tsv(iv_b, [{"rsid": "rs3", "variant_id": "2:100:C:T"}])
        exposure_manifest = work / "exposure.tsv"
        iv_manifest = work / "ivs.tsv"
        write_tsv(exposure_manifest, [
            {"trait": "Trait_B", "source_file": source_b},
            {"trait": "Trait_A", "source_file": source_a},
        ])
        write_tsv(iv_manifest, [
            {"trait": "Trait_B", "iv_file": iv_b},
            {"trait": "Trait_A", "iv_file": iv_a},
        ])

        output, qc_output = work / "variant_n.tsv", work / "variant_n_qc.tsv"
        command = [
            sys.executable, str(SCRIPT),
            "--exposure-manifest", str(exposure_manifest),
            "--iv-manifest", str(iv_manifest),
            "--out", str(output), "--qc-out", str(qc_output),
            "--workers", "2", "--max-n", "599249",
        ]
        subprocess.run(command, cwd=ROOT, check=True, capture_output=True, text=True)
        recovered, qc = read_tsv(output), read_tsv(qc_output)
        assert [(row["trait"], row["rsid"], row["n_exposure"]) for row in recovered] == [
            ("Trait_A", "rs1", "413897"),
            ("Trait_A", "rs2", "599249"),
            ("Trait_B", "rs3", "599249"),
        ]
        assert [row["trait"] for row in qc] == ["Trait_A", "Trait_B"]
        assert sum(int(row["target_iv_n"]) for row in qc) == 3
        assert sum(int(row["matched_iv_n"]) for row in qc) == 3
        assert sum(int(row["missing_iv_n"]) for row in qc) == 0
        assert qc[0]["n_values"] == "413897;599249"

        failed_command = command.copy()
        failed_command[failed_command.index("599249")] = "500000"
        failed = subprocess.run(
            failed_command, cwd=ROOT, capture_output=True, text=True,
        )
        assert failed.returncode != 0 and "out-of-range n=599249" in failed.stderr
    print("PASS: variant-specific exposure N recovered once from synthetic gzip sources")


if __name__ == "__main__":
    main()
