#!/usr/bin/env python3
"""Small deterministic test for exposure variant_id-to-rsID mapping."""
from __future__ import annotations

import gzip
import importlib.util
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "analysis" / "scripts" / "00_map_exposure_rsids.py"
SPEC = importlib.util.spec_from_file_location("exposure_rsid_mapping", SCRIPT)
assert SPEC and SPEC.loader
MAPPING = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MAPPING)


def main() -> None:
    with tempfile.TemporaryDirectory() as temporary:
        base = Path(temporary)
        source = base / "fixture.tsv.gz"
        header = [
            "variant_id", "effect_allele", "other_allele", "beta", "standard_error",
            "effect_allele_frequency", "neg_log_10_p_value", "n",
        ]
        rows = [
            ("1_100_A_G", "G", "A", "0.1", "0.01", "0.4", "9", "599249"),
            ("1_101_C_T", "T", "C", "0.1", "0.01", "0.4", "9", "599249"),
            ("1_102_G_A", "A", "G", "0.1", "0.01", "0.4", "9", "599249"),
            ("1_103_A_C", "C", "A", "0.1", "0.01", "0.4", "9", "599249"),
            ("1_104_T_C", "C", "T", "0.1", "0.01", "0.4", "9", "599249"),
            ("1_105_-_T", "T", "-", "0.1", "0.01", "0.4", "9", "599249"),
            ("1_106_A_G", "G", "A", "0.1", "0.01", "0.4", "2", "599249"),
            ("1_107_A_G", "G", "A", "0.1", "0.01", "0.0005", "10", "599249"),
        ]
        with gzip.open(source, "wt", encoding="utf-8", newline="") as handle:
            handle.write("\t".join(header) + "\n")
            for row in rows:
                handle.write("\t".join(row) + "\n")

        ordered, variants, qc = MAPPING.collect_candidates(source, 5e-8)
        assert ordered == [row[0] for row in rows[:6]] + [rows[7][0]]
        assert qc["source_rows_scanned"] == 8 and qc["candidate_rows"] == 7
        records = {
            (MAPPING.CONTIGS["1"], 100): [("A", "G", "rs100")],
            (MAPPING.CONTIGS["1"], 101): [("T", "C", "rs101")],
            (MAPPING.CONTIGS["1"], 102): [("G", "A", "rs102;rs202")],
            (MAPPING.CONTIGS["1"], 103): [("A", "T", "rs103")],
            (MAPPING.CONTIGS["1"], 107): [("A", "G", "rs107")],
        }
        mapped = MAPPING.assign_rsids(variants, records)
        assert mapped["1_100_A_G"] == ("rs100", "matched")
        assert mapped["1_101_C_T"] == ("rs101", "matched")
        assert mapped["1_102_G_A"] == ("rs102,rs202", "ambiguous")
        assert mapped["1_103_A_C"] == (".", "no_allele_match")
        assert mapped["1_104_T_C"] == (".", "no_dbsnp_entry")
        assert mapped["1_105_-_T"] == (".", "skipped")

        output = base / "maps"
        qc_rows = MAPPING.write_outputs(output, {"Fixture": ordered}, mapped, {"Fixture": qc})
        MAPPING.write_qc(output, qc_rows)
        assert qc_rows[0]["matched"] == 3 and qc_rows[0]["match_fraction"] == "0.428571428571"
        lines = (output / "Fixture_rsid_map.tsv").read_text(encoding="utf-8").splitlines()
        assert len(lines) == 8 and lines[0] == "variant_id\trsid\trsid_status"
        assert (output / "exposure_rsid_mapping_qc.tsv").read_text(encoding="utf-8").count("\n") == 2
    print("PASS: exposure rsID candidate selection, allele matching and status output")


if __name__ == "__main__":
    main()
