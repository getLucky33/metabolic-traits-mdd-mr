#!/usr/bin/env python3
"""Small dependency-free test for source preflight and manifest construction."""
from __future__ import annotations

import csv
import gzip
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def write_gzip(path: Path, columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(columns) + "\n")


def main() -> None:
    catalog = list(csv.DictReader((ROOT / "analysis/manifests/metabolic_traits_249.tsv").open(encoding="utf-8"), delimiter="\t"))
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp)
        metabolic = base / "metabolic"
        maps = base / "maps"
        metabolic_columns = [
            "variant_id", "effect_allele", "other_allele", "beta", "standard_error",
            "effect_allele_frequency", "neg_log_10_p_value",
        ]
        for index, row in enumerate(catalog):
            write_gzip(metabolic / row["accession"] / f"{row['accession']}.tsv.gz", metabolic_columns)
            maps.mkdir(parents=True, exist_ok=True)
            suffix = "_rsid_map.tsv" if index % 2 == 0 else "_forward_ivs_rsid_b157.tsv"
            (maps / f"{row['trait']}{suffix}").write_text("variant_id\trsid\trsid_status\n", encoding="utf-8", newline="\n")

        pgc_columns = ["CHR", "BP", "SNP", "A1", "A2", "OR", "SE", "P", "Nca", "Nco", "FRQ_A_test", "FRQ_U_test"]
        fg_columns = ["#chrom", "pos", "ref", "alt", "rsids", "nearest_genes", "pval", "mlogp", "beta", "sebeta", "af_alt", "af_alt_cases", "af_alt_controls"]
        write_gzip(base / "pgc_main.gz", pgc_columns)
        write_gzip(base / "pgc_noukbb.gz", pgc_columns)
        write_gzip(base / "finngen.gz", fg_columns)
        for name in ("plink", "bcftools", "dbsnp.gz", "dbsnp.gz.tbi", "EUR.bed", "EUR.bim", "EUR.fam"):
            (base / name).write_text("test\n", encoding="utf-8", newline="\n")

        config = base / "local.tsv"
        values = {
            "catalog": ROOT / "analysis/manifests/metabolic_traits_249.tsv",
            "metabolic_download_dir": metabolic,
            "rsid_map_dir": maps,
            "pgc_main": base / "pgc_main.gz",
            "pgc_noukbb": base / "pgc_noukbb.gz",
            "finngen": base / "finngen.gz",
            "ld_bfile": base / "EUR",
            "plink": base / "plink",
            "bcftools": base / "bcftools",
            "dbsnp": base / "dbsnp.gz",
            "dbsnp_index": base / "dbsnp.gz.tbi",
        }
        with config.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(("key", "path"))
            writer.writerows(values.items())

        command = [sys.executable, str(ROOT / "analysis/prepare_inputs.py")]
        subprocess.run(command + ["check", "--config", str(config), "--level", "sources"], cwd=ROOT, check=True)
        subprocess.run(command + ["check", "--config", str(config), "--level", "analysis"], cwd=ROOT, check=True)
        out_dir = base / "manifests"
        subprocess.run(command + ["build-manifests", "--config", str(config), "--out-dir", str(out_dir)], cwd=ROOT, check=True)
        exposures = list(csv.DictReader((out_dir / "exposure_manifest.local.tsv").open(encoding="utf-8"), delimiter="\t"))
        sources = list(csv.DictReader((out_dir / "source_manifest.local.tsv").open(encoding="utf-8"), delimiter="\t"))
        assert len(exposures) == 249 and len(sources) == 250
        assert sources[-1]["trait"] == "Major_depressive_disorder" and sources[-1]["source_type"] == "mdd"
    print("PASS: source preflight and full-manifest construction")


if __name__ == "__main__":
    main()
