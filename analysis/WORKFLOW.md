# Analysis workflow

Run all commands from the repository root. Paths are supplied only through command-line arguments or TSV manifests; the scripts contain no machine-specific data paths and never download or install packages at run time.

## 1. Environment and local-only data

Install R 4.5.1, PLINK 1.9 and bcftools 1.24. Restore `renv.lock`, or install the exact versions in `analysis/environment/package-versions.tsv`. Create `data-local/`; this directory is excluded by `.gitignore`.

Copy the examples in `analysis/manifests/` and replace their placeholder paths. Keep trait identifiers stable across all manifests. Full journal-facing names are used only in tables and figures; analysis identifiers must not be renamed after the run is frozen.

## 2. Instrument selection

`01_select_instruments.R` reads each GWAS once, applies the locked significance/MAF rules, filters to autosomal biallelic variants represented in the European LD panel, runs local PLINK clumping and removes variants with F<10.

```bash
Rscript analysis/scripts/01_select_instruments.R \
  --manifest analysis/manifests/source_manifest.local.tsv \
  --ld-bfile data-local/ld/1kg.v3/EUR \
  --plink data-local/bin/plink \
  --out-dir data-local/ivs
```

For metabolic-trait sources, `rsid_map_file` must contain `variant_id`, `rsid` and `rsid_status`. Only unique `matched` rows enter clumping. The stage writes per-trait MAF, F, approximate standardized-trait R², LD-panel missingness and source hashes.

## 3. Forward MR

Build `iv_manifest.local.tsv` from the final IV files, then harmonize against both PGC outcomes. Each large PGC file is read once per run.

```bash
Rscript analysis/scripts/02_harmonise_forward.R \
  --iv-manifest analysis/manifests/iv_manifest.local.tsv \
  --mdd-main data-local/mdd/daner_main.tsv.gz \
  --mdd-noukbb data-local/mdd/daner_no_ukbb.tsv.gz \
  --out-dir data-local/harmonised

Rscript analysis/scripts/03_run_forward_mr.R \
  --harm-dir data-local/harmonised \
  --out-dir results/analysis/forward \
  --n-tests 249 --presso-nb 1000 --presso-seed 20260815
```

The `pilot_fdr` column is descriptive for the executed subset and is not promoted to the locked full-family significance classification. Formal candidates must be rerun with at least 10,000 MR-PRESSO simulations.

## 4. Reverse MR

Select MDD instruments with stage 1 using `source_type=mdd`. Prepare the MDD-instrument × metabolic-trait association table. The Python extractor accepts a precomputed GRCh38 coordinate map or creates one by scanning dbSNP once with bcftools.

```bash
python analysis/scripts/00_prepare_reverse_associations.py \
  --mdd-iv data-local/ivs/Major_depressive_disorder_final_ivs.tsv \
  --coords data-local/maps/mdd_iv_grch38.tsv \
  --exposure-manifest analysis/manifests/exposure_manifest.local.tsv \
  --out data-local/reverse/mdd_iv_metabolite_assoc.tsv

Rscript analysis/scripts/04_run_reverse_mr.R \
  --assoc data-local/reverse/mdd_iv_metabolite_assoc.tsv \
  --out-dir results/analysis/reverse \
  --n-tests 249 --presso-nb 1000 --presso-seed 20260815
```

Repeat the two commands with the UK Biobank-excluded MDD instrument file and a separate output directory. Do not pool the main and UK Biobank-excluded instruments.

## 5. Robustness analyses

```bash
Rscript analysis/scripts/05_run_robustness.R \
  --harm-dir data-local/harmonised \
  --iv-manifest analysis/manifests/iv_manifest.local.tsv \
  --out-dir results/analysis/robustness \
  --traits Trait_1,Trait_2 \
  --presso-nb 10000 --presso-seed 20260815 \
  --prevalences 0.08,0.15,0.20
```

The trait list must be the locked formal-candidate set. Steiger results assess variance direction under the prevalence grid; they do not rule out clinical reverse causation.

## 6. FinnGen comparison

```bash
Rscript analysis/scripts/06_run_finngen.R \
  --finngen data-local/finngen/finngen_R13_F5_DEPRESSIO.gz \
  --iv-manifest analysis/manifests/iv_manifest.local.tsv \
  --pgc-summary results/analysis/forward/mr_main_summary.tsv \
  --out-dir results/analysis/finngen \
  --n-tests 15
```

The script reads FinnGen once, requires valid finite beta/SE/P/AF values, fails on duplicate matches and reuses the PGC instruments rather than selecting new FinnGen-specific instruments.

## 7. Colocalization

Regional exposure-window files must contain `rsid`, `beta`, `se`, `eaf`, `effect_allele` and `other_allele`. Prepare the local locus and outcome manifests from the examples.

```bash
Rscript analysis/scripts/07_run_coloc_abf.R \
  --locus-manifest analysis/manifests/locus_manifest.local.tsv \
  --outcome-manifest analysis/manifests/outcome_manifest.local.tsv \
  --out-dir results/analysis/coloc \
  --min-snps 50

Rscript analysis/scripts/08_classify_coloc.R \
  --locus-manifest analysis/manifests/locus_manifest.local.tsv \
  --results results/analysis/coloc/loci_coloc_results.tsv \
  --status results/analysis/coloc/loci_status.tsv \
  --final-snps results/analysis/coloc/loci_final_snps.tsv \
  --exposure-lead-qc results/analysis/coloc/exposure_lead_qc.tsv \
  --plink data-local/bin/plink \
  --ld-bfile data-local/ld/1kg.v3/EUR \
  --out-dir results/analysis/coloc
```

Use `--limit 1` or `--limit 2` for the required smoke test before a full regional run. The ABF stage records every locus-outcome failure rather than silently dropping it. Classification is a manifest-left-join operation and requires one final row per locus.

## 8. Verification

Before a full run, execute the synthetic smoke test and a one-trait/one-locus real-data smoke test. Compare the one-trait IVW estimate to the frozen result before launching the full family. Confirm that each large GWAS is read once, inspect the first output within 2–3 minutes and stop the run if observed timing materially exceeds the smoke-test estimate.

Every stage writes input hashes and parameter tables alongside its outputs. Do not edit an output after hashing; rerun its producing stage instead.
