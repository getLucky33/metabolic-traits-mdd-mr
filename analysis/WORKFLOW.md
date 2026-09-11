# Analysis workflow

Run all commands from the repository root. Paths are supplied only through command-line arguments or TSV manifests; the scripts contain no machine-specific data paths and never download or install packages at run time.

Release v0.2.13 adds an operational data checklist, the complete 249-accession map and fail-closed source/input preflight. It does not change the v0.2.12 analysis results or Figure 1 topology. Forward and reverse MR remain parallel direction-specific analyses with separately selected instruments. The 15-trait candidate set is selected only by the primary forward screen; matched reverse-MR estimates inform the directional assessment but do not enter the robustness or locus branches and are not interpreted as proof of reciprocal causality.

## 0. Establish the reproduction boundary

Read `DATA_SOURCES.md`, copy `analysis/config/local_paths.example.tsv` to an ignored `.local.tsv` file and run `analysis/prepare_inputs.py check --level sources`. This checks all 249 accession files, the two PGC outcomes, FinnGen, the PLINK LD prefix and the indexed dbSNP resource without reading the complete GWAS bodies.

The downloaded files are not yet analysis-ready. Instrument selection also requires one prepared `variant_id`-to-rsID map per metabolic trait. After those maps exist, run the `analysis` check and build the complete local manifests. The current repository consumes and validates the maps but does not regenerate the complete frozen mapping set.

Do not describe the integrated-evidence or exploratory SuSiE tables as independently executable outputs. They remain released aggregate audit records because the necessary regional association data, LD matrices and credible-set inputs are not public.

## 1. Environment and local-only data

Install R 4.5.1, PLINK 1.9 and bcftools 1.24. Restore the complete hard-dependency closure recorded in `renv.lock`; the GitHub packages are pinned to immutable commits. `analysis/environment/package-versions.tsv` records the principal analysis packages and command-line tools. Create `data-local/`; this directory is excluded by `.gitignore`.

Use `analysis/prepare_inputs.py build-manifests` to create the complete local source and exposure manifests after the analysis-level preflight passes. Build the IV manifest from stage 1 outputs, and copy only the remaining outcome, locus and reverse-tier examples when those stages are reached. Keep trait identifiers stable across all manifests. Full journal-facing names are used only in tables and figures; analysis identifiers must not be renamed after the run is frozen.

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

When no `--traits` subset is supplied, the script requires exactly 249 unique traits and writes `mr_main_summary_screen.tsv` with `p_fdr_bh`, fixed-family `p_bonf` and the four-level `screen_level`. For a subset run, `pilot_fdr` remains descriptive and is not promoted to the locked full-family classification. Formal candidates must be rerun with at least 10,000 MR-PRESSO simulations.

After instrument selection, generate the trait-level strength and NOME diagnostics:

```bash
Rscript analysis/scripts/09_summarise_instrument_qc.R \
  --iv-manifest analysis/manifests/iv_manifest.local.tsv \
  --screen results/analysis/forward/mr_main_summary_screen.tsv \
  --selection-qc data-local/ivs/instrument_selection_qc.tsv \
  --trait-labels data/derived/trait_display_dictionary.tsv \
  --out results/analysis/instrument_strength_qc.tsv
```

I²GX uses the inverse-variance-weighted SNP-exposure heterogeneity formula documented in the README. If the selection-stage LD-panel denominator is unavailable, the script writes `NA` and an explicit reason; it does not infer missingness from the final instruments.

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

Without `--traits`, each reverse run also requires exactly 249 traits and writes `reverse_mr_summary_screen.tsv` with the locked four-level classification.

## 5. Robustness analyses

```bash
# Required one-trait pilot before the formal run
Rscript analysis/scripts/05_run_robustness.R \
  --harm-dir data-local/harmonised \
  --out-dir results/analysis/robustness-pilot \
  --traits Trait_1 --pilot \
  --presso-nb 500 --presso-seed 20260815 \
  --prevalences 0.08,0.15,0.20

# Formal 15-trait run
Rscript analysis/scripts/05_run_robustness.R \
  --harm-dir data-local/harmonised \
  --iv-manifest analysis/manifests/iv_manifest.local.tsv \
  --out-dir results/analysis/robustness \
  --traits Trait_1,...,Trait_15 \
  --presso-nb 10000 --presso-seed 20260815 \
  --prevalences 0.08,0.15,0.20
```

Formal mode requires exactly 15 candidates, at least 10,000 MR-PRESSO simulations and a complete 249-trait IV manifest for the shared-instrument count. `--pilot` permits a smaller smoke test and prefixes its outputs with `pilot_`; a smaller simulation count never produces a file named `presso_rerun_10000.tsv`.

The Steiger calculation assigns every exposure SNP `N=599,249`, the maximum exposure meta-analysis sample size, because per-SNP effective N was not retained in the frozen harmonized inputs. This approximation tends to reduce exposure-side R² when the true SNP-specific N is smaller, but it can still change Steiger P values or near-boundary direction calls. The prevalence-grid results are sensitivity evidence and do not rule out clinical reverse causation.

## 6. FinnGen comparison

```bash
# Required pilot; pilot files do not contain a formal BH support label
Rscript analysis/scripts/06_run_finngen.R \
  --finngen data-local/finngen/finngen_R13_F5_DEPRESSIO.gz \
  --iv-manifest analysis/manifests/iv_manifest.local.tsv \
  --pgc-summary results/analysis/forward/mr_main_summary.tsv \
  --out-dir results/analysis/finngen-pilot \
  --traits Trait_1

# Formal run; use a manifest containing exactly the 15 candidates
Rscript analysis/scripts/06_run_finngen.R \
  --finngen data-local/finngen/finngen_R13_F5_DEPRESSIO.gz \
  --iv-manifest analysis/manifests/iv_manifest.candidates15.local.tsv \
  --pgc-summary results/analysis/forward/mr_main_summary.tsv \
  --out-dir results/analysis/finngen \
  --n-tests 15
```

The script reads FinnGen once, requires valid finite beta/SE/P/AF values, fails on duplicate matches and reuses the PGC instruments rather than selecting new FinnGen-specific instruments. Formal mode requires 15 unique candidate rows, one PGC row and a successful FinnGen MR result for each candidate before it computes the 15-test BH family. Any incomplete trait is written to the status/tool-loss files and aborts formal classification. A `--traits` subset is exploratory and writes only `pilot_*` files without `p_bh` or `cross_outcome_support_label`.

The primary cross-outcome field is `cross_outcome_support_label`, with values `FDR-supported`, `nominally-supported` and `no-nominal-support`. It depends on direction concordance and the FinnGen P value only. The two-study Q and I² fields are named `overlap_uncorrected_*` and are descriptive because the outcome datasets are not independent.

## 7. Colocalization

Build the locus manifest from the complete forward screen, final exposure IVs, the two GRCh38 MDD coordinate sets and the frozen reverse-tier table. The anchor rule joins a lead only when it is within 1 Mb of the first lead in the cluster, then adds a 500-kb flank.

```bash
Rscript analysis/scripts/10_build_locus_manifest.R \
  --forward-screen results/analysis/forward/mr_main_summary_screen.tsv \
  --iv-manifest analysis/manifests/iv_manifest.local.tsv \
  --mdd-main-coords data-local/reverse-main/mdd_iv_coords.tsv \
  --mdd-noukbb-coords data-local/reverse-noukbb/mdd_iv_coords.tsv \
  --r-tier analysis/manifests/r_tier.local.tsv \
  --out data-local/coloc/locus_manifest.tsv \
  --window-dir data-local/coloc/windows \
  --expected-loci 6434

Rscript analysis/scripts/11_extract_locus_windows.R \
  --locus-manifest data-local/coloc/locus_manifest.tsv \
  --source-manifest analysis/manifests/source_manifest.local.tsv \
  --out-dir data-local/coloc/windows \
  --traits GlycA
```

The extraction command above is the required one-trait smoke test. Remove `--traits` only after checking its row counts and timing. The extractor groups loci by trait, reads each large GWAS source once, joins the prepared unique rsID map and then writes the regional windows. Regional exposure-window files must contain `rsid`, `beta`, `se`, `eaf`, `effect_allele` and `other_allele`.

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

The classification can also be replayed from the released aggregate inputs without SNP-level data:

```bash
Rscript analysis/scripts/08_classify_coloc.R \
  --locus-manifest data/derived/coloc_locus_manifest.tsv \
  --results data/derived/coloc_abf_results.tsv \
  --status data/derived/coloc_status.tsv \
  --lead-qc data/derived/coloc_lead_qc.tsv \
  --out-dir results/analysis/coloc-replay
```

For the 6,434-row frozen manifest this command fails unless all five class counts equal `14/87/822/3430/2081`. MHC exclusion uses interval overlap between `window_start/window_end` and GRCh38 chr6:25–34 Mb, including boundary contact.

## 8. Verification

Before a full run, execute the synthetic smoke test and a one-trait/one-locus real-data smoke test. Compare the one-trait IVW estimate to the frozen result before launching the full family. Confirm that each large GWAS is read once, inspect the first output within 2–3 minutes and stop the run if observed timing materially exceeds the smoke-test estimate.

Every stage writes input hashes and parameter tables alongside its outputs. Do not edit an output after hashing; rerun its producing stage instead.

The repository releases aggregate integrated-evidence and exploratory SuSiE tables for audit, but not a clean-clone computation of those stages. The frozen computation used regional SNP data, LD matrices and author credible-set inputs that are not redistributed. Do not describe those aggregate snapshots as an independently executable SuSiE pipeline.
