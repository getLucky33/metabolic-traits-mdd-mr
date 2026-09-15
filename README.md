# Circulating metabolic traits and major depressive disorder

[![reproduce](https://github.com/getLucky33/metabolic-traits-mdd-mr/actions/workflows/reproduce.yml/badge.svg)](https://github.com/getLucky33/metabolic-traits-mdd-mr/actions/workflows/reproduce.yml)

Reproducibility code for the Scientific Reports manuscript by Zhouyi Wang and Qingmei Liu. Version 0.2.19 completes the 249-trait LD-panel selection-stage QC and propagates the observed missing counts, eligible counts and missing fractions into the instrument-strength table. MR estimates and inferential classifications are unchanged.

## Choose a reproduction target

| Target | What can be reproduced | Additional inputs |
| --- | --- | --- |
| Clean-clone verification | Figures 1–4, the complete 249-trait screen counts, the frozen 6,434-record colocalization classification and synthetic method checks | None beyond the recorded software environment |
| Prepared-input analysis | Exposure rsID mapping, instrument selection, forward and reverse MR, sensitivity analyses, the FinnGen alternative-outcome comparison, locus construction, ABF colocalization and classification | Third-party GWAS, LD and dbSNP resources plus the remaining prepared inputs listed in [`DATA_SOURCES.md`](DATA_SOURCES.md) |
| Aggregate-result audit only | The 101-row integrated-evidence table and exploratory SuSiE status, binding and stopping summaries | Independent recomputation is not possible from the released files alone |

This is not a one-command raw-data-to-paper repository. Exposure rsID maps can now be generated locally, but later stages still require selected MDD instruments, coordinate tables and regional inputs. The exact boundary is stated before the commands so that a successful figure build is not mistaken for a complete numerical rerun.

## Start here for a real-data rerun

1. Read [`DATA_SOURCES.md`](DATA_SOURCES.md) for exact accessions, filenames, genome builds, required columns and prepared-input dependencies.
2. Copy the local-path template and replace its paths. Files ending in `.local.tsv` are ignored by Git.
3. Run the source preflight before starting any long analysis.

```bash
cp analysis/config/local_paths.example.tsv analysis/config/local_paths.local.tsv
python analysis/prepare_inputs.py check \
  --config analysis/config/local_paths.local.tsv \
  --level sources
```

Generate one map first as the required real-data smoke test. The script reads that GWAS once, verifies the frozen dbSNP build-157 MD5, queries indexed positions and writes allele-matching QC. Check the observed row count, status counts and runtime before removing `--traits` for the complete 249-trait run.

```bash
python analysis/scripts/00_map_exposure_rsids.py \
  --config analysis/config/local_paths.local.tsv \
  --traits Cholesterol_in_very_large_HDL \
  --out-dir data-local/rsid-smoke

python analysis/scripts/00_map_exposure_rsids.py \
  --config analysis/config/local_paths.local.tsv \
  --out-dir data-local/maps
```

The mapping stage admits the genome-wide-significant (`P<5×10⁻⁸`) source rows as an upstream superset. The MAF-aware common/rare thresholds remain in instrument selection. Full runs process eight traits per indexed dbSNP query by default, so each source GWAS is streamed once without retaining all 249 traits in memory.

After the maps exist, validate the analysis-ready inputs and build the complete local manifests:

```bash
python analysis/prepare_inputs.py check \
  --config analysis/config/local_paths.local.tsv \
  --level analysis
python analysis/prepare_inputs.py build-manifests \
  --config analysis/config/local_paths.local.tsv \
  --out-dir analysis/manifests
```

Then follow [`analysis/WORKFLOW.md`](analysis/WORKFLOW.md). The workflow retains explicit one-trait and one-locus smoke-test stops; a single unattended runner is intentionally not provided because it would bypass those run-time and scientific checks.

## Clean-clone verification

Install R 4.5.1 and restore the complete environment recorded in `renv.lock`. `analysis/environment/package-versions.tsv` is a concise runtime summary, not a substitute for the lockfile.

```bash
Rscript scripts/make_figures.R
Rscript analysis/scripts/08_classify_coloc.R \
  --locus-manifest data/derived/coloc_locus_manifest.tsv \
  --results data/derived/coloc_abf_results.tsv \
  --status data/derived/coloc_status.tsv \
  --lead-qc data/derived/coloc_lead_qc.tsv \
  --out-dir results/analysis/coloc-replay
Rscript analysis/tests/smoke_test.R
```

`make_figures.R` verifies the aggregate-input hashes and regenerates PDF, SVG and 600-dpi PNG figures. `08_classify_coloc.R` must recover `14/87/822/3430/2081` records in the five frozen classes. `smoke_test.R` uses temporary synthetic data to exercise TwoSampleMR harmonization, random-effects IVW, MR-PRESSO parsing and the four-prior `coloc.abf` path. Synthetic results are never mixed with manuscript results.

Use `Rscript scripts/make_figures.R --out <directory>` to choose another output directory. The run also writes the exact 249-trait Figure 2 screen, the 45 candidate estimates, the Figure 3 evidence matrix, the Figure 4 classification, trait-distribution and integrated-disposition counts, and a checksum file beside the figures. Figure 1 uses neutral GWAS-source headings and presents forward and reverse MR as parallel direction-specific analyses with separately selected instruments. Only the primary forward screen selects the 15-trait candidate set; those candidates enter the robustness, locus and directional assessments, whereas matched reverse-MR results enter only the directional assessment. The layout does not treat reverse MR as validation or proof of reciprocal causality. Figure 2 first shows the complete primary forward screen and then the primary PGC, PGC excluding UK Biobank and FinnGen R13 estimates for each candidate. Figure 3 groups the recorded evidence by analytical role and labels its MR-PRESSO column as distortion-only. Figure 4 recomputes the five ABF classifications across 6,434 analysis-window records, explicitly separates the 101-record subset into 14 robust and 87 prior-sensitive records, and reports their four integrated dispositions and 0/101 endpoint.

On each push or pull request, the configured GitHub Actions workflow runs the repository audit, the synthetic input-preflight test, figure reproduction, code parsing, a secret/path scan and the synthetic MR/colocalization smoke test in a clean Linux environment. It does not download or analyse the third-party GWAS files.

## Prepared-input analysis coverage

The complete accession-to-trait mapping is released in `analysis/manifests/metabolic_traits_249.tsv`. Local source and exposure manifests can be built with `analysis/prepare_inputs.py` after the input checks pass.

The released code covers:

- MAF-aware metabolic-trait instrument selection and local PLINK clumping;
- deterministic GRCh38 exposure `variant_id`-to-rsID mapping against dbSNP build 157 using exact unordered allele-pair matching;
- MDD instrument selection for reverse MR;
- main and UK Biobank-excluded MDD harmonization with action 2, plus action 3 palindrome removal;
- random-effects IVW, MR-Egger, weighted-median and mode estimators;
- Cochran Q, MR-Egger intercept, leave-one-out analysis and MR-PRESSO;
- broad-pleiotropy removal and Steiger prevalence sensitivity analyses;
- reverse MR using main and UK Biobank-excluded MDD instrument sets;
- strict 249-row Benjamini–Hochberg/Bonferroni calculation and four-level screen assignment;
- FinnGen R13 `F5_DEPRESSIO` admission checks, harmonization, tool-loss accounting, MR and descriptive cross-study heterogeneity;
- neutral FinnGen support labels (`FDR-supported`, `nominally-supported`, `no-nominal-support`) that do not use the overlap-uncorrected Q statistic;
- anchor-based locus construction and one-read-per-source regional extraction;
- single-causal-variant `coloc.abf` under `p12={1e-6,5e-6,1e-5,1e-4}`; and
- fail-closed locus status reporting, lead/proxy QC, MHC window-overlap exclusion and the locked five-label colocalization classification.

`data/derived/instrument_strength_qc.tsv` reports all 249 traits with machine and display names, instrument counts, minimum/median F, cumulative approximate standardized-trait R², MAF quantiles, I²GX and selection-stage LD-panel coverage. I²GX is calculated as `max(0, (QX-(K-1))/QX)`, where `QX` is the inverse-variance-weighted heterogeneity statistic for SNP-exposure estimates. The `nome_i2gx_ge_0_90` field is a diagnostic flag, not proof that MR-Egger is unbiased. The companion `instrument_selection_ld_panel_qc_249.tsv` records the autosomal biallelic candidate count, count represented in the 1KG EUR panel and count absent from that panel for every trait. Across 15,105,492 eligible candidate rows, 961,597 (6.37%) were absent; trait-specific fractions ranged from 2.79% to 16.72% (median 6.36%). These counts were recomputed from the pre-clumping mapped candidate files and the registered panel, not inferred from final instruments.

Frozen parameters are listed in `analysis/config/parameters.tsv`. The mapping between the portable modules and the analysis-time source scripts, including source SHA-256 hashes, is in `analysis/provenance.tsv`.

The Steiger workflow now recovers the `n` value for each candidate instrument from its matching `meta_EUR` source row before harmonization. Recovery matched all 5,665 candidate IVs; the observed source values were 413,897 and 599,249. All 15 candidate estimates retained the exposure-to-outcome direction at outcome prevalences 0.08, 0.15 and 0.20. This remains a directionality sensitivity analysis rather than proof against reverse causation.

## Reproducibility boundary

The public repository does not include original GWAS rows, individual-level records, LD panels, dbSNP VCF/BCF files, regional SNP extracts, credentials or unpublished review material. A clean clone therefore reproduces the figures, screen/count checks, the frozen ABF classification layer and synthetic method checks. With the listed third-party data and all documented prepared inputs, the public scripts cover the main MR, FinnGen, locus and ABF-classification stages. Integrated-evidence and SuSiE computations remain audit-only because their regional data, LD matrices and credible-set inputs are not redistributed.

Released audit tables include:

- three 249-trait screens;
- the 15 candidate-level forward-MR estimates for the PGC outcome excluding UK Biobank used in Figure 2;
- the 15-trait leave-one-out summary, which reports aggregate stability without releasing SNP-level rows;
- the complete 15-trait MR-PRESSO summary and candidate-level Steiger/sample-size QC, without variant identifiers or local paths;
- the 6,434-row locus manifest and ABF classification inputs/outputs;
- the 101-row integrated evidence table and FinnGen tool-loss/cross-outcome summaries; and
- exploratory SuSiE status, binding and stopping summaries.

These tables permit auditing of the reported results. Recomputing the integrated-evidence and SuSiE stages still requires unreleased regional data, LD matrices and author credible-set inputs.

In `coloc_status.tsv`, `main` and `noUKBB` are the two outcome analyses; rows labelled `classification` preserve the frozen classification-admission record and are not a third GWAS outcome. The locus manifest releases lead rsID/position identifiers and window boundaries, but not regional SNP association statistics. Provider terms and the aggregate-only schema boundary are recorded in `data/ACCESS.md`. Version 0.2.19 adds only trait-level LD-panel selection counts and the corresponding completed instrument-QC fields; it does not add variant-level rows or change MR estimates.

The FinnGen analysis is a cross-outcome comparison, not an independent replication of the PGC major-depression phenotype. Directional concordance does not establish a shared causal mechanism. The single-variant ABF model is not a substitute for multi-signal fine-mapping.

## Frozen reporting and classification rules

For all 15 candidate reruns, the MR-PRESSO global test had zero simulated exceedances among 10,000 draws. The display field therefore reports `P<1e-04`, while the numeric field records the simulation resolution as `1e-4`; trait-level outlier counts range from 7 to 19. Thirteen distortion tests are nonsignificant and two are significant. The evidence matrix encodes distortion only. `data/derived/presso_sensitivity_15.tsv` contains the complete candidate-level global, outlier, distortion and regression fields, including standard errors retained directly from `Main MR results$Sd`, residual degrees of freedom and two-sided t-based 95% confidence intervals. No uncertainty value was reconstructed from a P value.

The candidate-level leave-one-out summary contains 4,993 single-variant deletions across the 15 traits. Every deletion retained the sign of the corresponding primary IVW estimate and remained nominally significant. The released table is trait-level and contains no SNP identifiers.

A candidate belongs to the priority reporting group only when its high-precision MR-PRESSO distortion P value is at least 0.05 and its IVW P value remains below 0.05 after removal of variants used as instruments for more than five metabolic traits. All other candidates are secondary. These are reporting categories, not new significance tiers or mechanism labels.

For ABF classification, a ratio pass requires `PP.H4>0.90` and either `PP.H3=0` with positive `PP.H4`, or `PP.H4/PP.H3>3`. Robust colocalization requires a ratio pass in both PGC datasets at `p12=5e-6` and `PP.H4>0.80` in both at `p12=1e-6`. If robust fails, prior-sensitive colocalization requires a ratio pass in the primary PGC dataset at `p12=5e-6` or `1e-5`. The exact distinct-signal, trait-specific/low-power and inconclusive branches are implemented in `analysis/scripts/08_classify_coloc.R`. Only robust, non-complex, non-hotspot records with the documented concordant forward and reverse regional exclusion results can meet the composite mechanism-support criterion; all 101 records fail at least one required condition.

## Upstream data-use boundary

This repository does not redistribute source GWAS summary statistics. It contains original analysis code and author-generated aggregate results only.

- The official Catalog records returned `agreedToCc0=true` for all 249 accessions, GCST90451106–GCST90451354; each per-accession v2 API record exposes [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/) as `terms_of_license`.
- [PGC MDD2025 Figshare v5](https://figshare.com/articles/dataset/GWAS_summary_statistics_for_major_depression_PGC_MDD2025_/27061255/5) is labelled CC BY 4.0, but its packaged [README](https://ndownloader.figshare.com/files/52419692) prohibits cross-posting the source data and limits use to scientific research unless the PGC Data Access Committee approves commercial use. No PGC source rows are included here.
- [FinnGen Data Freeze 13](https://www.finngen.fi/en/access_results) is a public release obtained through its access procedure. FinnGen requires acknowledgment and citation and asks anyone wishing to host its summary statistics to contact the project. No FinnGen source rows are included here.

The locus manifest contains lead variant identifiers, positions and author-defined windows. The ABF table contains locus-level prior-sensitivity summaries and posterior probabilities. Neither contains source alleles, effect estimates, standard errors, association P values, allele frequencies or regional SNP rows. See `data/ACCESS.md` for the recorded checks and release rules.

## Licenses and contact

Code is released under the MIT License. Author-created aggregate data are released under CC BY 4.0, subject to upstream rights and provider terms. Repository documentation is dual-licensed under the MIT License and CC BY 4.0; users may follow either license for those files. These licences grant no rights in third-party GWAS data and do not override provider-specific terms.

Zhouyi Wang: wangzhouyi1991@163.com
