# Circulating metabolic traits and major depressive disorder

[![reproduce](https://github.com/getLucky33/metabolic-traits-mdd-mr/actions/workflows/reproduce.yml/badge.svg)](https://github.com/getLucky33/metabolic-traits-mdd-mr/actions/workflows/reproduce.yml)

Reproducibility code for the Scientific Reports manuscript by Zhouyi Wang and Qingmei Liu. Version 0.2.13 is a local release candidate that reorganizes the reproduction entry points without changing the frozen analyses or results.

## Choose a reproduction target

| Target | What can be reproduced | Additional inputs |
| --- | --- | --- |
| Clean-clone verification | Figures 1–4, the complete 249-trait screen counts, the frozen 6,434-record colocalization classification and synthetic method checks | None beyond the recorded software environment |
| Prepared-input analysis | Instrument selection, forward and reverse MR, sensitivity analyses, the FinnGen alternative-outcome comparison, locus construction, ABF colocalization and classification | Third-party GWAS, LD and dbSNP resources plus the prepared inputs listed in [`DATA_SOURCES.md`](DATA_SOURCES.md) |
| Aggregate-result audit only | The 101-row integrated-evidence table and exploratory SuSiE status, binding and stopping summaries | Independent recomputation is not possible from the released files alone |

This is not a one-command raw-data-to-paper repository. Downloading the public GWAS files alone is insufficient: the prepared-input stages also require dbSNP-derived rsID maps, selected MDD instruments, coordinate tables and regional inputs. The exact boundary is stated before the commands so that a successful figure build is not mistaken for a complete numerical rerun.

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

After preparing all 249 dbSNP-derived rsID maps, validate the analysis-ready inputs and build the complete local manifests:

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

Use `Rscript scripts/make_figures.R --out <directory>` to choose another output directory. The run also writes the exact Figure 2 plotting values, Figure 3 matrix values, Figure 4 classification-to-disposition counts and a checksum file beside the figures. Figure 1 uses neutral GWAS-source headings and presents forward and reverse MR as parallel direction-specific analyses with separately selected instruments. Only the primary forward screen selects the 15-trait candidate set; those candidates enter the robustness, locus and directional assessments, whereas matched reverse-MR results enter only the directional assessment. The layout does not treat reverse MR as validation or proof of reciprocal causality. Figure 2 includes the primary PGC result, the PGC sensitivity analysis excluding UK Biobank and the FinnGen R13 alternative outcome for each candidate. Figure 4 recomputes the five ABF classifications across 6,434 analysis-window records and the disposition of the 101-record integrated-review subset from the released aggregate tables.

On each push or pull request, the configured GitHub Actions workflow runs the repository audit, the synthetic input-preflight test, figure reproduction, code parsing, a secret/path scan and the synthetic MR/colocalization smoke test in a clean Linux environment. It does not download or analyse the third-party GWAS files.

## Prepared-input analysis coverage

The complete accession-to-trait mapping is released in `analysis/manifests/metabolic_traits_249.tsv`. Local source and exposure manifests can be built with `analysis/prepare_inputs.py` after the input checks pass.

The released code covers:

- MAF-aware metabolic-trait instrument selection and local PLINK clumping;
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

`data/derived/instrument_strength_qc.tsv` reports all 249 traits with machine and display names, instrument counts, minimum/median F, cumulative approximate standardized-trait R², MAF quantiles and I²GX. I²GX is calculated as `max(0, (QX-(K-1))/QX)`, where `QX` is the inverse-variance-weighted heterogeneity statistic for SNP-exposure estimates. The `nome_i2gx_ge_0_90` field is a diagnostic flag, not proof that MR-Egger is unbiased. Selection-stage LD-panel denominators were not retained for all 249 frozen runs, so the corresponding missingness fields are `NA` with an explicit reason rather than reconstructed from final IVs.

Frozen parameters are listed in `analysis/config/parameters.tsv`. The mapping between the portable modules and the analysis-time source scripts, including source SHA-256 hashes, is in `analysis/provenance.tsv`.

The released Steiger code uses `599,249`, the maximum exposure meta-analysis sample size, for every exposure SNP because per-SNP effective sample sizes were not retained in the frozen harmonized inputs. This is an approximation: with beta and SE fixed, a larger substituted N tends to reduce the estimated exposure-side R², and it can change Steiger P values or near-boundary direction calls. The Steiger layer should therefore be read as a sensitivity analysis, not as proof against reverse causation.

## Reproducibility boundary

The public repository does not include original GWAS rows, individual-level records, LD panels, dbSNP VCF/BCF files, regional SNP extracts, credentials or unpublished review material. A clean clone therefore reproduces the figures, screen/count checks, the frozen ABF classification layer and synthetic method checks. With the listed third-party data and all documented prepared inputs, the public scripts cover the main MR, FinnGen, locus and ABF-classification stages. Integrated-evidence and SuSiE computations remain audit-only because their regional data, LD matrices and credible-set inputs are not redistributed.

Released audit tables include:

- three 249-trait screens;
- the 15 candidate-level forward-MR estimates for the PGC outcome excluding UK Biobank used in Figure 2;
- the 6,434-row locus manifest and ABF classification inputs/outputs;
- the 101-row integrated evidence table and FinnGen tool-loss/cross-outcome summaries; and
- exploratory SuSiE status, binding and stopping summaries.

These tables permit auditing of the reported results. Recomputing the integrated-evidence and SuSiE stages still requires unreleased regional data, LD matrices and author credible-set inputs.

In `coloc_status.tsv`, `main` and `noUKBB` are the two outcome analyses; rows labelled `classification` preserve the frozen classification-admission record and are not a third GWAS outcome. The locus manifest releases lead rsID/position identifiers and window boundaries, but not regional SNP association statistics. Provider terms and the current schema boundary were checked on 2026-09-07 and are recorded in `data/ACCESS.md`; responsible-author attestation for the exact current release schemas and publication boundary was confirmed by Zhouyi Wang on 2026-09-07.

The FinnGen analysis is a cross-outcome comparison, not an independent replication of the PGC major-depression phenotype. Directional concordance does not establish a shared causal mechanism. The single-variant ABF model is not a substitute for multi-signal fine-mapping.

## Upstream data-use boundary

This repository does not redistribute source GWAS summary statistics. It contains original analysis code and author-generated aggregate results only.

- The official Catalog records returned `agreedToCc0=true` for all 249 accessions, GCST90451106–GCST90451354; each per-accession v2 API record exposes [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/) as `terms_of_license`.
- [PGC MDD2025 Figshare v5](https://figshare.com/articles/dataset/GWAS_summary_statistics_for_major_depression_PGC_MDD2025_/27061255/5) is labelled CC BY 4.0, but its packaged [README](https://ndownloader.figshare.com/files/52419692) prohibits cross-posting the source data and limits use to scientific research unless the PGC Data Access Committee approves commercial use. No PGC source rows are included here.
- [FinnGen Data Freeze 13](https://www.finngen.fi/en/access_results) is a public release obtained through its access procedure. FinnGen requires acknowledgment and citation and asks anyone wishing to host its summary statistics to contact the project. No FinnGen source rows are included here.

The locus manifest contains lead variant identifiers, positions and author-defined windows. The ABF table contains locus-level prior-sensitivity summaries and posterior probabilities. Neither contains source alleles, effect estimates, standard errors, association P values, allele frequencies or regional SNP rows. See `data/ACCESS.md` for the recorded checks and release rules.

## Licenses and contact

Code is released under the MIT License. Author-created aggregate data are released under CC BY 4.0, subject to upstream rights and provider terms. Repository documentation is dual-licensed under the MIT License and CC BY 4.0; users may follow either license for those files. These licences grant no rights in third-party GWAS data and do not override provider-specific terms.

Zhouyi Wang: wangzhouyi1991@163.com
