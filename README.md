# Circulating metabolic traits and major depressive disorder

[![reproduce](https://github.com/getLucky33/metabolic-traits-mdd-mr/actions/workflows/reproduce.yml/badge.svg)](https://github.com/getLucky33/metabolic-traits-mdd-mr/actions/workflows/reproduce.yml)

Reproducibility code for the Scientific Reports manuscript by Zhouyi Wang and Qingmei Liu.

Version 0.2.5 has three reproducibility layers:

1. `scripts/make_figures.R` regenerates Figures 1–4 from the released aggregate tables.
2. Released aggregate inputs replay the complete 249-trait four-level screens and the 6,434-row colocalization classification, with fixed count assertions.
3. `analysis/` contains prepared-input code for instrument selection, forward and reverse Mendelian randomization (MR), sensitivity analyses, FinnGen cross-outcome comparison, locus-window construction, and ABF colocalization. These stages require the original third-party GWAS and reference files described in `data/ACCESS.md`; those files are not duplicated here.

## What can be run from a clean clone

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

Use `Rscript scripts/make_figures.R --out <directory>` to choose another output directory. The run also writes the exact Figure 2 plotting values, Figure 3 matrix values and a checksum file beside the figures. Figure 2 includes the primary PGC result, the PGC sensitivity analysis excluding UK Biobank and the FinnGen R13 alternative outcome for each of the 15 forward candidates.

On each push or pull request, the configured GitHub Actions workflow runs figure reproduction, code parsing, a secret/path scan and the synthetic analysis smoke test in a clean Linux environment.

## Upstream analysis pipeline

After obtaining the third-party data and a 1000 Genomes phase 3 European LD reference, copy the example manifests in `analysis/manifests/`, replace only their paths and follow `analysis/WORKFLOW.md`.

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

The public repository does not include original GWAS rows, individual-level records, LD panels, dbSNP VCF/BCF files, regional SNP extracts, credentials or unpublished review material. A clean clone therefore reproduces the figures, screen/count checks, the frozen ABF classification layer and synthetic method checks. It does not reproduce the complete MR, ABF, integrated-evidence or SuSiE computations without prepared upstream data and LD inputs.

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
