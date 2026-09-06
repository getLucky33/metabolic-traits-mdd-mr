# Circulating metabolic traits and major depressive disorder

[![reproduce](https://github.com/getLucky33/metabolic-traits-mdd-mr/actions/workflows/reproduce.yml/badge.svg)](https://github.com/getLucky33/metabolic-traits-mdd-mr/actions/workflows/reproduce.yml)

Reproducibility code for the Scientific Reports manuscript by Zhouyi Wang and Qingmei Liu.

The repository has two reproducibility layers:

1. `scripts/make_figures.R` regenerates Figures 1–4 from the released aggregate tables.
2. `analysis/` contains portable code for instrument selection, forward and reverse Mendelian randomization (MR), sensitivity analyses, FinnGen cross-outcome comparison, and ABF colocalization. These stages require the provider-controlled GWAS and reference files described in `data/ACCESS.md`; those files are not redistributed.

## What can be run from a clean clone

Install R 4.5.1 and restore the versions recorded in `renv.lock` or `analysis/environment/package-versions.tsv`.

```bash
Rscript scripts/make_figures.R
Rscript analysis/tests/smoke_test.R
```

The first command verifies the released aggregate-input hashes and regenerates PDF, SVG and 600-dpi PNG figures. The second creates synthetic data in a temporary directory and exercises the same TwoSampleMR harmonization, random-effects IVW and four-prior `coloc.abf` paths used by the analysis scripts. Synthetic results are tests only and are never mixed with manuscript results.

GitHub Actions runs figure reproduction, code parsing, a secret/path scan and the synthetic analysis smoke test in a clean Linux environment.

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
- FinnGen R13 `F5_DEPRESSIO` admission checks, harmonization, MR and cross-study heterogeneity;
- single-causal-variant `coloc.abf` under `p12={1e-6,5e-6,1e-5,1e-4}`; and
- fail-closed locus status reporting, lead/proxy QC and the locked five-label colocalization classification.

Frozen parameters are listed in `analysis/config/parameters.tsv`. The mapping between the portable modules and the analysis-time source scripts, including source SHA-256 hashes, is in `analysis/provenance.tsv`.

## Reproducibility boundary

The public repository does not include original GWAS rows, individual-level records, LD panels, dbSNP VCF/BCF files, licensed regional extracts, credentials or unpublished review material. Consequently, a clone alone reproduces the figures and synthetic checks; a full numerical rerun additionally requires lawful acquisition of the same upstream resources.

The FinnGen analysis is a cross-outcome comparison, not an independent replication of clinical MDD. Directional concordance does not establish a shared causal mechanism. The single-variant ABF model is not a substitute for multi-signal fine-mapping.

## Licenses and contact

Code is released under the MIT License. Author-created aggregate data and documentation are released under CC BY 4.0, subject to upstream rights and provider terms.

Zhouyi Wang: wangzhouyi1991@163.com
