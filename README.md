# Circulating metabolic traits and major depressive disorder: reproducibility code

This repository accompanies the Scientific Reports manuscript by Zhouyi Wang and Qingmei Liu. It is a minimal, self-contained reporting bundle: one deterministic R script regenerates Figures 1-4 from the frozen aggregate tables in `data/derived/` and checks the reported evidence counts.

## Reproduce

Requirements: R 4.3 or later and the packages `data.table`, `ggplot2`, `patchwork`, and `scales`.

```r
install.packages(c("data.table", "ggplot2", "patchwork", "scales"))
```

From the repository root:

```bash
Rscript scripts/make_figures.R
```

The script verifies the released-input hashes, writes PDF, SVG, and 600-dpi PNG files to `results/figures/`, and creates `results/figure_checksums.tsv`. It fails if an input hash changes, the frozen inputs are incomplete, trait labels are not unique/full-length, evidence-stage counts disagree, or an expected output is absent.

GitHub Actions runs the same command in a clean Linux environment on every push.

## Reproducibility scope

Included:

- aggregate estimates for the 15 Bonferroni-significant traits;
- primary and UK Biobank-excluded reverse-MR screen labels used in Figure 3;
- full journal-facing trait names;
- frozen study and colocalization counts;
- the manuscript GWAS-design table; and
- deterministic code for Figures 1-4.

Not included:

- original third-party GWAS summary statistics;
- participant-level data;
- linkage-disequilibrium or reference-genotype panels;
- dbSNP VCF files, provider-controlled regional extracts, or licensed credible sets; and
- access credentials, local paths, unpublished review files, or submission forms.

This repository reproduces the article's reporting layer from aggregate derived inputs. It does not claim to reproduce the upstream GWAS, MR, or colocalization analyses without the provider-controlled source data. See `data/ACCESS.md` for access boundaries.

## Integrity and interpretation

The evidence groups were locked before the FinnGen analysis. FinnGen R13 is a cross-outcome comparison, not an independent replication of clinical MDD. Directional concordance does not establish a shared causal mechanism. The included aggregate data contain no individual-level records.

## Licenses

Code is released under the MIT License. Author-created aggregate data and documentation are released under CC BY 4.0, subject to the rights and terms of the original data providers.

## Contact

Zhouyi Wang: wangzhouyi1991@163.com
