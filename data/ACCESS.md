# Data access and redistribution boundary

The files under `data/derived/` are compact, author-created aggregate reporting inputs. They do not contain participant-level data or original/regional GWAS association rows. `coloc_locus_manifest.tsv` contains only locus-level lead rsID/position identifiers and author-defined window boundaries. The provider terms and current release schemas were checked on 2026-09-07. Zhouyi Wang confirmed the responsible-author attestation below for the exact current release schemas and publication boundary.

## Verification status (checked 2026-09-07)

| Provider/source | Evidence status | Repository release rule |
|---|---|---|
| NHGRI–EBI GWAS Catalog, GCST90451106–GCST90451354 | `CC0_VERIFIED` (249/249) | Source data are not bundled; attribution and accessions are retained for provenance. |
| PGC MDD2025 Figshare v5 | `TERMS_VERIFIED__NO_SOURCE_REHOSTING` | Do not rehost source or raw-equivalent rows. Users who obtain or use the PGC source data must comply with its scientific-use, non-commercial-use and citation conditions. |
| FinnGen Data Freeze 13 | `PUBLIC_ACCESS_VERIFIED__NO_GENERAL_OPEN_LICENSE_FOUND` | Author-generated aggregate outputs only; no source or raw-equivalent rows. Acknowledge FinnGen and cite Kurki et al.; contact FinnGen before hosting source summary statistics. |

**Provider-terms review:** VERIFIED.

**Technical release boundary:** VERIFIED for the current schemas.

**Responsible-author attestation:** CONFIRMED by Zhouyi Wang on 2026-09-07 for the exact current release schemas and publication boundary stated in this file.

> 我确认本仓库仅公开作者生成的位点级与先验敏感性聚合结果，不含 PGC 或 FinnGen 原始、区域或可还原的 SNP 级汇总统计；本项目通过官方渠道获取并接受适用条款，用途为非商业科学研究，已履行不识别、引用与 FinnGen 致谢要求，并同意按 ACCESS.md 所列边界公开。

**Attested aggregate-manifest snapshot:** `data/derived/data_manifest.tsv`, SHA-256 `7aa9102539c535cee2c693842a204e0f1042095577a97511a01887f5e9b555c7`.

**Blanket upstream redistribution clearance:** NOT CLAIMED.

The 6,434-row locus manifest has no source allele, beta, standard-error, association-P-value, allele-frequency or regional-SNP fields. The 51,472-row ABF table contains eight author-generated prior-sensitivity summaries per locus (`p12`, `nsnp` and `PP.H0`–`PP.H4`) and no per-variant source statistics. These files cannot substitute for the official source datasets.

The original upstream files are not redistributed here, but the summary-statistics datasets used in the manuscript can be located as follows:

- Circulating metabolic traits: Tambets et al., *Nature* (2026), DOI `10.1038/s41586-026-10532-5`. The 249 `meta_EUR` datasets are the contiguous NHGRI–EBI GWAS Catalog accession range `GCST90451106`–`GCST90451354`: https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90451001-GCST90452000/.
- PGC major depressive disorder: PGC MDD2025 Figshare dataset version 5, DOI https://doi.org/10.6084/m9.figshare.27061255.v5. The primary file is `daner/daner_pgc_mdd_no23andMe_eur_hg19_v3.49.24.11.neff.gz` (https://ndownloader.figshare.com/files/52371878), and the UK Biobank-excluded sensitivity file is `daner/daner_pgc_mdd_no23andMe-noUKBB_eur_hg19_v3.49.24.11.neff.gz` (https://ndownloader.figshare.com/files/52371881). Both exclude 23andMe data.
- FinnGen: Data Freeze 13 endpoint `F5_DEPRESSIO`. Follow the official access procedure at https://www.finngen.fi/en/access_results. The endpoint-specific summary-statistics file is https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_F5_DEPRESSIO.gz, and the phenotype definition is at https://r13.risteys.finngen.fi/endpoints/F5_DEPRESSIO.
- LD and identifier references: the 1000 Genomes phase 3 European reference resources are available from https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/, and dbSNP build 157 resources are available from https://ftp.ncbi.nih.gov/snp/archive/b157/VCF/.

Public GWAS files are omitted to avoid redistribution and repository-size burdens. They remain available from the identifiers and links above, subject to provider terms. Individual-level data, restricted extracts, LD/reference genotypes and dbSNP files are not redistributed.

The metabolic-trait accessions record CC0 1.0 terms in the GWAS Catalog. PGC MDD2025 Figshare v5 is labelled CC BY 4.0, while its packaged README imposes the stricter no-cross-posting and scientific-use conditions followed here. FinnGen Data Freeze 13 is publicly released but no general Creative Commons/open-data licence was identified on its official pages; its access, acknowledgment, citation and no-rehosting boundary therefore controls.

The portable analysis code expects local paths supplied through manifests. Recommended local layout is `data-local/`, which is ignored by Git. The required schemas are documented by the examples in `analysis/manifests/` and by `analysis/WORKFLOW.md`.

Do not commit or redistribute original or substantially reconstructed GWAS rows, harmonized SNP-level tables unless provider terms permit it, LD/reference genotypes, dbSNP VCF/BCF files, restricted regional extracts, fine-mapping objects, credible sets, credentials, signed URLs, cookies or access tokens.

The analysis scripts record basenames and SHA-256 values for supplied inputs. Users should retain provider download receipts, accession/version identifiers and checksums in a private run record.
