# Data access and redistribution boundary

The files under `data/derived/` are compact, author-created aggregate reporting inputs. They do not contain participant-level data or original GWAS rows.

The original upstream files are not redistributed here, but the summary-statistics datasets used in the manuscript can be located as follows:

- Circulating metabolic traits: Tambets et al., *Nature* (2026), DOI `10.1038/s41586-026-10532-5`. The 249 `meta_EUR` datasets are the contiguous NHGRI–EBI GWAS Catalog accession range `GCST90451106`–`GCST90451354`: https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90451001-GCST90452000/.
- PGC major depressive disorder: PGC MDD2025 Figshare dataset version 5, DOI https://doi.org/10.6084/m9.figshare.27061255.v5. The primary file is `daner/daner_pgc_mdd_no23andMe_eur_hg19_v3.49.24.11.neff.gz` (https://ndownloader.figshare.com/files/52371878), and the UK Biobank-excluded sensitivity file is `daner/daner_pgc_mdd_no23andMe-noUKBB_eur_hg19_v3.49.24.11.neff.gz` (https://ndownloader.figshare.com/files/52371881). Both exclude 23andMe data.
- FinnGen: Data Freeze 13 endpoint `F5_DEPRESSIO`. Follow the official access procedure at https://www.finngen.fi/en/access_results. The endpoint-specific summary-statistics file is https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_F5_DEPRESSIO.gz, and the phenotype definition is at https://r13.risteys.finngen.fi/endpoints/F5_DEPRESSIO.
- LD and identifier references: the 1000 Genomes phase 3 European reference resources are available from https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/, and dbSNP build 157 resources are available from https://ftp.ncbi.nih.gov/snp/archive/b157/VCF/.

Not copying public GWAS files into this repository is a redistribution and repository-size choice; it is not a claim that those summary statistics are inaccessible. Users should retrieve them from the identifiers and links above and comply with the applicable provider terms. Individual-level data, restricted extracts, LD/reference genotypes and dbSNP files are not redistributed.

The portable analysis code expects local paths supplied through manifests. Recommended local layout is `data-local/`, which is ignored by Git. The required schemas are documented by the examples in `analysis/manifests/` and by `analysis/WORKFLOW.md`.

Do not commit or redistribute original or substantially reconstructed GWAS rows, harmonized SNP-level tables unless provider terms permit it, LD/reference genotypes, dbSNP VCF/BCF files, restricted regional extracts, fine-mapping objects, credible sets, credentials, signed URLs, cookies or access tokens.

The analysis scripts record basenames and SHA-256 values for supplied inputs. Users should retain provider download receipts, accession/version identifiers and checksums in a private run record.
