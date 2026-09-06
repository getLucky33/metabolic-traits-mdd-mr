# Data access and redistribution boundary

The files under `data/derived/` are compact, author-created aggregate reporting inputs. They do not contain participant-level data or original GWAS rows.

The upstream datasets remain under their providers' access and use conditions and are not redistributed here:

- Circulating metabolic traits: Tambets et al., *Nature* (2026), DOI `10.1038/s41586-026-10532-5`.
- PGC major depressive disorder: Major Depressive Disorder Working Group of the Psychiatric Genomics Consortium, *Cell* (2025), DOI `10.1016/j.cell.2024.12.002`.
- FinnGen: Kurki et al., *Nature* (2023), DOI `10.1038/s41586-022-05473-8`; the analysis used Data Freeze 13 endpoint `F5_DEPRESSIO`.
- LD reference data: 1000 Genomes phase 3 European reference panel; obtain from the original provider under its terms.

Access to upstream data must be arranged directly with the relevant provider. Users are responsible for complying with the applicable consent, institutional, contractual, and data-use requirements.

The portable analysis code expects local paths supplied through manifests. Recommended local layout is `data-local/`, which is ignored by Git. The required schemas are documented by the examples in `analysis/manifests/` and by `analysis/WORKFLOW.md`.

Do not commit or redistribute original or substantially reconstructed GWAS rows, harmonized SNP-level tables unless provider terms permit it, LD/reference genotypes, dbSNP VCF/BCF files, restricted regional extracts, fine-mapping objects, credible sets, credentials, signed URLs, cookies or access tokens.

The analysis scripts record basenames and SHA-256 values for supplied inputs. Users should retain provider download receipts, accession/version identifiers and checksums in a private run record.
