# Data required for real-data reruns

The repository does not download or redistribute third-party GWAS or reference files. This page is the operational checklist for a prepared-input rerun. Provider terms and the public-release boundary are recorded separately in [`data/ACCESS.md`](data/ACCESS.md).

## Downloaded sources

| Source | Exact version or identifier | Build | Suggested local location | Used for |
| --- | --- | --- | --- | --- |
| Circulating metabolic traits | Tambets `meta_EUR`, GCST90451106–GCST90451354; the complete accession-to-trait table is [`analysis/manifests/metabolic_traits_249.tsv`](analysis/manifests/metabolic_traits_249.tsv) | GRCh38 | `data-local/metabolite/`; retain each filename as `<accession>.tsv.gz` | Instrument selection, reverse MR outcomes and exposure-side colocalization |
| PGC major depression | MDD2025 Figshare v5, main file [`daner_pgc_mdd_no23andMe_eur_hg19_v3.49.24.11.neff.gz`](https://ndownloader.figshare.com/files/52371878) | GRCh37 | `data-local/mdd/` | Primary forward MR and primary reverse-MR instruments |
| PGC major depression excluding UK Biobank | MDD2025 Figshare v5, sensitivity file [`daner_pgc_mdd_no23andMe-noUKBB_eur_hg19_v3.49.24.11.neff.gz`](https://ndownloader.figshare.com/files/52371881) | GRCh37 | `data-local/mdd/` | Sample-overlap sensitivity analysis and second reverse-MR instrument set |
| FinnGen depression | Data Freeze 13, endpoint `F5_DEPRESSIO`, file [`finngen_R13_F5_DEPRESSIO.gz`](https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_F5_DEPRESSIO.gz) | GRCh38 | `data-local/finngen/` | Alternative-outcome comparison for the 15 forward candidates |
| Linkage-disequilibrium reference | [1000 Genomes phase 3](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) European samples | GRCh37 | PLINK prefix `data-local/ld/1kg.v3/EUR` with `.bed`, `.bim` and `.fam` files | Instrument clumping and lead/proxy checks |
| Variant identifier reference | dbSNP build 157 GRCh38 VCF/BCF plus index, from the [NCBI archive](https://ftp.ncbi.nih.gov/snp/archive/b157/VCF/) | GRCh38 | `data-local/dbsnp/` | Mapping metabolic `variant_id` values to rsIDs and preparing MDD coordinates |

The metabolic-trait files can be downloaded from the URLs in the 249-row accession table. PGC and FinnGen files remain subject to their provider-specific terms. The 1000 Genomes download is not itself the PLINK prefix expected by the scripts; prepare or obtain a verified European `.bed/.bim/.fam` set and record its checksum.

## Required source schemas

- Metabolic GWAS: `variant_id`, `effect_allele`, `other_allele`, `beta`, `standard_error`, `effect_allele_frequency` and `neg_log_10_p_value`.
- PGC files: `CHR`, `BP`, `SNP`, `A1`, `A2`, `OR`, `SE`, `P`, `Nca`, `Nco`, and dataset-specific `FRQ_A_*` and `FRQ_U_*` columns.
- FinnGen R13: `#chrom`, `pos`, `ref`, `alt`, `rsids`, `nearest_genes`, `pval`, `mlogp`, `beta`, `sebeta`, `af_alt`, `af_alt_cases` and `af_alt_controls`, in that order.
- Metabolic rsID maps: `variant_id`, `rsid` and `rsid_status`; only rows labelled `matched` are admitted.

## Prepared inputs that are not simple downloads

Downloading the six source groups above is necessary but not sufficient for every analysis stage.

1. Instrument selection requires one dbSNP-derived rsID map per metabolic trait in the configured `rsid_map_dir`. The preflight accepts `<trait>_rsid_map.tsv` or the frozen-project convention `<trait>_forward_ivs_rsid_b157.tsv`. The current public repository validates and consumes these maps but does not regenerate the complete frozen mapping set.
2. Reverse-MR extraction requires the selected MDD instruments and either a prepared GRCh38 coordinate map or `bcftools` plus the indexed dbSNP build 157 file.
3. Locus construction requires the complete forward screen, the 249-trait IV manifest, two MDD coordinate tables and the recorded reverse-tier table.
4. ABF colocalization requires the locus manifest, regional metabolic files, both regional PGC outcome files and their sample metadata. The scripts record failed or incomplete loci rather than silently dropping them.
5. The released integrated-evidence and exploratory SuSiE tables are audit snapshots. Independent recomputation requires regional association data, LD matrices and credible-set inputs that are not redistributed.

## Before running an analysis

Copy the local-path template and replace its paths:

```bash
cp analysis/config/local_paths.example.tsv analysis/config/local_paths.local.tsv
python analysis/prepare_inputs.py check --config analysis/config/local_paths.local.tsv --level sources
```

The `sources` check verifies the 249 accession set, file uniqueness, headers, the two PGC files, FinnGen, the PLINK LD prefix, dbSNP and its index. Once the 249 rsID maps have been prepared, run:

```bash
python analysis/prepare_inputs.py check --config analysis/config/local_paths.local.tsv --level analysis
python analysis/prepare_inputs.py build-manifests \
  --config analysis/config/local_paths.local.tsv \
  --out-dir analysis/manifests
```

This writes ignored local manifests for the 249 metabolic traits and the main PGC dataset. It does not download data, infer missing trait identities or bypass the one-trait smoke tests required by [`analysis/WORKFLOW.md`](analysis/WORKFLOW.md).
