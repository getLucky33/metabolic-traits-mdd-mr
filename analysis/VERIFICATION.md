# Verification record

Verification history: real-data one-trait checks were completed on 2026-09-06; v0.2.3, v0.2.4 and v0.2.5 local release checks were completed on 2026-09-07; v0.2.6, v0.2.7 and v0.2.8 local release checks were completed on 2026-09-08; v0.2.9, v0.2.10, v0.2.11 and v0.2.12 local release checks were completed on 2026-09-09. Versions 0.2.13 and 0.2.14 were checked locally and in clean GitHub Actions environments on 2026-09-11. Version 0.2.15 retains the same frozen analyses and clarifies the documented MR-PRESSO, reporting-group, ABF-classification and integrated-assessment rules. Runtime: R 4.5.1 on Windows 10 x64 using the recorded environment.

## v0.2.15 release verification (2026-09-12)

- The released MR-PRESSO table is checked for 15 global P values of `1e-4`, trait-level outlier counts of 7–19 and the frozen 13 nonsignificant/2 significant distortion-test split.
- Figure 1 and Figure 4 explicitly separate the 101-record H4-oriented subset into 14 robust and 87 prior-sensitive classifications. Figure 3 labels its MR-PRESSO column as distortion-only.
- The README and workflow state the Boolean rule for the nine priority and six secondary reporting groups, the exact robust/prior-sensitive ABF branches and the complete composite mechanism-support branch.
- The statistical results and all released aggregate inputs are unchanged. Local and clean-checkout verification are required for the release commit.

## v0.2.14 release verification (2026-09-11)

- Figure 2 adds a 249-row primary forward-screen panel and retains the 45 exact candidate estimates from the primary PGC, UK Biobank-excluded PGC and FinnGen sources.
- Figure 3 retains all 180 evidence cells and groups them into six analytical domains without changing their values.
- Figure 4 retains the five-class 6,434-record totals and 101-record integrated subset, adds the 30 trait-by-class counts, and reports the four mutually exclusive integrated dispositions and 0/101 endpoint.
- The plotting script validates every new display value against the released aggregate tables before rendering. Local and clean-checkout verification are required for the release commit.

## v0.2.13 release verification (2026-09-11)

- The README now starts with separate clean-clone, prepared-input and audit-only reproduction targets. It no longer presents figure generation as the implied entry point for every user.
- `DATA_SOURCES.md` records the exact 249-accession exposure series, two PGC files, FinnGen endpoint, reference builds, required schemas and non-downloadable prepared-input requirements.
- `metabolic_traits_249.tsv` maps the ordered GCST90451106–GCST90451354 range to 249 unique analysis identifiers and journal-facing display names (SHA-256 `4edabd80342a07dc0ab766f65d84d4cce6a1335a668ce10c77088a42021ea719`).
- The dependency-free input preflight checks file uniqueness, source headers, the PLINK LD prefix, indexed dbSNP resource and 249 rsID maps. Its synthetic test verifies both preflight levels and the generated 249-row exposure/250-row source manifests.
- The statistical code, released aggregate tables and figures are unchanged from v0.2.12. Local repository audit, Python compilation, R parsing, deterministic figure generation, the 6,434-row classification replay and synthetic MR/colocalization tests passed. The analysis and documentation changes also passed clean-environment GitHub Actions run 34602406692 on commit `e983c8db6ad10d7ad0cf7e91d321aa889f102603`.

## v0.2.12 local verification (2026-09-09)

- Figure 1 now uses neutral labels for the two GWAS inputs and shows forward and reverse MR as parallel direction-specific analyses with separately selected instruments.
- Only the primary forward screen selects the 15-trait candidate set. All 15 candidates enter parallel robustness, locus and directional assessments; matched reverse-MR estimates from the two complete 249-trait screens enter only the directional branch.
- The locus branch is centred and the directional branch is right-aligned. The rendered connectors do not cross, and the topology does not depict reverse MR as validation or proof of reciprocal causality.
- The 9/6 reporting-group counts, 6,434-record locus assessment, 101 records with shared-variant posterior support, 0/101 composite result and separate FinnGen comparison are unchanged.
- Local repository audit, R parsing, deterministic figure generation, colocalization-classification replay and the synthetic MR/colocalization smoke test passed. Clean-environment verification is recorded by GitHub Actions for the release commit.

## v0.2.11 local verification (2026-09-09)

- Figure 1 now uses the neutral title “Integrated evidence assessment” and removes the redundant note that the reporting groups were defined before FinnGen.
- The 9/6 reporting-group counts, 101-record locus assessment and 0/101 composite result are unchanged. FinnGen remains a separate downstream alternative-outcome comparison, and the figure still states that it did not change the reporting groups and is not an independent replication.
- Local repository audit, R parsing, deterministic figure generation and the synthetic MR/colocalization smoke test passed. The v0.2.11 tag and release must point to the verified commit whose clean-checkout GitHub Actions run passed.

## v0.2.10 local verification (2026-09-09)

- Figure 1 retains the verified v0.2.9 topology and removes no scientific content, but restructures the dense forward- and reverse-screen descriptions into shorter lines and enlarges the smallest workflow text for publication-size readability.
- The UK Biobank-excluded, locus-record and FinnGen summaries were reflowed without changing their registered counts or interpretation.
- Local repository audit, R parsing, deterministic figure generation and the synthetic MR/colocalization smoke test passed. The v0.2.10 tag and release must point to the verified commit whose clean-checkout GitHub Actions run passed.

## v0.2.9 local verification (2026-09-09)

- Figure 1 now separates the complete 249-trait forward and reverse MR screens. Only the forward screen produces the centered set of 15 Bonferroni-significant candidates.
- Robustness and sensitivity, candidate-linked directionality and locus evidence are shown as parallel assessments of the same 15 forward candidates. Colocalization is therefore positioned after candidate selection, while the corresponding subset of the two complete reverse screens is used only for directionality.
- The unnumbered stage headings and fill-only bands preserve the reading order without badge or border intersections. FinnGen remains a separate downstream alternative-outcome comparison and does not change the established reporting groups.
- Local repository audit, R parsing, deterministic figure generation and the synthetic MR/colocalization smoke test passed. The v0.2.9 tag and release must point to the verified commit whose clean-checkout GitHub Actions run passed.

## v0.2.8 local verification (2026-09-08)

- Figure 4 now combines a proportional overview of all 6,434 colocalization classifications with a magnified flow diagram for the 101 records entering integrated review and a neutral 0/101 composite-assessment endpoint.
- The plotted flow is recalculated from the released row-level aggregate tables: robust ABF records lead to 4 complex-region downgrades and 10 unassessable forward regional directions; prior-sensitive ABF records lead to 51 complex-region downgrades, 25 unassessable forward regional directions, 10 concordant but prior-sensitive records and 1 opposite regional direction.
- The plotting script writes `Figure_4_plot_data.tsv` and checks the five classification counts, four disposition counts, six cross-classification paths and zero mechanism-eligible records before rendering.
- Local repository audit, R parsing, deterministic figure generation and the synthetic MR/colocalization smoke test passed. The v0.2.8 tag and release must point to the verified commit whose clean-checkout GitHub Actions run passed.

## v0.2.7 local verification (2026-09-08)

- Figure 1 now places the separate FinnGen alternative-outcome analysis after the evidence-integration stage. The nine priority and six secondary reporting groups were not changed by the FinnGen results.
- FinnGen remains explicitly labelled as an alternative-outcome comparison rather than validation or independent replication.
- Local repository audit, R parsing, deterministic figure generation and the synthetic MR/colocalization smoke test passed. The v0.2.7 tag and release must point to the verified commit whose clean-checkout GitHub Actions run passed.

## v0.2.6 local verification (2026-09-08)

- Figure 1 retains the branching connectors but omits the redundant branch label.
- The FinnGen card states that the priority and secondary reporting groups were unchanged, so the FinnGen arrow denotes contribution to interpretation rather than group formation.
- Local repository audit, R parsing, deterministic figure generation and the synthetic MR/colocalization smoke test passed. The v0.2.6 tag and release must point to the verified commit whose clean-checkout GitHub Actions run passed.

## v0.2.5 local verification (2026-09-07)

- Figure 1 now distinguishes the primary PGC outcome from the UK Biobank-excluded sensitivity outcome and the separate FinnGen alternative outcome.
- Instrument construction, five-estimator and harmonization checks, Q/MR-Egger/MR-PRESSO diagnostics, exclusion analyses, Steiger and reverse MR, ABF/regional/complex-locus checks and exploratory SuSiE are visible in the workflow.
- Robustness, directionality, locus and FinnGen branches now feed evidence integration in parallel; none is drawn as a prerequisite for the FinnGen analysis.
- The locus-evidence text is line-wrapped inside its card. PDF, SVG and 600-dpi PNG inspection found no text outside the card boundaries.
- Local repository audit, R parsing, deterministic figure generation and the synthetic MR/colocalization smoke test passed. The v0.2.5 tag and release must point to the verified commit whose clean-checkout GitHub Actions run passed.

## v0.2.4 local verification (2026-09-07)

- Figure 2 now displays 45 estimates: primary PGC, PGC excluding UK Biobank and FinnGen R13 results for each of the 15 forward candidates.
- The 15-row `forward_noukbb_15.tsv` table matches the candidate subset of the frozen 249-trait UK Biobank-excluded forward-MR result. All 15 estimates are directionally concordant with the primary PGC estimates, all have P < 0.05 and 10 have P < 0.05/249.
- Figure 2 source rows retain dataset-specific instrument counts, beta, standard error and P value; plotted odds ratios and confidence limits are calculated as `exp(beta)` and `exp(beta ± 1.96 × SE)`.
- Local repository audit, R parsing, deterministic figure generation and the synthetic MR/colocalization smoke test passed. The v0.2.4 tag and release must point to the verified commit whose clean-checkout GitHub Actions run passed.

## v0.2.3 local verification (2026-09-07)

- The new path-parameterized locus builder reproduced all 6,434 frozen locus rows. All non-path columns and the generated relative window paths matched the released manifest exactly.
- Replaying `08_classify_coloc.R` from the released aggregate ABF/status/lead-QC inputs reproduced the frozen 6,434-row classification exactly: 14 robust, 87 prior-sensitive, 822 distinct-signal, 3,430 trait-specific/low-power and 2,081 inconclusive records.
- `Free_cholesterol_in_very_large_HDL__locus119` and `Triglycerides_to_total_lipids_ratio_in_very_small_VLDL__locus166` have leads outside chr6:25–34 Mb but windows that overlap it. Both now test as `inconclusive/mhc_excluded`.
- The three released 249-row screen tables have frozen four-level distributions of `15/53/35/146`, `122/49/9/69` and `89/67/9/84`.
- The final-IV QC table contains 249 uniquely mapped trait IDs/display names. The global minimum F is 29.7168760049588, and all 15 Bonferroni candidates have finite I²GX values.
- The 15-row broad-pleiotropy rerun retained SE and 95% CI. Its before/after beta and P columns match the prior frozen table exactly.
- A 20-instrument MR-PRESSO smoke test confirmed that `Main MR results$Sd` is parsed into raw/corrected SE and 95% CI fields. The earlier frozen 10,000-run summary did not retain `Sd`, so its released SE/CI fields remain `NA` with a reason and were not inferred from P values.
- The synthetic MR/colocalization smoke test, release-table audit, R parse check and deterministic Figure 1–4 regeneration passed in the local workspace under R 4.5.1. The v0.2.3 tag and release must point to a commit for which all jobs in `.github/workflows/reproduce.yml` have passed. The corresponding GitHub Actions run is the authoritative remote clean-checkout verification record.

## Real-data one-trait smoke tests

The portable scripts were run against the private frozen inputs for Acetate. The repository does not include original or regional SNP association-statistic rows; it releases locus-level lead identifiers and aggregate results for audit.

| Check | Instruments | Portable estimate | Frozen estimate | Result |
|---|---:|---:|---:|---|
| Forward random-effects IVW beta | 40 | -0.0432652800374973 | -0.0432652800374973 | exact match |
| Forward random-effects IVW SE | 40 | 0.0404536473106812 | 0.0404536473106812 | exact match |
| Forward random-effects IVW P | 40 | 0.284843268389097 | 0.284843268389097 | exact match |
| Reverse random-effects IVW beta | 198 | -0.0337492940727043 | -0.0337492940727043 | exact match |
| Reverse random-effects IVW SE | 198 | 0.0088452317405902 | 0.0088452317405902 | exact match |
| Reverse random-effects IVW P | 198 | 0.000135888019130649 | 0.000135888019130649 | exact match |

The hotspot-removal and Steiger outputs for Acetate also matched the frozen values, including exposure R² `0.00528761794062718` and outcome R² `0.0000505902508865417` at prevalence 0.15. Exposure R² used the fixed maximum exposure sample size of 599,249 rather than per-SNP effective N; this approximation can affect Steiger R², P values and near-boundary direction calls.

MR-PRESSO was disabled for the equality checks so the comparison isolates the deterministic estimator. The separate robustness stage fixes the MR-PRESSO seed at 20260815 and requires at least 10,000 simulations for formal candidates.

## Synthetic end-to-end smoke test

`analysis/tests/smoke_test.R` creates one 20-instrument MR trait and one 80-variant colocalization locus in a temporary directory. It verifies:

- action-2 and action-3 TwoSampleMR harmonization;
- random-effects IVW and the sensitivity-output structure;
- two case-control outcomes;
- all four p12 values for each outcome;
- finite, normalized ABF posterior probabilities;
- complete success-status rows and final SNP sets; and
- one manifest-preserving five-label classification row.

The test deletes its synthetic files on exit and does not alter the released manuscript inputs.
