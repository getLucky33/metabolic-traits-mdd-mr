# Verification record

Verification history: real-data one-trait checks were completed on 2026-09-06; v0.2.3, v0.2.4 and v0.2.5 local release checks were completed on 2026-09-07; v0.2.6 and v0.2.7 local release checks were completed on 2026-09-08. Runtime: R 4.5.1 on Windows 10 x64 using the recorded environment.

## v0.2.7 local verification (2026-09-08)

- Figure 1 now places the separate FinnGen alternative-outcome analysis after the evidence-integration stage. This matches the documented timing: the nine priority and six secondary reporting groups were defined before FinnGen and were not changed by its results.
- FinnGen remains explicitly labelled as an alternative-outcome comparison rather than validation or independent replication.
- Local repository audit, R parsing, deterministic figure generation and the synthetic MR/colocalization smoke test passed. The v0.2.7 tag and release must point to the verified commit whose clean-checkout GitHub Actions run passed.

## v0.2.6 local verification (2026-09-08)

- Figure 1 retains the branching connectors but omits the redundant branch label.
- The FinnGen card now states that the priority and secondary reporting groups were unchanged. The evidence-integration card identifies those groups as having been defined before FinnGen, so the FinnGen arrow denotes contribution to interpretation rather than group formation.
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
