# Verification record

Verification date: 2026-09-06. Runtime: R 4.5.1 on Windows 10 x64 using the versions in `analysis/environment/package-versions.tsv`.

## Real-data one-trait smoke tests

The portable scripts were run against the private frozen inputs for Acetate. No private input or SNP-level output is included in this repository.

| Check | Instruments | Portable estimate | Frozen estimate | Result |
|---|---:|---:|---:|---|
| Forward random-effects IVW beta | 40 | -0.0432652800374973 | -0.0432652800374973 | exact match |
| Forward random-effects IVW SE | 40 | 0.0404536473106812 | 0.0404536473106812 | exact match |
| Forward random-effects IVW P | 40 | 0.284843268389097 | 0.284843268389097 | exact match |
| Reverse random-effects IVW beta | 198 | -0.0337492940727043 | -0.0337492940727043 | exact match |
| Reverse random-effects IVW SE | 198 | 0.0088452317405902 | 0.0088452317405902 | exact match |
| Reverse random-effects IVW P | 198 | 0.000135888019130649 | 0.000135888019130649 | exact match |

The hotspot-removal and Steiger outputs for Acetate also matched the frozen values, including exposure R² `0.00528761794062718` and outcome R² `0.0000505902508865417` at prevalence 0.15.

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
