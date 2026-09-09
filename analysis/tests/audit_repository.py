#!/usr/bin/env python3
"""Fail on credentials, machine paths, runtime installers or restricted file types."""
from __future__ import annotations

import re
import subprocess
import csv
import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
listed = subprocess.run(
    ["git", "ls-files", "--cached", "--others", "--exclude-standard"],
    cwd=ROOT, check=True, capture_output=True, text=True,
).stdout.splitlines()

text_suffixes = {".R", ".r", ".py", ".md", ".tsv", ".txt", ".yml", ".yaml", ".json", ".cff", ""}
restricted_suffixes = {".rds", ".bed", ".bim", ".fam", ".bcf", ".vcf", ".gz", ".tbi"}
patterns = {
    "GitHub token": re.compile(r"\bgh[pousr]_[A-Za-z0-9]{20,}\b"),
    "OpenAI key": re.compile(r"\bsk-[A-Za-z0-9_-]{20,}\b"),
    "JWT assignment": re.compile(r"(?:OPENGWAS_JWT|JWT|TOKEN)\s*(?:=|<-|:)\s*['\"][^'\"]{12,}['\"]", re.I),
    "Windows absolute path": re.compile(r"(?<![A-Za-z0-9])[A-Za-z]:[\\/](?:Users|AI|Artical|developTools|data)[\\/]", re.I),
    "Unix home path": re.compile(r"/(?:home|Users)/[^/\s]+/"),
}
forbidden_installers = ("install.packages(", "install_github(", "pak::pkg_install(")
forbidden_public_phrases = {
    "deprecated PGC outcome label": "PGC clinical" + " MDD",
    "deprecated priority label": "priority_" + "replication",
    "deprecated replication result": "not_" + "replicated",
}

errors: list[str] = []
for rel in listed:
    path = ROOT / rel
    if not path.is_file():
        continue
    if path.suffix.lower() in restricted_suffixes:
        errors.append(f"restricted binary/source-data suffix: {rel}")
    if path.suffix not in text_suffixes and path.name not in {"DESCRIPTION", "LICENSE", ".gitignore", ".gitattributes"}:
        continue
    if b"\r\n" in path.read_bytes():
        errors.append(f"text file must use LF line endings: {rel}")
    try:
        text = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        continue
    for label, pattern in patterns.items():
        if pattern.search(text):
            errors.append(f"{label}: {rel}")
    for label, phrase in forbidden_public_phrases.items():
        if phrase in text:
            errors.append(f"{label}: {rel}")
    if rel.startswith("analysis/scripts/"):
        for marker in forbidden_installers:
            if marker in text:
                errors.append(f"runtime installer {marker}: {rel}")


def read_tsv(name: str) -> list[dict[str, str]]:
    with (ROOT / "data" / "derived" / name).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def counts(rows: list[dict[str, str]], field: str) -> dict[str, int]:
    out: dict[str, int] = {}
    for row in rows:
        out[row[field]] = out.get(row[field], 0) + 1
    return out


manifest_rows = read_tsv("data_manifest.tsv")
if b"\r\n" in (ROOT / "data" / "derived" / "data_manifest.tsv").read_bytes():
    errors.append("data_manifest.tsv must use LF line endings")
manifest_names = {row["file"] for row in manifest_rows}
derived_names = {
    path.name for path in (ROOT / "data" / "derived").glob("*.tsv")
    if path.name != "data_manifest.tsv"
}
if manifest_names != derived_names:
    errors.append("data_manifest.tsv does not enumerate every released TSV exactly once")
for row in manifest_rows:
    path = ROOT / "data" / "derived" / row["file"]
    if not path.is_file():
        errors.append(f"manifest file missing: {row['file']}")
        continue
    content = path.read_bytes()
    if b"\r\n" in content:
        errors.append(f"released TSV must use LF line endings: {row['file']}")
    if hashlib.md5(content).hexdigest() != row["md5"].lower():
        errors.append(f"MD5 mismatch: {row['file']}")
    if hashlib.sha256(content).hexdigest() != row["sha256"].lower():
        errors.append(f"SHA-256 mismatch: {row['file']}")
    with path.open(encoding="utf-8", newline="") as handle:
        n_rows = sum(1 for _ in csv.DictReader(handle, delimiter="\t"))
    if n_rows != int(row["rows"]):
        errors.append(f"row-count mismatch: {row['file']}")

description = (ROOT / "DESCRIPTION").read_text(encoding="utf-8")
citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
if "Version: 0.2.12" not in description or "version: 0.2.12" not in citation:
    errors.append("DESCRIPTION and CITATION.cff must both declare version 0.2.12")

readme = (ROOT / "README.md").read_text(encoding="utf-8")
workflow = (ROOT / "analysis" / "WORKFLOW.md").read_text(encoding="utf-8")
verification = (ROOT / "analysis" / "VERIFICATION.md").read_text(encoding="utf-8")
access = (ROOT / "data" / "ACCESS.md").read_text(encoding="utf-8")
figure_script = (ROOT / "scripts" / "make_figures.R").read_text(encoding="utf-8")
if ("Version 0.2.12" not in readme or
        "Release v0.2.12" not in workflow or
        "## v0.2.12 local verification" not in verification):
    errors.append("release-candidate version is not synchronized across repository documentation")
mermaid_start = figure_script.find("mermaid <- c(")
mermaid_end = figure_script.find("\n)\nstopifnot(", mermaid_start)
if mermaid_start < 0 or mermaid_end < 0:
    errors.append("Figure 1 Mermaid source block is missing")
    mermaid_block = ""
else:
    mermaid_block = figure_script[mermaid_start:mermaid_end]
for phrase in (
    "Circulating metabolic-trait GWAS",
    "PGC major-depression GWAS",
    "Parallel direction-specific MR analyses",
    "Exposure-specific instruments selected separately in each direction",
    "15 forward-selected traits",
    "Primary forward screen only",
    "Reverse MR for 15 forward-selected traits",
    "Both complete 249-trait screens",
    "Same 15 forward candidates",
    "Cochran's Q, MR-Egger and MR-PRESSO",
    "trait-specific analysis-window records",
    "reporting groups unchanged",
    '0.53, "Locus evidence"',
    '0.855, "Directional evidence assessment"',
    "connector(0.205, 0.530, 0.760, 0.530, end = FALSE)",
    "reverse_path_x <- c(0.72, 0.72, 0.90, 0.90)",
    "reverse_path_y <- c(0.603, 0.548, 0.548, 0.515)",
    'G[\\"Integrated evidence assessment<br/>',
):
    if phrase not in figure_script:
        errors.append(f"Figure 1 is missing a required scope marker: {phrase}")
for edge in (
    '"  E --> FW"',
    '"  O --> FW"',
    '"  E --> RV"',
    '"  O --> RV"',
    '"  FW --> C"',
    '"  RV --> RM"',
    '"  C --> R"',
    '"  C --> L"',
    '"  C --> D"',
    '"  RM --> D"',
    '"  R --> G"',
    '"  L --> G"',
    '"  D --> G"',
    '"  G -.-> F"',
):
    if edge not in mermaid_block:
        errors.append(f"Figure 1 is missing a required topology edge: {edge}")
if any(
    phrase in figure_script
    for phrase in ("Pre-FinnGen evidence synthesis", "Groups defined before FinnGen")
):
    errors.append("Figure 1 must use neutral evidence-synthesis wording")
for forbidden_edge in ('"  RV --> C"', '"  RM --> R"', '"  RM --> L"', '"  F --> G"'):
    if forbidden_edge in mermaid_block:
        errors.append(f"Figure 1 contains a forbidden topology edge: {forbidden_edge}")
if "stage_label" in figure_script or "grid.circle" in figure_script:
    errors.append("Figure 1 must use unnumbered stage headings without circular badges")
attestation = "**Responsible-author attestation:** CONFIRMED by Zhouyi Wang on 2026-09-07 for the exact current release schemas and publication boundary stated in this file."
if attestation not in access or "attestation remains pending" in readme or "attestation remains pending" in access:
    errors.append("responsible-author release-scope attestation is missing or still marked pending")

lock = json.loads((ROOT / "renv.lock").read_text(encoding="utf-8"))
locked = lock.get("Packages", {})
required_roots = {
    "data.table", "ggplot2", "patchwork", "scales", "MendelianRandomization",
    "ieugwasr", "coloc", "future", "future.apply", "TwoSampleMR", "MRPRESSO",
}
if len(locked) < 100 or not required_roots.issubset(locked):
    errors.append("renv.lock does not contain the verified hard-dependency closure")
for package, record in locked.items():
    if not record.get("Version") or not record.get("Source"):
        errors.append(f"renv.lock has incomplete metadata for {package}")
    if record.get("Source") == "GitHub" and not re.fullmatch(r"[0-9a-f]{40}", record.get("RemoteSha", "")):
        errors.append(f"renv.lock GitHub package is not pinned to an immutable SHA: {package}")
if locked.get("TwoSampleMR", {}).get("RemoteSha") != "3d119f20d6fc164b0c7f710f5590fee9580f2c7b":
    errors.append("TwoSampleMR v0.7.9 must be pinned to its immutable tag commit")

screen_expectations = {
    "forward_screen_249.tsv": {"bonferroni_hit": 15, "fdr_only": 53, "nominal": 35, "null": 146},
    "reverse_screen_main_249.tsv": {"bonferroni_hit": 122, "fdr_only": 49, "nominal": 9, "null": 69},
    "reverse_screen_noukbb_249.tsv": {"bonferroni_hit": 89, "fdr_only": 67, "nominal": 9, "null": 84},
}
for name, expected in screen_expectations.items():
    rows = read_tsv(name)
    if len(rows) != 249 or len({row["trait"] for row in rows}) != 249:
        errors.append(f"{name} must contain 249 unique traits")
    if counts(rows, "screen_level") != expected:
        errors.append(f"{name} four-level counts differ from the frozen values")

forward_main = read_tsv("forward_screen_249.tsv")
forward_noukbb = read_tsv("forward_noukbb_15.tsv")
forest = read_tsv("forest_estimates.tsv")
candidate_order = [row["trait"] for row in forest]
bonferroni_candidates = {
    row["trait"] for row in forward_main if row["screen_level"] == "bonferroni_hit"
}
if (len(forward_noukbb) != 15 or len({row["trait"] for row in forward_noukbb}) != 15 or
        [row["trait"] for row in forward_noukbb] != candidate_order or
        set(candidate_order) != bonferroni_candidates):
    errors.append("UK Biobank-excluded forward results must match the 15 primary candidates in frozen order")
forest_by_trait = {row["trait"]: row for row in forest}
noukbb_nominal = 0
noukbb_bonferroni = 0
for row in forward_noukbb:
    try:
        n_iv = int(row["n_iv"])
        beta = float(row["ivw_b"])
        se = float(row["ivw_se"])
        p_value = float(row["ivw_p"])
        expected_p = math.erfc(abs(beta / se) / math.sqrt(2))
        main_beta = float(forest_by_trait[row["trait"]]["b_pgc"])
    except (KeyError, TypeError, ValueError):
        errors.append(f"invalid UK Biobank-excluded forward row: {row.get('trait', '<missing>')}")
        continue
    if not (n_iv > 0 and math.isfinite(beta) and math.isfinite(se) and se > 0 and
            math.isfinite(p_value) and 0 < p_value <= 1 and abs(p_value - expected_p) < 1e-12):
        errors.append(f"invalid UK Biobank-excluded estimate or P value: {row['trait']}")
    if beta * main_beta <= 0:
        errors.append(f"UK Biobank-excluded direction differs from the primary estimate: {row['trait']}")
    noukbb_nominal += p_value < 0.05
    noukbb_bonferroni += p_value < 0.05 / 249
if noukbb_nominal != 15 or noukbb_bonferroni != 10:
    errors.append("UK Biobank-excluded sensitivity counts must remain 15 nominal and 10 at 0.05/249")

locus_manifest = read_tsv("coloc_locus_manifest.tsv")
classification = read_tsv("coloc_classification.tsv")
expected_classes = {
    "robust_coloc": 14, "prior_sensitive_coloc": 87, "distinct_signal": 822,
    "trait_specific_or_low_power": 3430, "inconclusive": 2081,
}
if len(locus_manifest) != 6434 or len(classification) != 6434:
    errors.append("colocalization manifest and classification must each contain 6,434 rows")
if counts(classification, "abf_class") != expected_classes:
    errors.append("five-class colocalization counts differ from 14/87/822/3430/2081")
manifest_by_key = {row["locus_key"]: row for row in locus_manifest}
class_by_key = {row["locus_key"]: row for row in classification}
boundary_keys = (
    "Free_cholesterol_in_very_large_HDL__locus119",
    "Triglycerides_to_total_lipids_ratio_in_very_small_VLDL__locus166",
)
for key in boundary_keys:
    m = manifest_by_key.get(key)
    c = class_by_key.get(key)
    if not m or not c:
        errors.append(f"MHC boundary record missing: {key}")
        continue
    lead = int(float(m["lead_pos"]))
    overlaps = int(m["chr"]) == 6 and float(m["window_end"]) >= 25_000_000 and float(m["window_start"]) <= 34_000_000
    lead_outside = lead < 25_000_000 or lead > 34_000_000
    if not (overlaps and lead_outside and c["abf_class"] == "inconclusive" and
            c["failure_reason"] == "mhc_excluded" and c["is_mhc"].upper() == "TRUE"):
        errors.append(f"MHC window-overlap boundary assertion failed: {key}")

integrated = read_tsv("integrated_mechanism_evidence.tsv")
if len(integrated) != 101 or any(row["mechanism_eligible"].upper() != "FALSE" for row in integrated):
    errors.append("integrated evidence must contain 101 rows and zero mechanism-eligible records")
expected_tiers = {
    "abf_label_only_complex_downgraded": 55,
    "abf_label_only_forward_unassessable": 35,
    "abf_prior_sensitive_forward_consistent": 10,
    "abf_label_only_forward_opposite": 1,
}
if counts(integrated, "evidence_tier") != expected_tiers:
    errors.append("integrated-evidence tier counts differ from 55/35/10/1")
expected_flows = {
    ("robust_coloc", "abf_label_only_complex_downgraded"): 4,
    ("robust_coloc", "abf_label_only_forward_unassessable"): 10,
    ("prior_sensitive_coloc", "abf_label_only_complex_downgraded"): 51,
    ("prior_sensitive_coloc", "abf_label_only_forward_unassessable"): 25,
    ("prior_sensitive_coloc", "abf_prior_sensitive_forward_consistent"): 10,
    ("prior_sensitive_coloc", "abf_label_only_forward_opposite"): 1,
}
observed_flows = {}
for row in integrated:
    key = (row["abf_class"], row["evidence_tier"])
    observed_flows[key] = observed_flows.get(key, 0) + 1
if observed_flows != expected_flows:
    errors.append("ABF-class to integrated-disposition counts differ from 4/10/51/25/10/1")

allowed_support = {"FDR-supported", "nominally-supported", "no-nominal-support"}
legacy_support_markers = (
    "priority_" + "replication", "replication_" + "label",
    "not_" + "replicated", "str" + "ong", "supp" + "ortive",
)
for name, field in (
    ("evidence_summary.tsv", "cross_outcome_support_label"),
    ("forest_estimates.tsv", "cross_outcome_support_label"),
    ("finngen_cross_outcome_15.tsv", "cross_outcome_support_label"),
):
    rows = read_tsv(name)
    if set(row[field] for row in rows) - allowed_support:
        errors.append(f"{name} contains a non-neutral cross-outcome support label")
    raw = (ROOT / "data" / "derived" / name).read_text(encoding="utf-8")
    if any(marker in raw for marker in legacy_support_markers):
        errors.append(f"legacy replication semantics remain in {name}")
if counts(read_tsv("finngen_cross_outcome_15.tsv"), "cross_outcome_support_label") != {
    "FDR-supported": 11, "nominally-supported": 1, "no-nominal-support": 3,
}:
    errors.append("FinnGen neutral support labels must have the frozen 11/1/3 distribution")

tool_loss = read_tsv("finngen_tool_loss.tsv")
expected_tool_totals = {
    "n_iv": 5665, "n_matched": 5516, "n_not_found": 149,
    "n_harmonise_used": 5355, "removed": 161,
}
if len(tool_loss) != 15:
    errors.append("FinnGen tool-loss table must contain 15 traits")
for field, expected in expected_tool_totals.items():
    if sum(int(row[field]) for row in tool_loss) != expected:
        errors.append(f"FinnGen {field} total differs from {expected}")

instrument_qc = read_tsv("instrument_strength_qc.tsv")
if len(instrument_qc) != 249 or len({row["trait_id"] for row in instrument_qc}) != 249 or \
        len({row["trait_display"] for row in instrument_qc}) != 249:
    errors.append("instrument-strength QC must map 249 unique machine and display names")
if min(float(row["min_f"]) for row in instrument_qc) < 10:
    errors.append("instrument-strength QC contains F below 10")
candidate_qc = [row for row in instrument_qc if row["screen_level"] == "bonferroni_hit"]
if len(candidate_qc) != 15 or any(not math.isfinite(float(row["i2gx"])) for row in candidate_qc):
    errors.append("all 15 forward candidates require finite I2GX values")
if any(row["ld_panel_missing_fraction"] != "NA" or not row["ld_panel_missingness_note"] for row in instrument_qc):
    errors.append("unavailable full-family LD-panel missingness must remain explicit NA")

pleio = read_tsv("broad_pleiotropy_sensitivity_15.tsv")
if len(pleio) != 15:
    errors.append("broad-pleiotropy sensitivity table must contain 15 traits")
for row in pleio:
    for prefix in ("full", "after"):
        b = float(row[f"{prefix}_b"])
        se = float(row[f"{prefix}_se"])
        if not (math.isfinite(se) and abs(float(row[f"{prefix}_ci_lo"]) - (b - 1.96 * se)) < 1e-12 and
                abs(float(row[f"{prefix}_ci_hi"]) - (b + 1.96 * se)) < 1e-12):
            errors.append(f"broad-pleiotropy SE/CI assertion failed: {row['trait_id']} {prefix}")

presso = read_tsv("presso_sensitivity_15.tsv")
if len(presso) != 15 or any(row["presso_raw_se"] != "NA" or row["presso_corrected_se"] != "NA"
                            for row in presso):
    errors.append("frozen MR-PRESSO SE fields must remain NA until the 10,000-run objects are rerun")
if any("not reconstructed from P values" not in row["presso_se_status"] for row in presso):
    errors.append("MR-PRESSO NA reason must prohibit reconstruction from P values")

for name, expected_rows in (
    ("susie_exploratory_summary.tsv", 44),
    ("susie_credible_set_binding.tsv", 334),
    ("susie_stop_rule_verdicts.tsv", 7),
):
    if len(read_tsv(name)) != expected_rows:
        errors.append(f"{name} row count differs from {expected_rows}")

if errors:
    raise SystemExit("Repository audit failed:\n- " + "\n- ".join(sorted(set(errors))))
print(f"PASS: audited {len(listed)} repository files")
