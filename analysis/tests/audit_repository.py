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
from statistics import NormalDist


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


def tsv_header(name: str) -> list[str]:
    with (ROOT / "data" / "derived" / name).open(encoding="utf-8", newline="") as handle:
        return next(csv.reader(handle, delimiter="\t"))


def counts(rows: list[dict[str, str]], field: str) -> dict[str, int]:
    out: dict[str, int] = {}
    for row in rows:
        out[row[field]] = out.get(row[field], 0) + 1
    return out


def student_t_critical_975(degrees_freedom: int) -> float:
    z = NormalDist().inv_cdf(0.975)
    inverse_df = 1.0 / degrees_freedom
    return (
        z
        + (z**3 + z) * inverse_df / 4
        + (5 * z**5 + 16 * z**3 + 3 * z) * inverse_df**2 / 96
        + (3 * z**7 + 19 * z**5 + 17 * z**3 - 15 * z) * inverse_df**3 / 384
        + (79 * z**9 + 776 * z**7 + 1482 * z**5 - 1920 * z**3 - 945 * z) * inverse_df**4 / 92160
    )


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
if "Version: 0.2.19" not in description or "version: 0.2.19" not in citation:
    errors.append("DESCRIPTION and CITATION.cff must both declare version 0.2.19")

readme = (ROOT / "README.md").read_text(encoding="utf-8")
workflow = (ROOT / "analysis" / "WORKFLOW.md").read_text(encoding="utf-8")
verification = (ROOT / "analysis" / "VERIFICATION.md").read_text(encoding="utf-8")
access = (ROOT / "data" / "ACCESS.md").read_text(encoding="utf-8")
figure_script = (ROOT / "scripts" / "make_figures.R").read_text(encoding="utf-8")
if ("Version 0.2.19" not in readme or
        "Release v0.2.19" not in workflow or
        "## v0.2.19 release verification" not in verification):
    errors.append("release version is not synchronized across repository documentation")
for rel in (
    "DATA_SOURCES.md",
    "analysis/config/local_paths.example.tsv",
    "analysis/manifests/metabolic_traits_249.tsv",
    "analysis/prepare_inputs.py",
    "analysis/scripts/00_map_exposure_rsids.py",
    "analysis/scripts/05a_recover_exposure_variant_n.py",
    "analysis/tests/test_exposure_rsid_mapping.py",
    "analysis/tests/test_exposure_variant_n.py",
):
    if not (ROOT / rel).is_file():
        errors.append(f"reproduction-entry file is missing: {rel}")
catalog: list[dict[str, str]] = []
catalog_file = ROOT / "analysis" / "manifests" / "metabolic_traits_249.tsv"
if catalog_file.is_file():
    if hashlib.sha256(catalog_file.read_bytes()).hexdigest() != "4edabd80342a07dc0ab766f65d84d4cce6a1335a668ce10c77088a42021ea719":
        errors.append("metabolic accession-to-trait mapping differs from the verified catalog")
    with catalog_file.open(encoding="utf-8", newline="") as handle:
        catalog = list(csv.DictReader(handle, delimiter="\t"))
    expected_accessions = [f"GCST{number}" for number in range(90451106, 90451355)]
    labels = {row["trait_id"]: row["trait_display"] for row in read_tsv("trait_display_dictionary.tsv")}
    if (len(catalog) != 249 or [row.get("accession") for row in catalog] != expected_accessions or
            len({row.get("trait") for row in catalog}) != 249 or
            any(row.get("genome_build") != "GRCh38" or row.get("sample_size") != "599249" for row in catalog)):
        errors.append("metabolic accession catalog must preserve the ordered 249-trait meta_EUR series")
    if {row.get("trait"): row.get("trait_display") for row in catalog} != labels:
        errors.append("metabolic accession catalog must match the frozen trait/display dictionary")
for marker in (
    "Choose a reproduction target",
    "Start here for a real-data rerun",
    "Exposure rsID maps can now be generated locally",
    "Aggregate-result audit only",
    "Generate one map first as the required real-data smoke test",
):
    if marker not in readme:
        errors.append(f"README is missing a reproduction-boundary marker: {marker}")
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
    "posterior support favoring H4 (",
    "robust; ",
    "prior-sensitive)",
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
attestation = "**Responsible-author attestation:** CONFIRMED by Zhouyi Wang on 2026-09-14 for the aggregate-only publication boundary stated in this file. Version 0.2.19 adds trait-level selection counts without variant identifiers and remains within that boundary."
if attestation not in access:
    errors.append("responsible-author aggregate release-scope attestation is missing")
manifest_sha256 = hashlib.sha256(
    (ROOT / "data" / "derived" / "data_manifest.tsv").read_bytes()
).hexdigest()
if f"`{manifest_sha256}`" not in access:
    errors.append("ACCESS.md attested manifest snapshot does not match data_manifest.tsv")
transition_markers = (
    "OPEN" + "/PENDING",
    "subject to renewed responsible-author approval",
    "Pending v0.2.19",
)
for name, document in (
        ("README.md", readme), ("analysis/WORKFLOW.md", workflow),
        ("analysis/VERIFICATION.md", verification), ("data/ACCESS.md", access)):
    if any(marker in document for marker in transition_markers):
        errors.append(f"transitional release-gate wording remains in {name}")

variant_n_script = (ROOT / "analysis" / "scripts" / "05a_recover_exposure_variant_n.py").read_text(encoding="utf-8")
for marker in ("--exposure-manifest", "--iv-manifest", "--out", "--qc-out",
               "--traits", "--workers", "--max-n", "gzip.open", "ProcessPoolExecutor"):
    if marker not in variant_n_script:
        errors.append(f"portable variant-N script is missing required behavior: {marker}")
if re.search(r"(?i)[A-Z]:[\\/]", variant_n_script):
    errors.append("portable variant-N script contains a machine-specific path")
prepare_script = (ROOT / "analysis" / "prepare_inputs.py").read_text(encoding="utf-8")
if not re.search(r'METABOLIC_COLUMNS\s*=\s*\{[^}]*["\']n["\']', prepare_script, re.S):
    errors.append("source preflight must require the metabolic GWAS n field")
for marker in ("validate_rsid_map", "RSID_STATUSES", "no valid uniquely matched rsID"):
    if marker not in prepare_script:
        errors.append(f"source preflight is missing rsID-map content validation: {marker}")
mapping_script = (ROOT / "analysis" / "scripts" / "00_map_exposure_rsids.py").read_text(encoding="utf-8")
for marker in (
    "EXPECTED_DBSNP_MD5", "--traits", "--batch-size", "--p-threshold", "allele_set_match",
    "no_allele_match", "no_dbsnp_entry", "source_decompressed_sha256", "map_sha256",
):
    if marker not in mapping_script:
        errors.append(f"exposure rsID mapper is missing required behavior: {marker}")
workflow_yaml = (ROOT / ".github" / "workflows" / "reproduce.yml").read_text(encoding="utf-8")
if "python analysis/tests/test_exposure_rsid_mapping.py" not in workflow_yaml:
    errors.append("GitHub Actions does not run the exposure rsID-mapping test")
robustness_script = (ROOT / "analysis" / "scripts" / "05_run_robustness.R").read_text(encoding="utf-8")
if ("Formal robustness mode requires --exposure-n-file" not in robustness_script or
        "source_variant_specific" not in robustness_script or
        re.search(r'arg_value\(args,\s*"exposure-n"', robustness_script)):
    errors.append("formal robustness code must require source-variant exposure N without a fixed-N fallback")

lock = json.loads((ROOT / "renv.lock").read_text(encoding="utf-8"))
locked = lock.get("Packages", {})
required_roots = {
    "data.table", "ggplot2", "patchwork", "scales", "MendelianRandomization",
    "ieugwasr", "coloc", "future", "future.apply", "TwoSampleMR", "MRPRESSO", "psych",
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
selection_qc = read_tsv("instrument_selection_ld_panel_qc_249.tsv")
if (tsv_header("instrument_selection_ld_panel_qc_249.tsv") !=
        ["trait", "autosomal_biallelic", "in_ld_panel", "ld_panel_missing"]):
    errors.append("LD-panel selection QC must preserve its four-column aggregate schema")
if len(selection_qc) != 249 or len({row["trait"] for row in selection_qc}) != 249:
    errors.append("LD-panel selection QC must contain 249 unique traits")
selection_by_trait = {row["trait"]: row for row in selection_qc}
selection_counts_valid = True
for row in selection_qc:
    try:
        eligible = int(row["autosomal_biallelic"])
        in_panel = int(row["in_ld_panel"])
        missing = int(row["ld_panel_missing"])
    except (KeyError, TypeError, ValueError):
        selection_counts_valid = False
        continue
    if eligible <= 0 or in_panel < 0 or missing < 0 or in_panel + missing != eligible:
        selection_counts_valid = False
if not selection_counts_valid:
    errors.append("LD-panel selection QC contains invalid or internally inconsistent counts")
if selection_counts_valid and (
        sum(int(row["autosomal_biallelic"]) for row in selection_qc) != 15105492 or
        sum(int(row["ld_panel_missing"]) for row in selection_qc) != 961597):
    errors.append("LD-panel selection QC aggregate counts differ from the verified 249-trait run")
if len(instrument_qc) != 249 or len({row["trait_id"] for row in instrument_qc}) != 249 or \
        len({row["trait_display"] for row in instrument_qc}) != 249:
    errors.append("instrument-strength QC must map 249 unique machine and display names")
if min(float(row["min_f"]) for row in instrument_qc) < 10:
    errors.append("instrument-strength QC contains F below 10")
candidate_qc = [row for row in instrument_qc if row["screen_level"] == "bonferroni_hit"]
if len(candidate_qc) != 15 or any(not math.isfinite(float(row["i2gx"])) for row in candidate_qc):
    errors.append("all 15 forward candidates require finite I2GX values")
for row in instrument_qc:
    source = selection_by_trait.get(row["trait_id"])
    try:
        missing = int(row["ld_panel_missing_n"])
        eligible = int(row["ld_panel_eligible_n"])
        fraction = float(row["ld_panel_missing_fraction"])
        source_missing = int(source["ld_panel_missing"]) if source else -1
        source_eligible = int(source["autosomal_biallelic"]) if source else -1
    except (KeyError, TypeError, ValueError):
        errors.append(f"instrument-strength LD-panel QC is nonnumeric: {row.get('trait_id', '<missing>')}")
        continue
    if (missing != source_missing or eligible != source_eligible or eligible <= 0 or
            not math.isfinite(fraction) or abs(fraction - missing / eligible) > 1e-14 or
            row["ld_panel_missingness_note"] !=
            "Computed from the complete 249-trait LD-panel selection QC table"):
        errors.append(f"instrument-strength LD-panel QC differs from its selection source: {row['trait_id']}")

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

loo = read_tsv("leave_one_out_summary_15.tsv")
labels_by_trait = {row["trait_id"]: row["trait_display"] for row in read_tsv("trait_labels.tsv")}
if (len(loo) != 15 or len({row["trait_id"] for row in loo}) != 15 or
        [row["trait_id"] for row in loo] != candidate_order):
    errors.append("leave-one-out summary must contain the 15 candidates in frozen order")
if len(loo) == 15:
    for row in loo:
        trait = row["trait_id"]
        try:
            n_estimates = int(row["n_leave_one_out_estimates"])
            min_beta = float(row["min_beta"])
            max_beta = float(row["max_beta"])
            max_p = float(row["max_p"])
            main_beta = float(next(item["ivw_b"] for item in forward_main if item["trait"] == trait))
            main_n = int(next(item["n_iv"] for item in forward_main if item["trait"] == trait))
        except (KeyError, TypeError, ValueError, StopIteration):
            errors.append(f"invalid leave-one-out summary row: {trait}")
            continue
        if (n_estimates != main_n or row["trait_display"] != labels_by_trait.get(trait) or
                not all(math.isfinite(value) for value in (min_beta, max_beta, max_p)) or
                min_beta > max_beta or max_p >= 0.05 or min_beta * main_beta <= 0 or
                max_beta * main_beta <= 0 or row["direction_change_n"] != "0" or
                row["nominal_loss_n"] != "0" or row["all_primary_direction"] != "TRUE" or
                row["all_nominal_p_lt_0_05"] != "TRUE"):
            errors.append(f"leave-one-out stability assertion failed: {trait}")
    if sum(int(row["n_leave_one_out_estimates"]) for row in loo) != 4993:
        errors.append("leave-one-out summary must contain 4,993 single-variant deletions")

expected_presso_header = [
    "trait_id", "trait_display", "n_iv", "presso_global_p", "presso_outlier_n",
    "presso_distortion_p", "presso_raw_b", "presso_raw_se", "presso_raw_df",
    "presso_raw_ci_lo", "presso_raw_ci_hi", "presso_raw_p", "presso_corrected_b",
    "presso_corrected_se", "presso_corrected_df", "presso_corrected_ci_lo",
    "presso_corrected_ci_hi", "presso_corrected_p", "presso_nb", "presso_seed",
    "presso_se_status", "presso_global_p_display", "presso_global_rssobs",
    "presso_global_exceedance_n", "presso_distortion_p_display",
    "presso_distortion_coefficient", "presso_raw_t_stat", "presso_corrected_t_stat",
    "presso_regression_p_sidedness", "presso_empirical_test_df",
]
presso = read_tsv("presso_sensitivity_15.tsv")
if tsv_header("presso_sensitivity_15.tsv") != expected_presso_header:
    errors.append("MR-PRESSO sensitivity table schema or column order has drifted")
if (len(presso) != 15 or len({row["trait_id"] for row in presso}) != 15 or
        {row["trait_id"] for row in presso} != set(candidate_order)):
    errors.append("MR-PRESSO sensitivity table must contain the 15 frozen candidates")
for row in presso:
    try:
        n_iv = int(row["n_iv"])
        outlier_n = int(row["presso_outlier_n"])
        if (row["trait_display"] != labels_by_trait.get(row["trait_id"]) or
                row["presso_nb"] != "10000" or row["presso_seed"] != "20260815"):
            raise ValueError("simulation metadata mismatch")
        for prefix, expected_df in (
                ("presso_raw", n_iv - 1),
                ("presso_corrected", n_iv - outlier_n - 1)):
            estimate = float(row[f"{prefix}_b"])
            se = float(row[f"{prefix}_se"])
            degrees_freedom = int(row[f"{prefix}_df"])
            lower = float(row[f"{prefix}_ci_lo"])
            upper = float(row[f"{prefix}_ci_hi"])
            if se <= 0 or degrees_freedom != expected_df:
                raise ValueError("invalid SE or residual degrees of freedom")
            critical = student_t_critical_975(degrees_freedom)
            tolerance = max(1e-10, abs(estimate) * 1e-8, abs(se) * 1e-8)
            if (abs(lower - (estimate - critical * se)) > tolerance or
                    abs(upper - (estimate + critical * se)) > tolerance):
                raise ValueError("confidence interval mismatch")
            if not math.isclose(float(row[f"{prefix}_t_stat"]), estimate / se,
                                rel_tol=1e-12, abs_tol=1e-12):
                raise ValueError("t statistic mismatch")
        if (row["presso_global_p_display"] != "<1e-04" or
                row["presso_global_exceedance_n"] != "0" or
                row["presso_distortion_p_display"] != row["presso_distortion_p"] or
                row["presso_regression_p_sidedness"] != "two-sided" or
                row["presso_empirical_test_df"] != "not_applicable" or
                not math.isfinite(float(row["presso_global_rssobs"])) or
                float(row["presso_global_rssobs"]) <= 0 or
                not math.isfinite(float(row["presso_distortion_coefficient"]))):
            raise ValueError("saved-object test metadata mismatch")
    except (KeyError, TypeError, ValueError):
        errors.append(f"MR-PRESSO uncertainty fields are invalid: {row.get('trait_id', '<missing>')}")
if any("directly" not in row["presso_se_status"].lower() or
       "not reconstructed from P values" not in row["presso_se_status"] for row in presso):
    errors.append("MR-PRESSO SE/CI provenance must identify direct Sd retention and prohibit P-value reconstruction")
if len(presso) == 15:
    if (sum(math.isclose(float(row["presso_global_p"]), 1e-4, rel_tol=0, abs_tol=1e-12) for row in presso) != 15 or
            min(int(row["presso_outlier_n"]) for row in presso) != 7 or
            max(int(row["presso_outlier_n"]) for row in presso) != 19 or
            sum(float(row["presso_distortion_p"]) >= 0.05 for row in presso) != 13 or
            sum(float(row["presso_distortion_p"]) < 0.05 for row in presso) != 2):
        errors.append("MR-PRESSO global P values, outlier counts or 13/2 distortion split differ from the frozen results")

steiger_header = [
    "trait_id", "trait_display", "source_accession", "n_iv_full",
    "n_iv_pleio_removed", "exposure_n_mode", "exposure_n_min", "exposure_n_max",
    "exposure_n_unique", "steiger_correct_008", "steiger_p_008", "steiger_z_008",
    "steiger_log10_p_008", "steiger_correct_015", "steiger_p_015", "steiger_z_015",
    "steiger_log10_p_015", "steiger_correct_020", "steiger_p_020", "steiger_z_020",
    "steiger_log10_p_020", "steiger_test_distribution", "steiger_p_sidedness",
    "steiger_df", "steiger_p_display", "r2_exposure", "r2_outcome",
]
qc_header = [
    "trait_id", "trait_display", "source_accession", "target_iv_n", "matched_iv_n",
    "missing_iv_n", "exposure_n_min", "exposure_n_max", "exposure_n_unique",
    "exposure_n_values", "source_rows_scanned", "compressed_bytes_read", "elapsed_seconds",
]
steiger = read_tsv("steiger_directionality_15.tsv")
sample_qc = read_tsv("exposure_sample_size_qc_15.tsv")
if tsv_header("steiger_directionality_15.tsv") != steiger_header:
    errors.append("Steiger directionality table schema or column order has drifted")
if tsv_header("exposure_sample_size_qc_15.tsv") != qc_header:
    errors.append("exposure sample-size QC table schema or column order has drifted")
for name, rows in (("steiger_directionality_15.tsv", steiger),
                   ("exposure_sample_size_qc_15.tsv", sample_qc)):
    forbidden_columns = {"rsid", "variant_id", "source_gwas_file", "source_file", "iv_file"}
    if forbidden_columns.intersection(tsv_header(name)):
        errors.append(f"candidate-level release contains a variant/path column: {name}")
    values = "\n".join(value for row in rows for value in row.values())
    if re.search(r"(?i)(?:[A-Z]:[\\/]|/(?:home|Users|mnt|tmp)/)", values):
        errors.append(f"candidate-level release contains a local path: {name}")
if ([row.get("trait_id") for row in steiger] != candidate_order or
        [row.get("trait_id") for row in sample_qc] != candidate_order):
    errors.append("Steiger and sample-size QC tables must contain the 15 candidates in frozen order")
accession_by_trait = {row["trait"]: row["accession"] for row in catalog}
presso_by_trait = {row["trait_id"]: row for row in presso}
steiger_by_trait = {row["trait_id"]: row for row in steiger}
for row in steiger:
    trait = row.get("trait_id", "<missing>")
    try:
        n_full = int(row["n_iv_full"])
        n_removed = int(row["n_iv_pleio_removed"])
        n_min, n_max = int(row["exposure_n_min"]), int(row["exposure_n_max"])
        n_unique = int(row["exposure_n_unique"])
        r2_exposure, r2_outcome = float(row["r2_exposure"]), float(row["r2_outcome"])
        if (row["trait_display"] != labels_by_trait.get(trait) or
                row["source_accession"] != accession_by_trait.get(trait) or
                n_full != int(presso_by_trait[trait]["n_iv"]) or
                not 0 < n_removed < n_full or
                row["exposure_n_mode"] != "source_variant_specific" or
                n_min not in {413897, 599249} or n_max != 599249 or
                n_unique != (1 if n_min == n_max else 2) or
                row["steiger_test_distribution"] != "standard_normal" or
                row["steiger_p_sidedness"] != "two-sided" or
                row["steiger_df"] != "not_applicable" or
                row["steiger_p_display"] != "P<1e-300 (double-precision underflow; see log10 P)" or
                not (math.isfinite(r2_exposure) and math.isfinite(r2_outcome) and
                     r2_exposure > r2_outcome > 0)):
            raise ValueError("metadata or aggregate mismatch")
        for suffix in ("008", "015", "020"):
            if (row[f"steiger_correct_{suffix}"] != "1" or
                    float(row[f"steiger_p_{suffix}"]) != 0 or
                    not math.isfinite(float(row[f"steiger_z_{suffix}"])) or
                    float(row[f"steiger_z_{suffix}"]) <= 0 or
                    not math.isfinite(float(row[f"steiger_log10_p_{suffix}"])) or
                    float(row[f"steiger_log10_p_{suffix}"]) >= -300):
                raise ValueError("directionality result mismatch")
    except (KeyError, TypeError, ValueError):
        errors.append(f"Steiger candidate-level assertion failed: {trait}")
for row in sample_qc:
    trait = row.get("trait_id", "<missing>")
    try:
        target, matched, missing = (int(row[name]) for name in
                                    ("target_iv_n", "matched_iv_n", "missing_iv_n"))
        if (row["trait_display"] != labels_by_trait.get(trait) or
                row["source_accession"] != accession_by_trait.get(trait) or
                target != matched or missing != 0 or
                target < int(steiger_by_trait[trait]["n_iv_full"]) or
                row["exposure_n_min"] != "413897" or
                row["exposure_n_max"] != "599249" or
                row["exposure_n_unique"] != "2" or
                row["exposure_n_values"] != "413897;599249" or
                int(row["source_rows_scanned"]) <= 0 or
                int(row["compressed_bytes_read"]) <= 0 or
                float(row["elapsed_seconds"]) <= 0):
            raise ValueError("sample-size recovery mismatch")
    except (KeyError, TypeError, ValueError):
        errors.append(f"exposure sample-size QC assertion failed: {trait}")
if (sum(int(row["target_iv_n"]) for row in sample_qc) != 5665 or
        sum(int(row["matched_iv_n"]) for row in sample_qc) != 5665 or
        sum(int(row["missing_iv_n"]) for row in sample_qc) != 0):
    errors.append("exposure sample-size recovery totals must be 5,665/5,665 with zero missing")

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
