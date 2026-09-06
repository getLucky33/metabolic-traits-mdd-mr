#!/usr/bin/env python3
"""Fail on credentials, machine paths, runtime installers or restricted file types."""
from __future__ import annotations

import re
import subprocess
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

errors: list[str] = []
for rel in listed:
    path = ROOT / rel
    if not path.is_file():
        continue
    if path.suffix.lower() in restricted_suffixes:
        errors.append(f"restricted binary/source-data suffix: {rel}")
    if path.suffix not in text_suffixes and path.name not in {"DESCRIPTION", "LICENSE", ".gitignore", ".gitattributes"}:
        continue
    try:
        text = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        continue
    for label, pattern in patterns.items():
        if pattern.search(text):
            errors.append(f"{label}: {rel}")
    if rel.startswith("analysis/scripts/"):
        for marker in forbidden_installers:
            if marker in text:
                errors.append(f"runtime installer {marker}: {rel}")

if errors:
    raise SystemExit("Repository audit failed:\n- " + "\n- ".join(sorted(set(errors))))
print(f"PASS: audited {len(listed)} repository files")
