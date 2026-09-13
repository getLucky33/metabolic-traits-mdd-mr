#!/usr/bin/env python3
"""Summarise candidate-level leave-one-out results without releasing SNP rows."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--leave-one-out", type=Path, required=True)
    parser.add_argument("--forward-screen", type=Path, default=ROOT / "data/derived/forward_screen_249.tsv")
    parser.add_argument("--evidence-summary", type=Path, default=ROOT / "data/derived/evidence_summary.tsv")
    parser.add_argument("--trait-labels", type=Path, default=ROOT / "data/derived/trait_labels.tsv")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    candidates = read_tsv(args.evidence_summary)
    order = [row["trait"] for row in candidates]
    if len(order) != 15 or len(set(order)) != 15:
        raise SystemExit("evidence-summary must contain 15 unique candidates")

    forward = {row["trait"]: row for row in read_tsv(args.forward_screen)}
    labels = {row["trait_id"]: row["trait_display"] for row in read_tsv(args.trait_labels)}
    if any(trait not in forward or trait not in labels for trait in order):
        raise SystemExit("candidate missing from forward screen or trait labels")

    rows_by_trait: dict[str, list[dict[str, str]]] = {trait: [] for trait in order}
    all_rows: dict[str, dict[str, str]] = {}
    with args.leave_one_out.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            trait = row.get("trait", "")
            if trait not in rows_by_trait:
                continue
            if row.get("SNP") == "All":
                if trait in all_rows:
                    raise SystemExit(f"duplicate All row: {trait}")
                all_rows[trait] = row
            else:
                rows_by_trait[trait].append(row)

    output: list[dict[str, object]] = []
    for trait in order:
        main_row = forward[trait]
        main_beta = float(main_row["ivw_b"])
        main_se = float(main_row["ivw_se"])
        main_p = float(main_row["ivw_p"])
        all_row = all_rows.get(trait)
        rows = rows_by_trait[trait]
        if all_row is None or len(rows) != int(main_row["n_iv"]):
            raise SystemExit(f"leave-one-out row count differs from n_iv: {trait}")
        if any(
            not math.isclose(float(all_row[field]), expected, rel_tol=1e-12, abs_tol=1e-15)
            for field, expected in (("b", main_beta), ("se", main_se), ("p", main_p))
        ):
            raise SystemExit(f"All row differs from the released forward estimate: {trait}")
        betas = [float(row["b"]) for row in rows]
        p_values = [float(row["p"]) for row in rows]
        if not all(math.isfinite(value) for value in betas + p_values):
            raise SystemExit(f"non-finite leave-one-out value: {trait}")
        direction_changes = sum(beta * main_beta <= 0 for beta in betas)
        nominal_losses = sum(p_value >= 0.05 for p_value in p_values)
        output.append(
            {
                "trait_id": trait,
                "trait_display": labels[trait],
                "n_leave_one_out_estimates": len(rows),
                "min_beta": format(min(betas), ".15g"),
                "max_beta": format(max(betas), ".15g"),
                "max_p": format(max(p_values), ".15g"),
                "direction_change_n": direction_changes,
                "nominal_loss_n": nominal_losses,
                "all_primary_direction": str(direction_changes == 0).upper(),
                "all_nominal_p_lt_0_05": str(nominal_losses == 0).upper(),
            }
        )

    if sum(int(row["n_leave_one_out_estimates"]) for row in output) != 4993:
        raise SystemExit("candidate leave-one-out total differs from the frozen 4,993 estimates")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)


if __name__ == "__main__":
    main()
