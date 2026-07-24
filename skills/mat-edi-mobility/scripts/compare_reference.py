"""Check a computed EDI mobility table against a reference within a tolerance.

Loads two ``prefix_transport.dat`` files (a freshly computed one and a reference,
e.g. the upstream ``reference/mos2_transport.dat``), matches rows by temperature,
and reports the relative deviation of each mobility component. Cross-toolchain
reproductions are not expected to match to the last digit, so the pass criterion
is a user-set relative tolerance (default 10%).

Usage:
    python compare_reference.py mos2_transport.dat reference_transport.dat \
        --tol 0.10 --output-dir results/

Requirements:
    - Conda environment: base-agent
    - Required packages: pyyaml (standard library otherwise)
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parse_transport import COLUMNS, parse_transport

MOBILITY_COLUMNS = [c for c in COLUMNS if c != "T_K"]


def compare_tables(
    computed: list[dict[str, float]],
    reference: list[dict[str, float]],
    tol: float,
) -> tuple[list[dict[str, object]], bool]:
    """Compare computed vs reference mobility rows by temperature.

    Args:
        computed: Rows from the computed transport.dat.
        reference: Rows from the reference transport.dat.
        tol: Relative tolerance (fraction) for the pass criterion.

    Returns:
        A tuple ``(rows, all_pass)`` where ``rows`` lists per-temperature,
        per-component relative deviations and pass flags, and ``all_pass`` is
        True only if every compared component is within ``tol``.
    """
    ref_by_temp = {round(r["T_K"], 3): r for r in reference}
    rows: list[dict[str, object]] = []
    all_pass = True
    for comp in computed:
        temp = round(comp["T_K"], 3)
        ref = ref_by_temp.get(temp)
        if ref is None:
            continue
        for col in MOBILITY_COLUMNS:
            ref_val = ref[col]
            rel = (
                abs(comp[col] - ref_val) / abs(ref_val)
                if ref_val != 0
                else float("inf")
            )
            passed = rel <= tol
            all_pass = all_pass and passed
            rows.append(
                {
                    "T_K": temp,
                    "component": col,
                    "computed": comp[col],
                    "reference": ref_val,
                    "rel_dev": rel,
                    "pass": passed,
                }
            )
    if not rows:
        raise ValueError("No matching temperatures between the two tables")
    return rows, all_pass


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare a computed EDI mobility table against a reference."
    )
    parser.add_argument("computed_dat", help="Computed prefix_transport.dat")
    parser.add_argument("reference_dat", help="Reference prefix_transport.dat")
    parser.add_argument(
        "--tol",
        type=float,
        default=0.10,
        help="Relative tolerance for the pass criterion (default 0.10 = 10%%)",
    )
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Directory for the JSON summary and input_configs.yaml",
    )
    args = parser.parse_args()

    computed, _ = parse_transport(args.computed_dat)
    reference, _ = parse_transport(args.reference_dat)
    rows, all_pass = compare_tables(computed, reference, args.tol)

    os.makedirs(args.output_dir, exist_ok=True)
    print(
        f"{'T(K)':>7} {'component':>14} {'computed':>13} {'reference':>13} {'rel_dev':>9} pass"
    )
    for row in rows:
        print(
            f"{row['T_K']:7.1f} {row['component']:>14} {row['computed']:13.4f} "
            f"{row['reference']:13.4f} {row['rel_dev']:9.4f} {row['pass']}"
        )
    print(f"\nOverall: {'PASS' if all_pass else 'FAIL'} at tol={args.tol}")

    with open(
        os.path.join(args.output_dir, "reference_comparison.json"), "w"
    ) as handle:
        json.dump(
            {"tol": args.tol, "all_pass": all_pass, "rows": rows}, handle, indent=2
        )

    configs = {
        "computed_file": os.path.abspath(args.computed_dat),
        "reference_file": os.path.abspath(args.reference_dat),
        "output_dir": os.path.abspath(args.output_dir),
        "tol": args.tol,
    }
    with open(os.path.join(args.output_dir, "input_configs.yaml"), "w") as handle:
        yaml.safe_dump(configs, handle, default_flow_style=False)


if __name__ == "__main__":
    main()
