"""Parse an EDI ``prefix_transport.dat`` file into a mobility-vs-temperature table.

EDI's transport pass (``do_transport = .true.``) writes ``prefix_transport.dat``
with three comment lines (title, grid, energy window) followed by one row per
temperature. Each data row is:

    T(K)  mu_SERTA_xx  mu_MRTA_xx  mu_SERTA_yy  mu_MRTA_yy   (cm^2/Vs)

as emitted by ``transport.f90`` (header at line 331, data at lines 337-650).
This script extracts those rows, prints a table, and writes a JSON summary plus
an ``input_configs.yaml`` capturing the run parameters.

Usage:
    python parse_transport.py mos2_transport.dat --output-dir results/

Requirements:
    - Conda environment: base-agent
    - Required packages: pyyaml (standard library otherwise)
"""

from __future__ import annotations

import argparse
import json
import os

import yaml

COLUMNS = ["T_K", "mu_SERTA_xx", "mu_MRTA_xx", "mu_SERTA_yy", "mu_MRTA_yy"]


def parse_transport(path: str) -> tuple[list[dict[str, float]], dict[str, str]]:
    """Parse an EDI transport.dat file.

    Args:
        path: Path to a ``prefix_transport.dat`` file.

    Returns:
        A tuple ``(rows, meta)`` where ``rows`` is a list of dicts keyed by
        COLUMNS (all cm^2/Vs except T_K in kelvin) and ``meta`` holds the
        comment-header strings (grid, window) for provenance.
    """
    rows: list[dict[str, float]] = []
    meta: dict[str, str] = {}
    with open(path) as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                if "Grid" in stripped:
                    meta["grid"] = stripped.lstrip("# ").strip()
                elif "Window" in stripped:
                    meta["window"] = stripped.lstrip("# ").strip()
                continue
            fields = stripped.split()
            if len(fields) < len(COLUMNS):
                continue
            rows.append({col: float(fields[i]) for i, col in enumerate(COLUMNS)})
    if not rows:
        raise ValueError(f"No data rows found in {path}")
    return rows, meta


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Parse an EDI transport.dat file into a mobility table."
    )
    parser.add_argument("transport_dat", help="Path to prefix_transport.dat")
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Directory for the JSON summary and input_configs.yaml",
    )
    args = parser.parse_args()

    rows, meta = parse_transport(args.transport_dat)
    os.makedirs(args.output_dir, exist_ok=True)

    header = "  ".join(f"{c:>13}" for c in COLUMNS)
    print(header)
    for row in rows:
        print("  ".join(f"{row[c]:13.4f}" for c in COLUMNS))

    summary = {"metadata": meta, "mobility_cm2_Vs": rows}
    with open(os.path.join(args.output_dir, "transport_summary.json"), "w") as handle:
        json.dump(summary, handle, indent=2)

    configs = {
        "input_file": os.path.abspath(args.transport_dat),
        "output_dir": os.path.abspath(args.output_dir),
        "columns": COLUMNS,
    }
    with open(os.path.join(args.output_dir, "input_configs.yaml"), "w") as handle:
        yaml.safe_dump(configs, handle, default_flow_style=False)


if __name__ == "__main__":
    main()
