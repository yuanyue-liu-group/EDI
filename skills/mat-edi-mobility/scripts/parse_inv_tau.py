"""Parse an EDI ``prefix_inv_tau.dat`` file of state-resolved scattering rates.

EDI's transport pass writes ``prefix_inv_tau.dat`` with one header comment line
and one row per (k-point, band) state that has a nonzero rate, as emitted by
``transport.f90`` (header at line 634, data at lines 637-643):

    ik  ibnd  E(eV)  inv_tau_SERTA(Ry)  inv_tau_MRTA(Ry)  tau_SERTA(fs)  tau_MRTA(fs)

This script loads those rows, reports simple statistics (state count, energy
range, mean SERTA/MRTA lifetimes), and writes a JSON summary plus an
``input_configs.yaml``.

Usage:
    python parse_inv_tau.py mos2_inv_tau.dat --output-dir results/

Requirements:
    - Conda environment: base-agent
    - Required packages: pyyaml (standard library otherwise)
"""

from __future__ import annotations

import argparse
import json
import os

import yaml

COLUMNS = [
    "ik",
    "ibnd",
    "E_eV",
    "inv_tau_SERTA_Ry",
    "inv_tau_MRTA_Ry",
    "tau_SERTA_fs",
    "tau_MRTA_fs",
]


def parse_inv_tau(path: str) -> list[dict[str, float]]:
    """Parse an EDI inv_tau.dat file.

    Args:
        path: Path to a ``prefix_inv_tau.dat`` file.

    Returns:
        A list of dicts keyed by COLUMNS. ``ik`` and ``ibnd`` are ints; the
        remaining fields are floats.
    """
    rows: list[dict[str, float]] = []
    with open(path) as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split()
            if len(fields) < len(COLUMNS):
                continue
            row: dict[str, float] = {}
            for i, col in enumerate(COLUMNS):
                row[col] = int(fields[i]) if col in ("ik", "ibnd") else float(fields[i])
            rows.append(row)
    if not rows:
        raise ValueError(f"No data rows found in {path}")
    return rows


def summarize(rows: list[dict[str, float]]) -> dict[str, float]:
    """Compute summary statistics over parsed inverse-lifetime rows.

    Args:
        rows: Output of :func:`parse_inv_tau`.

    Returns:
        Dictionary with keys: 'n_states', 'E_min_eV', 'E_max_eV',
        'mean_tau_SERTA_fs', 'mean_tau_MRTA_fs'.
    """
    n = len(rows)
    energies = [r["E_eV"] for r in rows]
    return {
        "n_states": n,
        "E_min_eV": min(energies),
        "E_max_eV": max(energies),
        "mean_tau_SERTA_fs": sum(r["tau_SERTA_fs"] for r in rows) / n,
        "mean_tau_MRTA_fs": sum(r["tau_MRTA_fs"] for r in rows) / n,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Parse an EDI inv_tau.dat file of state-resolved lifetimes."
    )
    parser.add_argument("inv_tau_dat", help="Path to prefix_inv_tau.dat")
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Directory for the JSON summary and input_configs.yaml",
    )
    args = parser.parse_args()

    rows = parse_inv_tau(args.inv_tau_dat)
    stats = summarize(rows)
    os.makedirs(args.output_dir, exist_ok=True)

    for key, value in stats.items():
        print(f"{key:>20}: {value}")

    with open(os.path.join(args.output_dir, "inv_tau_summary.json"), "w") as handle:
        json.dump({"statistics": stats, "states": rows}, handle, indent=2)

    configs = {
        "input_file": os.path.abspath(args.inv_tau_dat),
        "output_dir": os.path.abspath(args.output_dir),
        "columns": COLUMNS,
    }
    with open(os.path.join(args.output_dir, "input_configs.yaml"), "w") as handle:
        yaml.safe_dump(configs, handle, default_flow_style=False)


if __name__ == "__main__":
    main()
