"""Parse EDI electron-defect matrix-element output files ``|M(k_i, k_f)|^2``.

EDI writes several matrix-element text files that all share a fixed-width layout
with ``|M|^2  Re(M)  Im(M)`` as the last physical columns. This parser
auto-detects which of the three EDI layouts a file uses, from the number of
whitespace-separated fields in the first data row:

* 9 fields  -> single-k list (``prefix_edmatf.dat`` diagonal interpolation, or
  ``prefix_edmat_scatter.dat`` from-K interpolation):
  ``ik  kx ky kz  ibnd jbnd  |M|^2  Re(M)  Im(M)``
  (``ed_coarse.f90`` headers at lines 71 and 198).
* 13 fields -> Wannier-interpolated pair list (``prefix_edmat_interp.dat``,
  produced by ``edmat_interp_from_file``):
  ``iki ikf  kix kiy kiz  kfx kfy kfz  ibnd jbnd  |M|^2  Re(M)  Im(M)``
  (``ed_coarse.f90`` header at line 2526, format at line 2551).
* 15 fields -> direct / coarse Bloch pair list (``prefix_edmat_direct.dat`` from
  ``edmat_direct_from_file``, or ``prefix_edmat_bloch.dat``): the 13-field
  layout plus ``|M_loc|^2  |M_nl|^2``
  (``ed_coarse.f90`` headers at lines 2145 and 739, format at line 2405).

Each parsed record is normalised to keys: ``ikf`` (final-k index), ``kf`` (3
floats), ``ibnd``, ``jbnd``, ``m2`` (``|M|^2`` in Ry^2), ``reM``, ``imM``, and,
when present, ``m2_loc`` and ``m2_nl``.

Usage:
    python parse_edmat.py mos2_edmat_interp.dat --ibnd 2 --jbnd 2 --output-dir results/

Requirements:
    - Conda environment: base-agent
    - Required packages: pyyaml (standard library otherwise)
"""

from __future__ import annotations

import argparse
import json
import os

import yaml


def parse_edmat(path: str) -> tuple[list[dict[str, float]], int]:
    """Parse an EDI matrix-element file, auto-detecting its layout.

    Args:
        path: Path to an EDI ``*_edmat*.dat`` file.

    Returns:
        A tuple ``(records, ncol)`` where ``records`` is a list of normalised
        dicts and ``ncol`` is the detected field count (9, 13, or 15).

    Raises:
        ValueError: If the file has no data rows or an unrecognised layout.
    """
    records: list[dict[str, float]] = []
    ncol = 0
    with open(path) as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split()
            if ncol == 0:
                ncol = len(fields)
                if ncol not in (9, 13, 15):
                    raise ValueError(
                        f"Unrecognised EDI matrix-element layout: {ncol} columns in {path}"
                    )
            records.append(_parse_row(fields, ncol))
    if not records:
        raise ValueError(f"No data rows found in {path}")
    return records, ncol


def _parse_row(fields: list[str], ncol: int) -> dict[str, float]:
    """Normalise one data row to a common record schema.

    Args:
        fields: Whitespace-split tokens of a single data line.
        ncol: The layout width detected for the file (9, 13, or 15).

    Returns:
        Dict with keys ikf, kf, ibnd, jbnd, m2, reM, imM (+ m2_loc, m2_nl
        for the 15-column direct layout).
    """
    if ncol == 9:
        # ik  kx ky kz  ibnd jbnd  |M|^2  Re  Im
        return {
            "ikf": int(fields[0]),
            "kf": [float(fields[1]), float(fields[2]), float(fields[3])],
            "ibnd": int(fields[4]),
            "jbnd": int(fields[5]),
            "m2": float(fields[6]),
            "reM": float(fields[7]),
            "imM": float(fields[8]),
        }
    # 13- and 15-column layouts share the first 13 fields.
    record = {
        "ikf": int(fields[1]),
        "kf": [float(fields[5]), float(fields[6]), float(fields[7])],
        "ibnd": int(fields[8]),
        "jbnd": int(fields[9]),
        "m2": float(fields[10]),
        "reM": float(fields[11]),
        "imM": float(fields[12]),
    }
    if ncol == 15:
        record["m2_loc"] = float(fields[13])
        record["m2_nl"] = float(fields[14])
    return record


def select_band_pair(
    records: list[dict[str, float]], ibnd: int, jbnd: int
) -> list[dict[str, float]]:
    """Filter records to a single (ibnd, jbnd) band pair, ordered by ikf.

    Args:
        records: Output of :func:`parse_edmat`.
        ibnd: Initial-state band index to keep.
        jbnd: Final-state band index to keep.

    Returns:
        The matching records sorted by final-k index ``ikf``.
    """
    picked = [r for r in records if r["ibnd"] == ibnd and r["jbnd"] == jbnd]
    return sorted(picked, key=lambda r: r["ikf"])


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Parse an EDI electron-defect matrix-element file."
    )
    parser.add_argument("edmat_dat", help="Path to an EDI *_edmat*.dat file")
    parser.add_argument("--ibnd", type=int, default=None, help="Keep only this ibnd")
    parser.add_argument("--jbnd", type=int, default=None, help="Keep only this jbnd")
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Directory for the JSON summary and input_configs.yaml",
    )
    args = parser.parse_args()

    records, ncol = parse_edmat(args.edmat_dat)
    if args.ibnd is not None and args.jbnd is not None:
        records = select_band_pair(records, args.ibnd, args.jbnd)

    os.makedirs(args.output_dir, exist_ok=True)
    print(f"Detected {ncol}-column layout; {len(records)} records retained.")
    if records:
        m2_values = [r["m2"] for r in records]
        print(f"max |M|^2 = {max(m2_values):.6e} Ry^2")

    with open(os.path.join(args.output_dir, "edmat_records.json"), "w") as handle:
        json.dump({"ncol": ncol, "records": records}, handle, indent=2)

    configs = {
        "input_file": os.path.abspath(args.edmat_dat),
        "output_dir": os.path.abspath(args.output_dir),
        "ibnd": args.ibnd,
        "jbnd": args.jbnd,
    }
    with open(os.path.join(args.output_dir, "input_configs.yaml"), "w") as handle:
        yaml.safe_dump(configs, handle, default_flow_style=False)


if __name__ == "__main__":
    main()
