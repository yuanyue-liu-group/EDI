"""Compare directly computed vs Wannier-interpolated EDI matrix elements.

This reproduces EDI's Wannier-interpolation validation (README figure: direct
DFT dots vs interpolated lines along Gamma-K-M-Gamma). It reads a *direct*
matrix-element file (``prefix_edmat_direct.dat`` from ``edmat_direct_from_file``)
and an *interpolated* file (``prefix_edmat_interp.dat`` from
``edmat_interp_from_file``) computed on the SAME ``ki.dat`` / ``kf.dat`` k-path,
matches records by ``(ikf, ibnd, jbnd)``, and reports agreement statistics on
``|M|^2``: RMS deviation, max absolute deviation, and relative RMS.

With ``--band-sum`` it instead compares the gauge-invariant band-summed
``sum_{mn} |M_mn|^2`` per k-point over the shared band subspace. This is the
physically meaningful comparison near band crossings, where the arbitrary
Wannier gauge mixes individual band pairs (pairwise deviations are inflated
there even when the interpolation is correct).

Usage:
    python compare_edmat.py mos2_edmat_direct.dat mos2_edmat_interp.dat \
        --band-offset 12 --band-sum --output-dir results/

Requirements:
    - Conda environment: base-agent
    - Required packages: pyyaml (standard library otherwise)
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parse_edmat import parse_edmat


def _key(record: dict) -> tuple[int, int, int]:
    """Return the match key ``(ikf, ibnd, jbnd)`` for a record."""
    return (record["ikf"], record["ibnd"], record["jbnd"])


def compare(
    direct: list[dict],
    interp: list[dict],
    ibnd: int | None,
    jbnd: int | None,
    band_offset: int = 0,
) -> dict[str, float]:
    """Compute agreement statistics between direct and interpolated |M|^2.

    Args:
        direct: Records from the direct matrix-element file.
        interp: Records from the interpolated matrix-element file.
        ibnd: If set, restrict the comparison to this initial band
            (direct-file numbering).
        jbnd: If set, restrict the comparison to this final band
            (direct-file numbering).
        band_offset: Direct-file band index minus interpolated-file band
            index. The direct mode numbers bands absolutely (all NSCF
            bands) while the interpolated file numbers only the Wannier
            subspace, so with ``bands_skipped = 'exclude_bands = 1-12'``
            the offset is 12.

    Returns:
        Dictionary with keys 'n_matched', 'rms_abs', 'max_abs',
        'rms_rel', and 'mean_m2_direct' (all over |M|^2 in Ry^2).
    """
    interp_map = {_key(r): r for r in interp}
    diffs: list[float] = []
    direct_vals: list[float] = []
    for rec in direct:
        if ibnd is not None and rec["ibnd"] != ibnd:
            continue
        if jbnd is not None and rec["jbnd"] != jbnd:
            continue
        match = interp_map.get(
            (rec["ikf"], rec["ibnd"] - band_offset, rec["jbnd"] - band_offset)
        )
        if match is None:
            continue
        diffs.append(rec["m2"] - match["m2"])
        direct_vals.append(rec["m2"])

    n = len(diffs)
    if n == 0:
        raise ValueError("No overlapping (ikf, ibnd, jbnd) records to compare")
    rms_abs = math.sqrt(sum(d * d for d in diffs) / n)
    mean_direct = sum(direct_vals) / n
    denom = math.sqrt(sum(v * v for v in direct_vals) / n)
    return {
        "n_matched": n,
        "rms_abs": rms_abs,
        "max_abs": max(abs(d) for d in diffs),
        "rms_rel": rms_abs / denom if denom > 0 else float("nan"),
        "mean_m2_direct": mean_direct,
    }


def compare_band_sum(
    direct: list[dict],
    interp: list[dict],
    band_offset: int = 0,
) -> dict[str, float]:
    """Compare gauge-invariant band-summed |M|^2 per k-point.

    For each final k-point, sum ``|M|^2`` over every band pair present in the
    interpolated file (the Wannier subspace) and over the corresponding
    band-offset-shifted pairs of the direct file, then compare the two sums
    pointwise. The band sum is invariant under the arbitrary Wannier gauge,
    so it stays meaningful across band crossings where pairwise elements mix.

    Args:
        direct: Records from the direct matrix-element file.
        interp: Records from the interpolated matrix-element file.
        band_offset: Direct-file band index minus interpolated-file band
            index (see :func:`compare`).

    Returns:
        Dictionary with keys 'n_points', 'rms_rel' (RMS of the band-sum
        deviations normalized by the RMS of the direct band sums, the same
        convention as :func:`compare`), 'max_rel' and 'mean_rel' (pointwise
        relative deviations (direct - interp) / direct), and
        'mean_sum_direct' (mean band-summed |M|^2 in Ry^2).
    """
    subspace = {(r["ibnd"], r["jbnd"]) for r in interp}
    interp_sum: dict[int, float] = {}
    for rec in interp:
        interp_sum[rec["ikf"]] = interp_sum.get(rec["ikf"], 0.0) + rec["m2"]
    direct_sum: dict[int, float] = {}
    for rec in direct:
        pair = (rec["ibnd"] - band_offset, rec["jbnd"] - band_offset)
        if pair not in subspace:
            continue
        direct_sum[rec["ikf"]] = direct_sum.get(rec["ikf"], 0.0) + rec["m2"]

    common = sorted(set(direct_sum) & set(interp_sum))
    if not common:
        raise ValueError("No common k-points between direct and interp files")
    diffs = [direct_sum[ik] - interp_sum[ik] for ik in common]
    rel = [(direct_sum[ik] - interp_sum[ik]) / direct_sum[ik] for ik in common]
    n = len(rel)
    rms_abs = math.sqrt(sum(d * d for d in diffs) / n)
    denom = math.sqrt(sum(direct_sum[ik] ** 2 for ik in common) / n)
    return {
        "n_points": n,
        "rms_rel": rms_abs / denom if denom > 0 else float("nan"),
        "max_rel": max(abs(r) for r in rel),
        "mean_rel": sum(abs(r) for r in rel) / n,
        "mean_sum_direct": sum(direct_sum[ik] for ik in common) / n,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare direct vs Wannier-interpolated EDI matrix elements."
    )
    parser.add_argument("direct_dat", help="Direct M file (prefix_edmat_direct.dat)")
    parser.add_argument(
        "interp_dat", help="Interpolated M file (prefix_edmat_interp.dat)"
    )
    parser.add_argument("--ibnd", type=int, default=None, help="Restrict to this ibnd")
    parser.add_argument("--jbnd", type=int, default=None, help="Restrict to this jbnd")
    parser.add_argument(
        "--band-offset",
        type=int,
        default=0,
        help=(
            "Direct-file band index minus interp-file band index "
            "(= number of excluded bands, e.g. 12 for "
            "bands_skipped='exclude_bands = 1-12')"
        ),
    )
    parser.add_argument(
        "--band-sum",
        action="store_true",
        help=(
            "Compare gauge-invariant band-summed |M|^2 per k-point over the "
            "Wannier subspace instead of pairwise elements"
        ),
    )
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Directory for the JSON summary and input_configs.yaml",
    )
    args = parser.parse_args()
    if args.band_sum and (args.ibnd is not None or args.jbnd is not None):
        parser.error("--band-sum sums over all band pairs; drop --ibnd/--jbnd")

    direct, _ = parse_edmat(args.direct_dat)
    interp, _ = parse_edmat(args.interp_dat)
    if args.band_sum:
        stats = compare_band_sum(direct, interp, args.band_offset)
    else:
        stats = compare(direct, interp, args.ibnd, args.jbnd, args.band_offset)

    os.makedirs(args.output_dir, exist_ok=True)
    for key, value in stats.items():
        print(f"{key:>16}: {value}")

    with open(os.path.join(args.output_dir, "edmat_comparison.json"), "w") as handle:
        json.dump(stats, handle, indent=2)

    configs = {
        "direct_file": os.path.abspath(args.direct_dat),
        "interp_file": os.path.abspath(args.interp_dat),
        "output_dir": os.path.abspath(args.output_dir),
        "ibnd": args.ibnd,
        "jbnd": args.jbnd,
        "band_offset": args.band_offset,
        "band_sum": args.band_sum,
    }
    with open(os.path.join(args.output_dir, "input_configs.yaml"), "w") as handle:
        yaml.safe_dump(configs, handle, default_flow_style=False)


if __name__ == "__main__":
    main()
