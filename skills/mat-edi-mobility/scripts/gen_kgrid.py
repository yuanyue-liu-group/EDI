"""Generate an explicit uniform k-point list for a Quantum ESPRESSO NSCF deck.

EDI folds the primitive-cell Bloch states on an *explicit* uniform coarse
k-grid into Wannier functions. The grid dimensions written here MUST equal the
`coarse_nk1/nk2/nk3` values in the EDI input, or the Wannier-gauge folding is
inconsistent and the interpolation is meaningless. QE's `K_POINTS automatic`
compression cannot be used because EDI indexes the k-points directly.

The emitted block is a `K_POINTS crystal` card: a count line followed by one
`kx ky kz weight` row per k-point, in the row-major order
(i1 fastest, then i2, then i3) that EDI assumes.

Usage:
    python gen_kgrid.py 12 12 1
    python gen_kgrid.py 12 12 1 --output kpoints_block.txt

Requirements:
    - Conda environment: base-agent
    - Required packages: Python 3 standard library only
"""

from __future__ import annotations

import argparse
import sys


def build_kgrid(nk1: int, nk2: int, nk3: int) -> list[tuple[float, float, float]]:
    """Build a Gamma-centred uniform k-grid in crystal coordinates.

    Args:
        nk1: Number of subdivisions along reciprocal vector b1.
        nk2: Number of subdivisions along reciprocal vector b2.
        nk3: Number of subdivisions along reciprocal vector b3.

    Returns:
        List of (kx, ky, kz) fractional coordinates in the interval [0, 1),
        ordered with the first index fastest (EDI's expected ordering).
    """
    points: list[tuple[float, float, float]] = []
    for i3 in range(nk3):
        for i2 in range(nk2):
            for i1 in range(nk1):
                points.append((i1 / nk1, i2 / nk2, i3 / nk3))
    return points


def format_card(points: list[tuple[float, float, float]]) -> str:
    """Format the k-points as a QE `K_POINTS crystal` card.

    Args:
        points: List of (kx, ky, kz) fractional coordinates.

    Returns:
        The full card as a string (header + count + one row per k-point).
    """
    weight = 1.0 / len(points)
    lines = ["K_POINTS crystal", str(len(points))]
    for kx, ky, kz in points:
        lines.append(f"  {kx:.12f}  {ky:.12f}  {kz:.12f}  {weight:.12e}")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Generate an explicit uniform K_POINTS crystal card for an EDI "
            "primitive-cell NSCF calculation."
        )
    )
    parser.add_argument("nk1", type=int, help="Grid subdivisions along b1")
    parser.add_argument("nk2", type=int, help="Grid subdivisions along b2")
    parser.add_argument("nk3", type=int, help="Grid subdivisions along b3")
    parser.add_argument(
        "--output",
        default=None,
        help="Write the card to this file instead of stdout",
    )
    args = parser.parse_args()

    if min(args.nk1, args.nk2, args.nk3) < 1:
        raise ValueError("Grid dimensions must be positive integers")

    card = format_card(build_kgrid(args.nk1, args.nk2, args.nk3))

    if args.output:
        with open(args.output, "w") as handle:
            handle.write(card)
    else:
        sys.stdout.write(card)


if __name__ == "__main__":
    main()
