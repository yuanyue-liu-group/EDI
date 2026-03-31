#!/usr/bin/env python3
"""
EDI Supercell Generator

Generates all input files needed for an EDI electron-defect interaction
calculation from a primitive cell scf.in and a defect specification.

Usage examples:

  # S vacancy in a 3x3x1 MoS2 supercell:
  python gen_supercell.py --input scf.in --nx 3 --ny 3 --nz 1 \
      --defect vacancy --site-species S --nprocs 16 --output ./edi_run/

  # Substitute S with O:
  python gen_supercell.py --input scf.in --nx 4 --ny 4 --nz 1 \
      --defect substitution --site-species S \
      --new-species O --new-mass 15.999 \
      --new-pseudo O_ONCV_PBE_FR-1.0.upf --output ./edi_run/
"""

import argparse
import os
import numpy as np
import copy

from qe_input import parse_qe_input, write_qe_input, QEInput
from supercell import build_supercell_centered, DefectSpec, generate_nscf_kpoints
from gen_pp_inputs import write_edi_input, write_run_script


def main():
    parser = argparse.ArgumentParser(
        description="Generate supercell inputs for EDI calculation"
    )
    parser.add_argument("--input", required=True,
                        help="Path to primitive cell scf.in")

    parser.add_argument("--nx", type=int, required=True,
                        help="Supercell size along a1")
    parser.add_argument("--ny", type=int, required=True,
                        help="Supercell size along a2")
    parser.add_argument("--nz", type=int, default=1,
                        help="Supercell size along a3 (default: 1 for 2D)")

    parser.add_argument("--defect", required=True,
                        choices=["vacancy", "substitution", "interstitial"],
                        help="Type of defect")
    parser.add_argument("--site", type=int, default=None,
                        help="Atom index (0-based) in primitive cell "
                             "for vacancy/substitution")
    parser.add_argument("--site-species", default=None,
                        help="Species label to select defect site "
                             "(e.g. 'S' for S vacancy, 'Mo' for Mo vacancy). "
                             "Selects the first atom with this label. "
                             "Takes precedence over --site if both given.")
    parser.add_argument("--new-species", default="",
                        help="New species label for substitution/interstitial "
                             "(e.g. 'O' to substitute S with O)")
    parser.add_argument("--new-mass", type=float, default=0.0,
                        help="Mass of new species")
    parser.add_argument("--new-pseudo", default="",
                        help="Pseudopotential file for new species")
    parser.add_argument("--interstitial-pos", nargs=3, type=float, default=None,
                        help="Fractional position for interstitial (in supercell coords)")

    parser.add_argument("--nprocs", type=int, default=4,
                        help="Number of MPI processes for the run script")
    parser.add_argument("--nscf-kgrid", nargs=3, type=int, default=None,
                        help="Explicit NSCF k-grid for primitive cell "
                             "(overrides automatic)")
    parser.add_argument("--vac-idx", type=int, default=0,
                        help="Vacuum plane index for alignment (2D systems)")
    parser.add_argument("--output", default="./edi_run/",
                        help="Output directory")
    args = parser.parse_args()

    if args.defect in ("vacancy", "substitution"):
        if args.site is None and args.site_species is None:
            parser.error(
                f"--defect {args.defect} requires --site or --site-species"
            )

    prim = parse_qe_input(args.input)
    print(f"Parsed primitive cell: {prim.nat} atoms, {prim.ntyp} types")
    print(f"  Species: {[sp.label for sp in prim.species]}")
    print(f"  Atoms:   {[a.label for a in prim.atoms]}")
    print(f"  Noncollinear: {prim.noncolin}, SOC: {prim.lspinorb}")

    defect_spec = DefectSpec(
        defect_type=args.defect,
        site_index=args.site,
        site_species=args.site_species,
        new_species=args.new_species,
        new_mass=args.new_mass,
        new_pseudo=args.new_pseudo,
        position=np.array(args.interstitial_pos) if args.interstitial_pos else None
    )

    pristine, defect, aligned_prim, shift = build_supercell_centered(
        prim, args.nx, args.ny, args.nz, defect_spec
    )

    pristine.prefix = prim.prefix + "_pristine"
    pristine.outdir = "./dout/"
    defect.prefix = prim.prefix + "_defect"
    defect.outdir = "./dout/"
    aligned_prim.prefix = prim.prefix
    aligned_prim.outdir = "./dout/"

    print(f"\nPristine supercell: {pristine.nat} atoms, "
          f"{args.nx}x{args.ny}x{args.nz}")
    print(f"Defect supercell:   {defect.nat} atoms ({args.defect})")

    if args.nscf_kgrid:
        nscf_kgrid = args.nscf_kgrid
    elif prim.kpoints_mode == 'automatic' and prim.kpoints_grid:
        nscf_kgrid = prim.kpoints_grid[:3]
    else:
        nscf_kgrid = [6, 6, 1]

    nscf_kpoints = generate_nscf_kpoints(nscf_kgrid)
    print(f"Primitive NSCF k-points: {len(nscf_kpoints)} points "
          f"(grid {nscf_kgrid[0]}x{nscf_kgrid[1]}x{nscf_kgrid[2]})")

    base = os.path.abspath(args.output)
    prim_dir = os.path.join(base, "primitive")
    prist_dir = os.path.join(base, "pristine_super")
    defect_dir = os.path.join(base, "defect_super")
    edi_dir = os.path.join(base, "edi")
    for d in [prim_dir, prist_dir, defect_dir, edi_dir]:
        os.makedirs(d, exist_ok=True)

    write_qe_input(aligned_prim, os.path.join(prim_dir, "scf.in"))

    prim_nscf = copy.deepcopy(aligned_prim)
    prim_nscf.calculation = "bands"
    prim_nscf.kpoints_mode = "crystal"
    prim_nscf.kpoints_list = nscf_kpoints
    write_qe_input(prim_nscf, os.path.join(prim_dir, "nscf.in"))

    write_qe_input(pristine, os.path.join(prist_dir, "scf.in"))

    write_qe_input(defect, os.path.join(defect_dir, "scf.in"))

    write_edi_input(
        filepath=os.path.join(edi_dir, "edi.in"),
        prim=prim,
        prim_prefix=aligned_prim.prefix,
        prim_outdir=os.path.join(prim_dir, aligned_prim.outdir),
        pristine_prefix=pristine.prefix,
        pristine_outdir=os.path.join(prist_dir, pristine.outdir),
        defect_prefix=defect.prefix,
        defect_outdir=os.path.join(defect_dir, defect.outdir),
        is_noncolin=prim.noncolin,
        is_lspinorb=prim.lspinorb,
        lvacalign=(args.nz == 1),
        vac_idx=args.vac_idx
    )

    write_run_script(
        filepath=os.path.join(base, "run_edi.sh"),
        nprocs=args.nprocs,
        prim_dir=prim_dir,
        pristine_dir=prist_dir,
        defect_dir=defect_dir,
        is_noncolin=prim.noncolin,
        edi_dir=edi_dir
    )

    print(f"\n{'=' * 60}")
    print(f"Generated files in: {base}")
    print(f"{'=' * 60}")
    print(f"  primitive/scf.in          - Aligned primitive cell SCF")
    print(f"  primitive/nscf.in         - Primitive cell NSCF (wavefunctions)")
    print(f"  pristine_super/scf.in     - Pristine supercell SCF (potential)")
    print(f"  defect_super/scf.in       - Defect supercell SCF (potential)")
    print(f"  edi/edi.in                - EDI input file (read via stdin)")
    print(f"  run_edi.sh                - Full pipeline script")
    print(f"\nTo run: cd {base} && bash run_edi.sh")


if __name__ == "__main__":
    main()
