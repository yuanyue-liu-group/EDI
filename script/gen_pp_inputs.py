import os
from qe_input import QEInput


def write_edi_input(filepath: str, prim: QEInput,
                      prim_prefix: str, prim_outdir: str,
                      pristine_prefix: str, pristine_outdir: str,
                      defect_prefix: str, defect_outdir: str,
                      wt_filename: str = "scfwt.dat",
                      is_noncolin: bool = False,
                      is_lspinorb: bool = False,
                      lvacalign: bool = True,
                      vac_idx: int = 0):
    with open(filepath, 'w') as f:
        f.write("&edicontrol\n")

        f.write(f"  prefix='{prim_prefix}',\n")
        f.write(f"  outdir='{prim_outdir}'\n")

        f.write(f"  pristine_prefix='{pristine_prefix}',\n")
        f.write(f"  pristine_outdir='{pristine_outdir}'\n")
        f.write(f"  defect_prefix='{defect_prefix}',\n")
        f.write(f"  defect_outdir='{defect_outdir}'\n")

        if lvacalign:
            f.write("  lvacalign=.true.\n")
            f.write(f"  vac_idx={vac_idx}\n")
            f.write("  lcorealign=.false.\n")
        else:
            f.write("  lvacalign=.false.\n")
            f.write("  lcorealign=.true.\n")
            f.write("  core_v_d=0.0\n")
            f.write("  core_v_p=0.0\n")

        f.write(f"  wt_filename='{wt_filename}'\n")

        f.write(f"  noncolin ={'.true.' if is_noncolin else '.false.'}\n")
        f.write(f"  lspinorb ={'.true.' if is_lspinorb else '.false.'}\n")

        f.write("  calcmlocal = .true.\n")
        f.write("  calcmnonlocal = .true.\n")

        f.write("/\n")


def write_run_script(filepath: str, nprocs: int,
                     prim_dir: str, pristine_dir: str, defect_dir: str,
                     is_noncolin: bool, edi_dir: str):
    with open(filepath, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write("set -e\n\n")
        f.write(f"NPROCS={nprocs}\n\n")

        f.write("#" + "=" * 60 + "\n")
        f.write("# Step 1: Supercell SCF calculations (potentials only)\n")
        f.write("#" + "=" * 60 + "\n\n")

        f.write(f"echo '>>> Running pristine supercell SCF...'\n")
        f.write(f"cd {pristine_dir}\n")
        f.write(f"mpirun -np $NPROCS pw.x < scf.in > scf.out\n\n")

        f.write(f"echo '>>> Running defect supercell SCF...'\n")
        f.write(f"cd {defect_dir}\n")
        f.write(f"mpirun -np $NPROCS pw.x < scf.in > scf.out\n\n")

        f.write("#" + "=" * 60 + "\n")
        f.write("# Step 2: Primitive cell SCF + NSCF (wavefunctions)\n")
        f.write("#" + "=" * 60 + "\n\n")

        f.write(f"echo '>>> Running primitive cell SCF...'\n")
        f.write(f"cd {prim_dir}\n")
        f.write(f"mpirun -np $NPROCS pw.x < scf.in > scf.out\n\n")

        f.write(f"echo '>>> Running primitive cell NSCF...'\n")
        f.write(f"cd {prim_dir}\n")
        f.write(f"mpirun -np $NPROCS pw.x < nscf.in > nscf.out\n\n")

        f.write("#" + "=" * 60 + "\n")
        f.write("# Step 3: Run EDI (potentials extracted on-the-fly)\n")
        f.write("#" + "=" * 60 + "\n\n")

        f.write(f"echo '>>> Running EDI...'\n")
        f.write(f"cd {edi_dir}\n")
        f.write(f"mpirun -np $NPROCS edi.x < edi.in > edic.out\n\n")

        f.write("echo '>>> Done!'\n")

    os.chmod(filepath, 0o755)
