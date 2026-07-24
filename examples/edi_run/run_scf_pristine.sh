#!/bin/bash
#SBATCH -A your_allocation
#SBATCH -p shared
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -t 04:00:00
#SBATCH -J scf_prist
#SBATCH -o scf_prist.%j.out
#SBATCH -e scf_prist.%j.err

# Site-specific environment (Anvil example); adjust to your system
module reset
module load aocc/3.1.0 openmpi/4.1.6
module load amdblis/3.0 amdlibflame/3.0 amdlibm/3.0 fftw

QEDIR=${QEDIR:-/path/to/qe-7.5}   # QE root containing edi-code/
PW=$QEDIR/PW/src/pw.x
NPROC=128

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/pristine_super"
mkdir -p dout
srun -n $NPROC $PW < scf.in > scf.out 2>&1
grep "convergence has been achieved" scf.out && echo "Pristine SCF OK" || echo "Pristine SCF FAILED"
