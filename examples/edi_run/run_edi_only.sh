#!/bin/bash
#SBATCH -A your_allocation
#SBATCH -p shared
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -t 12:00:00
#SBATCH -J edi_only
#SBATCH -o edi_only.%j.out
#SBATCH -e edi_only.%j.err

# Site-specific environment (Anvil example); adjust to your system
module reset
module load aocc/3.1.0 openmpi/4.1.6
module load amdblis/3.0 amdlibflame/3.0 amdlibm/3.0 fftw

QEDIR=${QEDIR:-/path/to/qe-7.5}   # QE root containing edi-code/
PW=$QEDIR/PW/src/pw.x
EDI=$QEDIR/edi-code/src/edi.x
NPROC=108

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"

# Step 1: Primitive cell SCF + NSCF (skip if already done)
#if [ ! -f primitive/dout/mos2.save/data-file-schema.xml ]; then
#    echo '>>> Running primitive cell SCF...'
#    cd primitive
#    mkdir -p dout
#    srun -n $NPROC $PW < scf.in > scf.out 2>&1
#    grep "convergence has been achieved" scf.out && echo "Primitive SCF OK" || { echo "Primitive SCF FAILED"; exit 1; }
#
#    echo '>>> Running primitive cell NSCF...'
#    srun -n $NPROC $PW < nscf.in > nscf.out 2>&1
#    grep "JOB DONE" nscf.out && echo "Primitive NSCF OK" || { echo "Primitive NSCF FAILED"; exit 1; }
#    cd ..
#else
#    echo '>>> Primitive cell data already exists, skipping SCF/NSCF'
#fi

# Step 2: Run EDI
echo '>>> Running EDI (edi.x)...'
cd edi
srun -n 1 $QEDIR/edi-code/src/extract_pot.x < extract_pot.in > extract_pot.out

srun -n $NPROC $EDI -nk $NPROC -i edi.setup.in > edi.setup.out 
srun -n $NPROC $EDI -nk $NPROC -i edi.in > edi.out 
cd ..

echo '>>> Done!'
