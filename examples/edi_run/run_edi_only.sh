#!/bin/bash
#SBATCH -A che190065
#SBATCH -p shared
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -t 12:00:00
#SBATCH -J edi_only
#SBATCH -o edi_only.%j.out
#SBATCH -e edi_only.%j.err

module reset
module load aocc/3.1.0 openmpi/4.1.6
module load amdblis/3.0 amdlibflame/3.0 amdlibm/3.0 fftw

QEDIR=/anvil/projects/x-che190065/rjguo/qe-7.5
PW=$QEDIR/PW/src/pw.x
EDI=$QEDIR/edi-dev/src/edi.x
NPROC=108

cd /anvil/projects/x-che190065/rjguo/qe-7.5/edi-dev/test_edinterp/edi_run

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
srun -n 1 $QEDIR/edi-dev/extract_pot.x < extract_pot.in > extract_pot.out

srun -n $NPROC $EDI -nk $NPROC -i edi.setup.in > edi.setup.out 
srun -n $NPROC $EDI -nk $NPROC -i edi.in > edi.out 
cd ..

echo '>>> Done!'
