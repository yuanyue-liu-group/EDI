
mpirun_bin=ibrun
pw_bin=/home1/06975/zxiao/bin/pw.x

mpirun_bin=mpirun
pw_bin=/home/can/run/2-lu-defect-FFT-scattering/coding/q-e-M/bin/pw.x
pw_bin=/home/can/build/QE/q-e/bin/pw.x


prefix=mos2

for i in [1-9]x[1-9]
do
cd $i

ls
pwd

for j in pristine defect
do
cd $j
fp=${j:0:1}
$mpirun_bin -n 4 $pw_bin <$prefix.$i.$fp.in >$prefix.$i.$fp.out &
cd ../
done

cd ..

done
