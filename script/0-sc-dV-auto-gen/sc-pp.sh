cat >pp.in<<EOF
&inputpp
   prefix='mos2',
   outdir = 'dout/'
   filplot = 'vloc'
   plot_num= 1
/
&plot
   nfile = 1
   filepp(1) = 'vloc'
   weight(1) = 1.0
   iflag = 3
   output_format = 5
   fileout = 'vloc.xsf'
/
EOF

mpirun_bin=ibrun
pw_bin=/home1/06975/zxiao/bin/pw.x

mpirun_bin=mpirun
pw_bin=/home/can/run/2-lu-defect-FFT-scattering/coding/q-e-M/bin/pw.x
pw_bin=/home/can/build/QE/q-e/bin/pw.x
pp_bin=/home/can/build/QE/q-e/bin/pp.x




for i in [1-9]x[1-9]
 do
 echo $i
 cd $i
 for j in *
  do
  cd $j
  cp ../../pp.in .

  $mpirun_bin -n 4 $pp_bin <pp.in >pp.out &
  cd ..
 done
 cd ..
done
