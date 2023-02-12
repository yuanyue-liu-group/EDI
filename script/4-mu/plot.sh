cat >plot-ki.gplt<<EOF

set term png

###########################
#plot ki
###########################

Ef=0.6779762221
Ec=0.6779762221
###########################
## vs E

###### E
set xlabel 'ki'
set ylabel 'E(eV)'

set out 'Elinear.png'
plot 'ki.plt' u 4   w lp  title 'E' 

set xrange [0: 0.3]
set xlabel 'E-Ec(eV)'

set yrange [0: 5e12]
set ylabel 'γ(s^{-1})'
set out 'gammavsE.png'
plot 'ki.plt' u (\$4)-Ec:(\$3)  w p lt 5 title 'γ(s^{-1})' 

set yrange [0: 1.2e15]
set ylabel 'v^2((cm/s)^{2})'
set out 'v2vsE.png'
plot      'ki.plt' u (\$4)-Ec:((\$5)**2+(\$6)**2)   w p lt 5 title 'v^2'

unset yrange 
set ylabel 'sum_f δ(E_i-E_f)(eV^{-1})'
set out 'dosvsE.png'
plot      'ki.plt' u (\$4)-Ec:((\$9))   w p lt 5 title 'DOS'

set ylabel 'sum_f |M_{if}|^2(eV^2)'
set out 'summvsE.png'
plot      'ki.plt' u (\$4)-Ec:((\$10))   w p lt 5 title 'sum_f |M_{if}|^2'

set ylabel 'sum_f |M_{if}|^2(1-cosθ)(eV^2)'
set out 'summavsE.png'
plot      'ki.plt' u (\$4)-Ec:((\$11))   w p lt 5 title 'sum_f |M_{if}|^2(1-cosθ)'

unset yrange 
set ylabel 'mean_f δ(E_i-E_f)(eV^{-1})'
set out 'meanwtvsE.png'
plot      'ki.plt' u (\$4)-Ec:((\$9)/(\$14))   w p lt 5 title 'mean wt'

set ylabel 'mean_f |M_{if}|^2(eV^2)'
set out 'meanmvsE.png'
plot      'ki.plt' u (\$4)-Ec:((\$10)/(\$14))   w p lt 5 title 'mean_f |M_{if}|^2'

set ylabel 'mean_f |M_{if}|^2(1-cosθ)(eV^2)'
set out 'meanmavsE.png'
plot      'ki.plt' u (\$4)-Ec:((\$11)/(\$14))   w p lt 5 title 'mean_f |M_{if}|^2(1-cosθ)'

set ylabel 'rms_f |Ml_{if}|(eV)'
set out 'meanmlvsE.png'
plot      'ki.plt' u (\$4)-Ec:((\$12)/(\$14))**.5   w p lt 5 title 'rms_f |Ml_{if}|'

set ylabel 'rms_f |Mnl_{if}|(eV)'
set out 'meanmnlvsE.png'
plot      'ki.plt' u (\$4)-Ec:((\$13)/(\$14))**.5   w p lt 5 title 'rms_f |Mnl_{if}|^2'

set ylabel 'rms_f |M_{if}|(eV)'
set out 'rmsmlmnlvsE.png'
plot      'ki.plt' u (\$4)-Ec:((\$12)/(\$14))**.5   w p lt 5 title 'rms_f |Ml_{if}|',\
          'ki.plt' u (\$4)-Ec:((\$13)/(\$14))**.5   w p lt 7 title 'rms_f |Mnl_{if}|^2'



set ylabel 'Nkf'
set out 'nkfvsE.png'
plot      'ki.plt' u (\$4)-Ec:((\$12))   w p lt 5 title 'Nkf'



unset yrange
set logscale y
set format y '%2.0e'

set ylabel 'f'
set out 'fvsE.png'
plot      'ki.plt' u (\$4)-Ec:7   w p lt 5 title 'f'

set ylabel 'df/dE(ev^{-1})'
set out 'dfvsE.png'
plot      'ki.plt' u (\$4)-Ec:8   w p lt 5 title 'df'


#unset logscale y
#unset format y
#unset xrange
#unset yrange

reset
set size ratio -1
set xlabel 'kx(2pi/a)'
set ylabel 'ky(2pi/a)'

#set xrange [0.19: 0.30]
#set yrange [0.42: 0.51]

set out 'gamma.png'
 plot 'ki.plt' u (\$1)*0.866:(\$1)*0.5+(\$2):(\$3)*1e-12 w p lt 5 palette title 'γ(ps^{-1})'


set out 'E.png'
 plot 'ki.plt' u (\$1)*0.866:(\$1)*0.5+(\$2):4 w p lt 5 palette title 'E(eV)'

set out 'f.png'
 plot 'ki.plt' u (\$1)*0.866:(\$1)*0.5+(\$2):7 w p lt 5 palette title 'Fermi distribution'

set out 'df.png'
 plot 'ki.plt' u (\$1)*0.866:(\$1)*0.5+(\$2):8 w p lt 5 palette title 'df/dE (eV^{-1})'

set out 'vmod.png'
 plot 'ki.plt' u (\$1)*0.866:(\$1)*0.5+(\$2):((\$5)**2+(\$6)**2)**.5 w p lt 5 palette title '|v|(cm/s)'

set out 'vvec.png'
 plot 'ki.plt' u (\$1)*0.866:(\$1)*0.5+(\$2):(\$5)*1e-10:(\$6)*1e-10 with vectors head filled lt 2  title 'v'

EOF
cat >plot-kf.gplt<<EOF
set term png
###########################
#plot kf
###########################

Ef=0.6779762221
Ec=0.6779762221
###########################
## vs theta


set xlabel '(1-cosθ)'


set ylabel '|M|(eV)'
set out 'KFPLT-mvsA.png'
plot 'KFPLT' u 17:((((\$5)-(\$7))*((\$5)-(\$7)))+(((\$6)-(\$8))*((\$6)-(\$8))))**.5 w p lt 2  title '|M|(eV)',\
     'KFPLT' u 17:((\$9)**2+(\$10)**2)**.5 w p lt 3  title '|M|(eV)',\
     'KFPLT' u 17:13 w p lt 4  title '|M|(eV)'

set ylabel '|M_l|(eV)'
set out 'KFPLT-mlvsA.png'
plot 'KFPLT' u 17:((((\$5))*((\$5)))+(((\$6))*((\$6))))**.5 w p lt 2  title '|Ml|(eV)',\
     'KFPLT' u 17:11 w p lt 3  title '|Ml|(eV)'

set ylabel '|M_nl|(eV)'
set out 'KFPLT-mnlvsA.png'
plot 'KFPLT' u 17:((((\$7))*((\$7)))+(((\$8))*((\$8))))**.5 w p lt 2   title '|Mnl|(eV)',\
     'KFPLT' u 17:12 w p lt 3   title '|Mnl|(eV)'

set ylabel 'wt=δ(E_i-E_f)(eV^{-1})'
set out 'KFPLT-wtvsA.png'
plot 'KFPLT' u 17:4 w p lt 2  title 'weight from δ(Ei-Ef)/Nk (eV^{-1})'

set ylabel 'v(cm/s)'
set out 'KFPLT-v2vsA.png'
plot 'KFPLT' u 17:((\$14)**2+(\$15)**2)**.5  w p lt 2  title 'v(cm/s)'

#set ylabel 'v_i(cm/s)'
#set out 'KFPLT-vvecvsA.png'
#plot 'KFPLT' u 17:(\$14)*1e-10:(\$15)*1e-10 with vectors head filled lt 2  title 'v'


set ylabel 'E_{final}(eV)'
set out 'KFPLT-EvsA.png'
plot 'KFPLT' u 17:16 w p lt 5  title 'E_{final}(eV)'



###########################
## vs 2D k

set size ratio -1
set xlabel 'kx(2pi/a)'
set ylabel 'ky(2pi/a)'


set xrange [0.19: 0.40]
set yrange [0.40: 0.61]


set out 'KFPLT-m.png'
plot 'KFPLT' u (\$2)*0.866:(\$2)*0.5+(\$3):((((\$5)-(\$7))*((\$5)-(\$7)))+(((\$6)-(\$8))*((\$6)-(\$8))))**.5 w p lt 5 palette title '|M|(eV)',\
     'KFPLT' u (\$2)*0.866:(\$2)*0.5+(\$3):((\$9)**2+(\$10)**2)**.5 w p lt 5 palette title '|M|(eV)',\
     'KFPLT' u (\$2)*0.866:(\$2)*0.5+(\$3):13 w p lt 5 palette title '|M|(eV)'

set out 'KFPLT-ml.png'
plot 'KFPLT' u (\$2)*0.866:(\$2)*0.5+(\$3):((((\$5))*((\$5)))+(((\$6))*((\$6))))**.5 w p lt 5 palette title '|Ml|(eV)',\
     'KFPLT' u (\$2)*0.866:(\$2)*0.5+(\$3):11 w p lt 5 palette title '|Ml|(eV)'

set out 'KFPLT-mnl.png'
plot 'KFPLT' u (\$2)*0.866:(\$2)*0.5+(\$3):((((\$7))*((\$7)))+(((\$8))*((\$8))))**.5 w p lt 5 palette  title '|Mnl|(eV)',\
     'KFPLT' u (\$2)*0.866:(\$2)*0.5+(\$3):12 w p lt 5 palette  title '|Mnl|(eV)'

set out 'KFPLT-wt.png'
plot 'KFPLT' u (\$2)*0.866:(\$2)*0.5+(\$3):4 w p lt 5 palette title 'weight from δ(Ei-Ef)/Nk (eV^{-1})'

set out 'KFPLT-v2.png'
plot 'KFPLT' u (\$2)*0.866:(\$2)*0.5+(\$3):((\$14)**2+(\$15)**2)**.5  w p lt 5 palette title 'v(cm/s)'

set out 'KFPLT-vvec.png'
plot 'KFPLT' u (\$2)*0.866:(\$2)*0.5+(\$3):(\$14)*4e-9:(\$15)*4e-9 with vectors head filled lt 2  title 'v'


set out 'KFPLT-E.png'
plot 'KFPLT' u (\$2)*0.866:(\$2)*0.5+(\$3):16 w p lt 5 palette title 'Ef(eV)'



set xrange [0.19: 0.41]
set yrange [0.42: 0.61]

set xrange [0.5: 0.65]
set yrange [0.9: 1.1]

EOF

gnuplot <plot-ki.gplt
mkdir figs
mv *.png figs
cd kf-kis
for i in *.plt
do
cp ../plot-kf.gplt plot-kf$i.gplt
sed -i "s/KFPLT/$i/g"  plot-kf$i.gplt
gnuplot < plot-kf$i.gplt
done


