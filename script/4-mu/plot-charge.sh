cat >plot-ki.gplt<<EOF

set term png

###########################
#plot ki
###########################

Ef=0.6779762221
Ec=0.6779762221

############################################################################################################
##1D plot vs E
############################################################################################################

# 1  2 3               4           5        6 7        8  9       10  11         12   13     14                 15
#ki: kxyz(crystal), gamma (s^-1), E(eV) ,  vx vy(cm/s),vx vy(au), f, df(eV^-1), dos, sum m, sum m*angleterm,   Nkf

set xlabel 'ki'
set ylabel 'E(eV)'

set out 'Elinear.png'
plot 'ki.plt' u 1:5   w p  title 'E' 

set xrange [0: 0.3]
set xlabel 'E-Ec(eV)'

#set yrange [0: 5e12]
set ylabel 'γ(s^{-1})'
set out 'gammavsE.png'
plot 'ki.plt' u (\$5)-Ec:(\$4)  w p lt 5 title 'γ(s^{-1})' 

#set yrange [0: 1.2e15]
set ylabel 'v^2((cm/s)^{2})'
set out '1d-v2vsE.png'
plot      'ki.plt' u (\$5)-Ec:((\$6)**2+(\$7)**2)   w p lt 5 title 'v^2'

set logscale y
set format y '%2.0e'
set ylabel 'v^2/gamma((cm/s)^{2})'
set out '1d-v2bgammavsE.png'
plot      'ki.plt' u (\$5)-Ec:((\$6)**2+(\$7)**2)/(\$4)   w p lt 5 title 'v^2/gamma'

set ylabel 'v^2/gamma*df((cm/s)^{2})'
set out '1d-v2bgammadfvsE.png'
plot      'ki.plt' u (\$5)-Ec:((\$6)**2+(\$7)**2)/(\$4)*(\$11)   w p lt 5 title 'v^2/gamma*df'

unset logscale y
unset format y 


unset yrange 
set ylabel 'sum_f δ(E_i-E_f)(eV^{-1})'
set out '1d-dosvsE.png'
plot      'ki.plt' u (\$5)-Ec:((\$12))   w p lt 5 title 'DOS'

set ylabel 'sum_f |M_{if}|^2(eV^2)'
set out '1d-summvsE.png'
plot      'ki.plt' u (\$5)-Ec:((\$13))   w p lt 5 title 'sum_f |M_{if}|^2'

set ylabel 'sum_f |M_{if}|^2(1-cosθ)(eV^2)'
set out '1d-summavsE.png'
plot      'ki.plt' u (\$5)-Ec:((\$14))   w p lt 5 title 'sum_f |M_{if}|^2(1-cosθ)'

unset yrange 
set ylabel 'mean_f δ(E_i-E_f)(eV^{-1})'
set out '1d-meanwtvsE.png'
plot      'ki.plt' u (\$5)-Ec:((\$12)/(\$15))   w p lt 5 title 'mean wt'

set ylabel 'mean_f |M_{if}|^2(eV^2)'
set out '1d-meanmvsE.png'
plot      'ki.plt' u (\$5)-Ec:((\$13)/(\$15))   w p lt 5 title 'mean_f |M_{if}|^2'

set ylabel 'mean_f |M_{if}|^2(1-cosθ)(eV^2)'
set out '1d-meanmavsE.png'
plot      'ki.plt' u (\$5)-Ec:((\$14)/(\$15))   w p lt 5 title 'mean_f |M_{if}|^2(1-cosθ)'

set ylabel 'mean_f |M_{if}|^2(eV^2)'
set out '1d-meanmvsE.png'
plot      'ki.plt' u (\$5)-Ec:((\$13)/(\$15))**0.5   w p lt 5 title 'rms_f |M_{if}|'


#set ylabel 'rms_f |M_{if}|(eV)'
#set out '1d-rmsmlmnlvsE.png'
#plot      'ki.plt' u (\$5)-Ec:((\$12)/(\$15))**.5   w p lt 5 title 'rms_f |Ml_{if}|',\
#          'ki.plt' u (\$5)-Ec:((\$13)/(\$15))**.5   w p lt 7 title 'rms_f |Mnl_{if}|'


#set ylabel 'rms_f |Ml_{if}|(eV)'
#set out '1d-meanmlvsE.png'
#plot      'ki.plt' u (\$5)-Ec:((\$12)/(\$15))**.5   w p lt 5 title 'rms_f |Ml_{if}|'

#set ylabel 'rms_f |Mnl_{if}|(eV)'
#set out '1d-meanmnlvsE.png'
#plot      'ki.plt' u (\$5)-Ec:((\$13)/(\$15))**.5   w p lt 5 title 'rms_f |Mnl_{if}|'




set ylabel 'Nkf'
set out '1d-nkfvsE.png'
plot      'ki.plt' u (\$5)-Ec:((\$15))   w p lt 5 title 'Nkf'



unset yrange
set logscale y
set format y '%2.0e'

set ylabel 'f'
set out '1d-fvsE.png'
plot      'ki.plt' u (\$5)-Ec:10   w p lt 5 title 'f'

set ylabel 'df/dE(ev^{-1})'
set out '1d-dfvsE.png'
plot      'ki.plt' u (\$5)-Ec:11   w p lt 5 title 'df'


#unset logscale y
#unset format y
#unset xrange
#unset yrange
##################################################
#2D plot vs kx ky
##################################################
# 1  2 3               4           5        6 7        8  9       10  11         12   13     14                 15
#ki: kxyz(crystal), gamma (s^-1), E(eV) ,  vx vy(cm/s),vx vy(au), f, df(eV^-1), dos, sum m, sum m*angleterm,   Nkf
reset
set size ratio -1
set xlabel 'kx(2pi/a)'
set ylabel 'ky(2pi/a)'

#set xrange [0.19: 0.30]
#set yrange [0.42: 0.51]

set out '2d-gamma.png'
 plot 'ki.plt' u (\$2)*0.866:(\$2)*0.5+(\$3):(\$4)*1e-12 w p lt 5 palette title 'γ(ps^{-1})'

set out '2d-E.png'
 plot 'ki.plt' u (\$2)*0.866:(\$2)*0.5+(\$3):5 w p lt 5 palette title 'E(eV)'

set out '2d-vmod.png'
 plot 'ki.plt' u (\$2)*0.866:(\$2)*0.5+(\$3):((\$6)**2+(\$7)**2)**.5 w p lt 5 palette title '|v|(cm/s)'

set out '2d-vvec.png'
 plot 'ki.plt' u (\$2)*0.866:(\$2)*0.5+(\$3):(\$6)*1e-10:(\$7)*1e-10 with vectors head filled lt 2  title 'v'

set out '2d-f.png'
 plot 'ki.plt' u (\$2)*0.866:(\$2)*0.5+(\$3):10 w p lt 5 palette title 'Fermi distribution'

set out '2d-df.png'
 plot 'ki.plt' u (\$2)*0.866:(\$2)*0.5+(\$3):11 w p lt 5 palette title 'df/dE (eV^{-1})'

set out '2d-dos.png'
 plot 'ki.plt' u (\$2)*0.866:(\$2)*0.5+(\$3):12 w p lt 5 palette title 'df/dE (eV^{-1})'

set out '2d-summ.png'
 plot 'ki.plt' u (\$2)*0.866:(\$2)*0.5+(\$3):13 w p lt 5 palette title 'df/dE (eV^{-1})'

set out '2d-summa.png'
 plot 'ki.plt' u (\$2)*0.866:(\$2)*0.5+(\$3):14 w p lt 5 palette title 'df/dE (eV^{-1})'

set out '2d-Nkf.png'
 plot 'ki.plt' u (\$2)*0.866:(\$2)*0.5+(\$3):15 w p lt 5 palette title 'df/dE (eV^{-1})'

EOF
cat >plot-kf.gplt<<EOF
set term png
###########################
#plot kf
###########################

Ef=0.6779762221
Ec=0.6779762221

#set xlabel '(1-cosθ)'
#set ylabel 'arg(M)(eV)'
#set out '1d-KFPLT-mvsA.png'
#plot 'KFPLT' u 12:9 w p lt 4  title 'arg(M)'


##################################################
#1D plot vs theta
##################################################

# 1   2    3            4    5  6   7    8     9    10  11        12                       13     14
#kf: kfx, kfy(crystal), Wt, M1 M2, |M|, arg(M),vf1,vf2, E(eV), angleterm(1-cos(theta)), thetai, thetaf

set xlabel 'θ'

set ylabel '|M|(eV)'
set out '1d-KFPLT-mvsA.png'
plot 'KFPLT' u (\$14)-(\$13):7 w p lt 4  title '|M|(eV)'

set ylabel 'arg(M)(eV)'
set out '1d-KFPLT-argmvsA.png'
plot 'KFPLT' u (\$14)-(\$13):8 w p lt 4  title 'arg(M)'


#set ylabel '|M_l|(eV)'
#set out '1d-KFPLT-mlvsA.png'
#plot 'KFPLT' u (\$14)-(\$13):((((\$5))*((\$5)))+(((\$6))*((\$6))))**.5 w p lt 2  title '|Ml|(eV)',\
#     'KFPLT' u (\$14)-(\$13):11 w p lt 3  title '|Ml|(eV)'
#
#set ylabel '|M_nl|(eV)'
#set out '1d-KFPLT-mnlvsA.png'
#plot 'KFPLT' u (\$14)-(\$13):((((\$7))*((\$7)))+(((\$8))*((\$8))))**.5 w p lt 2   title '|Mnl|(eV)',\
#     'KFPLT' u (\$14)-(\$13):12 w p lt 3   title '|Mnl|(eV)'

set ylabel 'wt=δ(E_i-E_f)(eV^{-1})'
set out '1d-KFPLT-wtvsA.png'
plot 'KFPLT' u (\$14)-(\$13):4 w p lt 2  title 'weight from δ(Ei-Ef)/Nk (eV^{-1})'

set ylabel 'v(cm/s)'
set out '1d-KFPLT-v2vsA.png'
plot 'KFPLT' u (\$14)-(\$13):((\$9)**2+(\$10)**2)**.5  w p lt 2  title 'v(cm/s)'

#set ylabel 'v_i(cm/s)'
#set out '1d-KFPLT-vvecvsA.png'
#plot 'KFPLT' u (\$14)-(\$13):(\$9)*1e-10:(\$10)*1e-10 with vectors head filled lt 2  title 'v'


set ylabel 'E_{final}(eV)'
set out '1d-KFPLT-EvsA.png'
plot 'KFPLT' u (\$14)-(\$13):11 w p lt 5  title 'E_{final}(eV)'



##################################################
#2D plot vs kx ky
##################################################

set size ratio -1
set xlabel 'kx(2pi/a)'
set ylabel 'ky(2pi/a)'


set xrange [0.19: 0.40]
set yrange [0.40: 0.61]

# 1   2    3            4    5  6   7    8     9    10  11        12                       13     14
#kf: kfx, kfy(crystal), Wt, M1 M2, |M|, arg(M),vf1,vf2, E(eV), angleterm(1-cos(theta)), thetai, thetaf

#set out '2d-KFPLT-m.png'
#plot 'KFPLT' u (\$2)*0.866:(\$2)*0.5+(\$3):((((\$5)-(\$7))*((\$5)-(\$7)))+(((\$6)-(\$8))*((\$6)-(\$8))))**.5 w p lt 5 palette title '|M|(eV)',\
#     'KFPLT' u (\$2)*0.866:(\$2)*0.5+(\$3):((\$9)**2+(\$10)**2)**.5 w p lt 5 palette title '|M|(eV)',\
#     'KFPLT' u (\$2)*0.866:(\$2)*0.5+(\$3):13 w p lt 5 palette title '|M|(eV)'

set out '2d-KFPLT-m.png'
plot 'KFPLT' u (\$2)*0.866:(\$2)*0.5+(\$3):7 w p lt 5 palette title '|M|(eV)'
set out '2d-KFPLT-argm.png'
plot 'KFPLT' u (\$2)*0.866:(\$2)*0.5+(\$3):8 w p lt 5 palette title 'arg(M)(eV)'

#set out '2d-KFPLT-ml.png'
#plot 'KFPLT' u (\$2)*0.866:(\$2)*0.5+(\$3):((((\$5))*((\$5)))+(((\$6))*((\$6))))**.5 w p lt 5 palette title '|Ml|(eV)',\
#     'KFPLT' u (\$2)*0.866:(\$2)*0.5+(\$3):11 w p lt 5 palette title '|Ml|(eV)'
#
#set out '2d-KFPLT-mnl.png'
#plot 'KFPLT' u (\$2)*0.866:(\$2)*0.5+(\$3):((((\$7))*((\$7)))+(((\$8))*((\$8))))**.5 w p lt 5 palette  title '|Mnl|(eV)',\
#     'KFPLT' u (\$2)*0.866:(\$2)*0.5+(\$3):12 w p lt 5 palette  title '|Mnl|(eV)'

set out '2d-KFPLT-wt.png'
plot 'KFPLT' u (\$2)*0.866:(\$2)*0.5+(\$3):4 w p lt 5 palette title 'weight from δ(Ei-Ef)/Nk (eV^{-1})'

set out '2d-KFPLT-v2.png'
plot 'KFPLT' u (\$2)*0.866:(\$2)*0.5+(\$3):((\$9)**2+(\$10)**2)**.5  w p lt 5 palette title 'v(cm/s)'

set out '2d-KFPLT-vvec.png'
plot 'KFPLT' u (\$2)*0.866:(\$2)*0.5+(\$3):(\$9)*4e-9:(\$10)*4e-9 with vectors head filled lt 2  title 'v'


set out '2d-KFPLT-E.png'
plot 'KFPLT' u (\$2)*0.866:(\$2)*0.5+(\$3):11 w p lt 5 palette title 'Ef(eV)'

set out '2d-KFPLT-angleterm.png'
plot 'KFPLT' u (\$2)*0.866:(\$2)*0.5+(\$3):12 w p lt 5 palette title '1-cosθ'


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


