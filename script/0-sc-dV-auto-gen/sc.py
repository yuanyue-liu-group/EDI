from numpy import array
import os

inputheadlines='''
&control
   calculation = '%s'
   verbosity='high'
   restart_mode='from_scratch',
   prefix='mos2',
   pseudo_dir = '  /home/can/run/2-lu-defect-FFT-scattering/pseudo',
 
   outdir='dout/'
/
&system
    ibrav=  0, 
    A=3.09
    nat=  %d, ntyp= 5,
    ecutwfc =60.0,
    !nbnd=14
    occupations='smearing'
    smearing='gaussian'
    degauss=0.001
/
&electrons
    !electron_maxstep=2
    diagonalization='david'
    mixing_mode = 'plain'
    mixing_beta = 0.7
    conv_thr =  1.0d-6
/

&ions
/


&cell
cell_dofree='2Dxy'
/


&calcmcontrol

kpoint_initial=4
kpoint_final=10
bnd_initial=10
bnd_final=10
vperturb_filename='vloc.dat'
eps_filename='eps.dat'
/


ATOMIC_SPECIES
Mo    95.96   Mo.pw-mt_fhi.UPF
S     32.06    S.pw-mt_fhi.UPF
W    183.84    W.pw-mt_fhi.UPF
Se    78.96   Se.pw-mt_fhi.UPF
O     16.00    O.pw-mt_fhi.UPF

CELL_PARAMETERS alat
'''



kpline='''
K_POINTS gamma
'''
moc_init=['Mo',  0.0,   0.0, 0.5]
s1c_init=['S ',  1./3, 2./3, 0.5626550000]
s2c_init=['S ',  1./3, 2./3, 0.4373450000]
atpc=[moc_init,s1c_init,s2c_init]
iatm_vac=1# index of vacancy atom

# move vacancy to center
moc=['Mo',   0.0-atpc[iatm_vac][1],  0.0-atpc[iatm_vac][2], 0.5000000000-0.5]
s1c=['S ',  1./3-atpc[iatm_vac][1], 2./3-atpc[iatm_vac][2], 0.5626550000-0.5]
s2c=['S ',  1./3-atpc[iatm_vac][1], 2./3-atpc[iatm_vac][2], 0.4373450000-0.5]

atpc=[moc,s1c,s2c]

vlat_p1=array([1.0,         0.0,         0.0])
vlat_p2=array([-0.5,  3**0.5/2.0,        0.0])
vlat_p3=array([0.0,         0.0,         8.0])

Nscx=2
Nscy=2


for Nsc in range(1,10):
    Nscx=Nsc
    Nscy=Nsc
    inputheadlines_p=inputheadlines%('scf',len(atpc)*Nscx*Nscy)
    inputheadlines_d=inputheadlines%('relax',len(atpc)*Nscx*Nscy-1)
    
    inputheadlines_p=inputheadlines_p+'%18.10f %18.10f %18.10f\n'%tuple(vlat_p1*Nscx)
    inputheadlines_p=inputheadlines_p+'%18.10f %18.10f %18.10f\n'%tuple(vlat_p2*Nscy)
    inputheadlines_p=inputheadlines_p+'%18.10f %18.10f %18.10f\n'%tuple(vlat_p3)
    
    inputheadlines_d=inputheadlines_d+'%18.10f %18.10f %18.10f\n'%tuple(vlat_p1*Nscx)
    inputheadlines_d=inputheadlines_d+'%18.10f %18.10f %18.10f\n'%tuple(vlat_p2*Nscy)
    inputheadlines_d=inputheadlines_d+'%18.10f %18.10f %18.10f\n'%tuple(vlat_p3)

    atomsline_p='ATOMIC_POSITIONS crystal\n'
    atomsline_d='ATOMIC_POSITIONS crystal\n'
    
#    atom=atpc[iatm_vac]
#    atc_vac=(atom[0],(atom[1]+Nscx/2*1.0)/Nscx,(atom[2]+Nscy/2*1.0)/Nscy,atom[3])
    for i in range(Nscx):
     for j in range(Nscy):
      for atom in atpc:
        atc=(atom[0],(atom[1]+i*1.0)/Nscx,(atom[2]+j*1.0)/Nscy,atom[3])
        atomsline_p=atomsline_p+'%s     %18.10f %18.10f %18.10f\n'%atc
        if i==Nscx/2 and j == Nscy/2:
          if atom[3]<=0.5:
            atomsline_d=atomsline_d+'%s     %18.10f %18.10f %18.10f\n'%atc
        else:
          atomsline_d=atomsline_d+'%s     %18.10f %18.10f %18.10f\n'%atc
    
    dirn='%dx%d/'%(Nscx,Nscy)
    if not os.path.isdir(dirn):   os.mkdir(dirn)
    if not os.path.isdir(dirn+'/pristine'):  os.mkdir(dirn+'/pristine')
    if not os.path.isdir(dirn+'/defect'):  os.mkdir(dirn+'/defect')
    
    inputlines=inputheadlines_p+atomsline_p+kpline
    fn_input=dirn+'pristine/mos2.%dx%d.p.in'%(Nscx,Nscy)
    fh_input=open(fn_input,'w')
    fh_input.write(inputlines)
    
    inputlines=inputheadlines_d+atomsline_d+kpline
    fn_input=dirn+'defect/mos2.%dx%d.d.in'%(Nscx,Nscy)
    fh_input=open(fn_input,'w')
    fh_input.write(inputlines)
    
