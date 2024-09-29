#from numpy import *
import json
import numpy as np
import pickle
import os
import pickle
import cmath

def get_velocity_angle(vix,viy):
   if vix>eps :
       theta=np.arctan(viy/vix)
   elif vix<-eps :
       theta=np.arctan(viy/vix)+pi
   elif abs(vix)<=eps :
       if viy>0:
           theta=pi/2
       elif viy<0:
           theta=-pi/2
       elif viy==0:
           theta=0
   return theta
 
def calc_gamma(M,Wt,V,withangle):

    gamma={}
    angleterm={}
    theta={}
    summ={}
    summa={}
    Nkf={}
    dos={}
    dosspin={}

    for ki in M.keys():
        vix=V[ki][0]
        viy=V[ki][1]
        vimod=(vix**2+viy**2)**.5
        _thetai=get_velocity_angle(vix,viy)

        theta[ki]=_thetai
        Nkf[ki]=len(M[ki])
        gamma[ki]=0
        summ[ki]=0
        summa[ki]=0
        dos[ki]=0
        dosspin[ki]=0
        angleterm[ki]={}

        for kf in M[ki].keys():
            vfx=V[kf][0]
            vfy=V[kf][1]
            vfmod=(vfx**2+vfy**2)**.5
            _thetaf=get_velocity_angle(vfx,vfy)
            if theta.get(kf)!=None:
               if abs(theta[kf]-_thetaf)>1e-4:
                 print 'velocity inconsistency, kf, theta[kf],_thetaf', ki,kf, theta[kf],_thetaf
            else:
               theta[kf]=_thetaf

            _angleterm=1
            if withangle:
               _angleterm=1-(vix*vfx+viy*vfy)/vimod/vfmod
               _angletermp=1-np.cos(_thetai-_thetaf)
               if abs(_angleterm-_angletermp)>1e-2:
                   print 'angle term not fit: cos(dtheta)!=vi.vf/|vi||vf|, %f!=%f'%(_angletermp,_angleterm)
            angleterm[ki][kf]=_angleterm
            gamma[ki]=gamma[ki]+abs(M[ki][kf])**2*Wt[ki][kf]*angleterm[ki][kf] 
            summa[ki]=summa[ki]+abs(M[ki][kf])**2*angleterm[ki][kf] 
            summ[ki]=summ[ki]+abs(M[ki][kf])**2
            dos[ki]=dos[ki]+Wt[ki][kf]
            if ki[0]==kf[0]:
                dosspin[ki]=dosspin[ki]+Wt[ki][kf]

        gamma[ki]=gamma[ki]*2/hbar*Cd
        #gamma[ki]=gamma[ki]*2*pi/hbar*Cd
        # pi is in wt
    return [gamma, angleterm,theta,summ,summa,Nkf,dos,dosspin]
   
def calc_mu(gamma,V,f,df,Emin,Emax):
    mu=0
    Nc=0
    Ncounter=0
    for ki in gamma.keys():
      if E[ki]>Emin and E[ki]<Emax:
        if gamma[ki] >1e6:
          mu=mu+0.5*np.linalg.norm(V[ki], ord=2)**2/gamma[ki]*df[ki]
        if eh=='e':
          Nc=Nc+f[ki]
        if eh=='h':
          Nc=Nc+1-f[ki]
        Ncounter=Ncounter+1
    #    print E[ki],Nc/Ngrid**2,f[ki]
    if Nc>0:
      mu=mu/Nc
      Nc=Nc/Ngrid**2
    # mu in (eVA/hbar)^2*s*eV^-1=hbar_in_ev^-2 e_in_SI^-1 A_in_cm^2 cm^2 V^-1 s^-1=6.5821e-16^2 * 1.6e-19 *1e-16=0.00003693099 wrong
    # mu in (eVA/hbar)^2*s*eV^-1=hbar_in_ev^-2 e_in_SI^-1 A_in_cm^2 cm^2 V^-1 s^-1=6.5821e-16^2 * 1       *1e-16=2.3081873e+14 
    mu=mu*2.3081873e+14 
    Nc=Nc/(alat**2*0.866*1e-16)
    return [mu,Nc,Ncounter]


#############chenmu code
#  mob_w=v^2(1-f)f/kT/nktot
#  ncarrier=2 * sum_nk sum_nbnd f / nktot - nelec2 
#  selfen_k_im = sum_ibnd sum_jbnd        sum_iq w*g2*angle_term*modes_weight
#  selfen_k_im = sum_ibnd sum_jbnd sum_im sum_iq w*g2*angle_term*modes_weight
#  mob=sum_nktot sum_nbnd  1/gamma*mob_w= sum_nktot sum_nbnd  1/gamma*v^2(1-f)f/kT/nktot


########## constant
#####    pi hbar e nat Cd

########## data files
#ki file:  E v
#kif file: M wt

########## fomular
#####    wt=sum Ea*Ea/(Ea-Eb)/(Ea-Ec)/(Ea-Ec) /Nktot*pi
#####    tau^-1= 2pi/hbar nat*Cd/Nk' sum_nk' |M_nk,nk'|^2 delta(Enk-Enk')(1-cos theta)
#####    mu=2* e/Ncarrier sum_nk tau v^2(-df/dE)
 
##########    to print
#####    ki: E, tau, f, df/dE, v, 
#####    kf: dE, wt, Mif, theta 




##################################
##  constant
##################################

# units: eV, s, A 
eps=1e-16

kbt=0.0256 #eV
hbar=6.5821e-16 #eV s
pi=3.1415926
RytoeV=13.6
Ha2eV=27.2
evinv_a2cminv_s = 1.519267582e7
evinv_a2cminv_s = 4.25473778848e6 #2.18769126364e8/(27.2/.529)

Ngrid=180
Ngrid=144

alat=3.18
Cd=1e12*(alat**2*0.866*1e-16)
print('Cd',Cd)

Ec=-0.31132
Ef=Ec-0.04



withangle=True
restart=True
eh='e'

fppn='pp.dat'
fppnl=open(fppn,'r')
ml=fppnl.readlines()

E={} #eV
V={} #eV/A
f={}
df={}#eV^-1
Kxyz={}
M={} #eV
Wt={}
nkil=set(())
kil=set(())
nfl=set(())
kfl=set(())

for l in ml[1:]:
    data=['']+l.split()
    #print(data)
    ni=int(data[1])
    ki=int(data[2])
    nkil.add((ni,ki))

for i in nkil:
    M[i]={}
    Wt[i]={}

for l in ml[1:]:
    data=['']+l.split()
    ni=int(data[1])
    ki=int(data[2])
    kxyzi=np.array([float(data[3]),float(data[4]),float(data[5])])
    ei=float(data[6])
    vxyzi=np.array([float(data[7]),float(data[8]),float(data[9])])

    nf=int(data[10])
    kf=int(data[11])
    kxyzf=np.array([float(data[12]),float(data[13]),float(data[14])])
    ef=float(data[15])
    vxyzf=np.array([float(data[16]),float(data[17]),float(data[18])])

    wt=float(data[19])

    mt=data[20][1:-1].split(',')
    m= float(mt[0])+1j*float(mt[1])

    
    M[(ni,ki)][(nf,kf)]=m*Ha2eV
    M[(ni,ki)][(nf,kf)]=m*RytoeV
    Wt[(ni,ki)][(nf,kf)]=wt

    E[(ni,ki)]=ei
    E[(nf,kf)]=ef

    V[(ni,ki)]=vxyzi
    V[(nf,kf)]=vxyzf

    Kxyz[(ni,ki)]=kxyzi
    Kxyz[(nf,kf)]=kxyzf

    f[(ni,ki)]=1/(1+np.exp((ei-Ef)/kbt))
    df[(ni,ki)]=(1-f[(ni,ki)])*f[(ni,ki)]/kbt

    f[(nf,kf)]=1/(1+np.exp((ef-Ef)/kbt))
    df[(nf,kf)]=(1-f[(nf,kf)])*f[(nf,kf)]/kbt

[gamma, angleterm,theta,summ,summa,Nkf,dos,dosspin]=calc_gamma(M,Wt,V,withangle)
[mu,Nc,Ncounter]=calc_mu(gamma,V,f,df,Emin=-10,Emax=10)
###########################################################################################################
###    debug:           print ki property: 
###########################################################################################################

pl= '#mu: %f, sigma: %f, nc: %e\n'%(mu,mu*Nc,Nc)
pl= pl+'#ni ki: kxyz(crystal), gamma (s^-1), E(eV) ,  vx vy(cm/s),vx vy(au), f, df(eV^-1), dos, sum m, sum m*angleterm,   Nkf\n'
for ki in M.keys():
   lp=tuple([ki[0],ki[1],Kxyz[ki][0],Kxyz[ki][1],gamma[ki],E[ki],V[ki][0]*evinv_a2cminv_s,V[ki][1]*evinv_a2cminv_s,V[ki][0],V[ki][1],f[ki],df[ki],dos[ki],summ[ki],summa[ki],Nkf[ki],dosspin[ki]])
   pl=pl+'%d %d %10.14e %10.14e  %10.14e %10.14e %10.14e %10.14e %10.14e %10.14e  %10.14e   %10.14e  %10.14e     %10.14e   %10.14e   %d  %10.14e   \n'%lp


file_plotki_name='ki.plt'
file_plotki_handle=open(file_plotki_name,'w')
file_plotki_handle.write(pl)



###########################################################################################################
###    debug:           print ki-kf property: 
###########################################################################################################
kfdir='kf-kis/'
if not os.path.exists(kfdir):
   os.mkdir(kfdir)

for ki in M.keys():
    pl= '#ni, ki: kxyz(crystal), gamma (s^-1), E(eV) ,  vx vy(cm/s),vx vy(au), f, df(eV^-1), dos, sum m, sum m*angleterm,   Nkf\n'
    lp=tuple([ki[0],ki[1],Kxyz[ki][0],Kxyz[ki][1],gamma[ki],E[ki],V[ki][0]*evinv_a2cminv_s,V[ki][1]*evinv_a2cminv_s,V[ki][0],V[ki][1],f[ki],df[ki],dos[ki],summ[ki],summa[ki],Nkf[ki]])
    pl=pl+'# %d %d : %10.14e %10.14e  %10.14e %10.14e %10.14e %10.14e %10.14e %10.14e  %10.14e   %10.14e  %10.14e     %10.14e   %10.14e   %d \n'%lp
    pl=pl+'#nf kf: kfx, kfy(crystal), Wt, M1 M2, |M|, arg(M),vf1,vf2, E(eV), angleterm(1-cos(theta)), thetai, thetaf\n'
    for kf in M[ki].keys():
       lp=tuple([ kf[0],kf[1], Kxyz[kf][0],Kxyz[kf][1],Wt[ki][kf],M[ki][kf].real,M[ki][kf].imag,abs(M[ki][kf]),cmath.phase(M[ki][kf]), V[kf][0],V[kf][1],E[kf],angleterm[ki][kf],theta[ki],theta[kf]])
       pl=pl+'%d %d %10.14e %10.14e     %10.14e    %10.14e %10.14e  %10.14e %10.14e     %10.14e %10.14e     %10.14e   %10.14e %10.14e  %10.14e \n'%lp
       
    file_plotkf_name=kfdir+'/n-%d-k-%d.plt'%ki
    file_plotkf_handle=open(file_plotkf_name,'w')
    file_plotkf_handle.write(pl)

##########################################################################################################
##   END
##########################################################################################################
