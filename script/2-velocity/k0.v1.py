import sys 
import numpy as np 
import copy 
import scipy.sparse as sp 
import io_files 
import qirr_sym


def f(ei,ef,kbt):
    return 1/(np.exp((ei-ef)/kbt)+1)

if __name__ == '__main__':

 for nqf in [144]: 
  for [alat,sys] in zip([3.184,3.3197,3.1859,3.3197],['mos2','mose2','ws2','wse2']):
    out_folder = "./"+sys+'-%d/'%nqf
    velocity_fn="v.dat"
    kbt=0.025
    kbtha=0.025/27.2

    RangeE=0.2
    Ne=400
    efl=np.array(range(Ne+1))*RangeE/Ne*2-RangeE
    #efl=[-0.05,0,0.05]
    #efl=np.array(range(400))*0.001-0.2
    
    pseudo_freq = np.zeros((nqf*nqf, 1))


    _, bande1, velocity1 = io_files.reader_velocity(out_folder+velocity_fn , nelec=12, ibndlist=[1]) 
    _, bande2, velocity2 = io_files.reader_velocity(out_folder+velocity_fn , nelec=12, ibndlist=[2]) 
    _, bande3, velocity3 = io_files.reader_velocity(out_folder+velocity_fn , nelec=12, ibndlist=[3]) 
    _, bande4, velocity4 = io_files.reader_velocity(out_folder+velocity_fn , nelec=12, ibndlist=[4]) 

    print(np.max(bande1),np.max(bande2),np.min(bande3),np.min(bande4))
    [ev1,ev,ec,ec1]=[np.max(bande1),np.max(bande2),np.min(bande3),np.min(bande4)]
    print([ev1,ev,ec,ec1])
    #print(velocity1) 

    lne='#VBM-1 VBM CBM CBM+1 (eV): %f  %f %f %f\n'%(ev1,ev,ec,ec1)
    lne=lne+'#ef(eV),  ne(cm^-2) VBM-1 VBM CBM CBM+1 , hole, electron\n'
    lme='#VBM-1 VBM CBM CBM+1 (eV): %f  %f %f %f\n'%(ev1,ev,ec,ec1)
    lme=lme+'#ef(eV),  me(au) VBM-1 VBM CBM CBM+1 , hole, electron\n'
    lchi='#VBM-1 VBM CBM CBM+1 (eV): %f  %f %f %f\n'%(ev1,ev,ec,ec1)
    lchi=lchi+'#ef(eV),  chi(au) VBM-1 VBM CBM CBM+1 , hole, electron\n'
    lall=''
    omega=(alat**2*3**0.5/2.0*1e-16)
    omegaau=((alat/.529)**2*3**0.5/2.0)
    for ef in efl:
      ne1=0.0
      me1=0.0
      chi1=0.0
      for [ei ,vi] in zip( bande1,velocity1):
        fev=f(ei[0]-ev,ef,kbt)
        ne1=ne1+1-fev
        chi1=chi1+(1-fev)*fev/kbtha
        me1=me1+(1-fev)*fev*0.5*np.linalg.norm(np.array(vi[0]))**2/kbt
      ne1=ne1/nqf**2/omega
      me1=me1/nqf**2/omega*(14.6**-2*27.6)
      chi1=chi1/nqf**2/omegaau

      ne2=0.0
      me2=0.0
      chi2=0.0
      for [ei ,vi] in zip( bande2,velocity2):
        fev=f(ei[0]-ev,ef,kbt)
        ne2=ne2+1-fev
        chi2=chi2+(1-fev)*fev/kbtha
        me2=me2+(1-fev)*fev*0.5*np.linalg.norm(np.array(vi[0]))**2/kbt
      ne2=ne2/nqf**2/omega
      me2=me2/nqf**2/omega*(14.6**-2*27.6)
      chi2=chi2/nqf**2/omegaau

      ne3=0.0
      me3=0.0
      chi3=0.0
      for [ei ,vi] in zip( bande3,velocity3):
        fec=f(ei[0]-ec,ef,kbt)
        ne3=ne3+fec
        chi3=chi3+(1-fec)*fec/kbtha
        me3=me3+(1-fec)*fec*0.5*np.linalg.norm(np.array(vi[0]))**2/kbt
      ne3=ne3/nqf**2/omega
      me3=me3/nqf**2/omega*(14.6**-2*27.6)
      chi3=chi3/nqf**2/omegaau


      ne4=0.0
      me4=0.0
      chi4=0.0
      for [ei ,vi] in zip( bande4,velocity4):
        fec=f(ei[0]-ec,ef,kbt)
        ne4=ne4+fec
        chi4=chi4+(1-fec)*fec/kbtha
        me4=me4+(1-fec)*fec*0.5*np.linalg.norm(np.array(vi[0]))**2/kbt
      ne4=ne4/nqf**2/omega
      me4=me4/nqf**2/omega*(14.6**-2*27.6)
      chi4=chi4/nqf**2/omegaau

      lne  =lne+'%8.4f    %14.12e %14.12e %14.12e %14.12e %14.12e %14.12e \n'%(ef,ne1,ne2,ne3,ne4,(ne1+ne2),(ne3+ne4))
      lme  =lme+'%8.4f    %14.12e %14.12e %14.12e %14.12e %14.12e %14.12e \n'%(ef,(me1/ne1)**-1,(me2/ne2)**-1,(me3/ne3)**-1,(me4/ne4)**-1,((me1+me2)/(ne1+ne2))**-1,((me3+me4)/(ne3+ne4))**-1)
      lchi=lchi+'%8.4f    %14.12e %14.12e %14.12e %14.12e %14.12e %14.12e \n'%(ef,chi1,chi2,chi3,chi4,(chi1+chi2),(chi3+chi4))
      
      lall=lall+'%8.4f    %6.4e %6.4e %6.4e %6.4e %6.4e %6.4e   '%(ef,ne1,ne2,ne3,ne4,(ne1+ne2),(ne3+ne4))
      lall=lall+'m m m m  %6.4e %6.4e %6.4e %6.4e %6.4e %6.4e   '%(   (me1/ne1)**-1,(me2/ne2)**-1,(me3/ne3)**-1,(me4/ne4)**-1,((me1+me2)/(ne1+ne2))**-1,((me3+me4)/(ne3+ne4))**-1)
      lall=lall+'c c c c  %6.4e %6.4e %6.4e %6.4e %6.4e %6.4e \n'%(   chi1,chi2,chi3,chi4,(chi1+chi2),(chi3+chi4))

    with open ('2.n-vs-ef.'+sys+'-%d'%nqf+'.dat','w') as fh:
      fh.write(lne)
    with open ('2.m-vs-ef.'+sys+'-%d'%nqf+'.dat','w') as fh:
      fh.write(lme)
    with open ('2.chi-vs-ef.'+sys+'-%d'%nqf+'.dat','w') as fh:
      fh.write(lchi)
    with open ('2.all-vs-ef.'+sys+'-%d'%nqf+'.dat','w') as fh:
      fh.write(lall)

