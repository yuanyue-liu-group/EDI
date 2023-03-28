import numpy as np
#in unit a: lattice constant of fcc
rofuc_relative=np.array([[0,0,0],[0.25,0.25,0.25]])
avecofuc=np.array([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]])
Roffcc_relative=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
usecubic=True

Nsc=2
alat=5.4298 # in Angstrom
for isc in range(Nsc):
 for jsc in range(Nsc):
  for ksc in range(Nsc):
    for atomuc in rofuc_relative:
      if usecubic:
       for R in Roffcc_relative:
        r=np.dot(atomuc+R,avecofuc)+[isc,jsc,ksc]+Nsc/2
        [x,y,z]=r*alat
        print('%10.8f %10.8f %10.8f'%(x,y,z))
      else:
        r=np.dot(np.dot(atomuc,avecofuc))+[isc,jsc,ksc]
        [x,y,z]=r*alat
        print('%10.8f %10.8f %10.8f'%(x,y,z))
