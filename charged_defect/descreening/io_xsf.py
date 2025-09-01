import numpy as np

class Rho(object):

    def __init__(self,name_):
        self.name = name_
        self.file = open(self.name)
        self.f = self.file.readlines()

        self.a1 = self.f[2].split()
        self.a2 = self.f[3].split()
        self.a3 = self.f[4].split()
        self.Natom = self.f[6].split()
        self.Ngrids = self.f[10+int(self.Natom[0])].split()
        self.nx = int(self.Ngrids[0])
        self.ny = int(self.Ngrids[1])
        self.nz = int(self.Ngrids[2])
        self.Nhead = 15
        self.Nheader = self.Nhead + int(self.Natom[0])
        self.vl = ''.join(self.f[self.Nheader:-2])
        self.vn = [float(i) for i in self.vl.split()]
        self.__rho__ = np.array(self.vn).reshape([self.nz,self.ny,self.nx])
        self.rho = self.__rho__.transpose(2,1,0) 
        self.omega = (float(self.a1[0]))**2*np.sqrt(3)/2*float(self.a3[2])/0.52917721**3


class Pot(object):

    def __init__(self,name_):
        self.name = name_
        self.file = open(self.name)
        self.f = self.file.readlines()

        self.a1 = self.f[2].split()
        self.a2 = self.f[3].split()
        self.a3 = self.f[4].split()
        self.Natom = self.f[6].split()
        self.Ngrids = self.f[10+int(self.Natom[0])].split()
        self.nx = int(self.Ngrids[0])
        self.ny = int(self.Ngrids[1])
        self.nz = int(self.Ngrids[2])
        self.Nhead = 15
        self.Nheader = self.Nhead + int(self.Natom[0])
        self.vl = ''.join(self.f[self.Nheader:-2])
        self.vn = [float(i) for i in self.vl.split()]
        self.__pot__ = np.array(self.vn).reshape([self.nz,self.ny,self.nx])
        self.pot = self.__pot__.transpose(2,1,0) 
        self.omega = (float(self.a1[0]))**2*np.sqrt(3)/2*float(self.a3[2])/0.52917721**3


def write_xsf(Rho, xsf_type='rho', filedir='rho.xsf'):
    if xsf_type == 'rho':
        with open(filedir, 'w') as v_file:
            print('Crystal\n PRIMVEC\n', end=' ',file=v_file)
            print("%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n" 
                % (float(Rho.a1[0]), float(Rho.a1[1]), float(Rho.a1[2]),
                    float(Rho.a2[0]), float(Rho.a2[1]), float(Rho.a2[2]), 
                    float(Rho.a3[0]), float(Rho.a3[1]), float(Rho.a3[2])),
                    end='', file=v_file)
            print(' PRIMCOORD\n 0    1\n', end=' ',file=v_file)
            print(' BEGIN_BLOCK_DATAGRID_3D\n 3D_PWSCF\n BEGIN_DATAGRID_3D_UNKNOWN\n', end=' ',file=v_file)
            print("%8d %8d %8d\n" % (Rho.nx, Rho.ny, Rho.nz), end=' ', file=v_file)
            print(' 0.00 0.00 0.00\n', end=' ',file=v_file)
            print("%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n" 
                % (float(Rho.a1[0]), float(Rho.a1[1]), float(Rho.a1[2]),
                    float(Rho.a2[0]), float(Rho.a2[1]), float(Rho.a2[2]), 
                    float(Rho.a3[0]), float(Rho.a3[1]), float(Rho.a3[2])),
                    end='', file=v_file)

            count = 0
            for k in range(Rho.nz):
                for j in range(Rho.ny):
                    for i in range(Rho.nx):
                        print("%16.10f" % Rho.rho[i,j,k], end='', file=v_file)
                        count+=1
                        if count%5==0:
                            print(end='\n', file=v_file)
            if (Rho.nx*Rho.ny*Rho.nz)%5 != 0:
                print(end='\n', file=v_file)
            
            print(' END_DATAGRID_3D\n', end=' ',file=v_file)
            print(' END_BLOCK_DATAGRID_3D', end=' ',file=v_file)

    elif xsf_type == 'pot':
        with open(filedir, 'w') as v_file:
            print('Crystal\n PRIMVEC\n', end=' ',file=v_file)
            print("%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n" 
                % (float(Rho.a1[0]), float(Rho.a1[1]), float(Rho.a1[2]),
                    float(Rho.a2[0]), float(Rho.a2[1]), float(Rho.a2[2]), 
                    float(Rho.a3[0]), float(Rho.a3[1]), float(Rho.a3[2])),
                    end='', file=v_file)
            print(' PRIMCOORD\n 0    1\n', end=' ',file=v_file)
            print(' BEGIN_BLOCK_DATAGRID_3D\n 3D_PWSCF\n BEGIN_DATAGRID_3D_UNKNOWN\n', end=' ',file=v_file)
            print("%8d %8d %8d\n" % (Rho.nx, Rho.ny, Rho.nz), end=' ', file=v_file)
            print(' 0.00 0.00 0.00\n', end=' ',file=v_file)
            print("%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n" 
                % (float(Rho.a1[0]), float(Rho.a1[1]), float(Rho.a1[2]),
                    float(Rho.a2[0]), float(Rho.a2[1]), float(Rho.a2[2]), 
                    float(Rho.a3[0]), float(Rho.a3[1]), float(Rho.a3[2])),
                    end='', file=v_file)

            count = 0
            for k in range(Rho.nz):
                for j in range(Rho.ny):
                    for i in range(Rho.nx):
                        print("%16.10f" % Rho.pot[i,j,k], end='', file=v_file)
                        count+=1
                        if count%5==0:
                            print(end='\n', file=v_file)
            if (Rho.nx*Rho.ny*Rho.nz)%5 != 0:
                print(end='\n', file=v_file)
            
            print(' END_DATAGRID_3D\n', end=' ',file=v_file)
            print(' END_BLOCK_DATAGRID_3D', end=' ',file=v_file)

       