import os
import sys


import read_eps
import io_xsf
import lab
import const
import numpy as np
import math
import pickle
import psutil


if "OMPI_COMM_WORLD_SIZE" in os.environ:
    #print("It seems you're using 'mpirun'.")
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank==0:
        print("Using MPI")
else:
    print("Not using MPI")



def progress_bar(iteration, total, bar_length=50):
    percent = float(iteration) / float(total)
    arrow = '#' * int(round(percent * bar_length) - 1)
    spaces = '.' * (bar_length - len(arrow))
    sys.stdout.write(f'\rProgress: [{arrow}{spaces}] {int(percent * 100)}%')
    sys.stdout.flush()
    
def round_tuple(t, precision):
    return tuple(round(num, precision) if isinstance(num, float) else num for num in t)

def filter_G(ttkw, q_vec, Eps, ecut=5 ):
    #print('Staring filtering G vector and index for ', ttkw.prefix)
    #print(ttkw.lattpara_unit)
    #ttkw.fft_init()
    #ttkw.get_vcoul()
    J_to_Ry = 1 / (2.17987e-18)
    a_0 = 0.52917721067e-10
    Ang= 1e-10
    h_bar = 1.0545718e-34
    m_electron = 9.10938356e-31
    
    if ttkw.ibrav == 4:

        a = float(ttkw.lattpara_unit[0])*Ang
        c = float(ttkw.lattpara_unit[2])*Ang

        b1 = (2 * np.pi / a) * np.array([1, 1/np.sqrt(3), 0])
        b2 = (2 * np.pi / a) * np.array([0, 2/np.sqrt(3), 0])
        b3 = (2 * np.pi / c) * np.array([0, 0, 1])
    elif ttkw.ibrav == 0:
        A = ttkw.A*const.Bohr_R/ttkw.sc
       # print('Lattice Parameter A:',A,' angstrom')
        cell_r = ttkw.cell_r/ttkw.A*A*Ang
        cell_k = 2*np.pi*np.array([np.cross(cell_r[1],cell_r[2])/(np.dot(cell_r[0],np.cross(cell_r[1],cell_r[2]))),
                            np.cross(cell_r[2],cell_r[0])/(np.dot(cell_r[1],np.cross(cell_r[2],cell_r[0]))),
                            np.cross(cell_r[0],cell_r[1])/(np.dot(cell_r[2],np.cross(cell_r[0],cell_r[1])))])
        b1 = cell_k[0]
        b2 = cell_k[1]
        b3 = cell_k[2]
    E_k_limit = ecut
    filtered_G_vectors = []
    Gind_list = []
    maxi = Eps.gind_rho2eps.shape[1]
    for i in range(maxi):
        coeff = Eps.G_ind2vec[i]
        G_vector = (coeff[0]+q_vec[0])*b1 + (coeff[1]+q_vec[1])*b2 + (coeff[2]+q_vec[2])*b3
        E_k = (h_bar**2 * np.linalg.norm(G_vector)**2) / (2 * m_electron)
        E_k_Ry = E_k * J_to_Ry
        if E_k_Ry <= E_k_limit:
            filtered_G_vectors.append(coeff)
            Gind_list.append(i)

  #  with open(dir_name+'filtered_G_vec.pkl','wb') as file:
  #      pickle.dump(filtered_G_vectors, file)

  #  with open(dir_name+'Gind_list.pkl','wb') as file:
  #      pickle.dump(Gind_list, file)

    print('Filtering complete')
    print('Number of filtered G:', len(Gind_list))
    return Gind_list

class PotCorr():
    def __init__(self, cell):
        self.prefix = cell['prefix']
        self.folder = cell['folder']
        #self.Chi0_filename = cell['folder_chi']+'chi0mat.h5'
        #self.Chi1_filename = cell['folder_chi']+'chimat.h5'
        #self.Chi0_doped_filename = cell['folder_chi_doped']+'chi0mat.h5'
        #self.Chi1_doped_filename = cell['folder_chi_doped']+'chi0mat.h5'
        self.is2d = cell['is2d']
        self.ibrav = cell['ibrav']
        self.lattpara_unit = np.array(cell['lattpara_unit'])
        self.sc = np.array(cell['sc'])
        self.fft_g = np.array(cell['fft_g'])
        self.lattpara = self.lattpara_unit*self.sc/const.Bohr_R
        if self.ibrav == 4:
            self.omega = self.lattpara_unit[0]*self.sc[0] * \
                         self.lattpara_unit[1]*self.sc[1] * \
                         self.lattpara_unit[2]*self.sc[2] * \
                         np.sqrt(3)/2 / (const.Bohr_R)**3
        if self.ibrav == 1:
            self.omega = self.lattpara_unit[0]*self.sc[0] * \
                         self.lattpara_unit[1]*self.sc[1] * \
                         self.lattpara_unit[2]*self.sc[2]  \
                         /(const.Bohr_R)**3
        if self.ibrav == 0:
            self.A = cell['A']/const.Bohr_R
            self.cell_r = cell['cell_r']*self.A
            self.cell_k = 2*np.pi*np.array([np.cross(self.cell_r[1],self.cell_r[2])/(np.dot(self.cell_r[0],np.cross(self.cell_r[1],self.cell_r[2]))),
                            np.cross(self.cell_r[2],self.cell_r[0])/(np.dot(self.cell_r[1],np.cross(self.cell_r[2],self.cell_r[0]))),
                            np.cross(self.cell_r[0],self.cell_r[1])/(np.dot(self.cell_r[2],np.cross(self.cell_r[0],self.cell_r[1])))])
            self.omega = np.dot(self.cell_r[0],np.cross(self.cell_r[1],self.cell_r[2]))
            self.omega_k = np.dot(self.cell_k[0],np.cross(self.cell_k[1],self.cell_k[2]))
        try:
            self.sym_arr = cell['sym_arr']
        except(KeyError):
            print('No sym_arr')

    def mpi_display(self):   #only for mpi
        print("Rank {} Prefix:".format(comm.Get_rank()), self.prefix)
    
    
    def fft_init(self):
        print('Begin FFT grids Initialization')

        self.fft_nx=self.fft_g[0]
        self.fft_ny=self.fft_g[1]
        self.fft_nz=self.fft_g[2]

        self.fft_dx=self.lattpara[0]/self.fft_g[0]
        self.fft_dy=self.lattpara[1]/self.fft_g[1]
        self.fft_dz=self.lattpara[2]/self.fft_g[2]

        self.fft_x=np.linspace(0,self.fft_nx*self.fft_dx, self.fft_nx)
        self.fft_y=np.linspace(0,self.fft_ny*self.fft_dy, self.fft_ny)
        self.fft_z=np.linspace(0,self.fft_nz*self.fft_dz, self.fft_nz)

        self.fft_xx, self.fft_yy, self.fft_zz = np.meshgrid(self.fft_x, self.fft_y, self.fft_z, indexing='ij')

        self.fft_kx = 2*np.pi*np.fft.fftfreq(self.fft_nx, d=self.fft_dx)
        self.fft_ky = 2*np.pi*np.fft.fftfreq(self.fft_ny, d=self.fft_dy)
        self.fft_kz = 2*np.pi*np.fft.fftfreq(self.fft_nz, d=self.fft_dz)

        self.fft_kx_tt = self.fft_kx*self.lattpara_unit[0]/(2*np.pi)/const.Bohr_R
        self.fft_ky_tt = self.fft_ky*self.lattpara_unit[1]/(2*np.pi)/const.Bohr_R
        self.fft_kz_tt = self.fft_kz*self.lattpara_unit[2]/(2*np.pi)/const.Bohr_R
        
        self.fft_kx_ttt = np.fft.fftfreq(self.fft_nx, d=1/self.fft_nx)
        self.fft_ky_ttt = np.fft.fftfreq(self.fft_ny, d=1/self.fft_ny)
        self.fft_kz_ttt = np.fft.fftfreq(self.fft_nz, d=1/self.fft_nz)
        
        self.fft_kxx, self.fft_kyy, self.fft_kzz = np.meshgrid(self.fft_kx, self.fft_ky, self.fft_kz, indexing='ij')
        self.fft_kxx_tt, self.fft_kyy_tt, self.fft_kzz_tt = np.meshgrid(self.fft_kx_tt, self.fft_ky_tt, self.fft_kz_tt, indexing='ij')
        self.fft_kxx_ttt, self.fft_kyy_ttt, self.fft_kzz_ttt = np.meshgrid(self.fft_kx_ttt, self.fft_ky_ttt, self.fft_kz_ttt, indexing='ij')
        print('FFT grids:%6d %6d %6d' % (self.fft_nx, self.fft_ny, self.fft_nz))
        
        
    def fft_init_mpi(self):
        
        self.fft_nx=self.fft_g[0]
        self.fft_ny=self.fft_g[1]
        self.fft_nz=self.fft_g[2]
        
        self.fft_dx=self.lattpara[0]/self.fft_g[0]
        self.fft_dy=self.lattpara[1]/self.fft_g[1]
        self.fft_dz=self.lattpara[2]/self.fft_g[2]
        
        self.fft_x=np.linspace(0,self.fft_nx*self.fft_dx, self.fft_nx)
        self.fft_y=np.linspace(0,self.fft_ny*self.fft_dy, self.fft_ny)
        self.fft_z=np.linspace(0,self.fft_nz*self.fft_dz, self.fft_nz)
        
        self.fft_xx, self.fft_yy, self.fft_zz = np.meshgrid(self.fft_x, self.fft_y, self.fft_z, indexing='ij')
        
        self.fft_kx_ttt = np.fft.fftfreq(self.fft_nx, d=1/self.fft_nx)
        self.fft_ky_ttt = np.fft.fftfreq(self.fft_ny, d=1/self.fft_ny)
        self.fft_kz_ttt = np.fft.fftfreq(self.fft_nz, d=1/self.fft_nz)
        
        self.fft_kxx_ttt, self.fft_kyy_ttt, self.fft_kzz_ttt = np.meshgrid(self.fft_kx_ttt, self.fft_ky_ttt, self.fft_kz_ttt, indexing='ij')
        
    def fft_init_pert_mpi(self):
        
        self.fft_nx=self.fft_g[0]
        self.fft_ny=self.fft_g[1]
        self.fft_nz=self.fft_g[2]
        
        self.fft_dx=self.lattpara[0]/self.fft_g[0]
        self.fft_dy=self.lattpara[1]/self.fft_g[1]
        self.fft_dz=self.lattpara[2]/self.fft_g[2]
        
        self.fft_x=np.linspace(0,self.fft_nx*self.fft_dx, self.fft_nx)
        self.fft_y=np.linspace(0,self.fft_ny*self.fft_dy, self.fft_ny)
        self.fft_z=np.linspace(0,self.fft_nz*self.fft_dz, self.fft_nz)
        
        self.fft_xx, self.fft_yy, _ = np.meshgrid(self.fft_x, self.fft_y, self.fft_z, indexing='ij')
        
        self.fft_kx = 2*np.pi*np.fft.fftfreq(self.fft_nx, d=self.fft_dx)
        self.fft_ky = 2*np.pi*np.fft.fftfreq(self.fft_ny, d=self.fft_dy)
        
        del(_)
        
        #self.fft_kx_ttt = np.fft.fftfreq(self.fft_nx, d=1/self.fft_nx)
        #self.fft_ky_ttt = np.fft.fftfreq(self.fft_ny, d=1/self.fft_ny)
        #self.fft_kz_ttt = np.fft.fftfreq(self.fft_nz, d=1/self.fft_nz)
        
        #self.fft_kxx_ttt, self.fft_kyy_ttt, self.fft_kzz_ttt = np.meshgrid(self.fft_kx_ttt, self.fft_ky_ttt, self.fft_kz_ttt, indexing='ij')
        


    def get_vcoul(self):
        print('Constructing Coulomb kernel...')
        if self.ibrav == 1:
            print("ibrav=", self.ibrav)
            self.k_squared = self.fft_kxx**2 + self.fft_kyy**2 + self.fft_kzz**2
            self.v_coul = 8*np.pi/(self.k_squared)
            self.v_coul[0,0,0] = 0
            if self.is2d:
                self.kxy = np.sqrt(self.fft_kxx**2 + self.fft_kyy**2)
                self.v_coul2d = 8*np.pi/(self.k_squared)*(1-np.exp(-0.5*self.kxy*self.lattpara[2])*np.cos(0.5*self.fft_kzz*self.lattpara[2]))
                self.v_coul2d[0,0,0] = 0
        elif self.ibrav == 4:
            print("ibrav=", self.ibrav)
            self.k_squared = self.fft_kxx**2 + (np.sqrt(3)/3*self.fft_kxx+2*np.sqrt(3)/3*self.fft_kyy)**2 + self.fft_kzz**2
            self.v_coul = 8*np.pi/(self.k_squared)
            self.v_coul[0,0,0] = 0
            if self.is2d:
                self.kxy = np.sqrt(self.fft_kxx**2 + (np.sqrt(3)/3*self.fft_kxx+2*np.sqrt(3)/3*self.fft_kyy)**2)
                self.v_coul2d = 8*np.pi/(self.k_squared)*(1-np.exp(-0.5*self.kxy*self.lattpara[2])*np.cos(0.5*self.fft_kzz*self.lattpara[2]))
                self.v_coul2d[0,0,0] = 0
        elif self.ibrav == 0:
            print("ibrav=", self.ibrav)
            self.k_squared = (self.fft_kxx_ttt*self.cell_k[0,0]+self.fft_kyy_ttt*self.cell_k[1,0]+self.fft_kzz_ttt*self.cell_k[2,0])**2 + \
                            (self.fft_kxx_ttt*self.cell_k[0,1]+self.fft_kyy_ttt*self.cell_k[1,1]+self.fft_kzz_ttt*self.cell_k[2,1])**2  + \
                            (self.fft_kxx_ttt*self.cell_k[0,2]+self.fft_kyy_ttt*self.cell_k[1,2]+self.fft_kzz_ttt*self.cell_k[2,2])**2 
            self.v_coul = 8*np.pi/(self.k_squared)
            self.v_coul[0,0,0] = 0
        else:
            print('Lattice with this ibrav has not been implemented..')
        
        
    def get_vcoul_lr(self, alpha=1):
        print('Constructing Short-Range Coulomb kernel...')
        if self.ibrav == 1:
            print("ibrav=", self.ibrav)
            self.k_squared = self.fft_kxx**2 + self.fft_kyy**2 + self.fft_kzz**2
            self.v_coul_lr = 8*np.pi/(self.k_squared)*np.exp(-self.k_squared/(4*alpha**2))
            self.v_coul_lr[0,0,0] = 0
            if self.is2d:
                self.kxy = np.sqrt(self.fft_kxx**2 + self.fft_kyy**2)
                self.v_coul2d = 8*np.pi/(self.k_squared)*(1-np.exp(-0.5*self.kxy*self.lattpara[2])*np.cos(0.5*self.fft_kzz*self.lattpara[2]))
                self.v_coul2d[0,0,0] = 0
        elif self.ibrav == 4:
            print("ibrav=", self.ibrav)
            self.k_squared = self.fft_kxx**2 + (np.sqrt(3)/3*self.fft_kxx+2*np.sqrt(3)/3*self.fft_kyy)**2 + self.fft_kzz**2
            self.v_coul_lr = 8*np.pi/(self.k_squared)*np.exp(-self.k_squared/(4*alpha**2))
            self.v_coul_lr[0,0,0] = 0
            if self.is2d:
                self.kxy = np.sqrt(self.fft_kxx**2 + (np.sqrt(3)/3*self.fft_kxx+2*np.sqrt(3)/3*self.fft_kyy)**2)
                self.v_coul2d = 8*np.pi/(self.k_squared)*(1-np.exp(-0.5*self.kxy*self.lattpara[2])*np.cos(0.5*self.fft_kzz*self.lattpara[2]))
                self.v_coul2d[0,0,0] = 0
        else:
            print('Lattice with this ibrav has not been implemented..')
            
    def get_vcoul_sr(self, alpha=1):
        print('Constructing Short-Range Coulomb kernel...')
        if self.ibrav == 1:
            print("ibrav=", self.ibrav)
            self.k_squared = self.fft_kxx**2 + self.fft_kyy**2 + self.fft_kzz**2
            self.v_coul_sr = 8*np.pi/(self.k_squared)*(1-np.exp(-self.k_squared/(4*alpha**2)))
            self.v_coul_sr[0,0,0] = 2*np.pi/alpha**2
            if self.is2d:
                self.kxy = np.sqrt(self.fft_kxx**2 + self.fft_kyy**2)
                self.v_coul2d = 8*np.pi/(self.k_squared)*(1-np.exp(-0.5*self.kxy*self.lattpara[2])*np.cos(0.5*self.fft_kzz*self.lattpara[2]))
                self.v_coul2d[0,0,0] = 0
        elif self.ibrav == 4:
            print("ibrav=", self.ibrav)
            self.k_squared = self.fft_kxx**2 + (np.sqrt(3)/3*self.fft_kxx+2*np.sqrt(3)/3*self.fft_kyy)**2 + self.fft_kzz**2
            self.v_coul_sr = 8*np.pi/(self.k_squared)*(1-np.exp(-self.k_squared/(4*alpha**2)))
            self.v_coul_sr[0,0,0] = 2*np.pi/alpha**2
            if self.is2d:
                self.kxy = np.sqrt(self.fft_kxx**2 + (np.sqrt(3)/3*self.fft_kxx+2*np.sqrt(3)/3*self.fft_kyy)**2)
                self.v_coul2d = 8*np.pi/(self.k_squared)*(1-np.exp(-0.5*self.kxy*self.lattpara[2])*np.cos(0.5*self.fft_kzz*self.lattpara[2]))
                self.v_coul2d[0,0,0] = 0
        else:
            print('Lattice with this ibrav has not been implemented..')
    
    def get_greenfunc(self):
        print('Constructing Green Function...')
        x = np.linspace(-self.lattpara[0]/2, self.lattpara[0]/2, self.fft_nx)
        y = np.linspace(-self.lattpara[1]/2, self.lattpara[1]/2, self.fft_ny)
        z = np.linspace(-self.lattpara[2]/2, self.lattpara[2]/2, self.fft_nz)
        X, Y, Z = np.meshgrid(x, y, z)
        if self.ibrav == 1:
            print("ibrav=", self.ibrav)
            R = np.sqrt(X**2 + Y**2 + Z**2)
            
        elif self.ibrav == 4:
            print("ibrav=", self.ibrav)
            R = np.sqrt((X-0.5*Y)**2 + (np.sqrt(3)/2*Y)**2 + Z**2)
        epsilon = 1e-10  # to avoid divide by zero
        G = 1.0 / (4.0 * np.pi * (R + epsilon))
        V_r = np.roll(G, shift=-np.array(G.shape)//2, axis=(0, 1, 2))
        V_k = np.fft.fftn(V_r)
        self.greenfunc = G
        self.v_coulgfr = V_r
        self.v_coulgfk = V_k
        self.R = R
    
    def read_chi(self,chi1, chi0=None):
        print('Reading Chi matrix hdf5 files')
        if chi0 != None:
            self.Chi0 = read_eps.Epsmat(chi0)
        self.Chi1 = read_eps.Epsmat(chi1)
    
    def read_epsinv(self, eps1, eps0):
        print('Reading EpsInv matrix hdf5 files')
        if eps1 != None:
            self.Eps1 = read_eps.Epsmat(eps1)
            self.Eps1.read_epsh5()
        if eps0 != None:
            self.Eps0 = read_eps.Epsmat(eps0)
            self.Eps0.read_epsh5()
        
    def read_chi_doped(self):
        print('Read doped Chi matrix hdf5 files')
        self.Chi0 = read_eps.Epsmat(self.Chi0_doped_filename)
        self.Chi1 = read_eps.Epsmat(self.Chi1_doped_filename)

    def read_rho_bare(self, rho_bare_file):
        Rho = io_xsf.Rho(rho_bare_file)
        self.rho_bare_r = Rho.rho#*Rho.omega/Rho.nx/Rho.ny/Rho.nz
        self.rho_bare_k = np.fft.fftn(self.rho_bare_r)

    def read_rho_tot(self, rho_tot_file):
        Rho = io_xsf.Rho(rho_tot_file)
        self.rho_tot_r = Rho.rho#*Rho.omega/Rho.nx/Rho.ny/Rho.nz
        self.rho_tot_k = np.fft.fftn(self.rho_tot_r)  

    def read_pot_bare(self, pot_bare_file):
        Pot = io_xsf.Pot(pot_bare_file)
        self.pot_bare_r = Pot.pot
        self.pot_bare_k = np.fft.ifftn(self.pot_bare_r)

    def read_pot_tot(self, pot_tot_file):
        Pot = io_xsf.Pot(pot_tot_file)
        self.pot_tot_r = Pot.pot
        self.pot_tot_k = np.fft.ifftn(self.pot_tot_r)


    def get_pointc_pot_bare(self, position='center', kernel='3d'):
        print('Setting point-charge potential as bare potential')
        if position == 'center':
            phase_p=np.exp(1j*(self.fft_kxx*self.fft_nx*self.fft_dx*0.5+
                               self.fft_kyy*self.fft_ny*self.fft_dy*0.5+
                               self.fft_kzz*self.fft_nz*self.fft_dz*0.5))
        else:
            phase_p = np.exp(1j*(self.fft_kxx*self.fft_nx*self.fft_dx*position[0]+
                                 self.fft_kyy*self.fft_ny*self.fft_dy*position[1]+
                                 self.fft_kzz*self.fft_nz*self.fft_dz*position[2]))
        
        if kernel == '3d':
            v_coul_ = self.v_coul
        elif kernel == '2d':
            v_coul_ = self.v_coul2d  
        elif kernel == 'gf':
            v_coul_ = self.v_coulgfk  
        
        self.pot_bare_k =  v_coul_*phase_p
        self.pot_bare_r = np.fft.ifftn(self.pot_bare_k)
        
    def get_k_symmetry_map(self):
        k_rdc = self.Eps1.qlist
        k_symmetry_map = {}

        n = 0
        for idx, k in enumerate(k_rdc):
            #print('q-point id:',idx)
            k_tfm_set = set()
            for i, r in enumerate(self.sym_arr, start=0):

                k_transformed = np.dot(r.T, np.array(k))
                G_r = np.mod(k_transformed, 1)-k_transformed
                k_tfm = k_transformed + G_r
                k_tfm_rounded = round_tuple(tuple(k_tfm), 8)
              
             #   print('matrix:', r)
             #   print(k_transformed)
             #   print(k_tfm)
             #   print(G_r)
                # 如果这个变换后的点是新的，添加到集合和映射字典中
                if k_tfm_rounded not in k_tfm_set:
                    k_tfm_set.add(k_tfm_rounded)
                    # 保存原始的k点，变换后的k点和对称操作的编号
                    k_symmetry_map[k_tfm_rounded] = (idx, i, tuple(G_r))

            #print(f"Original k-point index {idx}: {k_tfm_set}")
            #print(f"Number of unique k-points: {len(k_tfm_set)}")
            n += len(k_tfm_set)
        print(f"Total number of unique k-points: {n}")
        
        fbz_q_vec2ind={}
        fbz_q_ind2vec={}
        q_ind_fbz2irrbz={}
        for ind_, vec_ in enumerate(sorted(k_symmetry_map.keys())):
            q_origin_id_, _, _ = k_symmetry_map[round_tuple(vec_, 8)]
            fbz_q_vec2ind[vec_]=ind_
            fbz_q_ind2vec[ind_]=vec_
            q_ind_fbz2irrbz[ind_] = q_origin_id_
        self.fbz_q_vec2ind=fbz_q_vec2ind
        self.fbz_q_ind2vec=fbz_q_ind2vec
        self.q_ind_fbz2irrbz=q_ind_fbz2irrbz

        return k_symmetry_map
    
    def get_k_symmetry_map_from_symdict(self,symdict):
        k_rdc = self.Eps1.qlist
        k_symmetry_map = {}

        n = 0
        for idx, k in enumerate(k_rdc):
            #print('q-point id:',idx)
            k_tfm_set = set()
            for i, r in enumerate(symdict.keys(), start=1):

                k_transformed = np.dot(symdict[r].T, np.array(k))
                G_r = np.mod(k_transformed, 1)-k_transformed
                k_tfm = k_transformed + G_r
                k_tfm_rounded = round_tuple(tuple(k_tfm), 8)
              
             #   print('matrix:', r)
             #   print(k_transformed)
             #   print(k_tfm)
             #   print(G_r)
                # 如果这个变换后的点是新的，添加到集合和映射字典中
                if k_tfm_rounded not in k_tfm_set:
                    k_tfm_set.add(k_tfm_rounded)
                    # 保存原始的k点，变换后的k点和对称操作的编号
                    k_symmetry_map[k_tfm_rounded] = (idx, i, tuple(G_r))

            #print(f"Original k-point index {idx}: {k_tfm_set}")
            #print(f"Number of unique k-points: {len(k_tfm_set)}")
            n += len(k_tfm_set)
        print(f"Total number of unique k-points: {n}")
        
        fbz_q_vec2ind={}
        fbz_q_ind2vec={}
        q_ind_fbz2irrbz={}
        for ind_, vec_ in enumerate(sorted(k_symmetry_map.keys())):
            q_origin_id_, _, _ = k_symmetry_map[round_tuple(vec_, 8)]
            fbz_q_vec2ind[vec_]=ind_
            fbz_q_ind2vec[ind_]=vec_
            q_ind_fbz2irrbz[ind_] = q_origin_id_
        self.fbz_q_vec2ind=fbz_q_vec2ind
        self.fbz_q_ind2vec=fbz_q_ind2vec
        self.q_ind_fbz2irrbz=q_ind_fbz2irrbz

        return k_symmetry_map

    def get_epsym_dict(self, k_symmetry_map, ecut=10):
        epsym_dict = {}
        for k_tfm_rounded, (original_idx, symm_op, G_r_vec) in sorted(k_symmetry_map.items()):
            #print(f"Transformed k-point: {k_tfm_rounded}, Original k-point index: {original_idx}, Symmetry operation: {symm_op}, G-vector: {G_r_vec}")
            print(f"k-point in first BZ: {k_tfm_rounded}")
            #print(f"original k-point index: {original_idx}")
            #print(f"original k-point in irrBZ: {Eps1.qlist[original_idx]}")
            #print(f"Symmetry operation: {symm_op}")
            #print(f"G-vector: {G_r_vec} \n")

            Gind_list = filter_G(self,q_vec=k_tfm_rounded, Eps=self.Eps1, ecut=ecut)
            #print(Gind_list)
            gin_rho2eps_ = {}
            for gin_rho in Gind_list:
                G1 = self.Eps1.G_ind2vec[gin_rho]
                #print(G1)
                symop_ = self.sym_arr[symm_op]
                G_ = np.dot(np.linalg.inv(symop_.T), np.array(G1)+np.array(G_r_vec))
                G_rhoind = self.Eps1.G_vec2ind[tuple(G_)]
                #print(G1)
                #print(G_rhoind)
                #print(G_)
                #print(Eps1.gind_rho2eps[4][G_rhoind])
                #gin_rho2eps_.append(Eps1.gind_rho2eps[original_idx][G_rhoind]-1)
                gin_rho2eps_[gin_rho]=self.Eps1.gind_rho2eps[original_idx][G_rhoind]-1
            #print(Gind_list)
            #print(gin_rho2eps_)
            epsym_dict[k_tfm_rounded] = gin_rho2eps_
        return epsym_dict
    
    def get_epsym_dict_from_symdict(self, k_symmetry_map,symdict, ecut=10):
        epsym_dict = {}
        for k_tfm_rounded, (original_idx, symm_op, G_r_vec) in sorted(k_symmetry_map.items()):
            #print(f"Transformed k-point: {k_tfm_rounded}, Original k-point index: {original_idx}, Symmetry operation: {symm_op}, G-vector: {G_r_vec}")
            print(f"k-point in first BZ: {k_tfm_rounded}")
            #print(f"original k-point index: {original_idx}")
            #print(f"original k-point in irrBZ: {Eps1.qlist[original_idx]}")
            #print(f"Symmetry operation: {symm_op}")
            #print(f"G-vector: {G_r_vec} \n")

            Gind_list = filter_G(self,q_vec=k_tfm_rounded, Eps=self.Eps1, ecut=ecut)
            #print(Gind_list)
            gin_rho2eps_ = {}
            for gin_rho in Gind_list:
                G1 = self.Eps1.G_ind2vec[gin_rho]
                #print(G1)
                if symm_op < 10:
                    symop_ = symdict['r0'+str(symm_op)]
                else:
                    symop_ = symdict['r'+str(symm_op)]
                G_ = np.dot(np.linalg.inv(symop_.T), np.array(G1)+np.array(G_r_vec))
                G_rhoind = self.Eps1.G_vec2ind[tuple(G_)]
                #print(G1)
                #print(G_rhoind)
                #print(G_)
                #print(Eps1.gind_rho2eps[4][G_rhoind])
                #gin_rho2eps_.append(Eps1.gind_rho2eps[original_idx][G_rhoind]-1)
                gin_rho2eps_[gin_rho]=self.Eps1.gind_rho2eps[original_idx][G_rhoind]-1
            #print(Gind_list)
            #print(gin_rho2eps_)
            epsym_dict[k_tfm_rounded] = gin_rho2eps_
        return epsym_dict
    
    def read_epsym_dict(self, epsym_file):
        with open(epsym_file, 'rb') as f:
            epsym_dict = pickle.load(f)
        return epsym_dict
    
    def read_epsmat_dict(self,epsmat_eps2rho_dict_file, epsmat_eps2eps_irrbz_dict_file, epsmat_dict_file):
        with open(epsmat_eps2rho_dict_file, 'rb') as f:
            epsmat_eps2rho_dict = pickle.load(f)
        with open(epsmat_eps2eps_irrbz_dict_file, 'rb') as f1:
            epsmat_eps2eps_irrbz_dict = pickle.load(f1)
        with open(epsmat_dict_file, 'rb') as f2:
            epsmat_dict = pickle.load(f2)
        
        return epsmat_eps2rho_dict, epsmat_eps2eps_irrbz_dict, epsmat_dict
    
    def read_wcoul_mat_dict(self,wcoul_mat_dict_file):
        with open(wcoul_mat_dict_file, 'rb') as f:
            wcoul_mat_dict = pickle.load(f)
        return wcoul_mat_dict
        
        return epsmat_eps2rho_dict, epsmat_eps2eps_irrbz_dict, epsmat_dict
    
    def gen_mat_q(self, q_vec, epsym_dict):
        nmtx_q = len(epsym_dict[q_vec])
        epsmat_q = np.zeros((nmtx_q,nmtx_q), dtype=complex)
        epsmat_eps2rho_q = np.array(sorted(epsym_dict[q_vec].keys()))
        q_id= self.fbz_q_vec2ind[q_vec]
        q_id_irrbz = self.q_ind_fbz2irrbz[q_id]
        epsmat_eps2eps_irrbz_q = []

        for i in range(nmtx_q):
            eps_id = i
            rho_id = epsmat_eps2rho_q[i]
            eps_id_irrbz = epsym_dict[q_vec][rho_id] 
            #print(eps_id, rho_id, eps_id_irrbz)
            epsmat_eps2eps_irrbz_q.append(eps_id_irrbz)
        epsmat_eps2eps_irrbz_q = np.array(epsmat_eps2eps_irrbz_q)

        for i in range(nmtx_q):
            for j in range(nmtx_q):
                rho_id_i = epsmat_eps2rho_q[i]
                rho_id_j = epsmat_eps2rho_q[j]
                eps_id_irrbz_i = epsym_dict[q_vec][rho_id_i]
                eps_id_irrbz_j = epsym_dict[q_vec][rho_id_j]

                epsmat_q[i,j] = self.Eps1.mat[q_id_irrbz, 0,0, eps_id_irrbz_j,eps_id_irrbz_i,0]+self.Eps1.mat[q_id_irrbz, 0,0, eps_id_irrbz_j,eps_id_irrbz_i,1]*1j

        return  epsmat_eps2rho_q, epsmat_eps2eps_irrbz_q, epsmat_q
    
    def gen_mat_dict(self, k_symmetry_map, epsym_dict):
        nq = len(sorted(k_symmetry_map.keys()))
        epsmat_dict = {}
        epsmat_eps2rho_dict = {}
        epsmat_eps2eps_irrbz_dict = {}
        for q_vec in (sorted(k_symmetry_map.keys())):
            print(f'Starting {q_vec} epsmat collection')
            
            if q_vec == (0.0,0.0,0.0):
                epsmat_dict[q_vec]=(self.Eps0.mat[0,0,0,:,:,0]+self.Eps0.mat[0,0,0,:,:,1]*1j).T
                #epsmat_dict[q_vec]=self.Eps0.mat[0,0,0,:,:,0]+self.Eps0.mat[0,0,0,:,:,1]*1j
                epsmat_eps2rho_dict[q_vec]=self.Eps0.gind_eps2rho[0,:self.Eps0.nmtx[0]]-1
                epsmat_eps2eps_irrbz_dict[q_vec]=np.arange(self.Eps0.nmtx[0])
            else:    
                epsmat_eps2rho_q, epsmat_eps2eps_irrbz_q, epsmat_q = self.gen_mat_q(q_vec, epsym_dict)
                epsmat_eps2rho_dict[q_vec]=epsmat_eps2rho_q
                epsmat_eps2eps_irrbz_dict[q_vec]=epsmat_eps2eps_irrbz_q
                epsmat_dict[q_vec] = epsmat_q
            
        
        return epsmat_eps2rho_dict, epsmat_eps2eps_irrbz_dict, epsmat_dict
    
    def vcoul2d0modify(self):
        q0_ = self.Eps0.qpts[0]
        q0xy_abs_ = np.sqrt((q0_[0])**2+ (np.sqrt(3)/3*q0_[0]+2*np.sqrt(3)/3*q0_[1])**2)*2*np.pi/self.lattpara_unit[0]*const.Bohr_R
        v0_coul2d = 8*np.pi/(q0xy_abs_)**2*(1-np.exp(-0.5*q0xy_abs_*self.lattpara[2]))
        print(v0_coul2d)
        self.v_coul2d[0,0,0]=v0_coul2d
        
    def vcoul0modify(self): # fixme
        q0_ = self.Eps0.qpts[0]
        q0xy_abs_ = np.sqrt((q0_[0])**2+ (np.sqrt(3)/3*q0_[0]+2*np.sqrt(3)/3*q0_[1])**2)*2*np.pi/self.lattpara_unit[0]*const.Bohr_R
        v0_coul = 8*np.pi/(q0xy_abs_)**2
        print(v0_coul)
        self.v_coul[0,0,0]=v0_coul

    def gen_wcoul_mat_q(self, q_vec, epsmat_eps2rho_dict, epsmat_eps2eps_irrbz_dict, epsmat_dict):
        wcoul_mat_q=np.zeros_like(epsmat_dict[q_vec])
        for eps_id, rho_id in enumerate(epsmat_eps2rho_dict[q_vec]):  
            g_vec =  self.Eps1.G_ind2vec[rho_id]
            #print(eps_id, rho_id, g_vec)
            x_ind = np.where(np.round(self.fft_kx_tt, 8) == np.round((q_vec[0]+g_vec[0]),8))[0][0]
            y_ind = np.where(np.round(self.fft_ky_tt, 8) == np.round((q_vec[1]+g_vec[1]),8))[0][0]
            z_ind = np.where(np.round(self.fft_kz_tt, 8) == np.round((q_vec[2]+g_vec[2]),8))[0][0]
            v_coul_tmp = self.v_coul2d[x_ind, y_ind, z_ind]
            wcoul_mat_q[:, eps_id]=epsmat_dict[q_vec][:,eps_id]*v_coul_tmp
        return wcoul_mat_q
    
    def gen_wcoul_mat_q_3d(self, q_vec, epsmat_eps2rho_dict, epsmat_eps2eps_irrbz_dict, epsmat_dict):
        wcoul_mat_q=np.zeros_like(epsmat_dict[q_vec])
        for eps_id, rho_id in enumerate(epsmat_eps2rho_dict[q_vec]):  
            g_vec =  self.Eps1.G_ind2vec[rho_id]
            #print(eps_id, rho_id, g_vec)
            x_ind = np.where(np.round(self.fft_kx_tt, 8) == np.round((q_vec[0]+g_vec[0]),8))[0][0]
            y_ind = np.where(np.round(self.fft_ky_tt, 8) == np.round((q_vec[1]+g_vec[1]),8))[0][0]
            z_ind = np.where(np.round(self.fft_kz_tt, 8) == np.round((q_vec[2]+g_vec[2]),8))[0][0]
            v_coul_tmp = self.v_coul[x_ind, y_ind, z_ind]
            wcoul_mat_q[:, eps_id]=epsmat_dict[q_vec][:,eps_id]*v_coul_tmp
        return wcoul_mat_q
    
    def gen_wcoul_mat_dict(self, k_symmetry_map, epsmat_eps2rho_dict, epsmat_eps2eps_irrbz_dict, epsmat_dict):
        nq = len(sorted(k_symmetry_map.keys()))
        wcoul_mat_dict = {}
        for q_vec in (sorted(k_symmetry_map.keys())):
            print(f'Starting {q_vec} woul_mat collection')
            wcoul_mat_q = self.gen_wcoul_mat_q(q_vec, epsmat_eps2rho_dict, epsmat_eps2eps_irrbz_dict, epsmat_dict)
            wcoul_mat_dict[q_vec]=wcoul_mat_q
        return wcoul_mat_dict
    
    def gen_wcoul_mat_dict_3d(self, k_symmetry_map, epsmat_eps2rho_dict, epsmat_eps2eps_irrbz_dict, epsmat_dict):
        nq = len(sorted(k_symmetry_map.keys()))
        wcoul_mat_dict = {}
        for q_vec in (sorted(k_symmetry_map.keys())):
            print(f'Starting {q_vec} woul_mat collection')
            wcoul_mat_q = self.gen_wcoul_mat_q_3d(q_vec, epsmat_eps2rho_dict, epsmat_eps2eps_irrbz_dict, epsmat_dict)
            wcoul_mat_dict[q_vec]=wcoul_mat_q
        return wcoul_mat_dict
    
    def gen_phi_G_q(self, rho_ext_k, q_vec, epsmat_eps2rho_dict, epsmat_eps2eps_irrbz_dict, wcoul_mat_dict):
        phi_G_q = np.zeros_like(wcoul_mat_dict[q_vec][:,0])
        for eps_id, rho_id in enumerate(epsmat_eps2rho_dict[q_vec]):  
            g_vec =  self.Eps1.G_ind2vec[rho_id]
            #print(eps_id, rho_id, g_vec)
            x_ind = np.where(np.round(self.fft_kx_tt, 8) == np.round((q_vec[0]+g_vec[0]),8))[0][0]
            y_ind = np.where(np.round(self.fft_ky_tt, 8) == np.round((q_vec[1]+g_vec[1]),8))[0][0]
            z_ind = np.where(np.round(self.fft_kz_tt, 8) == np.round((q_vec[2]+g_vec[2]),8))[0][0]
            rho_tmp = rho_ext_k[x_ind, y_ind, z_ind]
            phi_G_q += wcoul_mat_dict[q_vec][:,eps_id]*rho_tmp
        return phi_G_q
    
    
    def gen_phi_G_dict(self, rho_ext_k, k_symmetry_map, epsmat_eps2rho_dict, epsmat_eps2eps_irrbz_dict, wcoul_mat_dict):
        nq = len(sorted(k_symmetry_map.keys()))
        phi_G_dict = {}
        for q_vec in (sorted(k_symmetry_map.keys())):
            print(f'Starting {q_vec} phi_G collection')
            phi_G_q = self.gen_phi_G_q(rho_ext_k, q_vec, epsmat_eps2rho_dict, epsmat_eps2eps_irrbz_dict, wcoul_mat_dict)
            phi_G_dict[q_vec]=phi_G_q
        return phi_G_dict
    
    def map_phi_G(self, phi_G_dict, epsmat_eps2rho_dict):
        phi_k = np.zeros((self.fft_nx,self.fft_ny,self.fft_nz), dtype=complex)
        for q_vec in phi_G_dict.keys():
            print(f'mapping {q_vec} to phi_k')
            for eps_id, rho_id in enumerate(epsmat_eps2rho_dict[q_vec]):  
                g_vec =  self.Eps1.G_ind2vec[rho_id]
                #print(eps_id, rho_id, g_vec)
                x_ind = np.where(np.round(self.fft_kx_tt, 8) == np.round((q_vec[0]+g_vec[0]),8))[0][0]
                y_ind = np.where(np.round(self.fft_ky_tt, 8) == np.round((q_vec[1]+g_vec[1]),8))[0][0]
                z_ind = np.where(np.round(self.fft_kz_tt, 8) == np.round((q_vec[2]+g_vec[2]),8))[0][0]
                phi_k[x_ind, y_ind, z_ind]=phi_G_dict[q_vec][eps_id]
        return phi_k

    
    def epsmat_init(self):
        print('\n')
        print("Initalizing epsmat")
        try:
            Eps1=self.Chi1
            Eps0=self.Chi0
        except(AttributeError):
            Eps1=self.Eps1
            Eps0=self.Eps0
        nx_, ny_, nz_ = self.fft_nx, self.fft_ny, self.fft_nz
        g_tuple_tt = np.empty((nx_,ny_,nz_), dtype = tuple)
        g_tuple_floor_tt = np.empty((nx_,ny_,nz_), dtype = tuple)
        g_tuple_q_tt = np.empty((nx_,ny_,nz_), dtype = tuple)
        g_rhoind_tt = np.empty((nx_,ny_,nz_), dtype = int)
        q_ind_tt = np.empty((nx_,ny_,nz_), dtype = int)

        gg_eps = int(len(Eps1.G_vec2ind))
        for i in range(nx_):
            for j in range(ny_):
                for k in range(nz_):
                    g_tuple_tt[i,j,k] = (float(format(self.fft_kxx_tt[i,j,k],'.12f')), 
                                         float(format(self.fft_kyy_tt[i,j,k],'.12f')), 
                                         float(format(self.fft_kzz_tt[i,j,k],'.12f')))
                    g_tuple_floor_tt[i,j,k] = (math.floor(g_tuple_tt[i,j,k][0]), 
                                               math.floor(g_tuple_tt[i,j,k][1]), 
                                               math.floor(g_tuple_tt[i,j,k][2]))
                    g_tuple_q_tt[i,j,k] = (float(format(g_tuple_tt[i,j,k][0]-g_tuple_floor_tt[i,j,k][0],'.12f')),
                                           float(format(g_tuple_tt[i,j,k][1]-g_tuple_floor_tt[i,j,k][1],'.12f')),
                                           float(format(g_tuple_tt[i,j,k][2]-g_tuple_floor_tt[i,j,k][2],'.12f')))
                    if g_tuple_floor_tt[i,j,k] in Eps1.G_vec2ind:
                        g_rhoind_tt[i,j,k] = Eps1.G_vec2ind[g_tuple_floor_tt[i,j,k]]
                    else:
                        g_rhoind_tt[i,j,k] = gg_eps

                    if g_tuple_q_tt[i,j,k] in Eps1.q_vec2ind:
                        q_ind_tt[i,j,k] = Eps1.q_vec2ind[g_tuple_q_tt[i,j,k]]
                    else:
                        print("Error: No targeted q-points")
                        print(g_tuple_tt[i,j,k])
                        print(g_tuple_floor_tt[i,j,k])
                        print(g_tuple_q_tt[i,j,k])
                        return
            progress_bar(i + 1, nx_)


        Eps1.gind_rho2eps = np.hstack((Eps1.gind_rho2eps,200000*np.ones([len(Eps1.qpts),1],dtype = int)))
        Eps0.gind_rho2eps = np.hstack((Eps0.gind_rho2eps,200000*np.ones([len(Eps0.qpts),1],dtype = int)))
        g_epsind_tt = Eps1.gind_rho2eps[q_ind_tt, g_rhoind_tt]
        g_epsind0_tt = Eps0.gind_rho2eps[0, g_rhoind_tt]

        lmax = len(g_epsind_tt[g_epsind_tt==200000])
        g_epsind_tt[g_epsind_tt==200000] = 200000+np.arange(lmax)
        g_epsind0_tt[g_epsind0_tt==200000] = 200000+np.arange(lmax)

        self.q_ind_tt = q_ind_tt
        self.lmax = lmax
        self.g_epsind_tt = g_epsind_tt
        self.g_epsind0_tt = g_epsind0_tt

    def epsmat_irrbz_init(self,k_symmetry_map, epsym_dict):
        
        print('\n')
        print("Initalizing epsmat irrbz")
        try:
            Eps1=self.Chi1
            Eps0=self.Chi0
        except(AttributeError):
            Eps1=self.Eps1
            Eps0=self.Eps0
        nx_, ny_, nz_ = self.fft_nx, self.fft_ny, self.fft_nz
        g_tuple_tt = np.empty((nx_,ny_,nz_), dtype = tuple)
        g_tuple_floor_tt = np.empty((nx_,ny_,nz_), dtype = tuple)
        g_tuple_q_tt = np.empty((nx_,ny_,nz_), dtype = tuple)
        g_rhoind_tt = np.empty((nx_,ny_,nz_), dtype = int)
        q_ind_tt = np.empty((nx_,ny_,nz_), dtype = int)
        g_epsind_tt = np.empty((nx_,ny_,nz_), dtype = int)
        q_ind_irrbz_tt = np.empty((nx_,ny_,nz_), dtype = int)
        gg_eps = int(len(Eps1.G_vec2ind))
        for i in range(nx_):
            for j in range(ny_):
                for k in range(nz_):
                    g_tuple_tt[i,j,k] = (float(format(self.fft_kxx_tt[i,j,k],'.12f')), 
                                         float(format(self.fft_kyy_tt[i,j,k],'.12f')), 
                                         float(format(self.fft_kzz_tt[i,j,k],'.12f')))
                    g_tuple_floor_tt[i,j,k] = (math.floor(g_tuple_tt[i,j,k][0]), 
                                               math.floor(g_tuple_tt[i,j,k][1]), 
                                               math.floor(g_tuple_tt[i,j,k][2]))
                    g_tuple_q_tt[i,j,k] = (float(format(g_tuple_tt[i,j,k][0]-g_tuple_floor_tt[i,j,k][0],'.12f')),
                                           float(format(g_tuple_tt[i,j,k][1]-g_tuple_floor_tt[i,j,k][1],'.12f')),
                                           float(format(g_tuple_tt[i,j,k][2]-g_tuple_floor_tt[i,j,k][2],'.12f')))
                     
                    
                    
                    if g_tuple_floor_tt[i,j,k] in Eps1.G_vec2ind: 
                        g_rhoind_tt[i,j,k] = Eps1.G_vec2ind[g_tuple_floor_tt[i,j,k]]
                    else:
                        g_rhoind_tt[i,j,k] = gg_eps

                    if round_tuple(g_tuple_q_tt[i,j,k], 8) in k_symmetry_map.keys():
                        q_ind_tt[i,j,k] = self.fbz_q_vec2ind[round_tuple(g_tuple_q_tt[i,j,k], 8)]
                        q_ind_irrbz_tt[i,j,k] =  self.q_ind_fbz2irrbz[q_ind_tt[i,j,k]]
                    else:
                        print("Error: No targeted q-points")
                        print(g_tuple_tt[i,j,k])
                        print(g_tuple_floor_tt[i,j,k])
                        print(g_tuple_q_tt[i,j,k])
                        return
                   # q_origin_id_, _, _ = k_symmetry_map[round_tuple(g_tuple_q_tt[i,j,k], 8)]
                   # if g_rhoind_tt[i,j,k] in epsym_dict[round_tuple(g_tuple_q_tt[i,j,k], 8)].keys():
                   #     g_epsind_tt[i,j,k] = epsym_dict[round_tuple(g_tuple_q_tt[i,j,k], 8)][g_rhoind_tt[i,j,k]]
                   # else:
                   #     g_epsind_tt[i,j,k] = 200000
                    #g_epsind_tt[i,j,k] = Eps1.gind_rho2eps[q_origin_id_, ]
            progress_bar(i + 1, nx_)


        Eps1.gind_rho2eps = np.hstack((Eps1.gind_rho2eps,200000*np.ones([len(Eps1.qpts),1],dtype = int)))
        Eps0.gind_rho2eps = np.hstack((Eps0.gind_rho2eps,200000*np.ones([len(Eps0.qpts),1],dtype = int)))
        g_epsind_tt = Eps1.gind_rho2eps[q_ind_irrbz_tt, g_rhoind_tt]
        g_epsind0_tt = Eps0.gind_rho2eps[0, g_rhoind_tt]

        lmax = len(g_epsind_tt[g_epsind_tt==200000])
        g_epsind_tt[g_epsind_tt==200000] = 200000+np.arange(lmax)
        lmax = len(g_epsind_tt[g_epsind0_tt==200000])
        g_epsind0_tt[g_epsind0_tt==200000] = 200000+np.arange(lmax)

        self.q_ind_tt = q_ind_tt
        #self.lmax = lmax
        self.g_epsind_tt = g_epsind_tt
        self.g_epsind0_tt = g_epsind0_tt 
        self.q_ind_irrbz_tt = q_ind_irrbz_tt
       
        

    def epsmat_init_interp(self,G_vec2ind_dict, G_ind2vec_dict):
        with open(G_vec2ind_dict, "rb") as file_:
            G_vec2ind = pickle.load(file_)
        with open(G_ind2vec_dict, "rb") as file_:
            G_ind2vec = pickle.load(file_)
        nx_, ny_, nz_ = self.fft_nx, self.fft_ny, self.fft_nz
        g_tuple_tt = np.empty((nx_,ny_,nz_), dtype = tuple)
        g_tuple_floor_tt = np.empty((nx_,ny_,nz_), dtype = tuple)
        g_tuple_q_tt = np.empty((nx_,ny_,nz_), dtype = tuple)
        g_rhoind_tt = np.empty((nx_,ny_,nz_), dtype = int)
        q_ind_tt = np.empty((nx_,ny_,nz_), dtype = int)
        gg_eps = int(len(G_vec2ind))

        for i in range(nx_):
            for j in range(ny_):
                for k in range(nz_):
                    g_tuple_tt[i,j,k] = (float(format(self.fft_kxx_tt[i,j,k],'.8f')), 
                                         float(format(self.fft_kyy_tt[i,j,k],'.8f')), 
                                         float(format(self.fft_kzz_tt[i,j,k],'.8f')))
                    g_tuple_floor_tt[i,j,k] = (math.floor(g_tuple_tt[i,j,k][0]), 
                                               math.floor(g_tuple_tt[i,j,k][1]), 
                                               math.floor(g_tuple_tt[i,j,k][2]))
                    g_tuple_q_tt[i,j,k] = (float(format(g_tuple_tt[i,j,k][0]-g_tuple_floor_tt[i,j,k][0],'.8f')),
                                           float(format(g_tuple_tt[i,j,k][1]-g_tuple_floor_tt[i,j,k][1],'.8f')),
                                           float(format(g_tuple_tt[i,j,k][2]-g_tuple_floor_tt[i,j,k][2],'.8f')))
                    if g_tuple_floor_tt[i,j,k] in G_vec2ind:
                        g_rhoind_tt[i,j,k] = G_vec2ind[g_tuple_floor_tt[i,j,k]]
                    else:
                        g_rhoind_tt[i,j,k] = gg_eps

        
            print(i,"/",nx_)
        self.g_tuple_floor_tt = g_tuple_floor_tt
        self.g_rhoind_tt = g_rhoind_tt
        self.g_tuple_q_tt = g_tuple_q_tt
        self.G_ind2vec = G_ind2vec
        self.G_vec2ind = G_vec2ind


    def get_chimat(self, G_ind_cut = 1300):
        Eps1=self.Chi1
        chi_mat = np.zeros((len(Eps1.qpts),G_ind_cut, G_ind_cut), dtype=complex)
        #epsinv_mat = np.zeros((len(Eps1.qpts),G_ind_cut, G_ind_cut), dtype=complex)
        for  i in range(len(Eps1.qpts)):
            print("Constructing chi_mat of number %d/%d q-point" % (i+1, len(Eps1.qpts)))
            #G_ind_list_tmp = self.g_epsind_tt[self.q_ind_tt == i]
            #print(G_ind_list_tmp)

            for j in range(G_ind_cut):
               
                #print(v_coul_tmp)
                if i == 0:
                    chi_mat[i,j,:] = self.Chi0.mat[0,0,0,j,:G_ind_cut,0]+self.Chi0.mat[0,0,0,j,:G_ind_cut,1]*1.0j
                else:
                    chi_mat[i,j,:] = self.Chi1.mat[i,0,0,j,:G_ind_cut,0]+self.Chi1.mat[i,0,0,j,:G_ind_cut,1]*1.0j
               # for jp in range(G_ind_cut):
               #     if j == jp:
               #         chi_mat[i,j,jp] = chi_mat[i,j,jp]+1.0

        self.chi_mat = chi_mat

    def get_chimat_interp(self, q_ind, G_ind_cut):
        chi_mat = np.zeros((len(q_ind), G_ind_cut, G_ind_cut), dtype=complex)
        
        

    

    def get_epsmat_from_chi(self, G_ind_cut = 1300, kernel='3d', if_dope=False):
        print('\n')
        print("Getting epsmat")
        Eps1=self.Chi1
        eps_mat = np.zeros((len(Eps1.qpts),G_ind_cut, G_ind_cut), dtype=complex)
        #epsinv_mat = np.zeros((len(Eps1.qpts),G_ind_cut, G_ind_cut), dtype=complex)
        for  i in range(len(Eps1.qpts)):
            progress_bar(i + 1, len(Eps1.qpts))
            #print("Constructing eps_mat of number %d/%d q-point" % (i+1, len(Eps1.qpts)))
            G_ind_list_tmp = self.g_epsind_tt[self.q_ind_tt == i]
            #print(G_ind_list_tmp)
            if kernel == '3d':
                v_coul_list_tmp = self.v_coul[self.q_ind_tt == i]
                #print(v_coul_list_tmp)
            elif kernel == '2d':
                v_coul_list_tmp = self.v_coul2d[self.q_ind_tt == i]
            elif kernel == 'gf':
                v_coul_list_tmp = self.v_coulgfk[self.q_ind_tt == i]
                
            else:
                print('Undefined kernel type!')
                return

            for j in range(G_ind_cut):
                v_coul_tmp = v_coul_list_tmp[G_ind_list_tmp == j+1]
                #print(v_coul_tmp)
                if i == 0:
          #          eps_mat[i,j,:] = -1.0*v_coul_tmp[0]*(self.Chi0.mat[0,0,0,j,:G_ind_cut,0]+
          #                                               self.Chi0.mat[0,0,0,j,:G_ind_cut,1]*1.0j)
                    eps_mat[i,j,:] = -1.0*v_coul_tmp[0]*(self.Chi0.mat[0,0,0,:G_ind_cut,j,0]+
                                                         self.Chi0.mat[0,0,0,:G_ind_cut,j,1]*1.0j)
                else:
          #          eps_mat[i,j,:] = -1.0*v_coul_tmp[0]*(self.Chi1.mat[i,0,0,j,:G_ind_cut,0]+
          #                                               self.Chi1.mat[i,0,0,j,:G_ind_cut,1]*1.0j)
                    eps_mat[i,j,:] = -1.0*v_coul_tmp[0]*(self.Chi1.mat[i,0,0,:G_ind_cut,j,0]+
                                                         self.Chi1.mat[i,0,0,:G_ind_cut,j,1]*1.0j)
                for jp in range(G_ind_cut):
                    if j == jp:
                        eps_mat[i,j,jp] = eps_mat[i,j,jp]+1.0
                if if_dope:
                    eps_mat[0,0,0] = 0
        self.eps_mat = eps_mat

    

    def epsmat_inv(self, G_ind_cut = 1300):
        print('\n')
        print("Inversing epsmat")
        Eps1 = self.Chi1
        epsinv_mat = np.zeros((len(Eps1.qpts), G_ind_cut, G_ind_cut), dtype=complex)
        for i in range(len(Eps1.qpts)):
            progress_bar(i + 1, len(Eps1.qpts))
            #print("Inversing eps_mat of number %d/%d q-point" % (i+1, len(Eps1.qpts)))
            epsinv_mat[i,:,:] = np.linalg.inv(self.eps_mat[i,:,:])
        self.epsinv_mat = epsinv_mat
        
    def get_epsinv_mat(self, G_ind_cut=1300):
        print('\n')
        print("Getting EpsInv matrix")
        Eps1=self.Eps1
        Eps0=self.Eps0
        epsinv_mat = np.zeros((len(Eps1.qpts),G_ind_cut, G_ind_cut), dtype=complex)
        #epsinv_mat = np.zeros((len(Eps1.qpts),G_ind_cut, G_ind_cut), dtype=complex)
        for i in range(len(Eps1.qpts)):
            #print("Constructing chi_mat of number %d/%d q-point" % (i+1, len(Eps1.qpts)))
            #G_ind_list_tmp = self.g_epsind_tt[self.q_ind_tt == i]
            #print(G_ind_list_tmp)

            for j in range(G_ind_cut):   
                #print(v_coul_tmp)
                if i == 0:
                    epsinv_mat[i,j,:] = self.Eps0.mat[0,0,0,:G_ind_cut,j,0]+self.Eps0.mat[0,0,0,:G_ind_cut,j,1]*1.0j
                else:
                    epsinv_mat[i,j,:] = self.Eps1.mat[i,0,0,:G_ind_cut,j,0]+self.Eps1.mat[i,0,0,:G_ind_cut,j,1]*1.0j
               # for jp in range(G_ind_cut):
               #     if j == jp:
               #         chi_mat[i,j,jp] = chi_mat[i,j,jp]+1.0

        self.epsinv_mat = epsinv_mat
        
    def get_epsinv_mat_irrbz(self, k_symmetry_map, epsym_dict, G_ind_cut=1000):
        print('\n')
        print("Getting EpsInv matrix")
        Eps1=self.Eps1
        Eps0=self.Eps0
        epsinv_mat = np.zeros((len(k_symmetry_map.keys()),G_ind_cut, G_ind_cut), dtype=complex)
        #epsinv_mat = np.zeros((len(Eps1.qpts),G_ind_cut, G_ind_cut), dtype=complex)
        for i in range(len(k_symmetry_map.keys())):
            #print("Constructing chi_mat of number %d/%d q-point" % (i+1, len(Eps1.qpts)))
            #G_ind_list_tmp = self.g_epsind_tt[self.q_ind_tt == i]
            #print(G_ind_list_tmp)
            origin_i = self.q_ind_fbz2irrbz[i]
            for j in range(G_ind_cut):   
                #print(v_coul_tmp)
                if i == 0:
                    epsinv_mat[i,j,:] = self.Eps0.mat[0,0,0,:G_ind_cut,j,0]+self.Eps0.mat[0,0,0,:G_ind_cut,j,1]*1.0j
                else:
                    epsinv_mat[i,j,:] = self.Eps1.mat[origin_i,0,0,:G_ind_cut,j,0]+self.Eps1.mat[origin_i,0,0,:G_ind_cut,j,1]*1.0j
               # for jp in range(G_ind_cut):
               #     if j == jp:
               #         chi_mat[i,j,jp] = chi_mat[i,j,jp]+1.0

        self.epsinv_mat = epsinv_mat

        
    def get_wcoul_irrbz(self,k_symmetry_map,epsym_dict,G_ind_cut=1000, kernel='2d'):
        print('\n')
        print("Calculating screened Coulomb interaction W(q,G,G')")
        try:
            Eps1 = self.Chi1
            Eps0 = self.Chi0
        except(AttributeError):
            Eps1=self.Eps1
            Eps0=self.Eps0
        wcoul = np.zeros((len(k_symmetry_map.keys()), G_ind_cut, G_ind_cut), dtype=complex)
        for i  in range(len(k_symmetry_map.keys())):
           
            progress_bar(i + 1, len(k_symmetry_map.keys()))
            
            q_vec = self.fbz_q_ind2vec[i]
            q_orgid = self.q_ind_fbz2irrbz[i]
            
            G_ind_list_tmp = self.g_epsind_tt[self.q_ind_tt == i]
            if kernel == '3d':
                v_coul_list_tmp = self.v_coul[self.q_ind_tt == i]
                #print(v_coul_list_tmp)
            elif kernel == '2d':
                v_coul_list_tmp = self.v_coul2d[self.q_ind_tt == i]
            elif kernel == 'gf':
                v_coul_list_tmp = self.v_coulgfk[self.q_ind_tt == i]
            else:
                print('Undefined kernel type!')
                return
            
            for j in range(G_ind_cut):
                v_coul_tmp = v_coul_list_tmp[G_ind_list_tmp == j+1]
              
                
                
                if i==0 and j ==0 :
                    
                    print(self.Eps0.gind_eps2rho[i,j])
                    q0_ = self.Eps0.qpts[0]
                    qxy_abs_ = np.sqrt((q0_[0])**2+ (np.sqrt(3)/3*q0_[0]+2*np.sqrt(3)/3*q0_[1])**2)*2*np.pi/self.lattpara_unit[0]*const.Bohr_R
                    
                    #print(v_coul_tmp)
                    #gp_=[0,0,0]
                    #gp_vec_ = tuple(gp_)
                    #gind_rho_ = self.Eps0.G_vec2ind[gp_vec_]
                    #gind_eps_ = self.Eps0.gind_rho2eps[0, gind_rho_]
                    #print(gind_eps_-1)
                    if kernel == '3d':                        
                        v_coul_tmp0_ = 8*np.pi/(qxy_abs_)**2
                        wcoul[i,:,j] = v_coul_tmp0_*self.epsinv_mat[i,:G_ind_cut,j] 
                        #wcoul[i,:,j] = v_coul_tmp[0]*self.epsinv_mat[i,:,j] 
                    elif kernel == '2d':
                        v_coul_tmp0_ = 8*np.pi/(qxy_abs_)**2*(1-np.exp(-0.5*qxy_abs_*self.lattpara[2]))
                        wcoul[i,:,j] = v_coul_tmp0_*self.epsinv_mat[i,:G_ind_cut,j] 
                else:
                    j_rho_id = Eps1.gind_eps2rho[q_orgid, j]-1
                    if j_rho_id in epsym_dict[q_vec].keys():
                        j_eps_id = epsym_dict[q_vec][j_rho_id]
                    else:
                        #print('not in dict')
                        continue
                    if j_eps_id > G_ind_cut-1:
                        #print('out of range')
                        continue
                   # print(len(v_coul_tmp))
                    wcoul[i,:,j_eps_id] = v_coul_tmp[0]*self.epsinv_mat[i,:G_ind_cut,j_eps_id] 

        self.wcoul = wcoul
            
            
    def get_wcoul(self, G_ind_cut = 1300, kernel='2d'):
        print('\n')
        print("Calculating screened Coulomb interaction W(q,G,G')")
        try:
            Eps1 = self.Chi1
            Eps0 = self.Chi0
        except(AttributeError):
            Eps1=self.Eps1
            Eps0=self.Eps0
        wcoul = np.zeros((len(Eps1.qpts), G_ind_cut, G_ind_cut), dtype=complex)
        for i in range(len(Eps1.qpts)):
            progress_bar(i + 1, len(Eps1.qpts))
            G_ind_list_tmp = self.g_epsind_tt[self.q_ind_tt == i]
            #print(G_ind_list_tmp)
            if kernel == '3d':
                v_coul_list_tmp = self.v_coul[self.q_ind_tt == i]
                #print(v_coul_list_tmp)
            elif kernel == '2d':
                v_coul_list_tmp = self.v_coul2d[self.q_ind_tt == i]
            elif kernel == 'gf':
                v_coul_list_tmp = self.v_coulgfk[self.q_ind_tt == i]
                
            else:
                print('Undefined kernel type!')
                return
            
            for j in range(G_ind_cut):
                v_coul_tmp = v_coul_list_tmp[G_ind_list_tmp == j+1]
                if i==0 and j ==0 :
                    
                    q0_ = self.Eps0.qpts[0]
                    qxy_abs_ = np.sqrt((q0_[0])**2+ (np.sqrt(3)/3*q0_[0]+2*np.sqrt(3)/3*q0_[1])**2)*2*np.pi/self.lattpara_unit[0]*const.Bohr_R
                    #print(v_coul_tmp)
                    #gp_=[0,0,0]
                    #gp_vec_ = tuple(gp_)
                    #gind_rho_ = self.Eps0.G_vec2ind[gp_vec_]
                    #gind_eps_ = self.Eps0.gind_rho2eps[0, gind_rho_]
                    #print(gind_eps_-1)
                    if kernel == '3d':                        
                        v_coul_tmp0_ = 8*np.pi/(qxy_abs_)**2
                        wcoul[i,:,j] = v_coul_tmp0_*self.epsinv_mat[i,:,j] 
                        #wcoul[i,:,j] = v_coul_tmp[0]*self.epsinv_mat[i,:,j] 
                    elif kernel == '2d':
                        v_coul_tmp0_ = 8*np.pi/(qxy_abs_)**2*(1-np.exp(-0.5*qxy_abs_*self.lattpara[2]))
                        wcoul[i,:,j] = v_coul_tmp0_*self.epsinv_mat[i,:,j] 
                else:
                    wcoul[i,:,j] = v_coul_tmp[0]*self.epsinv_mat[i,:,j] 

        self.wcoul = wcoul

        
        
    def rho2pot_bare(self, kernel='3d'):
        if kernel == '3d':
            vcoul_tmp_ = self.v_coul
        elif kernel == '2d':
            vcoul_tmp_ = self.v_coul2d
        elif kernel == 'gf':
            vcoul_tmp_ = self.v_coulgfk
        else:
            print('Undefined kernel type!')
            return
        self.pot_bare_k = self.rho_bare_k*vcoul_tmp_
        self.pot_bare_r = np.fft.ifftn(self.pot_bare_k)
    
    def rho2pot_tot(self, kernel='3d'):
        if kernel == '3d':
            vcoul_tmp_ = self.v_coul
        elif kernel == '2d':
            vcoul_tmp_ = self.v_coul2d
        elif kernel == 'gf':
            vcoul_tmp_ = self.v_coulgfk
        else:
            print('Undefined kernel type!')
            return
        self.pot_tot_k = self.rho_tot_k*vcoul_tmp_
        self.pot_tot_r = np.fft.ifftn(self.pot_tot_k)

    def pot2rho_bare(self, kernel='3d', ncharge=0):
        if kernel == '3d':
            vcoul_tmp_ = self.v_coul
        elif kernel == '2d':
            vcoul_tmp_ = self.v_coul2d
        elif kernel == 'gf':
            vcoul_tmp_ = self.v_coulgfk
        else:
            print('Undefined kernel type!')
            return
        self.rho_bare_k = self.pot_bare_k/vcoul_tmp_
        self.rho_bare_k[0,0,0] = ncharge/self.omega*self.fft_nx*self.fft_ny*self.fft_nz
        self.rho_bare_r = np.fft.ifftn(self.rho_bare_k)


    def pot2rho_tot(self, kernel='3d', ncharge=0):
        if kernel == '3d':
            vcoul_tmp_ = self.v_coul
        elif kernel == '2d':
            vcoul_tmp_ = self.v_coul2d
        elif kernel == 'gf':
            vcoul_tmp_ = self.v_coulgfk
        else:
            print('Undefined kernel type!')
            return
        self.rho_tot_k = self.pot_tot_k/vcoul_tmp_
        self.rho_tot_k[0,0,0] = ncharge/self.omega*self.fft_nx*self.fft_ny*self.fft_nz
        self.rho_tot_r = np.fft.ifftn(self.rho_tot_k)

    def pot_bare2tot(self, kernel='3d', G_ind_cut=1300):
        if kernel == '3d':
            vcoul_tmp_ = self.v_coul
        elif kernel == '2d':
            vcoul_tmp_ = self.v_coul2d
        elif kernel == 'gf':
            vcoul_tmp_ = self.v_coulgfk
        else:
            print('Undefined kernel type!')
            return
        
        phi_k = np.zeros((self.fft_nx, self.fft_ny, self.fft_nz), dtype = complex)
        phi_k[self.g_epsind_tt > G_ind_cut] = self.pot_bare_k[self.g_epsind_tt > G_ind_cut]
        phik_list = []

        g_epsind_reduced = self.g_epsind_tt[self.g_epsind_tt<=G_ind_cut]
        g_epsind0_reduced = self.g_epsind_tt[self.g_epsind_tt<=G_ind_cut]
        q_ind_reduced = self.q_ind_tt[self.g_epsind_tt<=G_ind_cut]
        pot_reduced = self.pot_bare_k[self.g_epsind_tt<=G_ind_cut]
        v_coul_reduced = vcoul_tmp_[self.g_epsind_tt<=G_ind_cut]
        nnn = len(self.g_epsind_tt[self.g_epsind_tt<=G_ind_cut])
        for i in range(nnn):
            G_index = g_epsind_reduced[i]
            q_index = q_ind_reduced[i]
            pot0 = pot_reduced[q_ind_reduced==q_index]
            eps_tmp=self.epsinv_mat[q_ind_reduced[q_ind_reduced==q_index],
                               G_index-1,
                               g_epsind_reduced[q_ind_reduced==q_index]-1]
            phi_tmp = eps_tmp@pot0.T
            phik_list.append(phi_tmp)
            progress_bar(i + 1, nnn)
            #print(i)
        phi_k[self.g_epsind_tt<=G_ind_cut] = phik_list
        #phi_k[0,0,0] = 0
        if kernel == '3d' or kernel == '2d':
            phi_k[0,0,0] = 0
        self.pot_tot_k = phi_k
        self.pot_tot_r = np.fft.ifftn(self.pot_tot_k)
        
    def rho_ext2pot_tot(self, rho_ext_k, G_ind_cut=1300):
        phi_k = np.zeros((self.fft_nx, self.fft_ny, self.fft_nz), dtype = complex)
        phi_k[self.g_epsind_tt > G_ind_cut] = 0#self.pot_tot_k[self.g_epsind_tt > G_ind_cut]
        phik_list = []
        
        g_epsind_reduced = self.g_epsind_tt[self.g_epsind_tt<=G_ind_cut]
        g_epsind0_reduced = self.g_epsind0_tt[self.g_epsind_tt<=G_ind_cut]
        q_ind_reduced = self.q_ind_tt[self.g_epsind_tt<=G_ind_cut]
        rho_ext_reduced = rho_ext_k[self.g_epsind_tt<=G_ind_cut]
        #pot_reduced = self.pot_bare_k[self.g_epsind_tt<=G_ind_cut]
        #v_coul_reduced = vcoul_tmp_[self.g_epsind_tt<=G_ind_cut]
        nnn = len(self.g_epsind_tt[self.g_epsind_tt<=G_ind_cut])
        for i in range(nnn):
            #G_index = g_epsind_reduced[i]
            q_index = q_ind_reduced[i]
            rho_ext0 = rho_ext_reduced[q_ind_reduced==q_index]
            if q_index == 0:
                G_index = g_epsind0_reduced[i]
                wcoul_tmp = self.wcoul[q_ind_reduced[q_ind_reduced==q_index],
                                      G_index-1,
                                      g_epsind0_reduced[q_ind_reduced==q_index]-1]
            else:
                G_index = g_epsind_reduced[i]
                wcoul_tmp = self.wcoul[q_ind_reduced[q_ind_reduced==q_index],
                                      G_index-1,
                                      g_epsind_reduced[q_ind_reduced==q_index]-1]
                
            phi_tmp = wcoul_tmp@rho_ext0.T
            phik_list.append(phi_tmp)
            if i%(nnn//10)==0:
                print(i,'/',nnn)
        phi_k[self.g_epsind_tt<=G_ind_cut] = phik_list
        return phi_k
    
    
    def rho_ext2pot_tot_irrbz(self, epsym_dict ,rho_ext_k, G_ind_cut=1300):
        phi_k = np.zeros((self.fft_nx, self.fft_ny, self.fft_nz), dtype = complex)
        phi_k[self.g_epsind_tt > G_ind_cut] = 0#self.pot_tot_k[self.g_epsind_tt > G_ind_cut]
        phik_list = []
        
        g_epsind_reduced = self.g_epsind_tt[self.g_epsind_tt<=G_ind_cut]
        g_epsind0_reduced = self.g_epsind0_tt[self.g_epsind_tt<=G_ind_cut]
        q_ind_reduced = self.q_ind_tt[self.g_epsind_tt<=G_ind_cut]
        q_ind_origin_reduced = self.q_ind_irrbz_tt[self.g_epsind_tt<=G_ind_cut]
        rho_ext_reduced = rho_ext_k[self.g_epsind_tt<=G_ind_cut]
        #pot_reduced = self.pot_bare_k[self.g_epsind_tt<=G_ind_cut]
        #v_coul_reduced = vcoul_tmp_[self.g_epsind_tt<=G_ind_cut]
        nnn = len(self.g_epsind_tt[self.g_epsind_tt<=G_ind_cut])
        for i in range(nnn):
            #G_index = g_epsind_reduced[i]
            q_index = q_ind_reduced[i]
            q_index_origin = self.q_ind_fbz2irrbz[q_index]
            q_vec = self.fbz_q_ind2vec[q_index]
            rho_ext0 = rho_ext_reduced[q_ind_reduced==q_index]
            if q_index == 0:
                G_index = g_epsind0_reduced[i]
                wcoul_tmp = self.wcoul[q_ind_reduced[q_ind_reduced==q_index],
                                      G_index-1,
                                      g_epsind0_reduced[q_ind_reduced==q_index]-1]
            else:
                G_index = g_epsind_reduced[i]
                
                G_rho_id = self.Eps1.gind_eps2rho[q_index_origin, G_index-1]-1
                if G_rho_id in epsym_dict[q_vec].keys():
                    G_eps_id = epsym_dict[q_vec][G_rho_id]
                else:
                    #print('not in dict')
                    phik_list.append(0)
                    continue
                    
                if G_eps_id > G_ind_cut-1:
                    #print('out of range1, please reduce your cutoff')
                    phik_list.append(0)
                    continue
                
                G_eps_id_p_list=[]
                wcoul_tmp=[]
                for G_index_p in g_epsind_reduced[q_ind_reduced==q_index]:
                    G_rho_id_p = self.Eps1.gind_eps2rho[q_index_origin, G_index_p-1]-1
                    if G_rho_id_p in epsym_dict[q_vec].keys():
                        G_eps_id_p = epsym_dict[q_vec][G_rho_id_p]
                    else:
                        #print('not in dict')
                        wcoul_qGGp=0+0*1j
                        continue
                    if G_eps_id_p > G_ind_cut-1:
                        #print('out of range2, please reduce your cutoff')
                        wcoul_qGGp=0+0*1j
                        continue
                    wcoul_qGGp = self.wcoul[q_ind_reduced[q_ind_reduced==q_index],
                                      G_eps_id,
                                      G_eps_id_p]
                    #G_eps_id_p_list.append(G_eps_id_p)
                        
                        
                        
                    wcoul_tmp.append(wcoul_qGGp)
                
            phi_tmp = np.array(wcoul_tmp)@rho_ext0.T
            phik_list.append(phi_tmp)
            if i%(nnn//10)==0:
                print(i,'/',nnn)
        phi_k[self.g_epsind_tt<=G_ind_cut] = phik_list
        return phi_k
        


    def pot_tot2bare(self, kernel='3d', G_ind_cut=1300):
        if kernel == '3d':
            vcoul_tmp_ = self.v_coul
        elif kernel == '2d':
            vcoul_tmp_ = self.v_coul2d
        elif kernel == 'gf':
            vcoul_tmp_ = self.v_coulgfk
        else:
            print('Undefined kernel type!')
            return
        
        phi_k = np.zeros((self.fft_nx, self.fft_ny, self.fft_nz), dtype = complex)
        phi_k[self.g_epsind_tt > G_ind_cut] = self.pot_tot_k[self.g_epsind_tt > G_ind_cut]
        phik_list = []

        g_epsind_reduced = self.g_epsind_tt[self.g_epsind_tt<=G_ind_cut]
        g_epsind0_reduced = self.g_epsind_tt[self.g_epsind_tt<=G_ind_cut]
        q_ind_reduced = self.q_ind_tt[self.g_epsind_tt<=G_ind_cut]
        pot_reduced = self.pot_tot_k[self.g_epsind_tt<=G_ind_cut]
        v_coul_reduced = vcoul_tmp_[self.g_epsind_tt<=G_ind_cut]
        for i in range(len(self.g_epsind_tt[self.g_epsind_tt<=G_ind_cut])):
            G_index = g_epsind_reduced[i]
            q_index = q_ind_reduced[i]
            pot0 = pot_reduced[q_ind_reduced==q_index]
            eps_tmp=self.eps_mat[q_ind_reduced[q_ind_reduced==q_index],
                               G_index-1,
                               g_epsind_reduced[q_ind_reduced==q_index]-1]
            phi_tmp = eps_tmp@pot0.T
            phik_list.append(phi_tmp)
            progress_bar(i + 1, len(self.g_epsind_tt[self.g_epsind_tt<=G_ind_cut]))
            #print(i)
        phi_k[self.g_epsind_tt<=G_ind_cut] = phik_list
        #phi_k[0,0,0] = 0
        self.pot_bare_k = phi_k
        self.pot_bare_r = np.fft.ifftn(self.pot_bare_k)

    def write_xsf(self, ftype=None,filedir='./outfile.xsf'):


        if ftype == None:
            print("Please define ftype")
            return
        elif ftype == 'pot_bare':
            print("Writing pot_bare to xsf file")
            out_file = self.pot_bare_r

        elif ftype == 'pot_tot':
            print("Writing pot_tot to xsf file")
            out_file = self.pot_tot_r

        elif ftype == 'rho_bare':
            print("Writing rho_bare to xsf file")
            out_file = self.rho_bare_r

        elif ftype == 'rho_tot':
            print("Writing rho_tot to xsf file")
            out_file = self.rho_tot_r
        
        else:
            print('Undefined ftype')
            return

        with open(filedir, 'w') as v_file:
            print('Crystal\n PRIMVEC\n', end=' ',file=v_file)
            A = self.lattpara_unit[0]*self.sc[0]
            B = self.lattpara_unit[1]*self.sc[1]
            C = self.lattpara_unit[2]*self.sc[2]
            nx, ny, nz = self.fft_nx, self.fft_ny, self.fft_nz
            if self.ibrav == 4:
                print("%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n" % (A, 0, 0, -0.5*B, np.sqrt(3)/2*B, 0, 0, 0, C),end='', file=v_file)

            elif self.ibrav == 1:
                print("%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n" % (A, 0, 0, 0, B, 0, 0, 0, C),end='', file=v_file)

            print(' PRIMCOORD\n 0    1\n', end=' ',file=v_file)
            print(' BEGIN_BLOCK_DATAGRID_3D\n 3D_PWSCF\n BEGIN_DATAGRID_3D_UNKNOWN\n', end=' ',file=v_file)
            print("%8d %8d %8d\n" % (nx, ny, nz), end=' ', file=v_file)
            print(' 0.00 0.00 0.00\n', end=' ',file=v_file)

            if self.ibrav == 4:
                print("%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n" % (A, 0, 0, -0.5*B, np.sqrt(3)/2*B, 0, 0, 0, C),end='', file=v_file)

            elif self.ibrav == 1:
                print("%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n%16.10f %16.10f %16.10f\n" % (A, 0, 0, 0, B, 0, 0, 0, C),end='', file=v_file)
            
            count = 0

            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        print("%16.10f" % np.real(out_file[i,j,k]), end='', file=v_file)
                        count+=1
                        if count%5==0:
                            print(end='\n', file=v_file)
            if (nx*ny*nz)%5 != 0:
                print(end='\n', file=v_file)

            print(' END_DATAGRID_3D\n', end=' ',file=v_file)
            print(' END_BLOCK_DATAGRID_3D', end=' ',file=v_file)

            v_file.close


#    def eps_interp(self, inttype='dis'):
#        g_epsind_reduced = self.g_epsind_tt[self.g_epsind_tt<=G_ind_cut]
#        g_epsind0_reduced = self.g_epsind_tt[self.g_epsind_tt<=G_ind_cut]
#        q_ind_reduced = self.q_ind_tt[self.g_epsind_tt<=G_ind_cut]
#        pot_reduced = self.pot_bare_k[self.g_epsind_tt<=G_ind_cut]
#        v_coul_reduced = vcoul_tmp_[self.g_epsind_tt<=G_ind_cut]
        

        

        



if __name__ == '__main__':

    cell = lab.bn_12x12
   # cell = lab.mos2_tt
    #cell = lab.bn_6x6
    print('Starting job for '+cell['prefix']+'...')
    #Eps0 =  read_eps.Epsmat(cell['folder']+'chi0mat.h5')
    #Eps1 = read_eps.Epsmat(cell['folder']+'chimat.h5')
    #np.save(cell['folder']+'chi0mat.npy', Eps1.mat)

    ttkw=PotCorr(cell)
    print(ttkw.prefix)
    print(ttkw.lattpara_unit)
   # print(type(mos2.lattpara_unit))

    ttkw.fft_init()
    ttkw.get_vcoul()
    ttkw.get_greenfunc()
    ttkw.read_chi()
   # mos2.get_pointc_pot_bare(position='center', kernel='3d')
    #print(mos2.v_coul)
   # mos2.read_rho_bare(mos2.folder+'inp/rho_6pad12_coarse.xsf')
    ttkw.read_rho_tot(ttkw.folder+'drho.xsf')
    #print(np.shape(mos2.rho_bare_r)) 
    #mos2.rho_bare_r[:,:,:30] = 0
    #mos2.rho_bare_r[:,:,-30:] = 0
    #mos2.read_rho_bare(cell['rho_bare'])
    #mos2.read_rho_tot(cell['rho_tot'])
    
    ttkw.rho2pot_tot(kernel='3d')
    ttkw.epsmat_init()
    #print(mos2.q_ind_tt)
    ttkw.get_epsmat(G_ind_cut = 1000, kernel='3d')
    ttkw.epsmat_inv(G_ind_cut = 1000)
    del ttkw.Chi0
    del ttkw.Chi1
    ttkw.pot_tot2bare(kernel='3d', G_ind_cut=600)

    #print(mos2.pot_tot_r)
    ttkw.pot2rho_bare(kernel='3d', ncharge=1)

    ttkw.write_xsf(ftype='rho_bare', filedir=ttkw.folder+'rho_bare.xsf')
    #mos2.write_xsf(ftype='rho_bare', filedir=mos2.folder+'rho_bare.xsf')
    info = psutil.virtual_memory()
    print(u'Occupied memory：',psutil.Process(os.getpid()).memory_info().rss/1024/1024, 'MB')
    print(u'Total memory：',info.total/1024/1024, 'MB')


    
    print(u'Percent：',info.percent,'%')
    print(u'Number of Cup：',psutil.cpu_count())

