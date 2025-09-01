import sys
sys.path.append(r"/home/x-rg47749/code/potcorr")

import read_eps
import  matplotlib.pyplot as plt
import numpy as np
import lab
import potcorr
import const
import h5py
import pickle


#----- Constructing potcorr structure for unitcell -----
cell=lab.mos2_unit_tt
ttkw = potcorr.PotCorr(cell)
ttkw.fft_init()

#eps1 = '/anvil/scratch/x-rg47749/data/mos2/2DEG/cutoff_25/intrinsic/chimat.h5'
#eps0 = '/anvil/scratch/x-rg47749/data/mos2/2DEG/cutoff_25/intrinsic/chi0mat.h5'

eps1='/anvil/projects/x-che190065/rjguo/mos2/dielectric/chi_for_zz/12x12_cut30/chimat.h5'
eps0='/anvil/projects/x-che190065/rjguo/mos2/dielectric/chi_for_zz/12x12_cut30/chi0mat.h5'

ttkw.read_epsinv(eps1=eps1, eps0=eps0)
#------------------------------------------

#----- Lindhard parameters -----
E_f  = 2.87293035e-02 # eV
eff_m = 0.499 # electron mass
#-------------------------------


#----- Function for converting BerkeleyGW's chimat.h5 files -----
def BGW2obz(Chi=ttkw.Eps1, qind=1,Gy_= 0, nGz=71):
    nmtx_ = Chi.nmtx
    qind_ = qind
    #print(ttkw.Eps1.qlist[qind_])
    gin_eps_ = ttkw.Eps1.gind_eps2rho[qind_][:nmtx_[qind_]]

    Gz_list = np.fft.fftfreq(nGz,d=1/nGz)
    Gz_list2 = -np.fft.fftfreq(nGz,d=1/nGz)
    #print(Gz_list[-35:])

    Gind_list = []
    Gind_list2 = []
    for i in Gz_list:
        G_vec_ = tuple((0,Gy_,i))
        G_ind_rho = Chi.G_vec2ind[G_vec_]
        G_ind_eps = Chi.gind_rho2eps[qind_][G_ind_rho]
        #print(G_vec_)
        #print(G_ind_)
        Gind_list.append(G_ind_eps)

    for i in Gz_list2:
        G_vec_ = tuple((0,Gy_,i))
        G_ind_rho = Chi.G_vec2ind[G_vec_]
        G_ind_eps = Chi.gind_rho2eps[qind_][G_ind_rho]
        #print(G_vec_)
        #print(G_ind_)
        Gind_list2.append(G_ind_eps)

    chi_GzGz=np.zeros([len(Gind_list),len(Gind_list)], dtype=complex)
    for i,g0 in enumerate(Gind_list):
        for j,g1 in enumerate(Gind_list2):
            chi_real = Chi.mat[qind_,0,0,g1-1,g0-1,0]
            chi_imag = Chi.mat[qind_,0,0,g1-1,g0-1,1]
            chi_GzGz[i,j]=chi_real+1j*chi_imag
    return chi_GzGz
#-----------------------------------------------------------------

def get_kf(E_f, eff_m):
    return np.sqrt(2*eff_m*E_f/2/const.Ry2eV) # unit: bohr^(-1) 

def get_chi_2DEG(q, kf, eff_m, Ns, Nv):
    return np.where(q>2*kf, -Ns*Nv*eff_m/np.pi/4*(1-np.sqrt(1-4*kf**2/q**2)),-Ns*Nv*eff_m/np.pi/4 ) # energy unit: Ry


kf = get_kf(E_f, eff_m)
print(kf)
kf = 0.07
chi_rbf_path = '/anvil/projects/x-che190065/rjguo/mos2/dielectric/adler-wiser'
with open(chi_rbf_path+'/chi_200k_1e11_rbf.pkl', 'rb') as file:
    chi_1e11_rbf = pickle.load(file)

chi_rbf = chi_1e11_rbf

#----- Real space array along z direction -----
z = np.linspace(-12,12,225)/const.Bohr_R
dist_zz = np.abs(z[:, np.newaxis] - z[np.newaxis, :])
#----------------------------------------------


#----- Reading xi file -----
xi_Gz = np.load('/anvil/projects/x-che190065/rjguo/mos2/dielectric/xi/xi_Gz.npy')
xi = np.fft.ifftn(xi_Gz)/24*const.Bohr_R*225
#xi2 = np.zeros_like(xi)
#xi2[:225//2] = xi[225//2+1:]
#xi2[225//2:] = xi[:225//2+1]
#plt.plot(xi.real)
#xi_Gz2 = np.fft.fftn(xi2)

xizz = np.outer(xi.real,xi.real)
#---------------------------


#------ Calculating 2d dielectric functions and screening interactions -----
nGz=61
nGz_l = 225



epsf= ttkw.Eps0
wlist = []
qabs_list = []
for Gy in range(1):
    for qind in range(0,9):

        q1_ = epsf.qpts[:]
        q_abs_ = np.sqrt((q1_[qind,0])**2+ (np.sqrt(3)/3*(q1_[qind,0])+2*np.sqrt(3)/3*(q1_[qind,1]+Gy))**2)*2*np.pi/ttkw.lattpara_unit[0]*const.Bohr_R
        vcoul_2d= 4*np.pi * np.exp(-q_abs_*dist_zz) / q_abs_
       # print(q_abs_)
        if q_abs_ < 3.6:
            chi_GzGz = BGW2obz(Chi=epsf, qind=qind, Gy_=Gy, nGz=nGz)

            chi_GzGz_l = np.zeros([nGz_l,nGz_l], dtype=complex)
            chi_GzGz_l[:nGz//2+1,:nGz//2+1]=chi_GzGz[:nGz//2+1,:nGz//2+1]
            chi_GzGz_l[:nGz//2+1,-nGz//2+1:]=chi_GzGz[:nGz//2+1,-nGz//2+1:]
            chi_GzGz_l[-nGz//2+1:,:nGz//2+1]=chi_GzGz[-nGz//2+1:,:nGz//2+1]
            chi_GzGz_l[-nGz//2+1:,-nGz//2+1:]=chi_GzGz[-nGz//2+1:,-nGz//2+1:]

          #  chi_GzGz_s = np.zeros_like(chi_GzGz)
          #  chi_GzGz_s[:nGz//2+1,:nGz//2+1]=chi_GzGz[:nGz//2+1,:nGz//2+1]
          #  chi_GzGz_s[:nGz//2+1,-nGz//2+1:]=chi_GzGz[:nGz//2+1,-nGz//2+1:]
          #  chi_GzGz_s[-nGz//2+1:,:nGz//2+1]=chi_GzGz[-nGz//2+1:,:nGz//2+1]
          #  chi_GzGz_s[-nGz//2+1:,-nGz//2+1:]=chi_GzGz[-nGz//2+1:,-nGz//2+1:]

           # chi_zz_s=np.fft.ifftn(chi_GzGz_s)*71*71/(24/const.Bohr_R)**2
            chi_zz_l=np.fft.ifftn(chi_GzGz_l)*225*225/(24/const.Bohr_R)**2

            
            chi_zz_q2DEG = float(get_chi_2DEG(q_abs_, kf, eff_m, 2, 2)) * xizz
            chi_zz_ftq2DEG = float(0.5*chi_rbf([[q_abs_*const.Bohr_R]])) * xizz
            chi_zz_t =  chi_zz_l+chi_zz_ftq2DEG

            
            #epsilon = np.eye(225) - 2*np.pi * vcoul_2d@chi_zz_t*((24/const.Bohr_R)/225)
            epsilon = np.eye(225) - vcoul_2d@chi_zz_t*((24/const.Bohr_R)/225)
        else:
            epsilon = np.eye(225)
        epsilon_inv = np.linalg.inv(epsilon)
        w = epsilon_inv@vcoul_2d*((24/const.Bohr_R)/225)
        wlist.append(w)
        qabs_list.append(q_abs_)
        '''
        plt.imshow(w.real, cmap='hot', interpolation='nearest')
        plt.colorbar()  # 显示颜色条
        plt.title("Heatmap of chi_qxy(z, z')")
        plt.xlabel("z'")
        plt.ylabel("z")
        plt.gca().invert_yaxis()
        #plt.xlim(225//2-20,225//2+20)
        #plt.ylim(225//2-20,225//2+20)
        plt.show()
        '''
for Gy in range(24):
    epsf= ttkw.Eps1
    for qind in range(1,12):

        q1_ = epsf.qpts[:]
        #print(q1_)
        q_abs_ = np.sqrt((q1_[qind,0])**2+ (np.sqrt(3)/3*(q1_[qind,0])+2*np.sqrt(3)/3*(q1_[qind,1]+Gy))**2)*2*np.pi/ttkw.lattpara_unit[0]*const.Bohr_R

        vcoul_2d= 4*np.pi * np.exp(-q_abs_*dist_zz) / q_abs_
      #  print(q_abs_)
        
        if q_abs_ < 3.6:
            chi_GzGz = BGW2obz(Chi=epsf, qind=qind, Gy_=Gy, nGz=nGz)

            chi_GzGz_l = np.zeros([nGz_l,nGz_l], dtype=complex)
            chi_GzGz_l[:nGz//2+1,:nGz//2+1]=chi_GzGz[:nGz//2+1,:nGz//2+1]
            chi_GzGz_l[:nGz//2+1,-nGz//2+1:]=chi_GzGz[:nGz//2+1,-nGz//2+1:]
            chi_GzGz_l[-nGz//2+1:,:nGz//2+1]=chi_GzGz[-nGz//2+1:,:nGz//2+1]
            chi_GzGz_l[-nGz//2+1:,-nGz//2+1:]=chi_GzGz[-nGz//2+1:,-nGz//2+1:]

        #    chi_GzGz_s = np.zeros_like(chi_GzGz)
        #    chi_GzGz_s[:nGz//2+1,:nGz//2+1]=chi_GzGz[:nGz//2+1,:nGz//2+1]
        #    chi_GzGz_s[:nGz//2+1,-nGz//2+1:]=chi_GzGz[:nGz//2+1,-nGz//2+1:]
        #    chi_GzGz_s[-nGz//2+1:,:nGz//2+1]=chi_GzGz[-nGz//2+1:,:nGz//2+1]
        #    chi_GzGz_s[-nGz//2+1:,-nGz//2+1:]=chi_GzGz[-nGz//2+1:,-nGz//2+1:]

        #    chi_zz_s=np.fft.ifftn(chi_GzGz_s)*71*71/(24/const.Bohr_R)**2
            chi_zz_l=np.fft.ifftn(chi_GzGz_l)*225*225/(24/const.Bohr_R)**2


            chi_zz_q2DEG = float(get_chi_2DEG(q_abs_, kf, eff_m, 2, 2)) * xizz
            chi_zz_ftq2DEG = float(0.5*chi_rbf([[q_abs_*const.Bohr_R]])) * xizz
            chi_zz_t = chi_zz_l+chi_zz_ftq2DEG
           # chi_zz_t = chi_zz_q2DEG#+chi_zz_l
        
        #    epsilon = np.eye(225) - 2*np.pi * vcoul_2d@chi_zz_t#*((24/const.Bohr_R)/225)
            epsilon = np.eye(225) - vcoul_2d@chi_zz_t*((24/const.Bohr_R)/225)
        else:
            epsilon = np.eye(225)
        epsilon_inv = np.linalg.inv(epsilon)
        w = epsilon_inv@vcoul_2d*((24/const.Bohr_R)/225)
        wlist.append(w)
        qabs_list.append(q_abs_)
#---------------------------------------------------------------------------


#----- Interpolation for screened interactions -----
qabs_list = qabs_list
wlist = np.array(wlist)

from scipy.interpolate import interp1d
w_zz_interp = []
for i in range(225):
    print(i)
    w_z_interp = []
    for j in range(225):
        linear_interp = interp1d(np.array(qabs_list),np.array(wlist)[:,i,j])

        #q_ = np.linspace(1e-3,24,10000)
        #w_ = linear_interp(q_)
        w_z_interp.append(linear_interp)
    w_zz_interp.append(w_z_interp)
#---------------------------------------------------

#----- Constructing potcorr structure for supercell -----
cell=lab.mos2_12to48
ttkw2 = potcorr.PotCorr(cell)
ttkw2.fft_init()
#--------------------------------------------------------

rho_bare=np.load('/anvil/projects/x-che190065/rjguo/mos2/potential/dft/12x12/rho_dl.npy')

nx, ny, nz = rho_bare.shape

#----- Moving the defect center to (0,0)
rho_bare2 = np.zeros_like(rho_bare)
rho_bare2[:nx//2,:nx//2,:]=rho_bare[nx//2:,nx//2:,:]
rho_bare2[:nx//2,-nx//2:,:]=rho_bare[nx//2:,:nx//2,:]
rho_bare2[-nx//2:,:nx//2,:]=rho_bare[:nx//2,-nx//2:,:]
rho_bare2[-nx//2:,-nx//2:,:]=rho_bare[:nx//2,:nx//2,:]

fine_data = rho_bare2
f_transform = np.fft.fftn(fine_data)



new_f_transform = np.zeros((nx//2,ny//2,nz), dtype=complex)
nx, ny, nz = new_f_transform.shape

new_f_transform[:nx//2, :ny//2, :nz//2] = f_transform[:nx//2, :ny//2, :nz//2]
new_f_transform[-nx//2:, :ny//2, :nz//2] = f_transform[-nx//2:, :ny//2, :nz//2]
new_f_transform[:nx//2, -ny//2:, :nz//2] = f_transform[:nx//2, -ny//2:, :nz//2]
new_f_transform[-nx//2:, -ny//2:, :nz//2] = f_transform[-nx//2:, -ny//2:, :nz//2]
new_f_transform[:nx//2, :ny//2, -nz//2:] = f_transform[:nx//2, :ny//2, -nz//2:]
new_f_transform[-nx//2:, :ny//2, -nz//2:] = f_transform[-nx//2:, :ny//2, -nz//2:]
new_f_transform[:nx//2, -ny//2:, -nz//2:] = f_transform[:nx//2, -ny//2:, -nz//2:]
new_f_transform[-nx//2:, -ny//2:, -nz//2:] = f_transform[-nx//2:, -ny//2:, -nz//2:]


fine_data2 = np.fft.ifftn(new_f_transform).real/4 

rho_bare_l = np.zeros((720,720,225), dtype=complex )

rho_bare_l[:nx//2, :nx//2, :] = fine_data2[:nx//2, :nx//2, :]
rho_bare_l[:nx//2, -nx//2:, :] = fine_data2[:nx//2, -nx//2:, :]
rho_bare_l[-nx//2:, :nx//2, :] = fine_data2[-nx//2:, :nx//2, :]
rho_bare_l[-nx//2:, -nx//2:, :] = fine_data2[-nx//2:, -nx//2:, :]



fft_q_abs_ = np.sqrt(ttkw2.fft_kxx[:,:,ttkw.fft_nz//2]**2 + (ttkw2.fft_kyy[:,:,ttkw2.fft_nz//2]*2*np.sqrt(3)/3 +np.sqrt(3)/3* ttkw2.fft_kxx[:,:,ttkw2.fft_nz//2])**2)
fft_q_abs_[0,0] = 0.00057
fft_q_abs_[fft_q_abs_>24]=24

A_uc = ttkw2.lattpara_unit[0]*ttkw2.lattpara_unit[0]/const.Bohr_R/const.Bohr_R*np.sqrt(3)/2

pot_k_z = np.zeros([720,720,225], dtype=complex)
for i in range(225):
    print(i)
    for j in range(225//2-40,225//2+40):
        rho_k_l_ = np.fft.fftn(rho_bare_l[:,:,j])*ttkw2.lattpara[0]*ttkw2.lattpara[1]/720/720*np.sqrt(3)/2 / A_uc
        w_zz_ = w_zz_interp[i][j](fft_q_abs_)
        #w_zz_[0,0]=0
        #print(w_zz_[0,0])
        pot_k_z[:,:,i]= pot_k_z[:,:,i] + w_zz_ * rho_k_l_*ttkw2.lattpara[2]/225
        
np.save('/anvil/projects/x-che190065/rjguo/mos2/potential/standard/12x12/1e11_200K_dl/pot_k_z_2.npy',pot_k_z[:,:,225//2-40: 225//2+40] )
        
        
