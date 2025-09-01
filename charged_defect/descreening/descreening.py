import numpy as np
import potcorr
import utli
import lab

import matplotlib.pyplot as plt
import const

cell=lab.mos2_6x6_z24
ttkw = potcorr.PotCorr(cell)
ttkw.fft_init()
ttkw.get_vcoul()

rho_dn = utli.read_dat("/anvil/scratch/x-rg47749/mos2/z/24/Rho_dn.dat")
rho_dc = utli.read_dat("/anvil/scratch/x-rg47749/mos2/z/24-spin/Rho_d.dat")

rho_tot = rho_dc-rho_dn

eps1 = '/anvil/scratch/x-rg47749/mos2/bgw_results/6x6/20/chimat.h5'
eps0 = '/anvil/scratch/x-rg47749/mos2/bgw_results/6x6/20/chi0mat.h5'
ttkw.read_epsinv(eps1=eps1, eps0=eps0)

k_symmetry_map =  ttkw.get_k_symmetry_map()

epsym_dict = ttkw.get_epsym_dict(k_symmetry_map,ecut=60)
epsmat_eps2rho_dict, epsmat_eps2eps_irrbz_dict, epsmat_dict = ttkw.gen_mat_dict(k_symmetry_map, epsym_dict)
ttkw.vcoul2d0modify()
dpot_k = np.fft.fftn(rho_tot)*ttkw.v_coul2d
rho_induced_G_dict = ttkw.gen_phi_G_dict(dpot_k, k_symmetry_map, epsmat_eps2rho_dict, epsmat_eps2eps_irrbz_dict, epsmat_dict)
rho_induced_k = ttkw.map_phi_G(rho_induced_G_dict, epsmat_eps2rho_dict)
rho_ext = rho_tot-rho_induced_r
rho_ext = np.load('/anvil/scratch/x-rg47749/mos2/z/22/rho_bare.npy')