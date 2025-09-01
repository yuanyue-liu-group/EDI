 
import const
import numpy as np

tt = {'folder':'/anvil/scratch/x-rg47749/data/'}

#---mos2 sym_arr---
r1 = np.array([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
r2 = np.array([[-1, 1, 0], [-1, 0, 0], [0, 0, 1]])
r3 = np.array([[-1, 1, 0], [0, 1, 0], [0, 0, 1]])
r4 = np.array([[0, -1, 0], [-1, 0, 0], [0, 0, 1]])
r5 = np.array([[1, 0, 0], [1, -1, 0], [0, 0, 1]])
r6 = np.array([[0, -1, 0], [1, -1, 0], [0, 0, 1]])
sym_arr = np.array([r1, r2, r3, r4, r5, r6])
#-----

#----Si sym_arr----
r1 =   np.array([ [ 1,   0,   0] ,   [ 0,   1,   0] , [  -1,  -1,  -1] ])
r2 =   np.array([  [0,  -1,   0] ,  [ -1,   0,   0] ,  [  0,   0,  -1 ]])
r3 =   np.array([  [0,  -1,   0] ,  [ -1,   0,   0 ] , [  1,   1,   1]])
r4 =   np.array([  [1,   1,   1] ,  [  0,  -1,   0  ], [  0,   0,  -1]])
r5 =   np.array([  [0,   0,  -1] ,   [ 1,   1,   1  ], [ -1,   0,   0]])
r6 =   np.array([  [0,   0,  -1] ,    [1,   1,   1  ], [  0,  -1,   0]])
r7 =   np.array([ [-1,  -1,  -1] ,    [0,   0,   1  ], [  1,   0,   0]])
r8 =   np.array([ [-1,  -1,  -1] ,    [0,   0,   1  ], [  0,   1,   0]])
r9 =   np.array([ [ 1,   1,   1] ,    [0,   0,  -1  ],  [-1,   0,   0]])
r10 = np.array([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
si_sym_arr = np.array([r1, r2, r3, r4, r5, r6, r7, r8, r9,r10])

#------


gb = {
    'prefix': 'mos2',
    'folder': tt['folder']+'mos2/12x12/',
    'folder_chi': tt['folder']+'mos2_12x12/chi/',
    'folder_chi_doped': tt['folder']+'mos2_12x12/chi_fermi_0.3/',
    'folder_dft': tt['folder']+'mos2_12x12/dft/',
    'folder_inp': tt['folder']+'mos2_12x12/inp/',
    'folder_out': tt['folder']+'mos2_12x12/out/',
    'rho_bare': tt['folder']+'mos2_12x12/inp/rho_6pad12_coarse.xsf',
    'is2d': False,
    'ibrav': 0,
    'cell_r':np.array([[1,0,0],[0,1,0],[0,0,20]]),
    'lattpara_unit':[5,
                     5,
                     100.0],
    'sc': [1,1,1],
    'fft_g':[10,10,200],
    'A' : 5,
    'sym_arr': []
}

Si_888 = {
    'prefix': 'Si',
    'folder': tt['folder']+'mos2/12x12/',
    'folder_chi': tt['folder']+'mos2_12x12/chi/',
    'folder_chi_doped': tt['folder']+'mos2_12x12/chi_fermi_0.3/',
    'folder_dft': tt['folder']+'mos2_12x12/dft/',
    'folder_inp': tt['folder']+'mos2_12x12/inp/',
    'folder_out': tt['folder']+'mos2_12x12/out/',
    'rho_bare': tt['folder']+'mos2_12x12/inp/rho_6pad12_coarse.xsf',
    'is2d': False,
    'ibrav': 0,
    'cell_r':np.array([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]]),
    'lattpara_unit':[5.43,
                     5.43,
                     5.43],
    'sc': [8,8,8],
    'fft_g':[288,288,288],
    'A': 5.43*8,
    'sym_arr': si_sym_arr
}

Si_666 = {
    'prefix': 'Si',
    'folder': tt['folder']+'mos2/12x12/',
    'folder_chi': tt['folder']+'mos2_12x12/chi/',
    'folder_chi_doped': tt['folder']+'mos2_12x12/chi_fermi_0.3/',
    'folder_dft': tt['folder']+'mos2_12x12/dft/',
    'folder_inp': tt['folder']+'mos2_12x12/inp/',
    'folder_out': tt['folder']+'mos2_12x12/out/',
    'rho_bare': tt['folder']+'mos2_12x12/inp/rho_6pad12_coarse.xsf',
    'is2d': False,
    'ibrav': 0,
    'cell_r':np.array([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]]),
    'lattpara_unit':[5.43,
                     5.43,
                     5.43],
    'sc': [6,6,6],
    'fft_g':[216,216,216],
    'A': 5.43*6,
    'sym_arr': si_sym_arr
}

Si_444 = {
    'prefix': 'Si',
    'folder': tt['folder']+'mos2/12x12/',
    'folder_chi': tt['folder']+'mos2_12x12/chi/',
    'folder_chi_doped': tt['folder']+'mos2_12x12/chi_fermi_0.3/',
    'folder_dft': tt['folder']+'mos2_12x12/dft/',
    'folder_inp': tt['folder']+'mos2_12x12/inp/',
    'folder_out': tt['folder']+'mos2_12x12/out/',
    'rho_bare': tt['folder']+'mos2_12x12/inp/rho_6pad12_coarse.xsf',
    'is2d': False,
    'ibrav': 0,
    'cell_r':np.array([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]]),
    'lattpara_unit':[5.43,
                     5.43,
                     5.43],
    'sc': [4,4,4],
    'fft_g':[144,144,144],
    'A': 5.43*4,
    'sym_arr': si_sym_arr
}

Si_unit = {
    'prefix': 'Si',
    'folder': tt['folder']+'mos2/12x12/',
    'folder_chi': tt['folder']+'mos2_12x12/chi/',
    'folder_chi_doped': tt['folder']+'mos2_12x12/chi_fermi_0.3/',
    'folder_dft': tt['folder']+'mos2_12x12/dft/',
    'folder_inp': tt['folder']+'mos2_12x12/inp/',
    'folder_out': tt['folder']+'mos2_12x12/out/',
    'rho_bare': tt['folder']+'mos2_12x12/inp/rho_6pad12_coarse.xsf',
    'is2d': False,
    'ibrav': 0,
    'cell_r':np.array([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]]),
    'lattpara_unit':[5.43,
                     5.43,
                     5.43],
    'sc': [1,1,1],
    'fft_g':[32,32,32],
    'A': 5.43,
    'sym_arr': si_sym_arr
}


mos2_12to24_largez = {
    'prefix': 'mos2',
    'folder': tt['folder']+'mos2/12x12/',
    'folder_chi': tt['folder']+'mos2_12x12/chi/',
    'folder_chi_doped': tt['folder']+'mos2_12x12/chi_fermi_0.3/',
    'folder_dft': tt['folder']+'mos2_12x12/dft/',
    'folder_inp': tt['folder']+'mos2_12x12/inp/',
    'folder_out': tt['folder']+'mos2_12x12/out/',
    'rho_bare': tt['folder']+'mos2_12x12/inp/rho_6pad12_coarse.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[3.186,
                     3.186,
                     100.0],
    'sc': [24,24,1],
    'fft_g':[720,720,960],
    'sym_arr': sym_arr
}


mos2_12to24_irrbz = {
    'prefix': 'mos2',
    'folder': tt['folder']+'mos2/12x12/',
    'folder_chi': tt['folder']+'mos2_12x12/chi/',
    'folder_chi_doped': tt['folder']+'mos2_12x12/chi_fermi_0.3/',
    'folder_dft': tt['folder']+'mos2_12x12/dft/',
    'folder_inp': tt['folder']+'mos2_12x12/inp/',
    'folder_out': tt['folder']+'mos2_12x12/out/',
    'rho_bare': tt['folder']+'mos2_12x12/inp/rho_6pad12_coarse.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[3.186,
                     3.186,
                     24.0],
    'sc': [24,24,1],
    'fft_g':[180,180,60],
    'sym_arr': sym_arr
}


mos2_unit_largez = {
    'prefix': 'mos2',
    'folder': '/anvil/projects/x-che190065/rjguo/mos2/wfc/mos2.save',
    'folder_chi': tt['folder']+'mos2_tt/chi/',
    'folder_dft': tt['folder']+'mos2_tt/dft/',
    'folder_inp': tt['folder']+'mos2_tt/inp/',
    'folder_out': tt['folder']+'mos2_tt/out/',
    'rho_bare': tt['folder']+'mos2_tts/rho_test.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[3.186,
                     3.186,
                    100.0],
    'sc': [1,1,1],
    'fft_g':[30,30,960],
    'sym_arr': sym_arr
}


mos2_unit_tt = {
    'prefix': 'mos2',
    'folder': '/anvil/projects/x-che190065/rjguo/mos2/wfc/mos2.save',
    'folder_chi': tt['folder']+'mos2_tt/chi/',
    'folder_dft': tt['folder']+'mos2_tt/dft/',
    'folder_inp': tt['folder']+'mos2_tt/inp/',
    'folder_out': tt['folder']+'mos2_tt/out/',
    'rho_bare': tt['folder']+'mos2_tts/rho_test.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[3.186,
                     3.186,
                    24],
    'sc': [1,1,1],
    'fft_g':[30,30,225],
    'sym_arr': sym_arr
}
mos2_unit_tt_cut120 = {
    'prefix': 'mos2',
    'folder': '/anvil/projects/x-che190065/rjguo/mos2/wfc/mos2.save',
    'folder_chi': tt['folder']+'mos2_tt/chi/',
    'folder_dft': tt['folder']+'mos2_tt/dft/',
    'folder_inp': tt['folder']+'mos2_tt/inp/',
    'folder_out': tt['folder']+'mos2_tt/out/',
    'rho_bare': tt['folder']+'mos2_tts/rho_test.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[3.186,
                     3.186,
                    24],
    'sc': [1,1,1],
    'fft_g':[45,45,320],
    'sym_arr': sym_arr
}


bn_unit_new = {
    'prefix': 'bn',
    'folder': '/anvil/projects/x-che190065/rjguo/mos2/wfc/mos2.save',
    'folder_chi': tt['folder']+'mos2_tt/chi/',
    'folder_dft': tt['folder']+'mos2_tt/dft/',
    'folder_inp': tt['folder']+'mos2_tt/inp/',
    'folder_out': tt['folder']+'mos2_tt/out/',
    'rho_bare': tt['folder']+'mos2_tts/rho_test.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[2.51,
                     2.51,
                    24],
    'sc': [1,1,1],
    'fft_g':[24,24,225],
    'sym_arr': sym_arr
}

mos2_unit_new = {
    'prefix': 'mos2',
    'folder': '/anvil/projects/x-che190065/rjguo/mos2/wfc/mos2.save',
    'folder_chi': tt['folder']+'mos2_tt/chi/',
    'folder_dft': tt['folder']+'mos2_tt/dft/',
    'folder_inp': tt['folder']+'mos2_tt/inp/',
    'folder_out': tt['folder']+'mos2_tt/out/',
    'rho_bare': tt['folder']+'mos2_tts/rho_test.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[3.186,
                     3.186,
                    24],
    'sc': [1,1,1],
    'fft_g':[15,15,225],
    'sym_arr': sym_arr
}

mos2_18x18 = {
    'prefix': 'mos2',
    'folder': tt['folder']+'mos2/12x12/',
    'folder_chi': tt['folder']+'mos2_12x12/chi/',
    'folder_chi_doped': tt['folder']+'mos2_12x12/chi_fermi_0.3/',
    'folder_dft': tt['folder']+'mos2_12x12/dft/',
    'folder_inp': tt['folder']+'mos2_12x12/inp/',
    'folder_out': tt['folder']+'mos2_12x12/out/',
    'rho_bare': tt['folder']+'mos2_12x12/inp/rho_6pad12_coarse.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[3.186,
                     3.186,
                     24],
    'sc': [18,18,1],
    'fft_g':[540,540,225],
    'sym_arr': sym_arr
}

mos2_12x12_old = {
    'prefix': 'mos2',
    'folder': tt['folder']+'mos2_12x12/',
    'folder_chi': tt['folder']+'mos2_12x12/chi/',
    'folder_chi_doped': tt['folder']+'mos2_12x12/chi_fermi_0.3/',
    'folder_dft': tt['folder']+'mos2_12x12/dft/',
    'folder_inp': tt['folder']+'mos2_12x12/inp/',
    'folder_out': tt['folder']+'mos2_12x12/out/',
    'rho_bare': tt['folder']+'mos2_12x12/inp/rho_6pad12_coarse.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[3.186,
                     3.186,
                     3.186*8],
    'sc': [12,12,1],
    'fft_g':[180,180,120]
}

mos2_15x15 = {
    'prefix': 'mos2',
    'folder': tt['folder']+'mos2/12x12/',
    'folder_chi': tt['folder']+'mos2_12x12/chi/',
    'folder_chi_doped': tt['folder']+'mos2_12x12/chi_fermi_0.3/',
    'folder_dft': tt['folder']+'mos2_12x12/dft/',
    'folder_inp': tt['folder']+'mos2_12x12/inp/',
    'folder_out': tt['folder']+'mos2_12x12/out/',
    'rho_bare': tt['folder']+'mos2_12x12/inp/rho_6pad12_coarse.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[3.186,
                     3.186,
                     24.0],
    'sc': [15,15,1],
    'fft_g':[450,450,225],
    'sym_arr': sym_arr
}


mos2_12x12 = {
    'prefix': 'mos2',
    'folder': tt['folder']+'mos2/12x12/',
    'folder_chi': tt['folder']+'mos2_12x12/chi/',
    'folder_chi_doped': tt['folder']+'mos2_12x12/chi_fermi_0.3/',
    'folder_dft': tt['folder']+'mos2_12x12/dft/',
    'folder_inp': tt['folder']+'mos2_12x12/inp/',
    'folder_out': tt['folder']+'mos2_12x12/out/',
    'rho_bare': tt['folder']+'mos2_12x12/inp/rho_6pad12_coarse.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[3.186,
                     3.186,
                     24.0],
    'sc': [12,12,1],
    'fft_g':[360,360,225],
    'sym_arr': sym_arr
}

mos2_9x9 = {
    'prefix': 'mos2',
    'folder': tt['folder']+'mos2/12x12/',
    'folder_chi': tt['folder']+'mos2_12x12/chi/',
    'folder_chi_doped': tt['folder']+'mos2_12x12/chi_fermi_0.3/',
    'folder_dft': tt['folder']+'mos2_12x12/dft/',
    'folder_inp': tt['folder']+'mos2_12x12/inp/',
    'folder_out': tt['folder']+'mos2_12x12/out/',
    'rho_bare': tt['folder']+'mos2_12x12/inp/rho_6pad12_coarse.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[3.186,
                     3.186,
                     24.0],
    'sc': [9,9,1],
    'fft_g':[270,270,225],
    'sym_arr': sym_arr
}

mos2_6x6_z20 = {
    'prefix': 'mos2',
    'folder': tt['folder']+'mos2/12x12/',
    'folder_chi': tt['folder']+'mos2_12x12/chi/',
    'folder_chi_doped': tt['folder']+'mos2_12x12/chi_fermi_0.3/',
    'folder_dft': tt['folder']+'mos2_12x12/dft/',
    'folder_inp': tt['folder']+'mos2_12x12/inp/',
    'folder_out': tt['folder']+'mos2_12x12/out/',
    'rho_bare': tt['folder']+'mos2_12x12/inp/rho_6pad12_coarse.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[3.186,
                     3.186,
                     20.0],
    'sc': [6,6,1],
    'fft_g':[180,180,192],
    'sym_arr': sym_arr
}

mos2_6x6_z22 = {
    'prefix': 'mos2',
    'folder': tt['folder']+'mos2/12x12/',
    'folder_chi': tt['folder']+'mos2_12x12/chi/',
    'folder_chi_doped': tt['folder']+'mos2_12x12/chi_fermi_0.3/',
    'folder_dft': tt['folder']+'mos2_12x12/dft/',
    'folder_inp': tt['folder']+'mos2_12x12/inp/',
    'folder_out': tt['folder']+'mos2_12x12/out/',
    'rho_bare': tt['folder']+'mos2_12x12/inp/rho_6pad12_coarse.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[3.186,
                     3.186,
                     22.0],
    'sc': [6,6,1],
    'fft_g':[180,180,216],
    'sym_arr': sym_arr
}

mos2_6x6_z24 = {
    'prefix': 'mos2',
    'folder': tt['folder']+'mos2/12x12/',
    'folder_chi': tt['folder']+'mos2_12x12/chi/',
    'folder_chi_doped': tt['folder']+'mos2_12x12/chi_fermi_0.3/',
    'folder_dft': tt['folder']+'mos2_12x12/dft/',
    'folder_inp': tt['folder']+'mos2_12x12/inp/',
    'folder_out': tt['folder']+'mos2_12x12/out/',
    'rho_bare': tt['folder']+'mos2_12x12/inp/rho_6pad12_coarse.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[3.186,
                     3.186,
                     24.0],
    'sc': [6,6,1],
    'fft_g':[180,180,225],
    'sym_arr': sym_arr
}

mos2_6x6_z26 = {
    'prefix': 'mos2',
    'folder': tt['folder']+'mos2/12x12/',
    'folder_chi': tt['folder']+'mos2_12x12/chi/',
    'folder_chi_doped': tt['folder']+'mos2_12x12/chi_fermi_0.3/',
    'folder_dft': tt['folder']+'mos2_12x12/dft/',
    'folder_inp': tt['folder']+'mos2_12x12/inp/',
    'folder_out': tt['folder']+'mos2_12x12/out/',
    'rho_bare': tt['folder']+'mos2_12x12/inp/rho_6pad12_coarse.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[3.186,
                     3.186,
                     26.0],
    'sc': [6,6,1],
    'fft_g':[180,180,243],
    'sym_arr': sym_arr
}

mos2_6x6_z28 = {
    'prefix': 'mos2',
    'folder': tt['folder']+'mos2/12x12/',
    'folder_chi': tt['folder']+'mos2_12x12/chi/',
    'folder_chi_doped': tt['folder']+'mos2_12x12/chi_fermi_0.3/',
    'folder_dft': tt['folder']+'mos2_12x12/dft/',
    'folder_inp': tt['folder']+'mos2_12x12/inp/',
    'folder_out': tt['folder']+'mos2_12x12/out/',
    'rho_bare': tt['folder']+'mos2_12x12/inp/rho_6pad12_coarse.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[3.186,
                     3.186,
                     28.0],
    'sc': [6,6,1],
    'fft_g':[180,180,270],
    'sym_arr': sym_arr
}

mos2_6x6_z36 = {
    'prefix': 'mos2',
    'folder': tt['folder']+'mos2/12x12/',
    'folder_chi': tt['folder']+'mos2_12x12/chi/',
    'folder_chi_doped': tt['folder']+'mos2_12x12/chi_fermi_0.3/',
    'folder_dft': tt['folder']+'mos2_12x12/dft/',
    'folder_inp': tt['folder']+'mos2_12x12/inp/',
    'folder_out': tt['folder']+'mos2_12x12/out/',
    'rho_bare': tt['folder']+'mos2_12x12/inp/rho_6pad12_coarse.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[3.186,
                     3.186,
                     36.0],
    'sc': [6,6,1],
    'fft_g':[180,180,360],
    'sym_arr': sym_arr
}

mos2_6x6 = {
    'prefix': 'mos2',
    'folder': tt['folder']+'mos2/12x12/',
    'folder_chi': tt['folder']+'mos2_12x12/chi/',
    'folder_chi_doped': tt['folder']+'mos2_12x12/chi_fermi_0.3/',
    'folder_dft': tt['folder']+'mos2_12x12/dft/',
    'folder_inp': tt['folder']+'mos2_12x12/inp/',
    'folder_out': tt['folder']+'mos2_12x12/out/',
    'rho_bare': tt['folder']+'mos2_12x12/inp/rho_6pad12_coarse.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[3.186,
                     3.186,
                     24.0],
    'sc': [6,6,1],
    'fft_g':[180,180,225],
    'sym_arr': sym_arr
}

mos2_6x6_cut120 = {
    'prefix': 'mos2',
    'folder': tt['folder']+'mos2/12x12/',
    'folder_chi': tt['folder']+'mos2_12x12/chi/',
    'folder_chi_doped': tt['folder']+'mos2_12x12/chi_fermi_0.3/',
    'folder_dft': tt['folder']+'mos2_12x12/dft/',
    'folder_inp': tt['folder']+'mos2_12x12/inp/',
    'folder_out': tt['folder']+'mos2_12x12/out/',
    'rho_bare': tt['folder']+'mos2_12x12/inp/rho_6pad12_coarse.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[3.186,
                     3.186,
                     24.0],
    'sc': [6,6,1],
    'fft_g':[256,256,320],
    'sym_arr': sym_arr
}


mos2_12to48_new = {
    'prefix': 'mos2',
    'folder': tt['folder']+'mos2/12x12/',
    'folder_chi': tt['folder']+'mos2_12x12/chi/',
    'folder_chi_doped': tt['folder']+'mos2_12x12/chi_fermi_0.3/',
    'folder_dft': tt['folder']+'mos2_12x12/dft/',
    'folder_inp': tt['folder']+'mos2_12x12/inp/',
    'folder_out': tt['folder']+'mos2_12x12/out/',
    'rho_bare': tt['folder']+'mos2_12x12/inp/rho_6pad12_coarse.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[3.186,
                     3.186,
                     24.0],
    'sc': [48,48,1],
    'fft_g':[360,360,225],
    'sym_arr': sym_arr
}


mos2_12to48 = {
    'prefix': 'mos2',
    'folder': tt['folder']+'mos2/12x12/',
    'folder_chi': tt['folder']+'mos2_12x12/chi/',
    'folder_chi_doped': tt['folder']+'mos2_12x12/chi_fermi_0.3/',
    'folder_dft': tt['folder']+'mos2_12x12/dft/',
    'folder_inp': tt['folder']+'mos2_12x12/inp/',
    'folder_out': tt['folder']+'mos2_12x12/out/',
    'rho_bare': tt['folder']+'mos2_12x12/inp/rho_6pad12_coarse.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[3.186,
                     3.186,
                     24.0],
    'sc': [48,48,1],
    'fft_g':[720,720,225],
    'sym_arr': sym_arr
}


mos2_12to24 = {
    'prefix': 'mos2',
    'folder': tt['folder']+'mos2/12x12/',
    'folder_chi': tt['folder']+'mos2_12x12/chi/',
    'folder_chi_doped': tt['folder']+'mos2_12x12/chi_fermi_0.3/',
    'folder_dft': tt['folder']+'mos2_12x12/dft/',
    'folder_inp': tt['folder']+'mos2_12x12/inp/',
    'folder_out': tt['folder']+'mos2_12x12/out/',
    'rho_bare': tt['folder']+'mos2_12x12/inp/rho_6pad12_coarse.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[3.186,
                     3.186,
                     24.0],
    'sc': [24,24,1],
    'fft_g':[720,720,225],
    'sym_arr': sym_arr
}


mos2_12to24 = {
    'prefix': 'mos2',
    'folder': tt['folder']+'mos2/12x12/',
    'folder_chi': tt['folder']+'mos2_12x12/chi/',
    'folder_chi_doped': tt['folder']+'mos2_12x12/chi_fermi_0.3/',
    'folder_dft': tt['folder']+'mos2_12x12/dft/',
    'folder_inp': tt['folder']+'mos2_12x12/inp/',
    'folder_out': tt['folder']+'mos2_12x12/out/',
    'rho_bare': tt['folder']+'mos2_12x12/inp/rho_6pad12_coarse.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[3.186,
                     3.186,
                     24.0],
    'sc': [24,24,1],
    'fft_g':[720,720,225],
    'sym_arr': sym_arr
}


doped_chi_test = {
    'prefix': 'mos2',
    'folder': tt['folder']+'mos2/12x12/',
    'folder_chi': tt['folder']+'mos2_12x12/chi/',
    'folder_chi_doped': tt['folder']+'mos2_12x12/chi_fermi_0.3/',
    'folder_dft': tt['folder']+'mos2_12x12/dft/',
    'folder_inp': tt['folder']+'mos2_12x12/inp/',
    'folder_out': tt['folder']+'mos2_12x12/out/',
    'rho_bare': tt['folder']+'mos2_12x12/inp/rho_6pad12_coarse.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[3.186,
                     3.186,
                     24.0],
    'sc': [24,24,1],
    'fft_g':[180,180,120]
}


bn_3x3 = {
    'prefix': 'bn',
    'folder': tt['folder']+'bn/dft/3x3/',
    'folder_dft': tt['folder']+'bn/dft/3x3/',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[2.51,
                     2.51,
                     25],
    'sc': [3,3,1],
    'fft_g':[55,55,181]
}

bn_6x6 = {
    'prefix': 'bn',
    'folder': tt['folder']+'bn/6x6/',
    'folder_dft': tt['folder']+'bn/dft/6x6/',
    'folder_chi': tt['folder']+'bn/chi_6x6/',
    'folder_out': tt['folder']+'bn/out/6x6',
    'rho_tot': tt['folder']+'bn/dft/6x6/drho.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[2.51,
                     2.51,
                     25],
    'sc': [6,6,1],
   # 'fft_g':[109,109,181]
    'fft_g':[108,108,180]
}

bn_6x6_2d = {
    'prefix': 'bn_2d',
    'folder': tt['folder']+'bn/6x6_2d/',
    'folder_dft': tt['folder']+'bn/dft/6x6/',
    'folder_chi': tt['folder']+'bn/chi_6x6/',
    'folder_out': tt['folder']+'bn/out/6x6',
    'rho_tot': tt['folder']+'bn/dft/6x6_2d/drho.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[2.51,
                     2.51,
                     25],
    'sc': [6,6,1],
   # 'fft_g':[109,109,181]
    'fft_g':[108,108,180]
}



bn_7x7 = {
    'prefix': 'bn',
    'folder': tt['folder']+'bn/7x7/',
    'folder_dft': tt['folder']+'bn/dft/7x7/',
    'folder_chi': tt['folder']+'bn/chi_6x6/',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[2.51,
                     2.51,
                     25],
    'sc': [7,7,1],
    'fft_g':[128,128,180]
}

bn_8x8 = {
    'prefix': 'bn',
    'folder': tt['folder']+'bn/8x8/',
    'folder_dft': tt['folder']+'bn/dft/8x8/',
    'folder_chi': tt['folder']+'bn/chi_6x6/',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[2.51,
                     2.51,
                     25],
    'sc': [8,8,1],
    'fft_g':[144,144,180]
}

bn_9x9 = {
    'prefix': 'bn',
    'folder': tt['folder']+'bn/9x9/',
    'folder_dft': tt['folder']+'bn/dft/9x9/',
    'folder_chi': tt['folder']+'bn/chi_6x6/',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[2.51,
                     2.51,
                     25],
    'sc': [9,9,1],
    'fft_g':[162,162,180]
}

bn_10x10 = {
    'prefix': 'bn',
    'folder': tt['folder']+'bn/10x10/',
    'folder_dft': tt['folder']+'bn/dft/10x10/',
    'folder_chi': tt['folder']+'bn/chi_6x6/',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[2.51,
                     2.51,
                     25],
    'sc': [10,10,1],
    'fft_g':[180,180,180]
}

bn_11x11 = {
    'prefix': 'bn',
    'folder': tt['folder']+'bn/11x11/',
    'folder_dft': tt['folder']+'bn/dft/11x11/',
    'folder_chi': tt['folder']+'bn/chi_6x6/',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[2.51,
                     2.51,
                     25],
    'sc': [11,11,1],
    'fft_g':[200,200,180]
}

bn_12x12 = {
    'prefix': 'bn',
    'folder': tt['folder']+'bn/12x12/',
    'folder_dft': tt['folder']+'bn/dft/12x12/',
    'folder_chi': tt['folder']+'bn/chi_12x12/',
    'folder_out': tt['folder']+'bn/out/12x12',
    'rho_tot': tt['folder']+'bn/dft/12x12/drho.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[2.51,
                     2.51,
                     25],
    'sc': [12,12,1],
    'fft_g':[216,216,180],
    'folder_G_info':tt['folder']+'bn/G_info/',
    'folder_model':tt['folder']+'bn/model_interp/'
}
bn_13x13 = {
    'prefix': 'bn',
    'folder': tt['folder']+'bn/13x13/',
    'folder_dft': tt['folder']+'bn/dft/12x12/',
    'folder_chi': tt['folder']+'bn/chi_12x12/',
    'folder_out': tt['folder']+'bn/out/12x12',
    'rho_tot': tt['folder']+'bn/dft/12x12/drho.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[2.51,
                     2.51,
                     25],
    'sc': [13,13,1],
    'fft_g':[240,240,180]
}

bn_14x14 = {
    'prefix': 'bn',
    'folder': tt['folder']+'bn/14x14/',
    'folder_dft': tt['folder']+'bn/dft/12x12/',
    'folder_chi': tt['folder']+'bn/chi_12x12/',
    'folder_out': tt['folder']+'bn/out/12x12',
    'rho_tot': tt['folder']+'bn/dft/12x12/drho.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[2.51,
                     2.51,
                     25],
    'sc': [14,14,1],
    'fft_g':[256,256,180]
}

bn_15x15 = {
    'prefix': 'bn',
    'folder': tt['folder']+'bn/15x15/',
    'folder_dft': tt['folder']+'bn/dft/15x15/',
    'folder_chi': tt['folder']+'bn/chi_12x12/',
    'folder_out': tt['folder']+'bn/out/12x12',
    'rho_tot': tt['folder']+'bn/dft/12x12/drho.xsf',
    'is2d': True,
    'ibrav': 4,
    'lattpara_unit':[2.51,
                     2.51,
                     25],
    'sc': [15,15,1],
    'fft_g':[270,270,180]
}


Na_6x6x6 = {
    'prefix': 'Na',
    'folder': tt['folder']+'Na/6x6x6/',
    'folder_dft': tt['folder']+'bn/dft/15x15/',
    'folder_chi': tt['folder']+'bn/chi_12x12/',
    'folder_out': tt['folder']+'bn/out/12x12',
    'rho_tot': tt['folder']+'bn/dft/12x12/drho.xsf',
    'is2d': True,
    'ibrav': 1,
    'lattpara_unit':[3.5607451098,
                     3.5607451098,
                     3.5607451098],
    'sc': [6,6,6],
    'fft_g':[192,192,192],
    'sym_arr':sym_arr
        }

Na_6x6x6_fine = {
    'prefix': 'Na',
    'folder': tt['folder']+'Na/6x6x6/',
    'folder_dft': tt['folder']+'bn/dft/15x15/',
    'folder_chi': tt['folder']+'bn/chi_12x12/',
    'folder_out': tt['folder']+'bn/out/12x12',
    'rho_tot': tt['folder']+'bn/dft/12x12/drho.xsf',
    'is2d': True,
    'ibrav': 1,
    'lattpara_unit':[3.5607451098,
                     3.5607451098,
                     3.5607451098],
    'sc': [6,6,6],
    'fft_g':[384,384,384]

        }

Na_7x7x7 = {
    'prefix': 'Na',
    'folder': tt['folder']+'Na/7x7x7/',
    'folder_dft': tt['folder']+'bn/dft/15x15/',
    'folder_chi': tt['folder']+'bn/chi_12x12/',
    'folder_out': tt['folder']+'bn/out/12x12',
    'rho_tot': tt['folder']+'bn/dft/12x12/drho.xsf',
    'is2d': True,
    'ibrav': 1,
    'lattpara_unit':[3.5607451098,
                     3.5607451098,
                     3.5607451098],
    'sc': [7,7,7],
    'fft_g':[216,216,216]

        }
