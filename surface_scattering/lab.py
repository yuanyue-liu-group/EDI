tt = {'fs0': 12.4073, 'folder': '/work/pwcond/m2_cu/', 'area0': 242.1124}
cu001 = {'fs0': 12.4181, 'folder': '/work/pwcond/m3_cu001/', 'area0': 161.4083}
cu110 = {'fs0': 12.4181, 'folder': '/work/pwcond/m4_cu110/', 'area0': 161.4083}

dat = tt
fs0 = dat['fs0']
folder = dat['folder']
area = dat['area0']

cu111_final = {
    'ibrav':4,
    'fs0': tt['fs0'], 
    'fsthick': 0.15, 
    'assume_metal': True, 
    'isibz': True,
    'filinfo': tt['folder']+'scf.out', 
    'filpwc': None, 
    'bc': (1,1), 
    'area': tt['area0'],
    'spsym': 'noz',
    'nqg':None,
    'ngd': [30,30,72],
    'filvel': tt['folder']+'wan/tt_geninterp.dat.k30kz72',
    'fillw': tt['folder']+'inp/l_k30kz72_noz.dat',
    'filg2m': tt['folder']+'inp/fort.k30kz72_noz.709',
    }
cu110_final = {
    'ibrav':62,
    'fs0': cu110['fs0'], 
    'fsthick': 0.15, 
    'assume_metal': True, 
    'isibz': True,
    'filinfo': cu110['folder']+'scf.out', 
    'filpwc': cu110['folder']+'pwc/cond.kf604060.out.all', 
    'bc': (1,1), 
    'area': cu110['area0'],
    'spsym': 'noz',
    'nqg':None,
    'ngd': [60,40,60],
    'filvel': cu110['folder']+'wan/tt_geninterp.kc969_kf604060.dat',
    'fillw': cu110['folder']+'inp/l_qc969_qf604060.dat',
    'filg2m': cu110['folder']+'inp/fort.qc969_qf604060.709',
    }
cu001_final = {
    'ibrav':6,
    'fs0': cu001['fs0'], 
    'fsthick': 0.15, 
    'assume_metal': True, 
    'isibz': True,
    'filinfo': cu001['folder']+'scf.out', 
    'filpwc': None, 
    'bc': (1,1), 
    'area': cu001['area0'],
    'spsym': 'noz',
    'nqg':None,
    'ngd': [30,30,20],
    'filvel': cu001['folder']+'wan/tt_geninterp.kc12128.kf303020.dat',
    'fillw': cu001['folder']+'inp/l_qc664_kc12128.dat',
    'filg2m': cu001['folder']+'inp/fort.qc664_kc12128.709',
    }

cu210 = {'fs0': 12.4302, 'folder': '/work/pwcond/m5_cu210/', 'area0': 807.0414}
cu210_final = {
    'ibrav':704,
    'fs0': cu210['fs0'], 
    'fsthick': 0.15, 
    'assume_metal': True, 
    'isibz': True,
    'filinfo': cu210['folder']+'scf.out', 
    'filpwc': cu210['folder']+'pwc/cond.out.kz60.all', 
    'bc': (1,1), 
    'area': cu210['area0'],
    'spsym': 'noz',
    'nqg':None,
    'ngd': [40,30,60],
    'filvel': cu210['folder']+'wan/Cu_geninterp.kz60.dat',
    'fillw': cu210['folder']+'inp/l.kz60.dat',
    'filg2m': cu210['folder']+'inp/fort.709',
    }
cu001sup2 = {'fs0': 12.4149, 'folder': '/work/pwcond/m6_cu100_sup2/', 'area0': 322.8165}
cu001sup2_tt1 = {
    'ibrav':705,
    'fs0': cu001sup2['fs0'], 
    'fsthick': 0.15, 
    'assume_metal': True, 
    'isibz': True,
    'filinfo': cu001sup2['folder']+'scf.out', 
    'filpwc': cu001sup2['folder']+'pwc/cond.out.all', 
    'bc': (1,1), 
    'area': cu001sup2['area0'],
    'spsym': 'noz',
    'nqg':None,
    'ngd': [60,30,40],
    'filvel': cu001sup2['folder']+'wan/Cu_geninterp.dat',
    'fillw': cu001sup2['folder']+'inp/l.dat',
    'filg2m': cu001sup2['folder']+'inp/fort.709',
    }
cu001sup2_final = {
    'ibrav':705,
    'fs0': cu001sup2['fs0'], 
    'fsthick': 0.15, 
    'assume_metal': True, 
    'isibz': True,
    'filinfo': cu001sup2['folder']+'scf.out', 
    'filpwc': cu001sup2['folder']+'pwc/cond.out.kz80.all', 
    'bc': (1,1), 
    'area': cu001sup2['area0'],
    'spsym': 'noz',
    'nqg':None,
    'ngd': [60,30,80],
    'filvel': cu001sup2['folder']+'wan/Cu_geninterp.kz80.dat',
    'fillw': cu001sup2['folder']+'inp/l.kz80.dat',
    'filg2m': cu001sup2['folder']+'inp/fort.709',
    }
