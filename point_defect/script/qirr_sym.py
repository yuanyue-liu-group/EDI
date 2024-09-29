import numpy as np 
import time 
import scipy.sparse as sp 
import io_files
from math import floor, pi, exp 
from math import isinf, isnan 

def get_kindex1(kfrac_list, nqf1, nqf2):
    ik_list = np.zeros((len(kfrac_list),), dtype=int) 
    for ik in range(len(kfrac_list)):
        ikx = kfrac_list[ik][0] - floor(kfrac_list[ik][0])
        iky = kfrac_list[ik][1] - floor(kfrac_list[ik][1])
        ikx = int(floor(kfrac_list[ik][0] * nqf1))%nqf1
        iky = int(floor(kfrac_list[ik][1] * nqf2))%nqf2
        ik_list[ik] = iky + ikx * nqf2

    return ik_list

def get_kindex2(kfrac_list, nqf1, nqf2):
    ik_list = np.zeros((len(kfrac_list),), dtype=int)
    for ik in range(len(kfrac_list)):
        ikx = kfrac_list[ik][0] - floor(kfrac_list[ik][0])
        iky = kfrac_list[ik][1] - floor(kfrac_list[ik][1])
        ikx = int(round(kfrac_list[ik][0] * nqf1))%nqf1
        iky = int(round(kfrac_list[ik][1] * nqf2))%nqf2
        ik_list[ik] = iky + ikx * nqf2

    return ik_list

def kq_distance(kfx, qfx, reci_vec):
    qfx = np.array(qfx)
    dis = 1e10
    for ishift in [-1.0, 0.0, 1.0]:
        for jshift in [-1.0, 0.0, 1.0]:
            kfx_shift = np.array([kfx[0]+ishift, kfx[1]+jshift])
            kqdis = np.linalg.norm( kfrac2kcart([kfx_shift-qfx], reci_vec)[0] )
            if kqdis < dis:
                dis = kqdis 

    return dis 



def gen_qirr(nqf1, nqf2):
    sym_mat = io_files.reader_info()
    nsym = len(sym_mat) 

    ibz2bz = []
    bz_sym = np.ones((nqf1*nqf2), dtype=int)*-1
    bz2ibz = np.ones((nqf1*nqf2), dtype=int)*-1
    for iq in range(nqf1*nqf2):
        if bz_sym[iq] >= 0:
            continue 

        ix = iq//nqf2
        iy = iq%nqf2
        qx_bz = [1.0*ix/nqf1, 1.0*iy/nqf2]
        qxsym_fkx = np.array([np.dot(sym_mat[i], qx_bz) for i in range(nsym)])
        qxsym_ik = get_kindex2(qxsym_fkx, nqf1, nqf2)

        bz_sym[qxsym_ik] = np.arange(nsym)
        bz2ibz[qxsym_ik] = len(ibz2bz)
        ibz2bz.append(iq)

    ibz2bz = np.array(ibz2bz)

    if -1 in bz_sym:
        print("Problems here!!!") 
    if -1 in bz2ibz:
        print("Problems here!!!")
        
    return ibz2bz, bz2ibz, bz_sym 

def get_kfrac(ik_list, nqf1, nqf2):
    kcart_list = np.zeros((len(ik_list),2))
    for ik in range(len(ik_list)):
        kcart_list[ik][0] = float(ik_list[ik]//nqf2)/nqf1 
        kcart_list[ik][1] = float(ik_list[ik]%nqf2)/nqf2 

    return kcart_list

def kfrac2kcart(kfrac, reci_vec):
    kcart = np.array([np.dot(ikf, reci_vec) for ikf in kfrac])
    return kcart 

def get_kcart(ik_list, nqf1, nqf2, reci_vec):
    kfrac = get_kfrac(ik_list, nqf1, nqf2)
    kcart = kfrac2kcart(kfrac, reci_vec)

    return kcart 

def gen_qirr_wrap(nqf1, nqf2, reci_vec, cutoff=-1, iqbz=None):
    ibz2bz, _, _ = gen_qirr(nqf1, nqf2)

    if iqbz is None:
        iqibz = ibz2bz 
    else:
        iqibz = qbz2ibz(iqbz, nqf1, nqf2)
    
    kfrac = get_kfrac(iqibz, nqf1, nqf2)
    kcart = [np.dot(ifrac, reci_vec) for ifrac in kfrac]

    qirr = []
    qirr_fx = []
    qirr_cx = []
    if cutoff > 0.0:
        for i in range(len(kcart)):
            if (kcart[i][0]**2 + kcart[i][1]**2)**0.5 < cutoff:
                qirr.append(iqibz[i])
                qirr_fx.append(kfrac[i])
                qirr_cx.append(kcart[i])
    else:
        qirr = iqibz
        qirr_fx = kfrac
        qirr_cx = kcart

    qirr = np.array(qirr, dtype=int)
    qirr_fx = np.array(qirr_fx)
    qirr_cx = np.array(qirr_cx)

    return qirr, qirr_fx, qirr_cx 

def gen_kirr(nqf1, nqf2, reci_vec, bande, ecutoff=0.3, bz=False, efermi=0.0):
    if bz:
        ikibz = np.arange(nqf1*nqf2)
    else:
        ikibz, _, _ = gen_qirr(nqf1, nqf2)

    kfrac = get_kfrac(ikibz, nqf1, nqf2)
    kcart = [np.dot(ifrac, reci_vec) for ifrac in kfrac]

    kirr = []
    kirr_fx = []
    kirr_cx = []
    
    for i in range(len(ikibz)):
        if np.min(np.abs(bande[ikibz[i],:]-efermi)) < abs(ecutoff):
            kirr.append(ikibz[i])
            kirr_fx.append(kfrac[i])
            kirr_cx.append(kcart[i])
    
    kirr = np.array(kirr, dtype=int)
    kirr_fx = np.array(kirr_fx)
    kirr_cx = np.array(kirr_cx)

    return kirr, kirr_fx, kirr_cx 




def qbz2ibz(qbz, nqf1, nqf2):
    ibz2bz, bz2ibz, _ = gen_qirr(nqf1, nqf2)

    qibz = [ibz2bz[bz2ibz[iq]] for iq in qbz]
    qibz = list(set(qibz))

    return qibz 

def qibz2bz(qibz, nqf1, nqf2, output_mode=1):
    ibz2bz, bz2ibz, _ = gen_qirr(nqf1, nqf2)

    qibz = qbz2ibz(qibz, nqf1, nqf2)

    qbzs = []
    for i in range(len(qibz)):
        tmp = np.array([ibz2bz[bz2ibz[iq]] for iq in range(nqf1*nqf2)])
        qbzs.append(np.argwhere(tmp==qibz[i])[:,0])

    if output_mode == 1:
        return qbzs 
    elif output_mode == 2:
        qbz = []
        for qbz_slice in qbzs:
            qbz = qbz + list(qbz_slice)
        qbz = list(set(qbz))
        return np.sort(qbz)
    else:
        return qbzs 


def kindex_add(ik, iq, nqf1, nqf2=None):
    if nqf2 is None:
        nqf2 = nqf1 
    ikqx = (ik//nqf2 + iq//nqf2)%nqf1
    ikqy = (ik%nqf2 + iq%nqf2)%nqf2

    return ikqx*nqf2 + ikqy 

def get_kq_weight_mat(ikbz_list, freq, bande, nqf, efermi=-0.1, interp_w=False, za_qcut=0, reci_vec=None,triangular_wt=True):
    nmode = len(freq[0])
    nbnd = len(bande[0])
    nk = len(ikbz_list)

    kq_wmat = np.zeros((nbnd, nbnd, nmode), dtype=object)
    for ibnd in range(nbnd):
        for jbnd in range(nbnd):
            for im in range(nmode):
                kq_wmat[ibnd,jbnd,im] = sp.lil_matrix((nqf*nqf, nqf*nqf))
                
    for i in range(nk):
        ik = ikbz_list[i]
        for ibnd in range(nbnd):
            weight_iq = get_qweight_wrap(ik, ibnd, freq, bande, nqf, efermi, interp_w, za_qcut, reci_vec,triangular_wt)
            for jbnd in range(nbnd):
                for im in range(nmode):
                    kq_wmat[ibnd, jbnd, im][ik] = np.copy( weight_iq[:, im, jbnd] )
        fw = open("stdout", 'a+')
        fw.write("%i done!: %i / %i \n" %(ik, i, nk))
        fw.close() 
        print("%i done!: %i / %i " %(ik, i, nk))

    return kq_wmat 


def get_kq_weight_mat_mp(ikbz_list, freq, bande, nqf, efermi=-0.1, stdout=False, cores=1, interp_w=False, za_qcut=0, reci_vec=None,triangular_wt=True):
    import multiprocessing as mp 
    from functools import partial 

    if stdout:
        time_start = time.time() 
        print("kq_weight starts at: ", time_start) 

    nmode = len(freq[0])
    nbnd = len(bande[0])
    nk = len(ikbz_list)

    kq_wmat = np.zeros((nbnd, nbnd, nmode), dtype=object)
    for ibnd in range(nbnd):
        for jbnd in range(nbnd):
            for im in range(nmode):
                # kq_wmat[ibnd,jbnd,im] = sp.lil_matrix((nqf*nqf, nqf*nqf))
                kq_wmat[ibnd,jbnd,im] = sp.csc_matrix((nqf*nqf, nqf*nqf))

    # divide ikbz_list into slices
    slice_len = nk//cores + 1
    slices = []
    for i in range(cores):
        if i == cores - 1:
            slices.append(ikbz_list[i*slice_len:])
        else:
            slices.append(ikbz_list[i*slice_len:(i+1)*slice_len])

    # ftmp = lambda slice_ks: get_kq_weight_mat(slice_ks, freq, bande, nqf, efermi)
    ftmp = partial(get_kq_weight_mat, freq=freq, bande=bande, nqf=nqf, efermi=efermi, interp_w=interp_w, za_qcut=za_qcut, reci_vec=reci_vec)
    pool = mp.Pool(processes=cores)

    pools_count = 0
    for wmat_tmp in pool.imap_unordered(ftmp, slices):
        for ibnd in range(nbnd):
            for jbnd in range(nbnd):
                for im in range(nmode):
                    kq_wmat[ibnd,jbnd,im] = kq_wmat[ibnd,jbnd,im] + wmat_tmp[ibnd,jbnd,im].tocsr() 
        
        if stdout:
            print("Pools done: %i / %i" %(pools_count+1, cores))
        pools_count += 1
    
    for ibnd in range(nbnd):
        for jbnd in range(nbnd):
            for im in range(nmode):
                kq_wmat[ibnd,jbnd,im] = kq_wmat[ibnd,jbnd,im].tolil()

    if stdout:
        time_end = time.time() 
        print("kq_weight ends at: ", time_end) 
        print("kq_weight calculation time cost: ", time_end - time_start, 's') 

    return kq_wmat 


def get_qweight_wrap(ik, ibnd, freq, bande, nqf, efermi=0.0, interp_w=False, za_qcut=0, reci_vec=None,triangular_wt=True):
    # if efermi < 0.0:
    #     bande = bande - np.min(bande) - efermi 
    # else:
    #     bande = -1.0*(bande - np.max(bande)) + efermi 

    if isinstance(efermi, float):
        _bande = bande - efermi 

    nmode = len(freq[0])
    nqtot = len(freq)
    nbnd = len(bande[0])

    kT = 300 * 8.621738e-5
    ikibe = _bande[ik,ibnd]
    freq_cut = np.max(freq) + 100.0/nqf
    # here we apply the frequency cutoff to enhance efficiency

    weight_iq = np.zeros((nqtot, nmode, nbnd))
    # weight_iq = np.zeros((nmode, nbnd), dtype=object)
    # for im in range(nmode):
    #     for ibnd in range(nbnd):
    #         weight_iq[im,ibnd] = sp.lil_matrix((nqtot,1))
            
    for iq in range(nqtot):
        if interp_w:
            qfx = get_kfrac([iq], nqf, nqf)[0]
            distmp = kq_distance([0.0, 0.0], qfx, reci_vec)
            if distmp > za_qcut:
                continue

        for im in range(nmode):
            ikq = kindex_add(ik, iq, nqf)
            tmpfreq = freq[iq][im]
            if tmpfreq > 1e-50:
                wgq = 1./(exp(tmpfreq/kT)-1.0)
            else:
                wgq = 0.0 
                # print("cutoff at iq, im:", iq, im)
            for jbnd in range(nbnd):
              if triangular_wt:
                if abs(_bande[ikq,jbnd]-ikibe)>freq_cut:
                    continue 

                wgkq = 1./(exp(_bande[ikq, jbnd]/kT) + 1.0)
                F1 = 1.0 - wgkq + wgq 
                F2 = wgkq + wgq 
                weight_iq[iq][im][jbnd] = get_qweight(iq, ik, ibnd, jbnd, im, freq, _bande, nqf, F1, F2) 
              else:
                dE=bande[ik, ibnd] - bande[ikq, jbnd] 
                weight_iq[iq][im][jbnd] = exp(-(dE/degauss)**2/2)/((2*pi)**0.5*degauss)

    return weight_iq 

def get_vipq(ik, ibnd, freq, bande, nqf, efermi=-0.1):
    """ not used now """
    weight_iq = get_qweight_wrap(ik, ibnd, freq, bande, nqf, efermi)

    nqtot = len(freq)

    vipqw = np.array([np.sum(weight_iq[i,:,:]) for i in range(nqtot)])

    vipqi = []
    for iq in range(nqtot):
        if vipqw[iq] > 0.0:
            vipqi.append(iq)

    vipqi = np.array(vipqi, dtype=int) 
    vipkqi = np.array([kindex_add(ik, iq, nqf) for iq in vipqi])

    return vipqi, vipkqi

def get_vipq_wmat(wmat):
    nbnd = len(wmat)
    nmode = len(wmat[0,0])
    _, nqtot = wmat[0,0,0].get_shape()

    w_iqbz = np.zeros((nqtot,))
    for ibnd in range(nbnd):
        for jbnd in range(nbnd):
            for im in range(nmode):
                w_iqbz += np.array(wmat[ibnd,jbnd,im].sum(axis=0))[0]

    vipqi = []
    for iq in range(nqtot):
        if w_iqbz[iq] > 0.0:
            vipqi.append(iq)

    return vipqi


def get_qweight(iq, ik, ibnd, jbnd, imode, freq, bande, nqf, F1, F2):

    nest_ivec = np.array([(nqf-1)*nqf, (nqf-1)*nqf+1, 1, nqf, nqf*2-1, nqf-1], dtype=int)
    num_nest = len(nest_ivec) 
    iq_nest = np.zeros((num_nest, ), dtype=int)
    ikq_nest = np.zeros((num_nest, ), dtype=int)

    ikq = kindex_add(ik, iq, nqf)
    
    for i in range(num_nest):
        iq_nest[i] = kindex_add(iq, nest_ivec[i], nqf)
        ikq_nest[i] = kindex_add(ik, iq_nest[i], nqf)

    E0_f1 = bande[ik, ibnd] - bande[ikq, jbnd] - freq[iq, imode]
    E_nest_f1 = np.array([bande[ik, ibnd] - bande[ikq_nest[i], jbnd] - freq[iq_nest[i], imode] for i in range(num_nest)]) 

    E0_f2 = bande[ik, ibnd] - bande[ikq, jbnd] + freq[iq, imode]
    E_nest_f2 = np.array([bande[ik, ibnd] - bande[ikq_nest[i], jbnd] + freq[iq_nest[i], imode] for i in range(num_nest)]) 

    one_vertex = np.arange(num_nest)
    another_vertex = one_vertex + 1
    another_vertex[-1] = 0

    weight_f1 = 0.0
    if (min(E0_f1, np.min(E_nest_f1))<0.0 and max(E0_f1, np.max(E_nest_f1))>0.0):
        for i in range(num_nest):
            weight_f1 += get_triangular(E0_f1, E_nest_f1[one_vertex[i]], E_nest_f1[another_vertex[i]])
        weight_f1 = weight_f1*F1 

    weight_f2 = 0.0 
    if (min(E0_f2, np.min(E_nest_f2))<0.0 and max(E0_f2, np.max(E_nest_f2))>0.0):
        for i in range(num_nest):
            weight_f2 += get_triangular(E0_f2, E_nest_f2[one_vertex[i]], E_nest_f2[another_vertex[i]])
        weight_f2 = weight_f2*F2 

    weight = pi * (1./nqf**2) * (weight_f1 + weight_f2)

    return weight 

def get_triangular(ea, eb, ec):
    if ea >= 0.0 and eb >= 0.0 and ec >= 0.0:
        return 0.0
    elif ea <= 0.0 and eb <= 0.0 and ec <= 0.0:
        return 0.0 
    elif ea<eb and eb<ec:
        E0, E1, E2 = np.sort([ea, eb, ec])
        if eb>0.0:
            w_add = -1.*E0/(E1-E0)/(E2-E0)
            w_add = w_add*(2.0+E0/(E1-E0)+E0/(E2-E0))
        else:
            w_add = E2*E2/(E2-E1)/(E2-E0)/(E2-E0)
    elif ea<ec and ec<eb:
        E0, E1, E2 = np.sort([ea, eb, ec])
        if ec>0.0:
            w_add = -1.*E0/(E1-E0)/(E2-E0)
            w_add = w_add*(2.0+E0/(E1-E0)+E0/(E2-E0))
        else:
            w_add = E2*E2/(E2-E1)/(E2-E0)/(E2-E0)
    elif eb<ea and ea<ec:
        E0, E1, E2 = np.sort([ea, eb, ec])
        if ea>0.0:
            w_add = E0*E0/(E1-E0)/(E2-E0)/(E1-E0)
        else:
            w_add = E2*E2/(E2-E1)/(E2-E0)/(E2-E1)
    elif ec<ea and ea<eb:
        E0, E1, E2 = np.sort([ea, eb, ec])
        if ea>0.0:
            w_add = E0*E0/(E1-E0)/(E2-E0)/(E1-E0)
        else:
            w_add = E2*E2/(E2-E1)/(E2-E0)/(E2-E1)
    elif eb<ec and ec<ea:
        E0, E1, E2 = np.sort([ea, eb, ec])
        if ec>0.0:
            w_add = E0*E0/(E1-E0)/(E2-E0)/(E2-E0) 
        else:
            w_add = E2/(E2-E1)/(E2-E0)
            w_add = w_add*(2.0-E2/(E2-E1)-E2/(E2-E0))
    elif ec<eb and eb<ea:
        E0, E1, E2 = np.sort([ea, eb, ec])
        if eb>0.0:
            w_add = E0*E0/(E1-E0)/(E2-E0)/(E2-E0)
        else:
            w_add = E2/(E2-E1)/(E2-E0)
            w_add = w_add*(2.0-E2/(E2-E1)-E2/(E2-E0))
    else:
        return 0.0

    if isnan(w_add) or isinf(w_add):
        print("w_add error!!!")

    return w_add/2.0 
    

def elph_refolding_kirr(elphmat, nkf1, nkf2, nqf1, nqf2):
    sym_mat = io_files.reader_info()
    ibz2bz_k, _, _ = gen_qirr(nkf1, nkf2) 

    nbnd = len(elphmat)
    nmode = len(elphmat[0,0])
    nsym = len(sym_mat)

    elphmat_kirr = np.zeros((nbnd, nbnd, nmode), dtype=object)
    for ibnd in range(nbnd):
        for jbnd in range(nbnd):
            for im in range(nmode):
                elphmat_kirr[ibnd,jbnd,im] = sp.lil_matrix((nkf1*nkf2, nqf1*nqf2))

    for ibnd in range(nbnd):
        for jbnd in range(nbnd):
            for im in range(nmode):
                iks, iqs = elphmat[ibnd,jbnd,im].nonzero()
                for i in range(len(iks)):
                    ikbz = iks[i]
                    iqbz = iqs[i] 
                    gmat = elphmat[ibnd,jbnd,im][ikbz, iqbz]

                    ikbz_fx = get_kfrac([ikbz], nkf1, nkf2)[0]
                    ikbz_eq = get_kindex2( [np.dot(ism, ikbz_fx) for ism in sym_mat], nkf1, nkf2 )
                    # plot_kpoints(ikbz_eq, nkf1, nkf2, np.array([[0.304681, 0.175907],[0.000000, 0.351815]]))
                    for isym in range(nsym):
                        if ikbz_eq[isym] in ibz2bz_k:
                            iqbz_fx = get_kfrac([iqbz], nqf1, nqf2)[0]
                            iqbz_rot = get_kindex2( [np.dot(sym_mat[isym], iqbz_fx)], nqf1, nqf2 )[0]
                            elphmat_kirr[ibnd,jbnd,im][ikbz_eq[isym], iqbz_rot] = gmat 

    return elphmat_kirr 
