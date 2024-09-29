import numpy as np
import scipy.sparse as sp 
import constants 
import re 

folder_prefix='./'

def printer_matrix(m, output_file="matrix_output"):
    fw = open(folder_prefix+output_file, 'w')
    m = np.array(m) 
    if len(np.shape(m)) == 2:
        if len(m)<len(m[0]):
            m = np.transpose(m)
        for ik in range(len(m)):
            for ibnd in range(len(m[0])):
                fw.write("%20.12f "%(m[ik][ibnd]))

            fw.write("\n")
    elif len(np.shape(m)) == 1:
        for ik in range(len(m)):
            fw.write("%20.12f \n"%m[ik])

    fw.close() 

def reader_matrix():
    fr = open(folder_prefix+"matrix_output", "r")
    lines = fr.readlines()
    n1 = len(lines)
    n2 = len(lines[0].split())
    if n2 == 0:
        return []
    elif n2 == 1:
        m = np.zeros((n1, ))
        for i in range(n1):
            m[i] = float(lines[i])
    else:
        m = np.zeros((n1, n2))
        for i in range(n1):
            for j in range(n2):
                m[i][j] = float(lines[i].split()[j])
    
    return m 

def reader_info(FILENAME=folder_prefix+'info', ndim=2):
    fo = open(FILENAME, 'r')
    line = fo.readline()
    while line:
        line = line.split()
        if not line:
            line = fo.readline()
            continue

        if line[1] == 'symmetry':
            nsym = int(line[0])
            sym_matrix = np.zeros((nsym, 3, 3))
            for i in range(nsym//6):
                line = fo.readline()
                line = fo.readline().split()
                sym_matrix[i*6:(i+1)*6, 0] = np.reshape([int(s) for s in line], (6,3))
                line = fo.readline().split()
                sym_matrix[i*6:(i+1)*6, 1] = np.reshape([int(s) for s in line], (6,3))
                line = fo.readline().split()
                sym_matrix[i*6:(i+1)*6, 2] = np.reshape([int(s) for s in line], (6,3))

            if nsym%6 != 0:
                tmp = nsym%6
                line = fo.readline()
                line = fo.readline().split()
                sym_matrix[i*6:i*6+tmp, 0] = np.reshape([int(s) for s in line], (tmp,3))
                line = fo.readline().split()
                sym_matrix[i*6:i*6+tmp, 1] = np.reshape([int(s) for s in line], (tmp,3))
                line = fo.readline().split()
                sym_matrix[i*6:i*6+tmp, 2] = np.reshape([int(s) for s in line], (tmp,3))

        line = fo.readline()

    fo.close()

    if ndim == 2:
        sym_matrix = sym_matrix[:,0:2,0:2]

    return sym_matrix 

def reader_freq(FILENAME=folder_prefix+'freq.gp'):
    fo = open(FILENAME, 'r')
    lines = fo.readlines()

    nq = len(lines)
    nmode = len(lines[0].split()) - 1

    freq = np.zeros((nq, nmode))

    for iq in range(nq):
        for im in range(nmode):
            freq[iq][im] = abs(float(lines[iq].split()[im+1]))

    return freq*0.12398*0.001 
    # cm-1 to eV 

def reader_velocity(FILENAME=folder_prefix+'tt_geninterp.dat', nelec=-1, ibndlist=[0,1]):
    fo = open(FILENAME, 'r')
    lines = fo.readlines()
    
    nbnd = 0
    while int(lines[nbnd+3].split()[0]) == 1:
        nbnd += 1 

    if nelec < 0.0:
        nelec = nbnd 
        ibndlist = np.array([nelec-1], dtype=int)
    else:
        ibndlist = np.array(ibndlist, dtype=int) + nelec - 1 

    nbndout = len(ibndlist) 

    nk = int(lines[-1].split()[0])

    kvec = np.zeros((nk, 3))
    velocity = np.zeros((nk, nbndout, 3))
    bande = np.zeros((nk, nbndout))

    for ik in range(nk):
        for ibnd in range(nbnd):
            line = lines[ik*nbnd+ibnd+3]

            if ibnd in ibndlist:
                line = line.split()
                lowest_bnd = ibndlist[0]

                if (ibnd == 0):
                    kvec[ik][0] = float(line[1])
                    kvec[ik][1] = float(line[2])
                    kvec[ik][2] = float(line[3])

                bande[ik][ibnd-lowest_bnd] = float(line[4])

                velocity[ik][ibnd-lowest_bnd][0] = float(line[5])
                velocity[ik][ibnd-lowest_bnd][1] = float(line[6])
                velocity[ik][ibnd-lowest_bnd][2] = float(line[7])

    fo.close()

    return kvec, bande, velocity

def writer_phin(iqcx_list, folder=folder_prefix, ninput=1):
    nq = len(iqcx_list)
    for i in range(ninput):
        filename = folder+'ph_q'+str(i+1)+'.in'
        fo = open(filename, 'w')

        slice_len = nq//ninput + 1
        if i == 0:
            iqcx_slice = iqcx_list[0:slice_len] # exclude q==0 from dfpt g2 computation
        elif i == ninput - 1:
            iqcx_slice = iqcx_list[i*slice_len:]
        else:
            iqcx_slice = iqcx_list[i*slice_len:(i+1)*slice_len]
        slice_len = len(iqcx_slice)

        fo.write("--\n")
        fo.write("&inputph\n")
        fo.write("  prefix   = 'test'\n")
        fo.write("  tr2_ph   =  1.0d-23\n")
        fo.write("  alpha_mix(1) = 0.7\n")
        fo.write("  fildyn   = 'test.dyn'\n")
        fo.write("  fildvscf = 'dvscf'\n")
        if slice_len > 1:
            fo.write("  ldisp = .true.\n")
            fo.write("  qplot=.true.\n")
        fo.write(" /\n")

        if slice_len > 1:
            fo.write("%i\n"%(slice_len))
        for iqcx in iqcx_slice:
            fo.write("%20.12f %20.12f %20.12f 1 \n" %(iqcx[0], iqcx[1], 0.0))

def writer_elphin(iqcx_list, ikcx_list, folder=folder_prefix, ninput=1, nsubinput=1):
    nq = len(iqcx_list)
    nk = len(ikcx_list)
    for i in range(ninput):
        for j in range(nsubinput):
            slice_lenq = nq//ninput + 1
            if i == 0:
                iqcx_slice = iqcx_list[1:slice_lenq] # exclude q==0 from dfpt g2 computation
            elif i == ninput - 1:
                iqcx_slice = iqcx_list[i*slice_lenq:]
            else:
                iqcx_slice = iqcx_list[i*slice_lenq:(i+1)*slice_lenq]
            slice_lenq = len(iqcx_slice)

            slice_lenk = nk//nsubinput
            if j == nsubinput - 1:
                ikcx_slice = ikcx_list[j*slice_lenk:]
            else:
                ikcx_slice = ikcx_list[j*slice_lenk:(j+1)*slice_lenk]
            slice_lenk = len(ikcx_slice)

            filename = folder+'elph_q'+str(i+1)+'_s'+str(j+1)+'.in'
            fo = open(filename, 'w')
            
            fo.write("--\n")
            fo.write("&inputph\n")
            fo.write("  prefix   = 'test'\n")
            fo.write("  tr2_ph   =  1.0d-20\n")
            fo.write("  alpha_mix(1) = 0.7\n")
            fo.write("  fildyn   = 'test.dyn'\n")
            fo.write("  fildvscf = 'dvscf'\n")
            if slice_lenq > 1:
                fo.write("  ldisp = .true.\n")
                fo.write("  qplot=.true.\n")
            fo.write("  electron_phonon = 'chosen_ks'\n")
            fo.write("  trans = .false.\n")
            fo.write(" /\n")

            if slice_lenq > 1:
                fo.write("%i\n"%(slice_lenq))
            for iqcx in iqcx_slice:
                fo.write("%20.12f %20.12f %20.12f 1 \n" %(iqcx[0], iqcx[1], 0.0))

            fo.write("%i\n"%(slice_lenk))
            for ikcx in ikcx_slice:
                fo.write("%20.12f %20.12f %20.12f \n" %(ikcx[0], ikcx[1], 0.0))

            fo.close() 



def reader_elphmat(nkf1, nkf2, nqf1, nqf2, reci_vec, FILENAME=folder_prefix+"elphmat1", ibndlist=None):
    """ It's quite complicated for chosing the right ibndlist. e.g. if ibndlist =[0], only the lowest
        band in elphmat file will be considered. """

    from qirr_sym import get_kindex2 

    fo = open(FILENAME, 'r')
    
    line = fo.readline()
    line = fo.readline()
    nmode = int(line.split()[2])*3
    line = fo.readline().split()
    ibmin_f90 = int(line[3])
    ibmax_f90 = int(line[4])
    nbnd = ibmax_f90 - ibmin_f90 + 1

    if ibndlist is not None:
        if len(ibndlist) > nbnd:
            raise Exception("Requested band num is larger than bands num from elphmat files.")
        else:
            nbnd = len(ibndlist)
            ibmin_f90 += np.min(ibndlist)
            ibmax_f90 = ibmin_f90 + len(ibndlist) - 1

    elphmat = np.zeros((nbnd, nbnd, nmode), dtype=object)
    for ibnd in range(nbnd):
        for jbnd in range(nbnd):
            for im in range(nmode):
                elphmat[ibnd,jbnd,im] = sp.lil_matrix((nkf1*nkf2, nqf1*nqf2))

    reci_vec_inv = np.linalg.inv(reci_vec)

    fo.seek(0,0)
    line = fo.readline()
    while line:
        words = line.split()

        if len(words) >= 1 and words[0] == "q":
            qcx = [float(words[2]), float(words[3])]
            qfx = np.dot(qcx, reci_vec_inv)
            now_iq = get_kindex2([qfx], nqf1, nqf2)[0]

        if len(words) >= 2 and words[1] == "elph":
            kcx = [float(words[4]), float(words[5])]
            kfx = np.dot(kcx, reci_vec_inv)
            now_ik = get_kindex2([kfx], nkf1, nkf2)[0]
        
        if len(words) >= 2 and words[1] == "mode":
            now_mode = int(words[3])
            line = fo.readline()
            words = line.split()
            now_ibnd = int(words[0])
            now_jbnd = int(words[1])
            line = fo.readline()
            words = line.split()
            if now_ibnd >= ibmin_f90 and now_ibnd <= ibmax_f90 and now_jbnd >= ibmin_f90 and now_jbnd <= ibmax_f90:
                if float(words[2]) > 0:
                    elphmat[now_ibnd-ibmin_f90, now_jbnd-ibmin_f90, now_mode-1][now_ik, now_iq] = (float(words[2])*0.001)**2

        line = fo.readline() 
    return elphmat 
    # in eV*eV unit

def reader_modes(nqf1, nqf2, reci_vec, FILENAME=folder_prefix+"matdyn.modes", proj_mode=1, trgm=None):
    """ read the eigenvectors of the dynamic matrix, proj_mode = 1, use projection mode-resolved
        method; if proj_mode = 2, use continuity mode-resolved method. """

    from qirr_sym import get_kindex2

    reci_vec_inv = np.linalg.inv(reci_vec)
    trgm = np.array(trgm, dtype=complex)
    trgm = trgm/np.linalg.norm(trgm)

    with open(FILENAME, 'r') as fo:
        lines = fo.readlines()
    
    iline = 0
    freq_iline = []
    while len(freq_iline) < 2:
        if " freq " in lines[iline]:
            freq_iline.append(iline)
        iline += 1
    natm = freq_iline[1] - freq_iline[0] - 1
    nmode = 3*natm 

    qlist = []
    eignlist = []
    iline = 0
    while iline < len(lines):
        line = lines[iline]
        if 'q =' in line:
            q = [float(x) for x in re.findall('[- ][0-9].[0-9]*',line)]
            qcx = np.array(q[0:2]) 
            qfx = np.dot(qcx, reci_vec_inv) 
            qlist.append(qfx)  

            iline += 2
            eign = np.zeros((nmode, nmode), dtype=complex)
            for im in range(nmode):
                iline += 1
                for iatm in range(natm):
                    for j in [0,1,2]:
                        Re, Im = [float(x) for x in lines[iline].split()[2*j+1:2*j+3]]
                        eign[im,iatm*3+j] = complex(Re, Im)

                    iline += 1

            eignlist.append(np.copy(eign)) 

        iline += 1

    nqbz = len(qlist) 
    iqbz = get_kindex2(qlist, nqf1, nqf2)
    modes_weight = sp.lil_matrix((nqf1*nqf2, nmode)) 

    for i in range(nqbz):
        proj2 = get_proj2(eignlist[i], trgm, proj_mode)
        if abs(np.sum(proj2)-1.0)>1e-3:
            print("Somthing wrrong: sum of projection is larger than 1:")
        modes_weight[iqbz[i],:] = np.copy(proj2) 

    return modes_weight 

def get_proj2(eigns, trgm, proj_mode=1):
    if proj_mode == 1:
        proj2 = np.abs([np.dot(trgm, v) for v in eigns])**2
    elif proj_mode == 2:
        proj2 = np.zeros((len(trgm),))
        proj2[np.argmax(np.abs([np.dot(trgm, v) for v in eigns]))] = 1.0 

    return proj2 

def writer_kptdat(kfx_list, FILENAME=folder_prefix+"kpt.dat"):
    fo = open(FILENAME, 'w')

    fo.write("%5i  crystal\n" %(len(kfx_list)))
    for kfx in kfx_list:
        fo.write("%20.12f %20.12f %20.12f 1.0\n"%(kfx[0], kfx[1], 0.0))

    fo.close()

def reader_g2matrix(ik_list, nbnd, nmode, nkf1, nkf2, nqf1, nqf2, refer_w=None, FILENAME=folder_prefix+"fort.708"):
    # nqtot = nqf1*nqf2 

    fo = open(FILENAME, 'r')

    g2_mat = np.zeros((nbnd, nbnd, nmode), dtype=object)
    w_mat = np.zeros((nbnd, nbnd, nmode), dtype=object)
    for ibnd in range(nbnd):
        for jbnd in range(nbnd):
            for im in range(nmode):
                g2_mat[ibnd,jbnd,im] = sp.lil_matrix((nkf1*nkf2, nqf1*nqf2))
                w_mat[ibnd,jbnd,im] = sp.lil_matrix((nkf1*nkf2, nqf1*nqf2))

    iq = 0
    lines = fo.readlines() 
    for i in range(len(lines)):
        words = lines[i].split()
        if len(words) > 0:
            lastq = iq 
            iq = int(words[0]) - 1
            ikbz = ik_list[int(words[1]) - 1]
            ijbnd = int(words[2]) - 1
            ibnd = ijbnd//nbnd 
            jbnd = ijbnd%nbnd 
            im = int(words[3]) - 1

            if lastq != iq and (iq+1)%100==0:
                print("current iq:", iq+1, "/", nqf1*nqf2) 

            if (refer_w is not None) and (refer_w[ibnd,jbnd,im][ikbz,iq] <= 0.0):
                continue 

            g2_mat[ibnd,jbnd,im][ikbz,iq] = float(words[4])*unit2.ry2ev**2
            w_mat[ibnd,jbnd,im][ikbz,iq] = float(words[5])/unit2.ry2ev 

    fo.close()

    return g2_mat, w_mat 

def writer_scat_rate(scat, bande, kf_ibz, folder=folder_prefix):
    if sp.issparse(scat):
        nbnd = scat.get_shape()[0]

        fo = open(folder+"scattering.dat", 'w')

        for ibnd in range(nbnd):
            for indk, ikf in enumerate(kf_ibz):
                fo.write("%5i %5i %20.12f %20.12e \n" %(ibnd, indk, bande[ikf,ibnd], scat[ibnd,ikf]))

        fo.close()
    else:
        nmode = len(scat)
        nbnd = scat[0].get_shape()[0]

        fo = open(folder+"scattering_imode.dat", "w")

        for ibnd in range(nbnd):
            for indk, ikf in enumerate(kf_ibz):
                fo.write("%5i %5i %20.12f " %(ibnd, indk, bande[ikf,ibnd])) 
                for im in range(nmode):
                    fo.write("%20.12e " %(scat[im][ibnd,ikf]))
                fo.write("\n")

        fo.close() 

def writer_freq(freq, FILENAME=folder_prefix+'new_freq.gp'):
    _freq = freq /0.12398/0.001 # ev2cm-1
    nqtot = len(freq)
    nmode = len(freq[0])

    fo = open(FILENAME, "w")

    for iq in range(nqtot):
        fo.write("%14.8f " % (0.0))
        for im in range(nmode):
            fo.write("%20.12e" % (_freq[iq,im]))

        fo.write("\n")

    fo.close()
