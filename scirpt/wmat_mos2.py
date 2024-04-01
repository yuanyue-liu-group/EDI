import sys 
import numpy as np 
import copy 
import scipy.sparse as sp 
import io_files 
import qirr_sym

RECI_VEC = np.array([[0.313963, 0.181267],[0.000000, 0.362533]])

if __name__ == '__main__':
   #input set up
    plot_mode = True
    ecut = 0.100 # energy cutoff for k points in eV unit
    ef=-0.05
    nbnd_valence=14
    nbnd_valence_w90=8
    ibndlist=[1,2]
    nqf = 144 # fine q grid
    nkf = 144 # fine k grid
    out_folder = "./kq144/" # output files folder 
    velocity_fn="w90_mos2_geninterp.dat"
    ncore=4
    

   #read E, V 
    _, bande, velocity = io_files.reader_velocity(out_folder+velocity_fn , nbnd_valence_w90,ibndlist) 


   #kp set up 
    fermi_level = np.min(bande)+ef

    kf_ibz, _, _ = qirr_sym.gen_kirr(nkf, nkf, RECI_VEC, bande, ecut,   False, efermi=fermi_level) # generate irr k points within ecut 
    kf_bz = qirr_sym.qibz2bz(kf_ibz, nkf, nkf, output_mode=2)

    kf_ibz_fx = qirr_sym.get_kfrac(kf_ibz, nkf, nkf)
    kf_bz_fx = qirr_sym.get_kfrac(kf_bz, nkf, nkf)

    np.savez(out_folder+"kpoints.npz", kf_ibz=kf_ibz, kf_bz=kf_bz,    kf_ibz_fx=kf_ibz_fx, kf_bz_fx=kf_bz_fx)

   #calculate wt 
    pseudo_freq = np.zeros((nqf*nqf, 1))
    wmat = qirr_sym.get_kq_weight_mat_mp(kf_ibz, pseudo_freq, bande, nqf,efermi=fermi_level, stdout=True, cores=ncore, interp_w=False, za_qcut=0,     reci_vec=RECI_VEC,triangular_wt=True) # obtain weight matrix
    np.savez(out_folder+"weights.npz", wmat)


    print("bands min: ", np.min(bande), "bands max: ", np.max(bande))
    print("the fermi_level = ", fermi_level)
    print("E cut relative to Ef:",ecut)
    print("Number of initial k points in irrbz: ", len(kf_ibz))




    ## here we write the wmat in readable format
    fo = open(out_folder+"readable_wmat.dat", "w")
    fo1 = open(out_folder+"scfwt.dat", "w")
    fo2 = open(out_folder+"scfkidx.dat", "w")
    print("start to write wmat:")
    nbnd = len(bande[0])
    
    scfklist=[]
    for ib in range(nbnd):
        for jb in range(nbnd):
            iks, iqs = wmat[ib,jb,0].nonzero()
            for ik, iq in zip(iks, iqs):
                ikq = qirr_sym.kindex_add(ik, iq, nkf, nkf)
                if ik not in scfklist:
                    scfklist=scfklist+[ik]
                if ikq not in scfklist:
                    scfklist=scfklist+[ikq]

    for ik,ik_fx in zip(scfklist,qirr_sym.get_kfrac(scfklist, nkf, nkf)):
        fo2.write("%8d %8d %16.12f %16.12f %16.12f\n"  % (ik, scfklist.index(ik)+1,ik_fx[0], ik_fx[1], 0.0))
    fo2.close()

 
    fo1.write("#ibnd, ik, coord(crystal relative), E(eV),v(eV/A) -> fbnd, fk, fk coord,  E,v: wt \n")
    for ib in range(nbnd):
        for jb in range(nbnd):
            fo.write("\nibnd: %d -> jbnd: %d\n" %(ib, jb))
            iks, iqs = wmat[ib,jb,0].nonzero()
            for ik, iq in zip(iks, iqs):
                ikq = qirr_sym.kindex_add(ik, iq, nkf, nkf)
                [ik_fx, ikq_fx] = qirr_sym.get_kfrac([ik,ikq], nkf, nkf)
                fo.write("%8d ( %16.12f %16.12f ) ->%8d ( %16.12f %16.12f ): %20.12e\n"  % (ik, ik_fx[0], ik_fx[1], ikq, ikq_fx[0],  ikq_fx[1], wmat[ib,jb,0][ik,iq]))
                fo1l=''
                fo1l=fo1l+'%8d %8d '% (ib+nbnd_valence+1, scfklist.index(ik)+1) 
                fo1l=fo1l+'( %20.12e %20.12e  %20.12e ) '% (ik_fx[0],ik_fx[1],0.0)
                fo1l=fo1l+'%20.12e '% bande[ik,ib] 
                fo1l=fo1l+'( %20.12e %20.12e  %20.12e ) '% tuple(velocity[ik,ib])
                fo1l=fo1l+' ->  '
                fo1l=fo1l+'%8d %8d '% (jb+nbnd_valence+1, scfklist.index(ikq)+1) 
                fo1l=fo1l+'( %20.12e %20.12e  %20.12e ) '% (ikq_fx[0],ikq_fx[1],0.0)
                fo1l=fo1l+'%20.12e '% bande[ikq,jb] 
                fo1l=fo1l+'( %20.12e %20.12e  %20.12e ) '% tuple(velocity[ikq,jb])
                fo1l=fo1l+': %20.12e  \n'% wmat[ib,jb,0][ik,iq]
                fo1.write(fo1l)
#                fo1.write("%8d %8d ( %16.12f %16.12f  %16.12f )  %16.12f  ( %16.12f %16.12f  %16.12f )   -> %8d %8d ( %16.12f %16.12f  %16.12f )  %16.12f  ( %16.12f %16.12f  %16.12f )  : %20.12e\n"  % (ib+nbnd_valence+1, scfklist.index(ik)+1, ik_fx[0], ik_fx[1], ik_fx[2],bande[ik,ib],velocity jb+nbnd_valence+1, scfklist.index(ikq)+1, ikq_fx[0],   ikq_fx[1], ikq_fx[2], wmat[ib,jb,0][ik,iq]))
    fo.close()
    fo1.close()
    print("done!")


    # test and plot for special k point
    reci_vec_inv = np.linalg.inv(RECI_VEC)
    spk_fx = np.array([1./3., 1./3.])+np.dot([0,0.02],reci_vec_inv)
    spk_ik = np.array(qirr_sym.qbz2ibz(qirr_sym.get_kindex2([spk_fx], 
    nkf, nkf), nkf, nkf), dtype=int)
    print("energy of spk (fermi=0):", bande[spk_ik])
    print("num of final k points: ", len(wmat[0,0,0][spk_ik].nonzero()[0]))
    print("number of kf_ibz:", len(kf_ibz))


