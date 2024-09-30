import numpy as np 
import matplotlib.pyplot as plt
import scipy.sparse as sp 
from time import time
from slabcond import SlabCond
import os.path
import copy

class IterSlabCond(SlabCond):
    """ 
    Main driver for calculating slab conductivity using an iterative solution.

    Parameters:
    ----------
    ngd: list of int
        Number of grid points in each desired direction.
    filvel: str
        Filename containing the electronic group velocity data.
    fs0: float, optional
        Initial Fermi level value. Default is 0.
    fsthick: float, optional
        Thickness of the Fermi surface. Default is 0.15.
    assume_metal: bool, optional
        Whether to assume the material is metallic. Default is True.
    isibz: bool, optional
        Use the irreducible Brillouin Zone. Default is True.
    filinfo: str, optional
        Filename for QE standard output with lattice symmetry matrices. Default is 'hex.info'.
    fillw: str, optional
        Filename for electron self-energy due to e-ph scattering. Default is None.
    filpwc: str, optional
        Filename for PWCOND output. Default is None.
    slabwidth: float, optional
        Width of the slab, in units of Angstroms. Default is 1.
    bc: tuple of int (0 or 1), optional
        Boundary conditions for the film surfaces. Each element is 0 (without surface) or 1 (with surface). Default is (0, 0).
    area: float, optional
        Volume of the bulk unit cell in units of A^3 (cubic Angstroms). Default is 1.
    ibrav: int, optional
        Bravais lattice index. See `grids.py` for detailed lattice vectors. Default is 4 for hexagonal lattice.
    spsym: str or None, optional
        Special treatment of lattice symmetry. Options are 'nosym' to remove all lattice symmetry, 'noz' to remove symmetry related to the z-axis. Default is None.
    ----------
    filg2m: str, optional
        Filename containing e-ph coupling matrix elements (|g|^2 data). Default is 'fort.709'.
    nqg: list of int, optional
        Number of grid points in each desired direction for the phonon q-grid. Default is None.
    nz: int, optional
        Number of grid points in the z-direction for the slab. Default is 1.
    bcmod: str, optional
        Boundary condition mode for the iterative BTE solver. Options are 'periodic', 'simpleFS', or 'complexFS'. Default is 'periodic'.

    Notes:
    ------
    - **bcmod** specifies the type of boundary condition used in the calculation:
        - `'periodic'`: Periodic boundary conditions.
        - `'simpleFS'`: Simple Fermi surface boundary conditions.
        - `'complexFS'`: Complex Fermi surface boundary conditions.
    """

    def __init__(
        self, ngd, filvel, fs0=0, fsthick=0.15, assume_metal=True, isibz=True, filinfo='hex.info', fillw=None, filpwc=None, slabwidth=1, bc=(0,0), area=1, ibrav=4, spsym=None,
        filg2m='fort.709', nqg=None, nz=1, bcmod='periodic'
        ):
        """ bcmod: boundary condition, 'periodic', 'simpleFS' or 'complexFS'."""

        super().__init__(ngd=ngd, filvel=filvel, fs0=fs0, fsthick=fsthick, assume_metal=assume_metal, isibz=isibz, filinfo=filinfo, fillw=fillw, filpwc=filpwc, slabwidth=slabwidth, bc=bc, area=area, ibrav=ibrav, spsym=spsym)

        if os.path.isfile(filg2m+"_kk.npz"):
            self.mat_index = 'kk'
            tmp = np.load(filg2m+"_kk.npz", allow_pickle=True)
            self.g2xw = tmp['g2xw']
            self.trans = tmp['trans']
        elif os.path.isfile(filg2m+"_kq.npz"):
            self.mat_index = 'kq'
            tmp = np.load(filg2m+"_kq.npz", allow_pickle=True)
            self.g2xw = tmp['g2xw']
            self.trans = tmp['trans']
        else:
            if nqg == None:
                self.reader_g2m(filg2m=filg2m, index='kk')
                np.savez(filg2m+"_kk.npz", g2xw=self.g2xw, trans=self.trans)
            else:
                self.reader_g2m(filg2m=filg2m, nqg=nqg, index='kq')
                np.savez(filg2m+"_kq.npz", g2xw=self.g2xw, trans=self.trans)

        self.nz = nz
        self.bcmod = bcmod
        if bcmod == 'simpleFS':
            self.is_simplefs = True
        elif bcmod == 'complexFS':
            self.is_simplefs = False

    def reader_g2m(self, filg2m, nqg=None, index='kq'):
        """ here we don't force the nqgrid to be equal to nkgrid, 
        index = 'kq' for [ik, iq] and index = 'kk' for [ik, ik]. """

        fo = open(filg2m, 'r')

        #### here we initialize
        if nqg == None:
            nqgpts = self.ngpts
        else:
            nqgpts = int(nqg[0]*nqg[1]*nqg[2]) if len(nqg) == 3 else int(nqg[0]*nqg[1])
        g2xw = np.zeros((self.nbnd, self.nbnd), dtype=object)
        trans = np.zeros((self.nbnd, self.nbnd), dtype=object)

        for ibnd in range(self.nbnd):
            for jbnd in range(self.nbnd):
                if index == 'kq':
                    g2xw[ibnd,jbnd] = sp.lil_matrix((self.ngpts, nqgpts))
                    trans[ibnd,jbnd] = sp.lil_matrix((self.ngpts, nqgpts))
                elif index == 'kk':
                    assert nqgpts == self.ngpts
                    g2xw[ibnd,jbnd] = sp.lil_matrix((self.ngpts, self.ngpts))
                    trans[ibnd,jbnd] = sp.lil_matrix((self.ngpts, self.ngpts))

        #### read and check
        with open(filg2m, 'r') as fo:
            lines = fo.readlines()
        header = lines[:1000]
        guess_nbnd = max([int(l.split()[2]) for l in header])
        assert guess_nbnd == self.nbnd

        iq = 0
        ik_list = list(self.irrlockid.keys()) ## ordered values assumed
        if len(header[0].split()) > 6: ### transporb not implemented
            for i in range(len(lines)):
                words = lines[i].split()

                lastq = iq
                iq = int(words[0]) - 1
                ikbz = ik_list[int(words[1]) - 1]
                ibnd = int(words[2]) - 1
                jbnd = int(words[3]) - 1

                if lastq != iq and (iq+1)%100==0:
                    print("current iq:", iq+1, "/", nqgpts)
                
                if index == 'kq':
                    ind2 = iq
                elif index == 'kk':
                    ind2 = self._kindex_add(ikbz,iq)
                g2xw[ibnd,jbnd][ikbz,ind2] += float(words[5])*float(words[6])* self._ryd2ev
        else:
            for i in range(len(lines)):
                words = lines[i].split()

                lastq = iq
                iq = int(words[0]) - 1
                ikbz = ik_list[int(words[1]) - 1]
                ibnd = int(words[2]) - 1
                jbnd = int(words[3]) - 1

                if lastq != iq and (iq+1)%100==0:
                    print("current iq:", iq+1, "/", nqgpts)

                if index == 'kq':
                    ind2 = iq
                elif index == 'kk':
                    ind2 = self._kindex_add(ikbz,iq)
                g2xw[ibnd,jbnd][ikbz,ind2] += float(words[4])* self._ryd2ev
                trans[ibnd,jbnd][ikbz,ind2] += float(words[5])* self._ryd2ev

        self.mat_index = index
        self.g2xw = self.lil2csr_obj(g2xw)
        self.trans = self.lil2csr_obj(trans)
    
    def lil2csr_obj(self, input, inv=False):
        shape = np.shape(input)
        output = np.zeros(shape, dtype=object)
        if len(shape) == 2:
            for i, j in np.ndindex(shape):
                if inv:
                    output[i,j] = sp.lil_matrix(input[i,j])
                else:
                    output[i,j] = sp.csr_matrix(input[i,j])
        elif len(shape) == 3:
            for i, j, k in np.ndindex(shape):
                if inv:
                    output[i,j,k] = sp.lil_matrix(input[i,j,k])
                else:
                    output[i,j,k] = sp.csr_matrix(input[i,j,k])
        else:
            raise Exception("not implemented!")
        return output

    def set_scattering_rate(self, rta='se'):
        """ rta = 'se' for SERTA; 'm' for MRTA. """

        if rta == 'm':
            _vnrm = self.vel / np.linalg.norm(self.vel, axis=2)[:,:,np.newaxis]

        # _nk = self.irrnpts if self.has_irr else self.nk
        # _nbnd = self.nbnd
        # scat = np.zeros((_nbnd,_nk))

        self.scat = np.zeros(np.shape(self.scat))
        for ibnd in range(self.nbnd):
            for jbnd in range(self.nbnd):
                _g2xw = self.g2xw[ibnd,jbnd]

                if rta == 'm':
                    iks, iqs = _g2xw.nonzero()
                    for ik, iq in zip(iks, iqs):
                        ##### ikq for final state
                        if self.mat_index == 'kq':
                            ikq = self._kindex_add(ik,iq)
                        else:
                            ikq = iq

                        ##### angle term
                        try:
                            ikloc, ikqloc = self.b2i[ik], self.b2i[ikq]
                            angle_term = 1.0 - np.dot(_vnrm[ibnd,ikloc], _vnrm[jbnd,ikqloc])
                        except:
                            angle_term = 0
                        _g2xw[ik,iq] = _g2xw[ik,iq] * angle_term

                for ikbz in self.irrlockid.keys():
                    ikirrloc = self.irrlockid[ikbz]
                    self.scat[ibnd,ikirrloc] += np.array(_g2xw.sum(axis=1))[ikbz,0]
                
                print("scattering completed for ibnd, jbnd: ", ibnd, jbnd)
                        
        ##### FIXME: still need to update the self.scatfb and self.scatfb2


    def set_f1serta(self):
        try:
            self.scat
        except:
            self.set_scattering_rate()

        #### velocity
        _velirr = self.vel[:,self.i2locb,:]
        self.velirr = _velirr
        ## rotate velirr into f1bz, this should be symmetric and close to self.vel
        velbz_unrot = self.velirr[:,self.locb2i,:]
        velbz = np.zeros(np.shape(velbz_unrot))
        for (ib,ik), _ in np.ndenumerate(velbz[:,:,0]):
            _symmatf = self.rotmat[ik]
            _symmatc = self.f2c(_symmatf)
            velbz[ib,ik] = np.dot(_symmatc, velbz_unrot[ib,ik])
        self.velsym = velbz
        #### partial erg
        _mergirr = np.exp((self.fs-self.erg[:,self.i2locb])/self._kT)
        _mergirr = _mergirr / (_mergirr+1)**2 / self._kT
        #### tau
        _tau = np.zeros(np.shape(self.scat))
        _tau[self.scat>0] = 1./self.scat[self.scat>0]
        self.tauirr = _tau

        _f1serta = _velirr * _tau[:,:,np.newaxis] * _mergirr[:,:,np.newaxis]
        self.f1serta = np.repeat([_f1serta], repeats=self.nz, axis=0)
        ### in nz*nbnd*nk*ndim

    def gen_f1_empty(self):
        ##### initialize a empty f1
        f1ept = np.zeros((self.nz, self.nbnd, self.ndim), dtype=object)
        for iz, ib, ipol in np.ndindex(self.nz, self.nbnd, self.ndim):
                f1ept[iz,ib,ipol] = sp.csr_matrix((self.ngpts,1))
        return f1ept

    def iter_init(self):
        """ init variables for iter, should after set_f1serta"""

        ## here we also set the index for left-moving (vz<0) states and right-moving (vz>0) states
        _ikbzloc2ikbz = np.array(list(self.lockid.keys()), dtype=int)
        self.lmloc = [np.transpose(np.argwhere(self.velsym[ib,:,-1]<0))[0] for ib in range(self.nbnd)]
        self.rmloc = [np.transpose(np.argwhere(self.velsym[ib,:,-1]>=0))[0] for ib in range(self.nbnd)]
        self.lmind = [_ikbzloc2ikbz[self.lmloc[ib]] for ib in range(self.nbnd)]
        self.rmind = [_ikbzloc2ikbz[self.rmloc[ib]] for ib in range(self.nbnd)]

        ## calculate R related function
        if self.nz > 1:
            self.dz = self.a/(self.nz-1)
            _taubz = self.tauirr[:,self.locb2i]
            _R = self.velsym[:,:,-1]*self._eva2ms * _taubz/self._ev2inv_s / (self.dz*self._ang)
            _R = np.abs(_R)
            _eRinv = np.zeros(np.shape(_R))
            _eRinv[_R>0] = np.exp(-1./_R[_R>0])
            self.Rfz0 = _eRinv
            self.Rgz0 = _R - (1.0+_R) * _eRinv
            self.Rgz1 = 1 - _R + _R * _eRinv
            ## for boundary condition use
            _Ra = self.velsym[:,:,-1]*self._eva2ms * _taubz/self._ev2inv_s / (self.a*self._ang)
            _Ra = np.abs(_Ra)
            _eRainv = np.zeros(np.shape(_Ra))
            _eRainv[_Ra>0] = np.exp(-1./_Ra[_Ra>0])
            self.Lcoef1 = np.repeat(_eRainv[:,:,np.newaxis], repeats=self.ndim, axis=2)
            #### take this into tensor 
            if self.ndim == 3:
                self.Lcoef1fb = np.zeros((self.nbnd, self.ngd[0], self.ngd[1], self.ngd[2], self.ndim))
                for i in self.lockid.keys():
                    ix, iy, iz = self._get_ikil([i])[0]
                    self.Lcoef1fb[:,ix,iy,iz,:] = self.Lcoef1[:,self.lockid[i],:]
            elif self.ndim == 2:
                self.Lcoef1fb = np.zeros((self.nbnd, self.ngd[0], self.ngd[1], self.ndim))
                for i in self.lockid.keys():
                    ix, iy = self._get_ikil([i])[0]
                    self.Lcoef1fb[:,ix,iy,:] = self.Lcoef1[:,self.lockid[i],:]


    def get_slabcond_iter(self, f1=None):
        """calculate slab conductivity from f1 in ibz"""
        if f1 is None:
            try:
                self.f1serta
            except:
                self.set_f1serta()
            f1 = self.f1serta

        #### rotate f1 into f1bz
        f1bz_unrot = f1[:,:,self.locb2i,:]
        f1bz = np.zeros(np.shape(f1bz_unrot))
        for (iz,ib,ik), _ in np.ndenumerate(f1bz[:,:,:,0]):
            _symmatf = self.rotmat[ik]
            _symmatc = self.f2c(_symmatf)
            f1bz[iz,ib,ik] = np.dot(_symmatc, f1bz_unrot[iz,ib,ik])
            ### only fortest
            # if self.velsym[ib,ik,-1] < 0:
            #     f1bz[iz,ib,ik,:] = 0.0
        
        #### f1bz to slab conductivity
        condz = np.zeros((self.nz, self.ndim, self.ndim))
        # _velbz = self.vel # just for test, self.vel not symmetric
        _velbz = self.velsym
        for i in range(self.ndim):
            for j in range(self.ndim):
                condz[:,i,j] = np.sum(f1bz[:,:,:,i]*_velbz[np.newaxis,:,:,j], axis=(1,2))
        condz = condz * self._eva2ms**2 / self._ev2inv_s / self.ngpts
        condz = condz * 2 # without soc, spin degeneracy

        if self.ndim == 2:
            condz = condz * self._e / (self.area*self._bohr**2)
        else:
            condz = condz * self._e / (self.area*self._bohr**3)

        return condz

    def get_slabcond_sparse(self, f1, mod='+-'):
        """ f1 in sparse matrix in bz """

        f1bz = self.f1_array2sp(input=f1, inv=True, irrbz=False)

        #### f1bz to slab conductivity
        condz = np.zeros((self.nz, self.ndim, self.ndim))
        # _velbz = self.vel # just for test, self.vel not symmetric
        _velbz = self.velsym
        for i in range(self.ndim):
            for j in range(self.ndim):
                if mod == '+-':
                    condz[:,i,j] = np.sum(f1bz[:,:,:,i]*_velbz[np.newaxis,:,:,j], axis=(1,2))
                elif mod == '+':
                    for ib in range(self.nbnd):
                        _rmloc = self.rmloc[ib]
                        condz[:,i,j] += np.sum(f1bz[:,ib,_rmloc,i]*_velbz[np.newaxis,ib,_rmloc,j], axis=1)
                elif mod == '-':
                    for ib in range(self.nbnd):
                        _lmloc = self.lmloc[ib]
                        condz[:,i,j] += np.sum(f1bz[:,ib,_lmloc,i]*_velbz[np.newaxis,ib,_lmloc,j], axis=1)
        condz = condz * self._eva2ms**2 / self._ev2inv_s / self.ngpts
        condz = condz * 2 # without soc, spin degeneracy

        if self.ndim == 2:
            condz = condz * self._e / (self.area*self._bohr**2)
        else:
            condz = condz * self._e / (self.area*self._bohr**3)

        return condz

    def f1_array2sp(self, input, inv=False, irrbz=True):
        if not inv and irrbz:
            f1 = np.zeros((self.nz, self.nbnd, self.ndim), dtype=object)
            for iz, ib, ipol in np.ndindex(self.nz, self.nbnd, self.ndim):
                f1[iz,ib,ipol] = sp.lil_matrix((self.ngpts,1))
                for ikbz in self.irrlockid.keys():
                    ikirrloc = self.irrlockid[ikbz]
                    f1[iz,ib,ipol][ikbz,0] = input[iz,ib,ikirrloc,ipol]
            f1 = self.lil2csr_obj(f1)
        elif not inv and (not irrbz):
            f1 = np.zeros((self.nz, self.nbnd, self.ndim), dtype=object)
            for iz, ib, ipol in np.ndindex(self.nz, self.nbnd, self.ndim):
                f1[iz,ib,ipol] = sp.lil_matrix((self.ngpts,1))
                for ikbz in self.lockid.keys():
                    ikbzloc = self.lockid[ikbz]
                    f1[iz,ib,ipol][ikbz,0] = input[iz,ib,ikbzloc,ipol]
            f1 = self.lil2csr_obj(f1)
        elif inv and irrbz:
            nk = self.irrnpts if self.has_irr else self.nk
            f1 = np.zeros((self.nz, self.nbnd, nk, self.ndim))
            for iz, ib, ipol in np.ndindex(self.nz, self.nbnd, self.ndim):
                for ikbz in self.irrlockid.keys():
                    ikirrloc = self.irrlockid[ikbz]
                    f1[iz,ib,ikirrloc,ipol] = input[iz,ib,ipol][ikbz,0]
        elif inv and (not irrbz):
            f1 = np.zeros((self.nz, self.nbnd, self.nk, self.ndim))
            for iz, ib, ipol in np.ndindex(self.nz, self.nbnd, self.ndim):
                for ikbz in self.lockid.keys():
                    ikbzloc = self.lockid[ikbz]
                    f1[iz,ib,ikbzloc,ipol] = input[iz,ib,ipol][ikbz,0]
        else:
            raise Exception('Not implemented!')
        return f1
        # return sp.csr_matrix(f1)

    def f1_ibz2bz(self, input, dtype='sparse', inv=False):
        """ dtype for sparse matrix """
        if dtype == 'sparse' and inv == False:
            output = np.zeros((self.nz, self.nbnd, self.ndim), dtype=object)
            for iz, ib, ipol in np.ndindex(self.nz, self.nbnd, self.ndim):
                output[iz,ib,ipol] = sp.lil_matrix((self.ngpts,1))
            for ikbz in self.lockid.keys():
                ikbzloc = self.lockid[ikbz]
                ikirrloc = self.b2i[ikbz]
                ikbz2 = list(self.irrlockid.keys())[ikirrloc]
                _symmatc = self.f2c(self.rotmat[ikbzloc])
                for iz in range(self.nz):
                    for ib in range(self.nbnd):
                        if self.ndim == 3:
                            _f1 = [input[iz,ib,0][ikbz2,0], input[iz,ib,1][ikbz2,0], input[iz,ib,2][ikbz2,0]]
                            _f1 = np.dot(_symmatc, _f1)
                            output[iz,ib,0][ikbz,0] = _f1[0]
                            output[iz,ib,1][ikbz,0] = _f1[1]
                            output[iz,ib,2][ikbz,0] = _f1[2]
                        elif self.ndim == 2:
                            _f1 = [input[iz,ib,0][ikbz2,0], input[iz,ib,1][ikbz2,0]]
                            _f1 = np.dot(_symmatc, _f1)
                            output[iz,ib,0][ikbz,0] = _f1[0]
                            output[iz,ib,1][ikbz,0] = _f1[1]
            output = self.lil2csr_obj(output)
        elif dtype == 'sparse' and inv == True:
            ikbzirr = list(self.irrlockid.keys())
            output = np.zeros((self.nz, self.nbnd, self.ndim), dtype=object)
            for iz, ib, ipol in np.ndindex(self.nz, self.nbnd, self.ndim):
                output[iz,ib,ipol] = sp.lil_matrix((self.ngpts,1))
                output[iz,ib,ipol][ikbzirr,0] = input[iz,ib,ipol][ikbzirr,0].toarray()
            output = self.lil2csr_obj(output)
        else:
            raise Exception("not implemented. ")
        
        return output

    def set_f1_boundary_pd(self, f1, f2):
        """ f1 -> f2, in bz"""
        lz, rz = 0, self.nz-1
        for ib in range(self.nbnd):
            _lmind, _rmind = self.lmind[ib], self.rmind[ib]
            for ipol in range(self.ndim):
                f2[lz,ib,ipol][_rmind,0] = f1[rz,ib,ipol][_rmind,0].toarray()
                f2[rz,ib,ipol][_lmind,0] = f1[lz,ib,ipol][_lmind,0].toarray()
        return f2

    def set_f1_boundary_pd2(self, f1, f2):
        lz, rz = 0, self.nz-1
        for ib in range(self.nbnd):
            _lmloc, _rmloc = self.lmloc[ib], self.rmloc[ib]
            f2[lz,ib,_rmloc,:] = f1[rz,ib,_rmloc,:]
            f2[rz,ib,_lmloc,:] = f1[lz,ib,_lmloc,:]
        return f2

    def set_f1_boundary_zero(self, f1):
        """ set right-moving states at lz and left-moving states at rz zeros"""
        lz, rz = 0, self.nz-1
        for ib in range(self.nbnd):
            _lmind, _rmind = self.lmind[ib], self.rmind[ib]
            for ipol in range(self.ndim):
                f1[lz,ib,ipol][_rmind,0] = 0.0
                f1[rz,ib,ipol][_lmind,0] = 0.0
        return f1

    def do_drift(self, f1, g1, f2):
        """ f1 -> f2, do drift in bz"""

        if self.bcmod == 'periodic':
            f2 = self.set_f1_boundary_pd(g1, f2)
        else:
            raise Exception("Not implemented!")

        # condz = self.get_slabcond_iter(self.f1_array2sp(f2, inv=True))
        ## fortest


        #### from boundary conditions to full real space
        ## first the right moving states
        for iz in range(1, self.nz):
            for ib in range(self.nbnd):
                _rmind, _rmloc = self.rmind[ib], self.rmloc[ib]
                for ipol in range(self.ndim):
                    _fz0, _Rfz0 = f2[iz-1,ib,ipol][_rmind,0].toarray(), self.Rfz0[ib,_rmloc][:,np.newaxis]
                    _gz0, _Rgz0 = g1[iz-1,ib,ipol][_rmind,0].toarray(), self.Rgz0[ib,_rmloc][:,np.newaxis]
                    _gz1, _Rgz1 = g1[iz,ib,ipol][_rmind,0].toarray(), self.Rgz1[ib,_rmloc][:,np.newaxis]
                    f2[iz,ib,ipol][_rmind,0] = _fz0*_Rfz0 + _gz0*_Rgz0 + _gz1*_Rgz1
        ## then the left moving states
        for iz in list(range(self.nz-1))[::-1]:
            for ib in range(self.nbnd):
                _lmind, _lmloc = self.lmind[ib], self.lmloc[ib]
                for ipol in range(self.ndim):
                    _fz0, _Rfz0 = f2[iz+1,ib,ipol][_lmind,0].toarray(), self.Rfz0[ib,_lmloc][:,np.newaxis]
                    _gz0, _Rgz0 = g1[iz+1,ib,ipol][_lmind,0].toarray(), self.Rgz0[ib,_lmloc][:,np.newaxis]
                    _gz1, _Rgz1 = g1[iz,ib,ipol][_lmind,0].toarray(), self.Rgz1[ib,_lmloc][:,np.newaxis]
                    f2[iz,ib,ipol][_lmind,0] = _fz0*_Rfz0 + _gz0*_Rgz0 + _gz1*_Rgz1

        return f2

    def do_drift2(self, g1):
        """similar to do_drift but use ndarray to speedup"""
        f2 = np.zeros(np.shape(g1))

        if self.bcmod == 'periodic':
            f2 = self.set_f1_boundary_pd2(g1, f2)
        elif self.bcmod == 'simpleFS':
            f2 = self.set_f1_boundary_simpleFS2(g1, f2)
        elif self.bcmod == 'complexFS':
            f2 = self.set_f1_boundary_simpleFS2(g1, f2)
        else:
            raise Exception("Not implemented!")

        ##### from boundary conditions to full real space
        ## first the right moving states
        for iz in range(1, self.nz):
            for ib in range(self.nbnd):
                _rmloc = self.rmloc[ib]
                _fz0, _Rfz0 = f2[iz-1,ib,_rmloc,:], self.Rfz0[ib,_rmloc,np.newaxis]
                _gz0, _Rgz0 = g1[iz-1,ib,_rmloc,:], self.Rgz0[ib,_rmloc,np.newaxis]
                _gz1, _Rgz1 = g1[iz,ib,_rmloc,:], self.Rgz1[ib,_rmloc,np.newaxis]
                f2[iz,ib,_rmloc,:] = _fz0*_Rfz0 + _gz0*_Rgz0 + _gz1*_Rgz1
        ## then the left moving states
        for iz in list(range(self.nz-1))[::-1]:
            for ib in range(self.nbnd):
                _lmloc = self.lmloc[ib]
                _fz0, _Rfz0 = f2[iz+1,ib,_lmloc,:], self.Rfz0[ib,_lmloc,np.newaxis]
                _gz0, _Rgz0 = g1[iz+1,ib,_lmloc,:], self.Rgz0[ib,_lmloc,np.newaxis]
                _gz1, _Rgz1 = g1[iz,ib,_lmloc,:], self.Rgz1[ib,_lmloc,np.newaxis]
                f2[iz,ib,_lmloc,:] = _fz0*_Rfz0 + _gz0*_Rgz0 + _gz1*_Rgz1

        return f2

    def iter_solver(self, maxnstep=0):
        """ iteratively solve the BTE. """

        ###### set up for f10
        try:
            self.f1serta
        except:
            self.set_f1serta()
        f10 = self.f1_array2sp(self.f1serta, inv=False)

        stt = time()
        ##### set up G matrix: G=tau*\tilde T
        # _G = copy.deepcopy(self.trans)
        # _tau = self.tauirr
        # for ikbz in self.irrlockid.keys():
        #     ikirrloc = self.irrlockid[ikbz]
        #     for i in range(self.nbnd):
        #         for j in range(self.nbnd):
        #             _G[i,j][ikbz,:] = self.trans[i,j][ikbz,:]*_tau[i,ikirrloc]
        #             _G[i,j][ikbz,ikbz] = 0.0 ## enforce off-diagonal
        # fit = time()
        # print(f"G matrix set up, time: {fit-stt:10.3f} s.")
        ##### set up G matrix: G=tau*\tilde T
        _G = self.lil2csr_obj(self.trans, inv=True)
        _trans = self.lil2csr_obj(self.trans, inv=True)
        _tau = self.tauirr
        for ikbz in self.irrlockid.keys():
            ikirrloc = self.irrlockid[ikbz]
            for i in range(self.nbnd):
                for j in range(self.nbnd):
                    _G[i,j][ikbz,:] = _trans[i,j][ikbz,:]*_tau[i,ikirrloc]
                    _G[i,j][ikbz,ikbz] = 0.0 ## enforce off-diagonal
        _G = self.lil2csr_obj(_G)
        fit = time()
        print(f"G matrix set up, time: {fit-stt:10.3f} s.")

        ##### initialize a empty f1
        f1ept = np.zeros((self.nz, self.nbnd, self.ndim), dtype=object)
        for iz, ib, ipol in np.ndindex(self.nz, self.nbnd, self.ndim):
                f1ept[iz,ib,ipol] = sp.csr_matrix((self.ngpts,1))

        ##### other init
        self.iter_init()
        print("iteration initialized!")


        ##### do iteration
        f1tot = copy.deepcopy(f1ept)
        istep = 0
        f_in = copy.deepcopy(f10)
        f_in_bz = self.f1_ibz2bz(f_in, dtype='sparse')
        sttot = time()
        while (istep < maxnstep):

            #### k-space
            stt = time()
            if istep > 0:
                g_in = copy.deepcopy(f1ept)
                for iz, ib, ipol in np.ndindex(self.nz, self.nbnd, self.ndim):
                    for jb in range(self.nbnd):
                        g_in[iz,ib,ipol] = g_in[iz,ib,ipol] + _G[ib,jb].dot(f_in_bz[iz,jb,ipol])
            elif istep == 0:
                g_in = copy.deepcopy(f10)
            g_in_bz = self.f1_ibz2bz(g_in, dtype='sparse')
            fit = time()
            print(f"k-space step: {istep:5d}, time: {fit-stt:10.3f}/{fit-sttot:10.3f} s.")

            #### r-space
            stt2 = time()
            if self.nz == 1:
                f_out_bz = g_in_bz
            else:
                #### do drift in sparse matrix
                # f_out_bz = copy.deepcopy(f1ept)
                # f_out_bz = self.do_drift(f1=f_in_bz, g1=g_in_bz, f2=f_out_bz)
                #### or do drift2 in ndarray
                ginarray = self.f1_array2sp(g_in_bz, inv=True, irrbz=False)
                foutarray = self.do_drift2(g1=ginarray)
                f_out_bz = self.f1_array2sp(foutarray, inv=False, irrbz=False)

            fit2 = time()
            print(f"r-space step: {istep:5d}, time: {fit2-stt2:10.3f}/{fit2-sttot:10.3f} s.")


            # df1 = self.f1_array2sp(f_out, inv=True)
            f1tot = f1tot + f_out_bz
            if istep > 0: old_condz = np.copy(condz)
            condz = self.get_slabcond_sparse(f1=f1tot, mod='+-')

            ### simple report and check
            print(f"current step: {istep:5d}")
            # self.print_condz(condz, avg=False)

            # print(condz)
            if istep > 0 and np.allclose(old_condz[:,0,0], condz[:,0,0]):
                print(f"coverged in {istep+1:5d} steps!")
                self.print_condz(condz, avg=False)
                break

            f_in_bz = f_out_bz
            istep += 1
        
        if self.nz > 1:
            print("final average conductivity:")
            self.print_condz(condz, avg=True)
            
        ##### 
        return f1tot


    def print_condz(self, condz, avg=False):
        if avg == False:
            if self.ndim == 3:
                for icondz in condz: print(f"{icondz[0,0]:10.4e}, {icondz[1,1]:10.4e}, {icondz[2,2]:10.4e}")
            elif self.ndim == 2:
                for icondz in condz: print(f"{icondz[0,0]:10.4e}, {icondz[1,1]:10.4e}")
        elif avg == True:
            avgcond = [(np.sum(condz[0:self.nz-1,ipol,ipol])+np.sum(condz[1:self.nz,ipol,ipol]))/2/(self.nz-1) for ipol in range(self.ndim)]
            if self.ndim == 3:
                print(f"{avgcond[0]:10.4e}, {avgcond[1]:10.4e}, {avgcond[2]:10.4e}")
            elif self.ndim == 2:
                print(f"{avgcond[0]:10.4e}, {avgcond[1]:10.4e}")


    def set_Lcoef0(self, g1):
        """set up 0-order boundary connection function"""
        Lcoef0 = np.zeros((self.nbnd, self.nk, self.ndim))

        ### first calculate L2
        for ib in range(self.nbnd):
            _rmloc = self.rmloc[ib]
            _coef0 = np.zeros((len(_rmloc),self.ndim))
            _Rfz0 = self.Rfz0[ib,_rmloc,np.newaxis]
            _Rgz0 = self.Rgz0[ib,_rmloc,np.newaxis]
            _Rgz1 = self.Rgz1[ib,_rmloc,np.newaxis]
            for iz in range(1,self.nz):
                _coef0 = _coef0*_Rfz0 + g1[iz-1,ib,_rmloc,:]*_Rgz0 + g1[iz,ib,_rmloc,:]*_Rgz1
            Lcoef0[ib,_rmloc,:] = _coef0
        ### then calcualte L1
        for ib in range(self.nbnd):
            _lmloc = self.lmloc[ib]
            _coef0 = np.zeros((len(_lmloc),self.ndim))
            _Rfz0 = self.Rfz0[ib,_lmloc,np.newaxis]
            _Rgz0 = self.Rgz0[ib,_lmloc,np.newaxis]
            _Rgz1 = self.Rgz1[ib,_lmloc,np.newaxis]
            for iz in list(range(self.nz-1))[::-1]:
                _coef0 = _coef0*_Rfz0 + g1[iz+1,ib,_lmloc,:]*_Rgz0 + g1[iz,ib,_lmloc,:]*_Rgz1
            Lcoef0[ib,_lmloc,:] = _coef0
        self.Lcoef0 = Lcoef0
        #### take this into tensor
        if self.ndim == 3:
            self.Lcoef0fb = np.zeros((self.nbnd, self.ngd[0], self.ngd[1], self.ngd[2], self.ndim))
            for i in self.lockid.keys():
                ix, iy, iz = self._get_ikil([i])[0]
                self.Lcoef0fb[:,ix,iy,iz,:] = self.Lcoef0[:,self.lockid[i],:]
        elif self.ndim == 2:
            self.Lcoef0fb = np.zeros((self.nbnd, self.ngd[0], self.ngd[1], self.ndim))
            for i in self.lockid.keys():
                ix, iy = self._get_ikil([i])[0]
                self.Lcoef0fb[:,ix,iy,:] = self.Lcoef0[:,self.lockid[i],:]


    def get_f1ik_boundary(self, ib_ikbz):
        """ calculate f1 for ib, ikbz"""
        ib, ikbz = ib_ikbz
        
        e0 = self.get_dsik(ib_ikbz, sv='erg')
        z0 = self._get_kfl([ikbz])[0][-1]

        ###### find the roots
        ib_iz = [(ib, z0)]
        for ibf in range(self.nbnd):
            ibf_ikbz = (ibf, ikbz)
            if ib == ibf:
                fzs = self.find_roots_xxinz(ibf_ikbz, z0=e0, sv='erg', exclude=True)
            else:
                fzs = self.find_roots_xxinz(ibf_ikbz, z0=e0, sv='erg', exclude=False)
            for ifz in fzs: ib_iz.append((ibf,ifz))
        
        # print(f"ib: {ib:5d}, ikbz: {ikbz:5d}, nstates: {len(ib_iz):5d}")

        ##### resrot by: 0 current, other by their kz
        new_ib_iz = sorted(ib_iz[1:], key=lambda x: x[1])
        new_ib_iz.insert(0, ib_iz[0])

        L0, L1 = [], []
        kzl = []
        ##### get the interpolated final kz' and L1 and L2
        for iib, iiz in new_ib_iz:
            l0list = []
            l0list.append(self.interpinz((iib,ikbz), sv='L0x')(iiz))
            l0list.append(self.interpinz((iib,ikbz), sv='L0y')(iiz))
            if self.ndim == 3:
                l0list.append(self.interpinz((iib,ikbz), sv='L0z')(iiz))
            L0.append(l0list)
            l1list = []
            l1list.append(self.interpinz((iib,ikbz), sv='L1x')(iiz))
            l1list.append(self.interpinz((iib,ikbz), sv='L1y')(iiz))
            if self.ndim == 3:
                l1list.append(self.interpinz((iib,ikbz), sv='L1z')(iiz))
            L1.append(l1list)
            kzl.append(iiz)
        L0, L1 = np.array(L0), np.array(L1)
        kzl = np.array(kzl)

        ##### now slove the linalg equation
        rmat = self.get_rmatik(ib_ikbz, kzl) ## just to use the length of L0
        f1ik = np.zeros((self.ndim,))
        for ipol in range(self.ndim):
            #### fortest  FIXME: fixed for initial and final states
            # LinMatA = np.diag([1]*len(L0)) - np.dot(rmat, np.diag(L1[:,ipol]))
            # LinMatB = np.dot(rmat, L0[:,ipol])
            LinMatA = np.diag([1]*len(L0)) - np.dot(np.diag(L1[:,ipol]), rmat)
            LinMatB = np.dot(L0[:,ipol], rmat)
            x = np.linalg.solve(LinMatA, LinMatB)
            f1ik[ipol] = x[0]
        
        return f1ik


    def get_lcoefz(self, ib_ikbz, order=0, dir=2):
        ib, ikbz = ib_ikbz
        if order == 0:
            if self.ndim == 3:
                ikx, iky, _ = self._get_ikil([ikbz])[0]
                return self.Lcoef0fb[ib,ikx,iky,:,dir]
            elif self.ndim == 2:
                ikx, _ = self._get_ikil([ikbz])[0]
                return self.Lcoef0fb[ib,ikx,:,dir]
        elif order == 1:
            if self.ndim == 3:
                ikx, iky, _ = self._get_ikil([ikbz])[0]
                return self.Lcoef1fb[ib,ikx,iky,:,dir]
            elif self.ndim == 2:
                ikx, _ = self._get_ikil([ikbz])[0]
                return self.Lcoef1fb[ib,ikx,:,dir]

    def get_xxinz(self, ib_ikbz, sv='erg'):
        if sv in ['selfen', 'vx', 'vy', 'vz', 'erg']:
            cvals = super().get_xxinz(ib_ikbz, sv)
        elif sv == 'L0x':
            cvals = self.get_lcoefz(ib_ikbz, order=0, dir=0)
        elif sv == 'L0y':
            cvals = self.get_lcoefz(ib_ikbz, order=0, dir=1)
        elif sv == 'L0z':
            cvals = self.get_lcoefz(ib_ikbz, order=0, dir=2)
        elif sv == 'L1x':
            cvals = self.get_lcoefz(ib_ikbz, order=1, dir=0)
        elif sv == 'L1y':
            cvals = self.get_lcoefz(ib_ikbz, order=1, dir=1)
        elif sv == 'L1z':
            cvals = self.get_lcoefz(ib_ikbz, order=1, dir=2)
        return cvals
            

    def set_f1_boundary_simpleFS2(self, g1, f2):
        """ for current g(z), calculate the L1 and L2 in Lcoef"""

        ##### set up g-dependent Lcoef0
        self.set_Lcoef0(g1)

        ##### now do boundary conditions
        lz, rz = 0, self.nz-1
        for ib in range(self.nbnd):
            for ikloc, ikbz in zip(self.rmloc[ib], self.rmind[ib]):
                ib_ikbz = (ib, ikbz)
                erg = self.get_dsik(ib_ikbz, sv='erg')

                if np.abs(erg-self.fs) < self.ecut:
                    rmlzf1 = self.get_f1ik_boundary(ib_ikbz)
                    f2[lz,ib,ikloc,:] = rmlzf1
        ## then left-moving at rz
        for ib in range(self.nbnd):
            for ikloc, ikbz in zip(self.lmloc[ib], self.lmind[ib]):
                ib_ikbz = (ib, ikbz)
                erg = self.get_dsik(ib_ikbz, sv='erg')

                if np.abs(erg-self.fs) < self.ecut:
                    lmrzf1 = self.get_f1ik_boundary(ib_ikbz)
                    f2[rz,ib,ikloc,:] = lmrzf1

        return f2


if __name__ == '__main__':
    import lab
    ttkw = lab.tt1_kw
    folder = lab.tt['folder']

    ##### test basics
    fig, ax = plt.subplots()

    # for ISBW in [40,80,160,400,800,1200,2000,3000,4000]:
    for ISBW in [800]:
        for INZ in [13]:
            print(f"slabwidth:{ISBW} A; nz:{INZ}: ")
            metal1 = IterSlabCond(
                **ttkw, 
                slabwidth=ISBW,
                spsym='noz',
                ngd=[30,30,12], 
                filvel=folder+'wan/tt_geninterp.dat.k30kz12',
                fillw=folder+'inp/l_k30kz12_noz.dat',
                filg2m=folder+'inp/fort.k30kz12_noz.709', 
                nqg=None, 
                nz=INZ,
                bcmod='simpleFS',
                )

            metal1.set_scattering_rate(rta='se')
            # conref, *_ = metal1.get_conductivity(plot=False, area=lab.tt['area0'])
            mfpref = metal1.get_mean_free_path(outmod='vftauz')
            print("MPF:", mfpref)
            condz = metal1.get_slabcond_iter()
            print("SERTA conductivity:")
            metal1.print_condz(condz, avg=True)

            metal1.iter_solver(maxnstep=20)

