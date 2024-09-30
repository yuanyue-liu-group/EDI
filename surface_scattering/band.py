import numpy as np 
import matplotlib.pyplot as plt 
from grids import GridPts

class BandStruct(GridPts):
    """ class of band structure, containing:
        energy, occupation, fermi surface, etc. """
    
    _kT = 300 * 8.617333262e-5 
    _e = 1.60218e-19
    _me = 9.10938356e-31
    _bohr = 0.5291772109e-10
    _ang = 1e-10
    _hbar = 1.054571817E-34
    _eva2ms = _e * 1e-10 / _hbar
    _eff2invme = _eva2ms**2 / _e * _me 
    _ev2inv_s = _e * 2 / _hbar
    _ryd2ev = 13.605662285137

    def set_fermi(self, fs, fsthick=0.3, is_elec=True, assume_metal=False):

        self.fs = fs
        self.ecut = fsthick
        self.is_elec = is_elec 
        _nvb_all = np.sum( np.all(self.erg<self.fs, axis=1) )
        _nvb_out = np.sum( np.all(self.erg<(self.fs-self.ecut), axis=1) )
        self.nbnd_below_fs = _nvb_all - _nvb_out
        _ncb_all = np.sum( np.all(self.erg>self.fs, axis=1) )
        _ncb_out = np.sum( np.all(self.erg>(self.fs+self.ecut), axis=1) )
        self.nbnd_over_fs = _ncb_all - _ncb_out
        _nbnd_old = self.nbnd 
        self.nbnd = _nbnd_old - _nvb_out - _ncb_out

        if assume_metal:
            fsbndlist = list(range(_nvb_out, _nbnd_old-_ncb_out))
            self.cb_index = list(range(self.nbnd))
            self.vb_index = list(range(self.nbnd))
            self.erg = self.erg[fsbndlist,:]
            self.vel = self.vel[fsbndlist,:,:]
            print(_nvb_out, 'lower band disregared;', _ncb_out, 'higher band disregared.')
            bndlist = fsbndlist
        else:
            vb_index_old = [i for i in range(_nvb_out,_nvb_all)]
            cb_index_old = \
                [i for i in range(_nbnd_old - _ncb_all, _nbnd_old - _ncb_out)]
            self.vb_index = [i for i in range(0, self.nbnd_below_fs)]
            self.cb_index = \
                [i for i in range(self.nbnd - self.nbnd_over_fs, self.nbnd)]

            if is_elec:
                bndlist = cb_index_old
            else:
                bndlist = vb_index_old
            self.erg = self.erg[bndlist,:]
            self.vel = self.vel[bndlist,:,:]

        self.is_semi = not (assume_metal)

        self.cden = None 
        self.effmass = None 
        self.dosavg = None

        return bndlist

    def get_conductivity(self, plot=False, tfs=None, lmode=None, area=1.0):
        """ area: unit cell volumn in a.u.^3"""

        fs = self.fs if tfs is None else tfs 

        _bind = self.cb_index if self.is_elec else self.vb_index
        _vel = self.vel[_bind,:,:]
        _wt = 1.0/self.ngpts
        _merg = np.exp( (fs - self.erg[_bind,:]) / self._kT)
        _distr = _merg / (_merg + 1)**2 / self._kT

        if lmode is None:
            _scatfb = 1./(np.array(self.scat[:,self.locb2i])) \
                if self.has_irr else 1./self.scat
        else:
            dim1, dim2 = np.array(self.locb2i)[:,None], np.array(lmode)
            _scatfb = 1./(np.sum(self.scat_im[:,dim1,dim2], axis=2)) \
                if self.has_irr else 1./(np.sum(self.scat_im[:,:,dim2],axis=2))
        _scatfb[np.isinf(_scatfb)] = 0.0
        # _scatfb[:,:] = 1.0 #### FIXME: to obtain average tau
        
        self.cond = np.zeros((self.ndim,self.ndim))
        for i in range(self.ndim):
            for j in range(self.ndim):
                self.cond[i,j] = np.sum(_scatfb*_vel[:,:,i]*_vel[:,:,j]
                *_wt*_distr)

        _cond = self.cond / self._ev2inv_s * self._eva2ms**2 # to SI
        _cond = _cond * 2 # without soc, spin degeneracy
        if self.ndim == 2:
            _cond = _cond * self._e / (area*self._bohr**2)
        else:
            _cond = _cond * self._e / (area*self._bohr**3)

        print(_cond[0,0])
        if plot:
            print("Conductivity: (1/ohm/m)")
            print(_cond)
            plt.scatter(self.erg[_bind,:]-fs, 1./_scatfb)
            plt.ylabel("self-energy (eV)")
            plt.xlabel("energy (eV)")
            plt.xlim(-0.25, 0.25)
            plt.ylim(0,0.025)
            plt.show()
            plt.close()
            _tmpv = (_vel[:,:,0]**2 + _vel[:,:,1]**2)/2
            if self.ndim == 3: _tmpv = (_tmpv*2 + _vel[:,:,2]**2)/3
            valin = 1./_scatfb; valin[np.isinf(valin)] = 0
            # valin = _tmpv*_distr
            valin = valin[np.abs(self.erg[_bind,:]-fs)<0.05]
            print(np.average(valin))
            plt.scatter(self.erg[_bind,:]-fs, _tmpv*_distr)
            plt.ylabel(r'$v^2 \partial f/ \partial E$')
            # plt.scatter(self.erg[_bind,:]-fs, _tmpv)
            # plt.ylabel(r'$v^2$')
            plt.xlabel("energy (eV)")
            plt.xlim(-0.25, 0.25)
            plt.show()
            plt.close()

        # return _cond, self.erg[_bind,:], 1./_scatfb

        k1 = {et:i for i,et in enumerate(self.locb2i)}
        k2 = [k1[i] for i in range(self.irrnpts)]
        _tmpv = (_vel[:,:,0]**2 + _vel[:,:,1]**2 + _vel[:,:,2]**2)**0.5
        return _cond, self.erg[:,k2], self.scat, _tmpv[:,k2]*self._eva2ms
        # return _mob, self.erg[_bind,:], _tmpv*_wt*_distr/_cden*self._eff2invme

    def get_cden(self, tfs=None, area=None):

        if self.is_semi and tfs == None:
            _bind = self.cb_index if self.is_elec else self.vb_index
            if self.cden is None:
                _erg = (self.erg[_bind,:] - self.fs)/self._kT
                if not self.is_elec: _erg = -1.0*_erg
                self.cden = np.sum(1./(np.exp(_erg)+1))/self.ngpts 
            _cden = self.cden
        elif self.is_semi and tfs != None:
            _bind = self.cb_index if self.is_elec else self.vb_index
            _erg = (self.erg[_bind,:] - tfs)/self._kT
            if not self.is_elec: _erg = -1.0*_erg
            _cden = np.sum(1./(np.exp(_erg)+1))/self.ngpts
        elif not self.is_semi and tfs == None:
            _bind = self.cb_index if self.is_elec else self.vb_index
            _erg = (self.erg[_bind,:] - self.fs)/self._kT
            _felec = 1./(np.exp(_erg)+1.0)
            _cden = np.sum(_felec)/self.ngpts
            # print(_cden); _cden = 6 ## for test
        elif not self.is_semi and tfs != None:
            _bind = self.cb_index if self.is_elec else self.vb_index
            _erg = (self.erg[_bind,:] - tfs)/self._kT
            _felec = 1./(np.exp(_erg)+1.0)
            _cden = np.sum(_felec)/self.ngpts

        if area != None:
            """ output in m-3 unit if input with area, 2 for spin"""
            _cden = _cden*2 / (area*self._bohr**self.ndim)

        return _cden

    def get_eff_mass(self, plot=False):

        if not self.is_semi: raise Exception("No metal")

        if self.effmass is not None:
            _res = self.effmass * self._eff2invme
            return 1./_res 

        _bind = self.cb_index if self.is_elec else self.vb_index
        _vel = self.vel[_bind,:,:]
        _wt = 1.0/self.ngpts
        _merg = np.exp( (self.fs - self.erg[_bind,:]) / self._kT )
        _distr = _merg / (_merg + 1)**2 / self._kT

        if plot:
            _plt_tmp = _vel[:,:,0]*_vel[:,:,0] + _vel[:,:,1]*_vel[:,:,1] 
            if self.ndim == 3:
                _plt_tmp += _vel[:,:,2]*_vel[:,:,2]
            plt.scatter(self.erg[_bind,:], _plt_tmp*_distr*_wt)
            plt.xlim(self.fs-self.ecut, self.fs+self.ecut)
            plt.show()

        self.effmass = np.zeros((self.ndim,self.ndim))
        for i in range(self.ndim):
            for j in range(self.ndim):
                self.effmass[i,j] = np.sum(_vel[:,:,i] * _vel[:,:,j] * 
                    _wt * _distr )
        self.effmass = self.effmass / self.get_cden() 

        _res = self.effmass * self._eff2invme
        return 1./_res 

    def get_mean_free_path(self, tfs=None, lmode=None, outmod='L2', kT=None):
        """ area: unit cell volumn in a.u.^3"""

        fs = self.fs if tfs is None else tfs 

        if kT != None: 
            _kT = kT/300 * self._kT
        else:
            _kT = self._kT

        _bind = self.cb_index if self.is_elec else self.vb_index
        _vel = self.vel[_bind,:,:]
        _wt = 1.0/self.ngpts
        _merg = np.exp( (fs - self.erg[_bind,:]) / _kT)
        _distr = _merg / (_merg + 1)**2 / _kT

        if lmode is None:
            _scatfb = 1./(np.array(self.scat[:,self.locb2i])) \
                if self.has_irr else 1./self.scat
        else:
            dim1, dim2 = np.array(self.locb2i)[:,None], np.array(lmode)
            _scatfb = 1./(np.sum(self.scat_im[:,dim1,dim2], axis=2)) \
                if self.has_irr else 1./(np.sum(self.scat_im[:,:,dim2],axis=2))
        _scatfb[np.isinf(_scatfb)] = 0.0
        # _scatfb[:,:] = 1.0 #### FIXME: to obtain average tau
        
        if outmod == 'L2':
            self.cond = np.zeros((self.ndim,self.ndim))
            for i in range(self.ndim):
                for j in range(self.ndim):
                    self.cond[i,j] = np.sum(_scatfb*_scatfb*_vel[:,:,i]*_vel[:,:,j] *_wt*_distr) / np.sum(_wt*_distr)

            _cond = self.cond / self._ev2inv_s**2 * self._eva2ms**2 # to SI
            _cond = np.abs(_cond) ** 0.5 / 1e-10 # to A
        elif outmod == 'vftau2':
            _v2 = (_vel[:,:,0]**2 + _vel[:,:,1]**2)
            vf = np.sum(_v2**0.5 *_wt*_distr) / np.sum(_wt*_distr)
            tau = np.sum(_scatfb *_wt*_distr) / np.sum(_wt*_distr)
            _cond = vf*tau / self._ev2inv_s * self._eva2ms / 1e-10
        elif outmod == 'vftau3':
            _v2 = (_vel[:,:,0]**2 + _vel[:,:,1]**2 + _vel[:,:,2]**2)
            vf = np.sum(_v2**0.5 *_wt*_distr) / np.sum(_wt*_distr)
            tau = np.sum(_scatfb *_wt*_distr) / np.sum(_wt*_distr)
            _cond = vf*tau / self._ev2inv_s * self._eva2ms / 1e-10
        elif outmod == 'vftauz':
            _v2 = _vel[:,:,2]**2
            vf = np.sum(_v2**0.5 *_wt*_distr) / np.sum(_wt*_distr)
            tau = np.sum(_scatfb *_wt*_distr) / np.sum(_wt*_distr)
            _cond = vf*tau / self._ev2inv_s * self._eva2ms / 1e-10
        elif outmod == 'vf2':
            _v2 = (_vel[:,:,0]**2 + _vel[:,:,1]**2)
            vf = np.sum(_v2**0.5 *_wt*_distr) / np.sum(_wt*_distr)
            _cond = vf * self._eva2ms
        elif outmod == 'vf3':
            _v2 = (_vel[:,:,0]**2 + _vel[:,:,1]**2 + _vel[:,:,2]**2)
            vf = np.sum(_v2**0.5 *_wt*_distr) / np.sum(_wt*_distr)
            _cond = vf * self._eva2ms
        elif outmod == 'vfz':
            _v2 = _vel[:,:,0]**2
            vf = np.sum(_v2**0.5 *_wt*_distr) / np.sum(_wt*_distr)
            _cond = vf * self._eva2ms

        return _cond

    def check_info(self):
        print("nk = ", self.nk, "; nbnd = ", self.nbnd)
        print("Is semiconductor: ", self.is_semi)
        if self.is_semi:
            print("band index for valence bands: ", self.vb_index)
            print("band index for conduction bands: ", self.cb_index)
            print("Is electron quantities: ", self.is_elec)
            print("density of carriers: ", self.get_cden())
        else:
            print("band index for fermi-surface: ", self.cb_index)

    def reader_interp_dat(self, filename="tt_geninterp.dat", irrbz=False):

        with open(filename, 'r') as handle:
            data=handle.readlines()
        
        # guess nk and nbnd
        self.nk = int(data[-1].split()[0])
        if (not irrbz and self.nk != self.npts) or (irrbz and self.nk != self.irrnpts):
            raise Exception("num of pts from init and from doc inconsistent!")
        self.nbnd = [x.split()[0] == "1" \
            for x in data[:min(100,len(data))]].count(True)
        print("reader_interp_dat: nk = ", self.nk, "; nbnd = ", self.nbnd)

        self.vel = np.zeros((self.nbnd, self.npts, self.ndim)) 
        self.erg = np.zeros((self.nbnd, self.npts))

        if irrbz:
            _irrvel = np.zeros((self.nbnd, self.irrnpts, self.ndim))
            _irrerg = np.zeros((self.nbnd, self.irrnpts))
            for ik in range(self.irrnpts):
                for ibnd in range(self.nbnd):
                    fl = [float(x) \
                        for x in data[3+ik*self.nbnd+ibnd].strip('\n').split()]
                    _irrerg[ibnd,ik] = fl[4]
                    _irrvel[ibnd,ik] = fl[5:8] if self.ndim==3 else fl[5:7]
            ##### here assume we don't have ecut; not applicable when followed by ecut process
            for ikbz in range(self.npts):
                ikirrloc = self.b2i[ikbz]
                _symmatc = self.f2c(self.rotmat[ikbz])
                for ib in range(self.nbnd):
                    self.vel[ib,ikbz] = np.dot(_symmatc, _irrvel[ib,ikirrloc])
                    self.erg[ib,ikbz] = _irrerg[ib,ikirrloc]
        else:
            for ik in range(self.nk):
                for ibnd in range(self.nbnd):
                    fl = [float(x) \
                        for x in data[3+ik*self.nbnd+ibnd].strip('\n').split()]
                    self.erg[ibnd,ik] = fl[4]
                    self.vel[ibnd,ik] = fl[5:8] if self.ndim==3 else fl[5:7]

    def reader_selfen(self, filename='linewidth.elself.300.000K'):

        with open(filename, 'r') as handle:
            data=handle.readlines()

        #### iverbosity = 2 or 3
        has_mode = True if len(data[-1].split()) > 4 else False

        # guess nk and nbnd
        _nk = int(data[-1].split()[0])
        _nbnd = len(set([int(l.split()[1]) for l in data[2:]]))
        if has_mode: 
            _nmode = len(set([int(l.split()[3]) for l in data[2:]]))
        else:
            _nmode = 1
        # print("reader_selfen: nk = ", _nk, "; nbnd = ", _nbnd, 
        #     "; nmode = ", _nmode)
        if _nk != (self.irrnpts if self.has_irr else self.nk) :
            raise Exception("reader_selfen: inconsistent k points")
        if _nbnd != self.nbnd:
            raise Exception("reader_selfen: inconsistent nbnd")
        
        self.scat = np.zeros((_nbnd, _nk))
        self.scat_im = np.zeros((_nbnd, _nk, _nmode))

        _nbndmin = min([int(l.split()[1]) for l in data[2:]])
        for line in data[2:]:
            entry = line.split()
            ik = int(entry[0]) - 1
            ibnd = int(entry[1]) - _nbndmin
            if has_mode:
                imode = int(entry[3]) - 1
                self.scat_im[ibnd,ik,imode] = float(entry[4]) /1000
            else:
                self.scat_im[ibnd,ik,0] = float(entry[3]) /1000

        self.scat = np.sum(self.scat_im, axis=2)

    def writer_interp_kpt(self, kfx, filename="tt_geninterp.kpt"):

        with open(filename, 'w') as handle:
            handle.write("The .rst line is a comment (its maximum allowed     length is 500 characters).i\n")
            handle.write("crystal\n")
            handle.write("%i\n" %(len(kfx)))

            for i, kx in enumerate(kfx):
                handle.write("%i %20.12f %20.12f %20.12f \n" %(i+1,
                    kx[0], kx[1], kx[2] if self.ndim==3 else 0.0))
    
    def writer_kpt_dat(self, filename="kpt.dat", irrbz=False, spmod=''):

        if irrbz and (not self.has_irr):
            raise Exception("output irrbz but has no irrbz info.")

        _okfx = self.irrkfx if irrbz else self.kfx
        _onk = self.irrnpts if irrbz else self.nk

        if 'wan' in spmod:
            self.writer_interp_kpt(kfx=_okfx, filename=filename)
            # with open(filename, 'w') as handle:
            #     handle.write("The .rst line is a comment (its maximum allowed     length is 500 characters).i\n")
            #     handle.write("crystal\n")
            #     handle.write("%i\n" %(_onk))
            #     for i, kx in enumerate(_okfx):
            #         handle.write("%i %20.12f %20.12f %20.12f \n" %(i+1, kx[0], kx[1], kx[2] if self.ndim==3 else 0.0))
        else:
            with open(filename, 'w') as handle:
                handle.write("%d crystal\n" % _onk)
                for kx in _okfx:
                    handle.write("%20.12f %20.12f %20.12f 1.0\n" %(kx[0], kx[1],
                        kx[2] if self.ndim==3 else 0.0))

    def get_fine_kfx(self, afac=[2,2,2]):

        if not self.is_uniform: raise Exception("Need uniform grid.")

        # build up finer grid point
        _sub_fkl = []
        for x in [float(i)/afac[0] for i in range(afac[0])]:
            for y in [float(i)/afac[1] for i in range(afac[1])]:
                if self.ndim == 2:
                    _sub_fkl.append([x/self.ngd[0],y/self.ngd[1]])
                    continue
                for z in [float(i)/afac[2] for i in range(afac[2])]:
                    _sub_fkl.append([x/self.ngd[0],y/self.ngd[1],z/self.ngd[2]])
        _sub_fkl = np.array(_sub_fkl) 

        fsthick = self.ecut 
        _argl = np.where(np.any(np.abs(self.erg-self.fs)<fsthick,axis=0))[0]
        _ckfxl = self._get_kfl(_argl)

        fkl = []
        for ickfx in _ckfxl:
            for isub in _sub_fkl:
                fkl.append(ickfx+isub)
        fkl = np.array(fkl) 
        
        return fkl 

    def apply_ecut(self):
        """ apply energy cutoff for k points. """

        _irrkid = []
        for i in self.lockid.keys():
            kid = self.lockid[i]
            if np.any(np.abs(self.erg[:,kid]-self.fs)<self.ecut):
                _irrkid.append(i)

        _beforenk = self.nk

        self.nk = self.npts = len(_irrkid)
        _savekid = [self.lockid[i] for i in _irrkid]
        self.erg = self.erg[:,_savekid]
        self.vel = self.vel[:,_savekid,:]
        self.lockid = dict(zip(_irrkid, list(range(len(_irrkid)))))
        self.kfx = self.kfx[_savekid,:]

        print("apply energy cutoff: %d -> %d"%(_beforenk, self.nk))

    def set_vel_renorm(self):
        self.vel_rnm = np.array([[ ivel/np.linalg.norm(ivel) 
            for ivel in ibvel] for ibvel in self.vel])
















if __name__ == "__main__":
    ## test part
    gaas_e = BandStruct([50,50,50])
    gaas_e.reader_interp_dat("gaas/wan/tt_nk50.dat")
    gaas_e.set_fermi(fs=9.731, fsthick=0.30, is_elec=False)

    gaas_e.check_info()
    print(gaas_e.get_cden(), gaas_e.get_eff_mass(True))

    for iaf in [1,2,3,4]:
        fk_e = gaas_e.get_fine_kfx([iaf,iaf,iaf])
        print("iaf: ", iaf, np.shape(fk_e))
        gaas_af_e = BandStruct([50*iaf]*3, kfxl=fk_e)
        # gaas_af_e.writer_interp_kpt(kfx=fk_e, 
        #     filename="gaas/wan/tt_geninterp.kpt.h"+str(iaf))
        gaas_af_e.reader_interp_dat("gaas/wan/tt_geninterp.dat.h"+str(iaf))
        gaas_af_e.set_fermi(fs=9.731, fsthick=0.30, is_elec=False)
        print(gaas_af_e.get_cden(), gaas_af_e.get_eff_mass(True))

    # gaas_af8e = BandStruct([400,400,400], kfxl=fk_e)
    # gaas_af8e.reader_interp_dat("tt_nk50_af8_e.dat")
    # gaas_af8e.set_fermi(fs=10.3124, is_elec=True)

    # for i, ik in enumerate(gaas_af8e.lockid.keys()):
    #     print(i, "-th: ", ik, "->", gaas_af8e.lockid[ik])

    # gaas_af8e.check_info()
    # print(gaas_af8e.get_cden(), gaas_af8e.get_eff_mass(True))
