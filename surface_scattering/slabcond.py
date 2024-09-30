import numpy as np 
import matplotlib.pyplot as plt
from band import BandStruct
from grids import GridPts

Cyes = 0
Cno = 0

class SlabCond(BandStruct):
    """ 
    Main driver for calculating slab conductivity with relaxation approximation.

    Parameters:
    ----------
    ngd: list of int 
        Number of grid points in each desired direction.
    filvel: str
        Filename containing the electronic group velocity data.
    fs0: float, optional
        Initial Fermi surface value. Default is 0.
    fsthick: float, optional
        Thickness of the Fermi surface. Default is 0.15.
    assume_metal: bool, optional
        Whether to assume the material is metallic. Default is True.
    isibz: bool, optional
        Use irreducible Brillouin Zone. Default is True.
    filinfo: str, optional
        Filename for QE standard output with lattice symmetry matrics. Default is 'hex.info'.
    fillw: str, optional
        Filename for electron self-energy due to e-ph scattering. Default is None.
    filpwc: str, optional
        Filename for PWCOND output. Default is None.
    slabwidth: float, optional
        Width of the slab, in unit of Angstroms. Default is 1.
    bc: tuple of 0 or 1, optional
        Boundary conditions for the film surfaces. Each element is 0 (without surface) or 1 (with surface). Default is (0, 0).
    area: float, optional
        Volume of the bulk unit cell in units of A^3 (cubic Angstroms). Default is 1.
    ibrav: int, optional
        Bravais lattice index. See `grids.py` for detailed lattice vectors. Default is 4 for hexagonal lattice.
    spsym: bool, optional
        Special treatment of lattice symmetry. Options are 'nosym' to remove all lattice symmetry, 'noz' to remove symmetry related to the z-axis. Default is None.
    fixedp: float, optional
        Fixed specular parameter for e-s scattering. Mainly for test purpose. Default is None.
    fixedtau: float, optional 
        Fixed relaxation times. Mainly for test purpose. Default is None.
    fixrmat: bool, optional
        Whether to correct the reflection coefficient using velocity. Default is True.
    """

    def __init__(self, ngd, filvel, fs0=0, fsthick=0.15, assume_metal=True, isibz=True, filinfo='hex.info', fillw=None, filpwc=None, slabwidth=1, bc=(0,0), area=1, ibrav=4, spsym=None, fixedp=None, fixedtau=None, fixrmat=True):
        GridPts.__init__(self, ngd=ngd, ibrav=ibrav, spsym=spsym)

        self.reader_interp_dat(filename=filvel)

        self.set_fermi(fs=fs0, fsthick=fsthick, assume_metal=assume_metal)
        self.check_info()

        self.a = slabwidth # in angstrom
        self.bc = bc
        self.area = area
        self.fixedp = fixedp
        self.fixedtau = fixedtau

        self.fixrmat = fixrmat

        if isibz:
            self.apply_ecut()
            self.set_kbz2ibz(infofile=filinfo)
    
        if fillw is not None:
            self.reader_selfen(filename=fillw)

        if filpwc == None:
            self.is_simplefs = True
        else:
            self.is_simplefs = False
            self.reader_pwc(filpwc=filpwc)

        self.cz = np.linspace(0,1,num=self.ngd[-1], endpoint=False)
        self.has_Fbz = False



    def reader_pwc(self, filpwc='cond.out', spmat='noT'):
        ### init in writer_pwc_in
        try:
            self.pwc_ke
        except:
            self.writer_pwc_in(filpwcin=False)

        ### do read
        with open(filpwc, 'r') as fo:
            lines = fo.readlines()

        i = 0
        pwc_out, pwc_out_re = [], []
        while i < len(lines):
            if "---  ie =     " not in lines[i]:
                i += 1
                continue
            
            ii = int(lines[i].split()[3])-1; i+=1
            nch = int(lines[i].split()[6]); i+=3
            kzl = []
            for _ in range(nch):
                iz = float(lines[i].split()[0]); i+=1
                if iz < 0: iz = iz + 1.0
                kzl.append(iz)
            i += 2
            for _ in range(nch):
                iz = float(lines[i].split()[0]); i+=1
                if iz < 0: iz = iz + 1.0
                kzl.append(iz)
            
            if nch > 0:
                while not ("Band j to band i" in lines[i]):
                    i += 1

                i += 1
                rmat = np.zeros((nch*2, nch*2))
                for jch in range(nch):
                    i += 2
                    for ich in range(nch):
                        t = float(lines[i].split()[3])
                        r = float(lines[i].split()[4])
                        if spmat == 'noT':
                            r = r + t; t = 0.0
                        
                        rmat[jch, ich] = t
                        rmat[jch+nch, ich+nch] = t
                        rmat[jch, ich+nch] = r
                        rmat[jch+nch, ich] = r
                        i += 1

                rmat_reord = np.zeros((nch*2, nch*2))
                kz0 = self.pwc_ke[ii][2]
                argkzl_reord = np.argsort(kzl)
                kzl_reord = np.sort(kzl)
                ikz0_reord = np.argmin(np.abs(np.array(kzl_reord)-kz0))
                arg2 = list(argkzl_reord)
                ikz0 = arg2[ikz0_reord]
                arg2.pop(ikz0_reord)
                arg2.insert(0, ikz0)
                # kzl_reord = kzl[arg2]
                kzl_reord = [kzl[iarg] for iarg in arg2]

                for jch in range(nch*2):
                    for ich in range(nch*2):
                        jch_, ich_ = arg2[jch], arg2[ich]
                        rmat_reord[jch, ich] = rmat[jch_, ich_]

                pwc_out.append([kzl, np.copy(rmat)])
                pwc_out_re.append([kzl_reord, np.copy(rmat_reord)])
            else:
                pwc_out.append([[], np.zeros((0,0))])
                pwc_out_re.append([[], np.zeros((0,0))])

        self.pwc_out = pwc_out
        self.pwc_out_re = pwc_out_re


    def writer_pwc_in(self, filpwcin='kpt.dat'):
        _nk = self.irrnpts if self.has_irr else self.nk
        pwc_loc = np.ones((self.nbnd, _nk), dtype=int)
        pwc_loc = pwc_loc*-1

        pwc_ke = []
        ii = 0
        for ikbz in self.irrlockid.keys():
            ikloc = self.irrlockid[ikbz]
            for ib in range(self.nbnd):
                ierg = self.get_dsik((ib,ikbz), sv='erg') - self.fs

                if np.abs(ierg) < self.ecut:
                    if self.ndim == 3:
                        xx, yy, zz = self._get_kfl([ikbz])[0]
                    elif self.ndim == 2:
                        xx, zz = self._get_kfl([ikbz])[0]
                        yy = 0.0
                    pwc_ke.append([xx, yy, zz, ierg])
                    pwc_loc[ib, ikloc] = ii
                    ii += 1
        
        if filpwcin:
            with open(filpwcin, 'w') as fo:
                fo.write(f'{len(pwc_ke):5d}\n')
                for xx, yy, _, _ in pwc_ke:
                    fo.write(f'{xx:16.8f}{yy:16.8f}{1.0:10.4f}\n')
                fo.write(f'{len(pwc_ke):5d}\n')
                for _, _, _, ierg in pwc_ke:
                    fo.write(f'{ierg:16.8f}\n')

        self.pwc_loc = pwc_loc
        self.pwc_ke = pwc_ke

    def set_fermi(self, fs, fsthick=0.3, is_elec=True, assume_metal=False):
        bndlist = super().set_fermi(fs, fsthick, is_elec, assume_metal)
        self.vel2 = self.vel2[bndlist]
        self.erg2 = self.erg2[bndlist]
        
    def reader_interp_dat(self, filename="tt_geninterp.dat"):
        super().reader_interp_dat(filename=filename)
        if self.ndim == 3:
            self.vel2 = np.reshape(self.vel, (self.nbnd, self.ngd[0], self.ngd[1], self.ngd[2], self.ndim))
            self.erg2 = np.reshape(self.erg, (self.nbnd, self.ngd[0], self.ngd[1], self.ngd[2]))
        elif self.ndim == 2:
            self.vel2 = np.reshape(self.vel, (self.nbnd, self.ngd[0], self.ngd[1], self.ndim))
            self.erg2 = np.reshape(self.erg, (self.nbnd, self.ngd[0], self.ngd[1]))
    
    def reader_selfen(self, filename='linewidth.elself.300.000K'):
        """ mode resolved disabled"""
        super().reader_selfen(filename=filename)

        # self.scat = self.scat*2 #FIXME

        if self.fixedtau:
            self.scat[:,:] = self.fixedtau

            # self.scat[0,:] = 0.006
            # self.scat[1,:] = 0.0006
            # self.scat[2,:] = 0.06

        self.scatfb = self.scat[:,self.locb2i] if self.has_irr else self.scat
        # self.scat_imfb = self.scat_im[:,self.locb2i,:] if self.has_irr else self.scat_im

        if self.ndim == 3:
            self.scatfb2 = np.zeros((self.nbnd, self.ngd[0], self.ngd[1], self.ngd[2]))
            for i in self.lockid.keys():
                ix,iy,iz = self._get_ikil([i])[0]
                self.scatfb2[:,ix,iy,iz] = self.scatfb[:,self.lockid[i]]
        elif self.ndim == 2:
            self.scatfb2 = np.zeros((self.nbnd, self.ngd[0], self.ngd[1]))
            for i in self.lockid.keys():
                ix,iy = self._get_ikil([i])[0]
                self.scatfb2[:,ix,iy] = self.scatfb[:,self.lockid[i]]

    def get_velinz(self, ib_ikbz, dir=2):
        ib, ikbz = ib_ikbz
        if self.ndim == 3:
            ikx, iky, _ = self._get_ikil([ikbz])[0]
            return self.vel2[ib,ikx,iky,:,dir]
        elif self.ndim == 2:
            ikx, _ = self._get_ikil([ikbz])[0]
            return self.vel2[ib,ikx,:,dir]

    def get_erginz(self, ib_ikbz):
        ib, ikbz = ib_ikbz
        if self.ndim == 3:
            ikx, iky, _ = self._get_ikil([ikbz])[0]
            return self.erg2[ib,ikx,iky,:]
        elif self.ndim == 2:
            ikx, _ = self._get_ikil([ikbz])[0]
            return self.erg2[ib,ikx,:]

    def get_selfeninz(self, ib_ikbz):
        ib, ikbz = ib_ikbz
        if self.ndim == 3:
            ikx, iky, _ = self._get_ikil([ikbz])[0]
            return self.scatfb2[ib,ikx,iky,:]
        elif self.ndim == 2:
            ikx, _ = self._get_ikil([ikbz])[0]
            return self.scatfb2[ib,ikx,:]
    
    def get_xxinz(self, ib_ikbz, sv='erg'):
        if sv == 'selfen':
            cvals = self.get_selfeninz(ib_ikbz)
        elif sv == 'vx':
            cvals = self.get_velinz(ib_ikbz, dir=0)
        elif sv == 'vy':
            cvals = self.get_velinz(ib_ikbz, dir=1)
        elif sv == 'vz':
            cvals = self.get_velinz(ib_ikbz, dir=2)
        elif sv == 'erg':
            cvals = self.get_erginz(ib_ikbz)
        return cvals

    def find_roots_xxinz(self, ib_ikbz, z0=0, sv='erg', exclude=False):
        ### if exclude, we exclude the root == iz(ikbz)
        nz = self.ngd[-1]
        cz = np.linspace(0,1,num=nz+1,endpoint=True)
        cvals = self.get_xxinz(ib_ikbz, sv) - z0

        exroot = self._get_kfl([ib_ikbz[1]])[0][-1]

        roots = []
        #### if roots on points
        for i in np.where(cvals==0.0)[0]:
            iroot = cz[i]
            if exclude and np.abs(iroot-exroot)<1e-10:
                pass
            else:
                roots.append(iroot)

        #### if roots on interpolation
        ncvals = np.append(cvals, cvals[0])
        signs = ncvals[:-1]*ncvals[1:]
        for i in np.where(signs<-1e-8)[0]:
            iroot = cz[i]+(-ncvals[i])/(ncvals[i+1]-ncvals[i])*(cz[i+1]-cz[i])
            roots.append(iroot)

        return np.array(roots)

    def interpinz(self, ib_ikbz, sv='erg'):
        from scipy.interpolate import interp1d
        nz = self.ngd[-1]
        cz = np.linspace(0,1,num=nz+1,endpoint=True)
        cvals = self.get_xxinz(ib_ikbz, sv)
        ncvals = np.append(cvals, cvals[0])
        # ffun = interp1d(cz, ncvals, kind='quadratic')
        ffun = interp1d(cz, ncvals, kind='linear')
        return ffun

    def plot_fvals(self, ax, ib_ikbzs=[(0,0)], sv='erg', kw={}):
        for ib_ikbz in ib_ikbzs:
            cvals = self.get_xxinz(ib_ikbz, sv)
            if sv == 'erg': cvals = cvals - self.fs
            ax.scatter(self.cz, cvals, label=str(ib_ikbz))

            ffun = self.interpinz(ib_ikbz, sv)
            fz = np.linspace(0,1,num=3000)
            fvals = ffun(fz)
            if sv == 'erg': fvals = fvals - self.fs
            ax.plot(fz, fvals, **kw)

    def plot_fvals_all(self, ax, sv='erg', ibz=True, kw={}):
        if ibz:
            ib_ikbzs = []
            for ib in range(self.nbnd):
                for ikbz in self.irrlockid.keys():
                    ib_ikbzs.append((ib,ikbz))
        else:
            ib_ikbzs = []
            for ib in range(self.nbnd):
                for ikbz in self.lockid.keys():
                    ib_ikbzs.append((ib,ikbz))

        self.plot_fvals(ax=ax, ib_ikbzs=ib_ikbzs, sv=sv, kw=kw)


    ######## search in fullbz after ecut for random read
    def get_dsik(self, ib_ikbz, sv='erg'):
        ib, ikbz = ib_ikbz
        ikecut = self.lockid.get(ikbz)
        if ikecut is None: return 0.0

        if sv == 'selfen':
            return self.scatfb[ib, ikecut]
        elif sv == 'vx':
            return self.vel[ib, ikecut, 0]
        elif sv == 'vy':
            return self.vel[ib, ikecut, 1]
        elif sv == 'vz':
            return self.vel[ib, ikecut, 2]
        elif sv == 'erg':
            return self.erg[ib, ikecut]

    def get_mfpik(self, ib_ikbz, sv='mfp'):
        """ mfp in A """
        ib, ikbz = ib_ikbz
        ikecut = self.lockid.get(ikbz)
        if ikecut is None: return 0.0

        if sv == 'mfpz':
            vabs = np.abs(self.vel[ib,ikecut,2])
            selfen = self.get_dsik(ib_ikbz, sv='selfen')
            mfp = vabs/selfen / self._ev2inv_s * self._eva2ms / 1e-10
        elif sv == 'mfp':
            vabs = np.linalg.norm(self.vel[ib,ikecut])
            selfen = self.get_dsik(ib_ikbz, sv='selfen')
            mfp = vabs/selfen / self._ev2inv_s * self._eva2ms / 1e-10
        elif sv == 'Fx':
            if not self.has_Fbz: self.get_Fbz()
            mfp = np.abs(self.Fbz[ib,ikecut,0])
        elif sv == 'Fy':
            if not self.has_Fbz: self.get_Fbz()
            mfp = np.abs(self.Fbz[ib,ikecut,1])

        return mfp

    def get_dvik(self, ib_ikbz, sv='vx_'):
        """ generate the dv from rmat and velocity """
        ib, ikbz = ib_ikbz

        e0 = self.get_dsik(ib_ikbz, sv='erg')
        z0 = self._get_kfl([ikbz])[0][-1]

        #### store the initial state at 0
        ib_iz = [(ib, z0)]
        for ibf in range(self.nbnd):
            ibf_ikbz = (ibf, ikbz)
            if ib == ibf:
                fzs = self.find_roots_xxinz(ibf_ikbz, z0=e0, sv='erg', exclude=True)
            else:
                fzs = self.find_roots_xxinz(ibf_ikbz, z0=e0, sv='erg', exclude=False)
            for ifz in fzs: ib_iz.append((ibf,ifz))
            # if ikbz == 8137: ## fortest
            #     print("stop")

        ### fortest
        # print(f"ib: {ib:5d}, ikbz: {ikbz:5d}, nstates: {len(ib_iz):5d}")
        
        ##### resort by: 0 current, other by their kz
        new_ib_iz = sorted(ib_iz[1:], key=lambda x: x[1])
        new_ib_iz.insert(0, ib_iz[0])

        vxs, vys, kzl = [], [], []
        for iib, iiz in new_ib_iz:
            vx = self.interpinz((iib,ikbz), sv='vx')(iiz)
            vy = self.interpinz((iib,ikbz), sv='vy')(iiz)
            vxs.append(vx)
            vys.append(vy)
            kzl.append(iiz)

        ###### combine with rmat
        rmat = self.get_rmatik(ib_ikbz, kzl)
        initp = np.zeros((len(kzl),))
        initp[0] = 1.0
        finp = np.dot(rmat, initp)
        vxf = np.dot(finp, vxs)
        vyf = np.dot(finp, vys)
        vxi = np.dot(initp, vxs)
        vyi = np.dot(initp, vys)

        if 'vx_' in sv:
            dvik = vxf - vxi
        elif 'vy_' in sv:
            dvik = vyf - vyi
        elif 'vxy_' in sv:
            dvik = ((vxf-vxi)**2 + (vyf-vyi)**2)**0.5
        elif 'vxy2_' in sv:
            dvik = ((vxf-vxi)**2 + (vyf-vyi)**2)
        elif 'rrx_' in sv:
            dvik = (vxf - vxi) / np.abs(vxi)
        elif 'rry_' in sv:
            dvik = (vyf - vyi) / np.abs(vyi)
        elif 'rrxy_' in sv:
            dvik = ((vxf-vxi)**2 + (vyf-vyi)**2)**0.5 / np.linalg.norm([vxi, vyi])
        elif 'avgx_' in sv:
            dvik = np.average(vxs)
        elif 'avgy_' in sv:
            dvik = np.average(vys)
        elif 'avgxy_' in sv:
            dvik = np.linalg.norm([np.average(vxs), np.average(vys)])
        elif 'avgrx_' in sv:
            dvik = np.average(vxs)/np.max(np.abs(vxs))
        elif 'avgry_' in sv:
            dvik = np.average(vys)/np.max(np.abs(vys))
        elif 'avgrxy_' in sv:
            vnorms = [np.linalg.norm([vx,vy]) for vx, vy in zip(vxs, vys)]
            dvik = np.linalg.norm([np.average(vxs), np.average(vys)])/np.max(vnorms)
        elif 'varx_' in sv:
            dvik = np.var(vxs)
        elif 'vary_' in sv:
            dvik = np.var(vys)
        elif 'varxy_' in sv:
            vs = np.array([[vx,vy] for vx,vy in zip(vxs,vys)])
            dvik = np.mean(np.linalg.norm(vs-np.mean(vs,axis=0), axis=0)**2)
        elif 'stdx_' in sv:
            dvik = np.var(vxs)**0.5
        elif 'stdy_' in sv:
            dvik = np.var(vys)**0.5
        elif 'stdxy_' in sv:
            vs = np.array([[vx,vy] for vx,vy in zip(vxs,vys)])
            dvik = np.mean(np.linalg.norm(vs-np.mean(vs,axis=0), axis=0)**2)**0.5

        if 'abs_' in sv:
            dvik = np.abs(dvik)

        return dvik

    def get_dv_all(self, sv='vx_', fsth=0.15, filout='dv.dat'):
        
        dvall = []
        for ib in range(self.nbnd):
            for ikbz in self.lockid.keys():
                ib_ikbz = (ib, ikbz)

                erg0 = self.get_dsik(ib_ikbz, sv='erg')
                erg0 = erg0 - self.fs

                if np.abs(erg0) < fsth:
                    ikbzloc = self.lockid[ikbz]
                    kfx = self.kfx[ikbzloc]

                    if ('vx_' in sv) or ('vy_' in sv) or ('vxy_' in sv) or ('vxy2_' in sv) or ('rrx_' in sv) or ('rry_' in sv) or ('rrxy_' in sv) or ('avgx_' in sv) or ('avgy_' in sv) or ('avgxy_' in sv) or ('avgrx_' in sv) or ('avgry_' in sv) or ('avgrxy_' in sv) or ('var' in sv) or ('std' in sv):
                        dvik = self.get_dvik(ib_ikbz, sv=sv)
                    elif 'scat_' in sv:
                        dvik = self.get_dsik(ib_ikbz, sv='selfen')
                    elif 'mfp_' in sv:
                        dvik = self.get_mfpik(ib_ikbz, sv='mfp')
                    elif 'mfpz_' in sv:
                        dvik = self.get_mfpik(ib_ikbz, sv='mfpz')
                    elif 'Fx_' in sv:
                        dvik = self.get_mfpik(ib_ikbz, sv='Fx')
                    elif 'Fy_' in sv:
                        dvik = self.get_mfpik(ib_ikbz, sv='Fy')
                    else:
                        raise Exception("Not implemented in get_dv_all()")

                    if 'half_' in sv and kfx[2]>0.5:
                        pass
                    else:
                        dvall.append([kfx[:2], dvik])


        dvall = sorted(dvall, key=lambda x: x[0][0]*1e6+x[0][1])

        ##### do average to avoid overlap
        if not 'full_' in sv:
            new_dvall = []
            save_dvkf = dvall[0][0]
            new_dv = 0.0
            new_iv = 0
            for dvds in dvall:
                dvkf, dv = dvds[0], dvds[1]
                if np.allclose(dvkf, save_dvkf):
                    new_iv += 1
                    new_dv += dv
                else:
                    avg_dv = new_dv/new_iv
                    new_dvall.append([save_dvkf, avg_dv])

                    save_dvkf = dvkf
                    new_dv = dv
                    new_iv = 1

            print(f"{len(dvall)} dv points --> {len(new_dvall)}")
            dvall = new_dvall

        if filout:
            fo = open(filout, 'w')
            rv = np.transpose(self.bg) # reciprocal basis vectors
            for i in [0,1,2]:
                fo.write(f"{rv[i,0]:12.8f}{rv[i,1]:12.8f}{rv[i,2]:12.8f}\n")

            for dv in dvall:
                fo.write(f"{dv[0][0]:20.12f}{dv[0][1]:20.12f}{dv[1]:20.10e}\n")
            fo.close()

        return dvall

    def plot_dv_all(self, ax, sv='vx_', fsth=0.15, doabs=False, renorm=1.0, shift=False, kw={}):
        dvall = self.get_dv_all(sv=sv, fsth=fsth, filout=False)

        rv = np.transpose(self.bg[0:2,0:2])
        kfx = np.array([np.dot(dv_[0],rv) for dv_ in dvall])
        dv = np.array([dv_[1] for dv_ in dvall])
        dv = dv*renorm

        if doabs:
            dv = np.abs(dv)

        x = kfx[:,0]
        y = kfx[:,1]

        if shift:
            halfx = np.max(x)/2.0
            halfy = np.max(y)/2.0
            x[x>halfx] = x[x>halfx] - 2*halfx
            y[y>halfy] = y[y>halfy] - 2*halfy

        sc = ax.scatter(x, y, c=dv, **kw)
        ax.set_aspect('equal', adjustable='box')

        ax.set_xlim(np.min(kfx[:,0]), np.max(kfx[:,0]))
        ax.set_ylim(np.min(kfx[:,1]), np.max(kfx[:,1]))

        return sc

    def get_Aik(self, ib_ikbz):
        ib, ikbz = ib_ikbz

        e0 = self.get_dsik(ib_ikbz, sv='erg')
        z0 = self._get_kfl([ikbz])[0][-1]

        #### store the initial state at 0
        ib_iz = [(ib, z0)]
        for ibf in range(self.nbnd):
            ibf_ikbz = (ibf, ikbz)
            if ib == ibf:
                fzs = self.find_roots_xxinz(ibf_ikbz, z0=e0, sv='erg', exclude=True)
            else:
                fzs = self.find_roots_xxinz(ibf_ikbz, z0=e0, sv='erg', exclude=False)
            for ifz in fzs: ib_iz.append((ibf,ifz))
            # if ikbz == 8137: ## fortest
            #     print("stop")

        ### fortest
        # print(f"ib: {ib:5d}, ikbz: {ikbz:5d}, nstates: {len(ib_iz):5d}")
        
        ##### resort by: 0 current, other by their kz
        new_ib_iz = sorted(ib_iz[1:], key=lambda x: x[1])
        new_ib_iz.insert(0, ib_iz[0])

        ##### for test
        # if ib == 2 and ikbz == 1782:
        #     print('stop')

        Aik, expik, tauvz, kzl = [], [], [], []
        for iib, iiz in new_ib_iz:
            vx = self.interpinz((iib,ikbz), sv='vx')(iiz)
            vy = self.interpinz((iib,ikbz), sv='vy')(iiz)
            vz = self.interpinz((iib,ikbz), sv='vz')(iiz)
            scat = self.interpinz((iib,ikbz), sv='selfen')(iiz)
            erg = self.interpinz((iib,ikbz), sv='erg')(iiz) ## should be e0

            # assert scat > 1e-8 ### FIXME: errors when nkz too small
            tau = 1./scat if scat > 1e-8 else 0
            merg = np.exp((self.fs - erg)/self._kT)
            df = merg / (merg+1)**2 / self._kT
            Aik.append(np.array([vx,vy])*self._eva2ms * tau/self._ev2inv_s * df)
            expik.append(np.exp(-self.a*self._ang/np.abs(tau/self._ev2inv_s * vz*self._eva2ms)))
            tauvz.append(tau/self._ev2inv_s * vz*self._eva2ms / self._ang) # in unit of ang
            kzl.append(iiz)
        
        return np.array(Aik), np.array(expik), np.array(tauvz), np.array(kzl)

    def get_vzik(self, ib_ikbz):
        '''similar to get_Aik but only return np.abs(vz)'''
        ib, ikbz = ib_ikbz

        e0 = self.get_dsik(ib_ikbz, sv='erg')
        z0 = self._get_kfl([ikbz])[0][-1]

        #### store the initial state at 0
        ib_iz = [(ib, z0)]
        for ibf in range(self.nbnd):
            ibf_ikbz = (ibf, ikbz)
            if ib == ibf:
                fzs = self.find_roots_xxinz(ibf_ikbz, z0=e0, sv='erg', exclude=True)
            else:
                fzs = self.find_roots_xxinz(ibf_ikbz, z0=e0, sv='erg', exclude=False)
            for ifz in fzs: ib_iz.append((ibf,ifz))

        ##### resort by: 0 current, other by their kz
        new_ib_iz = sorted(ib_iz[1:], key=lambda x: x[1])
        new_ib_iz.insert(0, ib_iz[0])

        vzl = []
        for iib, iiz in new_ib_iz:
            vz = self.interpinz((iib,ikbz), sv='vz')(iiz)
            vzl.append(np.abs(vz))
        
        return np.array(vzl)

    def get_dfrmat(self, kzl): # created by genius ChatGPT
        # Sort kzl while keeping track of the original indices
        sorted_indices = sorted(range(len(kzl)), key=lambda k: kzl[k])
        
        # Initialize an empty matrix with zeros
        n = len(kzl)
        rmat = [[0] * n for _ in range(n)]
        
        # Create the pairs of lowest-highest, second-lowest-second-highest, etc.
        pairs = [(sorted_indices[i], sorted_indices[~i]) for i in range(n//2)]
        
        # If the length of kzl is odd, add the middle element to the pairs
        if n % 2 == 1:
            middle_index = sorted_indices[n // 2]
            pairs.append((middle_index, middle_index))
        
        # Check the pairs to determine which matrix elements should be set to 1
        for pair in pairs:
            rmat[pair[0]][pair[1]] = 1
            rmat[pair[1]][pair[0]] = 1
                    
        return rmat


    def get_rmatik(self, ib_ikbz, kzl):
        """ should reorder the ib, 0 is the current band, other band ordered by kz, from 0 to 1. """
        from scipy.linalg import circulant
        c = np.zeros((len(kzl),))
        c[-1] = 1
        # rmat_default =  circulant(c)
        rmat_default =  self.get_dfrmat(kzl)

        if self.is_simplefs == True:
            return rmat_default
        elif self.is_simplefs == False:
            ###### here we do several data alignment
            ### first we check if calculated from pwcond
            ib, ikbz = ib_ikbz
            ikirrloc = self.b2i[ikbz]
            ii = self.pwc_loc[ib, ikirrloc]
            if ii < 0:
                # print(f'not calculated in pwcond: ib:{ib}, ikbz:{ikbz}')
                return rmat_default
            kzl_, rmat_ = self.pwc_out_re[ii]
            nch2_ = len(kzl_)
            ### then check if nch matches
            nch2 = len(kzl)
            if nch2_ != nch2:
                # print(f'mismatch in nchannels*2: nch2:{nch2}, nch2_:{nch2_}')
                # print(self._get_kfl([ikbz]), ikbz)
                global Cno
                Cno += 1
                ## should be minority
                return rmat_default
            else:
                global Cyes
                Cyes += 1

            # assert np.allclose(kzl_, kzl, atol=0.2) #FIXME
            # if not np.allclose(kzl_, kzl, atol=0.2):
            #     print("error!!!")
            if self.fixrmat:
                vzl = self.get_vzik(ib_ikbz)
                scaled_rmat = rmat_ / vzl[np.newaxis,:]
                rmat_ = scaled_rmat / np.sum(scaled_rmat, axis=1)[:, np.newaxis]

            return rmat_

    def get_Fik(self, ib_ikbz):
        Aik, expik, tauvz, kzl = self.get_Aik(ib_ikbz)
        vz = self.get_dsik(ib_ikbz, sv='vz')

        if not (self.fixedp == None):
            assert (self.bc[0] == 1 and self.bc[1] == 1)
            p = self.fixedp
            f = -(1-p + p*(1-p)*expik[0])/(1-p**2*expik[0]**2)
            Fik = np.array([f,f])
            return Fik, Aik[0], tauvz[0]

        ### FIXED rmat intial and fianl states
        if self.bc[0] == 0 and self.bc[1] == 0:
            Fik = 0.0
        elif self.bc[0] == 1 and self.bc[1] == 0:
            if vz < 0: 
                Fik = 0.0
            else:
                rmat = self.get_rmatik(ib_ikbz, kzl)
                Fik = np.sum(Aik*rmat[:,0,np.newaxis], axis=0)/Aik[0] - 1
        elif self.bc[0] == 0 and self.bc[1] == 1:
            if vz > 0:
                Fik = 0.0
            else:
                rmat = self.get_rmatik(ib_ikbz, kzl)
                Fik = np.sum(Aik*rmat[:,0,np.newaxis], axis=0)/Aik[0] - 1
        elif self.bc[0] == 1 and self.bc[1] == 1:
            rmat = self.get_rmatik(ib_ikbz, kzl)
            Fik = np.zeros((2,))
            ##### for x
            LinMatA = np.diag(Aik[:,0]) - np.dot(np.transpose(rmat), np.diag(Aik[:,0]*expik))
            LinMatB = -Aik[:,0] + np.dot(Aik[:,0], rmat)
            try:
                x = np.linalg.solve(LinMatA, LinMatB)
                Fik[0] = x[0]
            except:
                Fik[0] = 0.0
            ##### for y
            LinMatA = np.diag(Aik[:,1]) - np.dot(np.transpose(rmat), np.diag(Aik[:,1]*expik))
            LinMatB = -Aik[:,1] + np.dot(Aik[:,1], rmat)
            try:
                x = np.linalg.solve(LinMatA, LinMatB)
                Fik[1] = x[0]
            except:
                Fik[1] = 0.0

            # if Fik[1] > -1000:
            #     ib, ikbz = ib_ikbz
            #     kfx = self._get_kfl([ikbz])[0]
            #     if abs(kfx[0]-0.766666667) < 1e-4:
            #         print(f'stop: kfx:{kfx}')

        return Fik, Aik[0], tauvz[0]

    def get_Fbz(self):
        nbnd = self.nbnd
        nkecut = self.nk
        # nkirr = self.irrnpts if self.has_irr else self.nk

        self.Fbz = np.zeros((nbnd, nkecut, 2))
        self.Abz = np.zeros((nbnd, nkecut, 2))
        self.tvzbz = np.zeros((nbnd, nkecut))
        for ib in range(nbnd):
            for ikbz in self.lockid.keys():
                ib_ikbz = (ib, ikbz)
                erg = self.get_dsik(ib_ikbz, sv='erg')

                if np.abs(erg-self.fs) < self.ecut:
                    ikbzecut = self.lockid[ikbz]
                    Fik, Aik, tvzik = self.get_Fik(ib_ikbz)
                    self.Fbz[ib, ikbzecut] = Fik
                    self.Abz[ib, ikbzecut] = Aik
                    self.tvzbz[ib, ikbzecut] = tvzik

        self.has_Fbz = True
        return self.Fbz

    def get_slabcondz(self, rz):
        if not self.has_Fbz: self.get_Fbz()

        expz = np.zeros((self.nbnd, self.nk))
        expz[self.tvzbz>0] = np.exp(-rz/self.tvzbz[self.tvzbz>0])
        expz[self.tvzbz<0] = np.exp((self.a-rz)/self.tvzbz[self.tvzbz<0])

        f1bz = self.Abz * (1+self.Fbz*expz[:,:,np.newaxis]) # in nbnd, nk, 2
        vxy = self.vel[:,:,0:2]
        condz = np.zeros((2,2))
        condz[0,0] = np.sum(vxy[:,:,0]*f1bz[:,:,0])
        condz[0,1] = np.sum(vxy[:,:,0]*f1bz[:,:,1])
        condz[1,0] = np.sum(vxy[:,:,1]*f1bz[:,:,0])
        condz[1,1] = np.sum(vxy[:,:,1]*f1bz[:,:,1])
        condz = condz * self._eva2ms / self.ngpts
        condz = condz * 2 # without soc, spin degeneracy

        if self.ndim == 2:
            condz = condz * self._e / (self.area*self._bohr**2)
        else:
            condz = condz * self._e / (self.area*self._bohr**3)

        return condz

    def get_slabcond(self, nrz=10, edir=None, ax=None, normz=False, normc=False,kw={}):
        if not self.has_Fbz: self.get_Fbz()

        self.rz = np.linspace(0, self.a, num=nrz, endpoint=True)
        self.condz = np.array([self.get_slabcondz(irz) for irz in self.rz])

        self.cond = np.sum(self.condz, axis=0)

        if edir is None: 
            edir = np.array([1,0])
            # edir = np.array([0,1])
        else:
            edir = np.array(edir)/np.linalg.norm(edir)
        condz1dir = np.array([np.dot(edir, np.dot(icz, edir)) for icz in self.condz])

        if ax is not None:
            _x = self.rz/np.max(self.rz) if normz else self.rz
            _y = condz1dir/normc if normc else condz1dir
            ax.plot(_x, _y, **kw)

        return self.rz, self.condz




if __name__ == '__main__':
    import lab

    ttkw = lab.cu111_final; NORMC=5.506e7
    # ttkw = lab.cu110_final; NORMC=5.044e7; NORMCY=5.044e7
    # ttkw = lab.cu001_final; NORMC=3.702e7
    # ttkw = lab.cu210_final; NORMC=4.9849e7; NORMCY=5.1537e7
    # ttkw = lab.cu001sup2_final; NORMC=4.415e7; NORMCY=4.164e7

    ########## test the conductivity solver
    fig, ax = plt.subplots()
    fo = open("output.dat", 'w') 
    for ISBW in [40,80,160,400,800,1200,2000,3000,4000]:
        metal1 = SlabCond(ngd=ttkw['ngd'], filvel=ttkw['filvel'], fs0=ttkw['fs0'], fsthick=ttkw['fsthick'], assume_metal=ttkw['assume_metal'], isibz=ttkw['isibz'], filinfo=ttkw['filinfo'], fillw=ttkw['fillw'], filpwc=ttkw['filpwc'], bc=ttkw['bc'], area=ttkw['area'], ibrav=ttkw['ibrav'], spsym=ttkw['spsym'],
        # fixedp=0.954801,
        # fixedtau=0.006,
        slabwidth=ISBW,
        )

        EDIR = [1,0]
        EDIR = [0,1]
        # EDIR = [1,1]

        svtag = 'vx_'
        svtag = 'vx_abs_'
        svtag = 'vy_abs_'
        svtag = 'vxy_abs_'
        # dv = metal1.get_dv_all(sv=svtag, fsth=0.15, filout='d'+svtag+'.out')
        # sc = metal1.plot_dv_all(ax, sv=svtag, fsth=0.15, kw={})

        # metal1.writer_pwc_in('/work/pwcond/m5_cu210/pwc/cond.ke.kz60')
        cond, _, _, _ = metal1.get_conductivity(plot=False, area=ttkw['area'])
        print(cond)

        rz, cond2 = metal1.get_slabcond(nrz=51, edir=EDIR, ax=ax, normz=False, normc=NORMC, kw={'label': f'{ISBW/10:.2e} A'})
        # print(cond2)
        _cond = np.mean(cond2, axis=0)
        # print(ISBW, ":\n", (_cond[0,1]+_cond[1,0])/2, ' S/m')
        print(ISBW, ":\n", (_cond[0,0]+_cond[1,1])/2, ' S/m')
        print(ISBW, ":\n", f'x: {_cond[0,0]:10.4e}  y: {_cond[1,1]:10.4e} S/m')
        # print(ISBW, ":\n", (_cond[0,0]+_cond[1,1])/2, ' S/m')
        print("boundary:\n", f'x: {(cond2[0,0,0]+cond2[-1,0,0])/2:10.4e}  y: {(cond2[0,1,1]+cond2[-1,1,1])/2:10.4e} S/m')


        ###### output condz
        diagcond = (cond2[:,0,0]+cond2[:,1,1])/2
        for iz, icx, icy in zip(rz, cond2[:,0,0], cond2[:,1,1]):
            fo.write(f"{iz:20.10f}{icx:20.10e}{icy:20.10e}\n")
            # fo.write(f"{icx:20.10e}\n")
        fo.write("\n")


    # plt.xlim(0, 400)
    # plt.ylim(0, 7e7)
    # plt.ylim(0.3, 1.1)
    plt.ylim(0, 1.4)
    plt.xlabel('slab position (A)')
    plt.ylabel('conductivity (1/ohm/m)')

    plt.legend()
    # plt.colorbar(sc)
    plt.show()
    plt.close()
