import numpy as np 

class GridPts():

    def _ikf_2d(self, ik):
        nx, ny = self.ngd[0:2]
        return [float(ik//ny)/nx, float(ik%ny)/ny]
    def _ikf_3d(self, ik):
        nx, ny, nz = self.ngd[0:3]
        return [float(ik//ny//nz)/nx, float(ik//nz%ny)/ny, float(ik%nz)/nz]
    def _get_kfl(self, ikl):
        if self.ndim == 2:
            return np.array([self._ikf_2d(ik) for ik in ikl])
        elif self.ndim == 3:
            return np.array([self._ikf_3d(ik) for ik in ikl])

    def _kindex_add(self, ik, iq):
        iki, iqi = self._get_ikil([ik,iq])
        ikqi = [(iki[i]+iqi[i])%self.ngd[i] for i in range(self.ndim)]
        if self.ndim == 2:
            nx, ny = self.ngd[0:2]
            ikq = ikqi[0]*ny + ikqi[1]
        elif self.ndim == 3:
            nx, ny, nz = self.ngd[0:3]
            ikq = ikqi[0]*ny*nz + ikqi[1]*nz + ikqi[2]
        return ikq

    def _get_kil(self, kfxl):
        if self.ndim == 2:
            nx, ny = self.ngd[0:2]
            return np.around(kfxl[:,0]*nx)*ny + np.around(kfxl[:,1]*ny)
        elif self.ndim == 3:
            nx, ny, nz = self.ngd[0:3]
            return np.around(kfxl[:,0]*nx)*ny*nz + np.around(kfxl[:,1]*ny)*nz +\
                np.around(kfxl[:,2]*nz) 

    def _get_ikil(self, ikl):
        if self.ndim == 2:
            nx, ny = self.ngd[0:2]
            return [[ik//ny, ik%ny] for ik in ikl]
        elif self.ndim == 3:
            nx, ny, nz = self.ngd[0:3]
            return [[ik//ny//nz, ik//nz%ny, ik%nz] for ik in ikl]

    def _set_uniform_grid(self):
        self.lockid = {}
        for i in range(self.ngpts):
            self.lockid[i] = i

        self.kfx = np.zeros((self.ngpts, self.ndim))
        if self.ndim == 2:
            counter = 0
            for x in [float(i)/self.ngd[0] for i in range(self.ngd[0])]:
                for y in [float(i)/self.ngd[1] for i in range(self.ngd[1])]:
                    self.kfx[counter] = [x,y]
                    counter += 1
        elif self.ndim == 3:
            counter = 0
            for x in [float(i)/self.ngd[0] for i in range(self.ngd[0])]:
                for y in [float(i)/self.ngd[1] for i in range(self.ngd[1])]:
                    for z in [float(i)/self.ngd[2] for i in range(self.ngd[2])]:
                        self.kfx[counter] = [x,y,z]
                        counter += 1

    def _set_irr_grid(self, infofile="info"):

        self.has_irr = True
        self.reader_info(infofile)

        _symm = self.sym_mat
        _lockid = np.ones((self.ngpts,), dtype=int)*-1
        _irrnpts = 0
        _irrkil = []

        for iq in range(self.ngpts):
            if _lockid[iq] >= 0:
                continue

            qx = self._get_kfl([iq])[0]
            sym_fkx = np.array([np.dot(isym, qx) for isym in _symm])
            sym_fkx = sym_fkx - np.floor(sym_fkx+1e-8)
            sym_ik = np.array(self._get_kil(sym_fkx), dtype=int)

            _lockid[sym_ik] = _irrnpts
            _irrkil.append(iq)
            _irrnpts += 1

        _irrkil = np.array(_irrkil)
        self.lockid = self.irrlockid = dict(zip(range(self.ngpts), _lockid))
        self.npts = self.irrnpts = _irrnpts
        self.kfx = self.irrkfx = self._get_kfl(_irrkil)

    def f2c(self, v):
        bg, bgi = self.bg, np.linalg.inv(self.bg)
        if len(np.shape(v)) == 1:
            vout = np.dot(bg, v)
        elif len(np.shape(v)) == 2:
            vout = np.dot(bg, np.dot(v, bgi))
        return vout

    def c2f(self, v):
        bg, bgi = self.bg, np.linalg.inv(self.bg)
        if len(np.shape(v)) == 1:
            vout = np.dot(bgi, v)
        elif len(np.shape(v)) == 2:
            vout = np.dot(bgi, np.dot(v, bg))
        return vout


    def __init__(self, ngd, ibrav=4, kfxl=None, kindl=None, spsym=None):
        self.ndim = len(ngd)
        self.ngd = ngd 
        self.ngpts = 1
        for i in range(self.ndim):
            self.ngpts *= ngd[i]

        self.lockid = {}
        if kfxl is not None:
            self.kfx = np.array(kfxl)
            for i, ik in enumerate(self._get_kil(kfxl)):
                self.lockid[ik] = i 
        elif kindl is not None:
            self.kfx = np.zeros((len(kindl, self.ndim)))
            for i, kx in enumerate(self._get_kfl(kindl)):
                self.lockid[kindl[i]] = i 
                self.kfx[i] = kx 
        else:
            self._set_uniform_grid()
            
        self.npts = len(self.lockid)
        self.is_uniform = True if self.npts == self.ngpts else False 
        self.has_irr = None

        self.ibrav = ibrav
        if ibrav == 4:
            self.at = np.array([[1.0,0,0], [-0.5,3**0.5/2,0], [0,0,5.]])
        elif ibrav == 6:
            self.at = np.array([[1.0,0,0], [0.0,1.0,0.0], [0,0,5.0]])
        elif ibrav == 62:
            self.at = np.array([[1.0,0,0], [0.0,2**0.5,0.0], [0,0,1.0]])
        elif ibrav == 2:
            self.at = np.array([[-0.5,0.0,0.5],[0.0,0.5,0.5],[-0.5,0.5,0.0]])
        elif ibrav == 704: ### predefined ibrav for cu210
            self.at = np.array([[1.0,0.0,0.0],[0.5,1.25**0.5,0.0],[0.0,0.0,5**0.5]])
        elif ibrav == 705: ### predefined ibrav for cu100 1x2x1 supercell
            self.at = np.array([[1.0,0.0,0.0],[0.0,2.0,0.0],[0.0,0.0,2**0.5]])
        else:
            raise Exception("ibrav not implemented! ")
        if self.ndim == 2: self.at = self.at[:2,:2] 
        self.bg = np.linalg.inv(self.at) 
        ## at[i,:] is the i-th direct basis vec
        ## bg[:,i] is the i-th reciprocal vec
        self.metric_g = np.array([[np.dot(self.at[i], self.at[j]) for j in range(self.ndim)] for i in range(self.ndim)])
        self.metric_gstar = np.linalg.inv(self.metric_g)

        self.spsym = spsym ### special tag for symmetry prohit


    def get_kindex(self):
        return np.array([i for i in self.lockid.keys()], dtype=int) 
    def get_kfrac(self):
        return np.array(self.kfx) 

    def reader_info_brav(self, FILENAME="info"):
        fo = open(FILENAME, 'r')
        line = fo.readline()
        while line:
            line = line.split()
            if not line:
                line = fo.readline()
                continue

            if line[1] == 'symmetry':
                self.nsym = int(line[0])
                self.sym_mat = np.zeros((self.nsym, 3, 3))
                for i in range(self.nsym//6):
                    line = fo.readline()
                    line = fo.readline().split()
                    self.sym_mat[i*6:(i+1)*6, 0] = np.reshape([int(s) 
                        for s in line], (6,3))
                    line = fo.readline().split()
                    self.sym_mat[i*6:(i+1)*6, 1] = np.reshape([int(s) 
                        for s in line], (6,3))
                    line = fo.readline().split()
                    self.sym_mat[i*6:(i+1)*6, 2] = np.reshape([int(s) 
                        for s in line], (6,3))

                if self.nsym%6 != 0:
                    tmp = self.nsym%6
                    line = fo.readline()
                    line = fo.readline().split()
                    self.sym_mat[i*6:i*6+tmp, 0] = np.reshape([int(s) 
                        for s in line], (tmp,3))
                    line = fo.readline().split()
                    self.sym_mat[i*6:i*6+tmp, 1] = np.reshape([int(s) 
                        for s in line], (tmp,3))
                    line = fo.readline().split()
                    self.sym_mat[i*6:i*6+tmp, 2] = np.reshape([int(s) 
                        for s in line], (tmp,3))

            line = fo.readline()

        fo.close()

        if self.ndim == 2:
            self.sym_mat = self.sym_mat[:,0:2,0:2]

    def reader_info(self, FILENAME='scf.out'):
        def search_patt(p, l):
            return [(p in et) for et in l].index(True)
        def sl2f(sl):
            return [float(et) for et in sl]

        if self.spsym == 'nosym':
            nsym = 1
            sym_mat = np.array([[[1,0,0],[0,1,0],[0,0,1]]])
            if self.ndim == 2: sym_mat = sym_mat[:,:2,:2]
            print(f"{nsym} sym. enforced.")
        else:    
            with open(FILENAME, 'r') as fo:
                lines = fo.readlines()

            try:
                ind0 = search_patt('Sym. Ops.', lines)
                has_inv = ('with inversion' in lines[ind0])
                nsym_ = int(lines[ind0].split()[0])
            except:
                ind0 = search_patt('No symmetry found', lines)
                has_inv = False
                nsym_ = 1

            nsym = nsym_ if has_inv else 2*nsym_
            sym_mat = np.zeros((nsym,3,3))

            isym = -1
            nstep = 7
            ind = ind0
            while (isym < nsym_ - 1):
                ind = search_patt('isym = ', lines[ind:ind+nstep])+ind
                isym = int(lines[ind].split()[2]) - 1
                sym_mat[isym,:,0] = sl2f(lines[ind+2][20:50].split())
                sym_mat[isym,:,1] = sl2f(lines[ind+3][20:50].split())
                sym_mat[isym,:,2] = sl2f(lines[ind+4][20:50].split())
                ind += nstep

            if not has_inv:
                for isym in range(nsym_):
                    sym_mat[isym+nsym_] = np.dot(np.diag([-1,-1,-1]), sym_mat[isym])

            if self.ndim == 2: sym_mat = sym_mat[:,:2,:2]

            ##### change the basis of sym_mat from direct to reciprocal
            for isym in range(nsym):
                sym_mat[isym] = np.dot(self.metric_g, np.dot(sym_mat[isym], self.metric_gstar))
            
            print(f"{nsym} sym. found. with inv.")

            if self.spsym == 'noz':
                z1ind = [i for i in range(nsym) if np.abs(sym_mat[i,-1,-1]-1)<1e-6]
                nsym = len(z1ind)
                sym_mat = sym_mat[z1ind]
                print(f"{nsym} sym. found. with z -> z")

        self.nsym = nsym
        self.sym_mat = sym_mat

    def set_kbz2ibz(self, infofile="info"):

        self.has_irr = True
        self.reader_info(infofile)

        _irrkid = []
        _savekid = []
        _tag = dict.fromkeys(self.lockid.keys(), -1)
        for ik in self.lockid.keys():
            i = self.lockid[ik]
            if self.npts > 100 and i%(self.npts//100) == 0:
                print("\ri: ", i, "/", self.npts, end='')

            if _tag[ik] >= 0:
                continue

            _symm = self.sym_mat
            ik = int(ik)
            kx = self.kfx[i]

            sym_fkx = np.array([np.dot(_symm[i],kx) 
                for i in range(self.nsym)])
            sym_fkx = sym_fkx - np.floor(sym_fkx+1e-8)
            sym_ik = np.array(self._get_kil(sym_fkx), dtype=int)

            if sym_ik[0] != ik or ik == 2704:
                print("?")

            locirrk = [_irrkid.index(i) for i in sym_ik if i in _irrkid]
            if locirrk:
                assert len(set(locirrk)) == 1
                for j in sym_ik:
                    if j in _tag.keys(): _tag[j] = locirrk[0]
            else:
                _irrkid.append(ik)
                _savekid.append(i)
                for j in sym_ik:
                    if j in _tag.keys(): _tag[j] = len(_irrkid)-1

        lstrng = list(range(len(_irrkid)))
        self.irrnpts = len(_irrkid)
        self.irrlockid = dict(zip(_irrkid, lstrng))
        self.irrkfx = self.kfx[_savekid]
        self.b2i = _tag
        self.locb2i = list(_tag.values())

        _i2locb = []
        for ikbz in self.irrlockid.keys():
            _i2locb.append(self.lockid[ikbz])
        self.i2locb = _i2locb

        self.set_rotmat()
        print("\nnumber of k:", self.npts, "->", self.irrnpts)
    
    def set_rotmat(self):

        _irot = [] ## rotation matrix index, irrkfx to bckfx, for each bckfx
        _symm = self.sym_mat

        for ikbzloc, ikbz in enumerate(list(self.b2i.keys())):
            ikirrloc = self.b2i[ikbz]
            # ikirrloc = self.locb2i[ikbzloc]

            # kbcfx = self.kfx[ikbzloc]
            kirrfx = self.irrkfx[ikirrloc]

            sym_fkx = np.array([np.dot(_symm[i], kirrfx) for i in range(self.nsym)])
            sym_fkx = sym_fkx - np.floor(sym_fkx+1e-8)
            sym_ik = list(np.array(self._get_kil(sym_fkx), dtype=int))

            _irot.append(sym_ik.index(ikbz))

        self.irot = _irot
        self.rotmat = _symm[_irot]



if __name__ == '__main__':
    tt = GridPts([6,6,3], ibrav=4)
    # tt.set_kbz2ibz(infofile='/work/pwcond/m2_cu/symm_test/scf.out')
    tt.set_kbz2ibz(infofile='/work/pwcond/m2_cu/symm_test/test2/scf.out')
    print(tt.sym_mat[1])
