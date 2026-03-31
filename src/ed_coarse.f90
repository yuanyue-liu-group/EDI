MODULE ed_coarse
  USE kinds, ONLY : dp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: load_supercell_pot, load_pot_from_file, read_filukk_edi
  PUBLIC :: ed_fine_interp, ed_fine_interp_offdiag, read_edmatw_file
  PUBLIC :: read_edmatw_2d_file
  PUBLIC :: ed_interp_from_file
  PUBLIC :: ed_coarse_full_q
  PUBLIC :: ed_fine_interp_2d
  PUBLIC :: ed_direct_from_files, generate_nscf_input
  PUBLIC :: build_vcolin_aligned, build_vcolin_corealign, write_vcolin_cube

CONTAINS

  ! [ed_coarse_calc and compute_edmatkq_fast removed — dead code,
  !  all matrix elements now go through ed_coarse_full_q (double-FT)]

  SUBROUTINE ed_fine_interp(nbndsub, nrr, irvec, ndegen, chw, edmatw, &
                             nk1f, nk2f, nk3f, prefix_in)
    !-----------------------------------------------------------------------
    ! Part B: Wannier → Bloch interpolation on fine k-grid.
    !
    ! For each k_fine on the fine grid:
    !   1. cfac = exp(ik_fine·R) for all R
    !   2. Diagonalize H(k_fine) → eigenvalues ε(k_fine), eigenvectors U(k_fine)
    !   3. Interpolate M: M_W(k_fine) = Σ_R exp(ik·R)/ndegen(R) · M(R)
    !   4. Rotate to Bloch: M_B = U† · M_W · U
    !   5. Output |M_B(n,m,k)|² for transport
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode, stdout
    USE constants, ONLY : rytoev
    USE wan2bloch_edi, ONLY : get_cfac, hamwan2bloch_with_evec, edmatwan2bloch
    USE ep_constants, ONLY : czero, cone
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nbndsub, nrr
    INTEGER, INTENT(IN) :: irvec(3, nrr), ndegen(nrr)
    COMPLEX(dp), INTENT(IN) :: chw(nbndsub, nbndsub, nrr)    ! H(R)
    COMPLEX(dp), INTENT(IN) :: edmatw(nbndsub, nbndsub, nrr)  ! M(R)
    INTEGER, INTENT(IN) :: nk1f, nk2f, nk3f
    CHARACTER(LEN=*), INTENT(IN) :: prefix_in

    INTEGER :: ik, i1, i2, i3, nktotf, ibnd, jbnd, iunit
    REAL(dp) :: xk_cryst(3)
    REAL(dp) :: eig(nbndsub)
    COMPLEX(dp) :: evec(nbndsub, nbndsub)
    COMPLEX(dp) :: cfac(nrr)
    COMPLEX(dp) :: edmatf_w(nbndsub, nbndsub)   ! M in Wannier basis at k
    COMPLEX(dp) :: edmatf_b(nbndsub, nbndsub)   ! M in Bloch basis at k
    COMPLEX(dp) :: tmp(nbndsub, nbndsub)
    CHARACTER(LEN=256) :: fname

    nktotf = nk1f * nk2f * nk3f

    IF (ionode) THEN
       WRITE(stdout, '(/,5X,A)') REPEAT('=', 60)
       WRITE(stdout, '(5X,A)')   'Part B: Wannier interpolation on fine grid'
       WRITE(stdout, '(5X,A)')   REPEAT('=', 60)
       WRITE(stdout, '(5X,A,3I4,A,I8)') 'Fine k-grid: ', nk1f, nk2f, nk3f, &
            '  total: ', nktotf
       FLUSH(stdout)
    ENDIF

    ! Open output file for interpolated matrix elements
    iunit = 86
    IF (ionode) THEN
       fname = TRIM(prefix_in) // '_edmatf.dat'
       OPEN(iunit, FILE=TRIM(fname), FORM='formatted')
       WRITE(iunit, '(A)') '# Interpolated e-d matrix elements |M(n,m,k)|^2'
       WRITE(iunit, '(A)') '# ik  kx  ky  kz  ibnd  jbnd  |M|^2(Ry^2)  Re(M)  Im(M)'
    ENDIF

    ik = 0
    DO i3 = 0, nk3f - 1
       DO i2 = 0, nk2f - 1
          DO i1 = 0, nk1f - 1
             ik = ik + 1
             xk_cryst(1) = DBLE(i1) / DBLE(nk1f)
             xk_cryst(2) = DBLE(i2) / DBLE(nk2f)
             xk_cryst(3) = DBLE(i3) / DBLE(nk3f)

             ! Phase factors exp(ik·R)
             CALL get_cfac(nrr, irvec, xk_cryst, cfac)

             ! Diagonalize H(k) → eigenvalues + eigenvectors U(k)
             CALL hamwan2bloch_with_evec(nbndsub, nrr, ndegen, eig, evec, chw, cfac)

             ! Interpolate M(R) → M_W(k) in Wannier basis
             CALL edmatwan2bloch(nbndsub, nrr, ndegen, edmatw, cfac, edmatf_w)

             ! Rotate to Bloch basis: M_B = U† · M_W · U
             ! tmp = M_W · U
             CALL ZGEMM('N', 'N', nbndsub, nbndsub, nbndsub, &
                         cone, edmatf_w, nbndsub, evec, nbndsub, czero, tmp, nbndsub)
             ! M_B = U† · tmp
             CALL ZGEMM('C', 'N', nbndsub, nbndsub, nbndsub, &
                         cone, evec, nbndsub, tmp, nbndsub, czero, edmatf_b, nbndsub)

             ! Write to file
             IF (ionode) THEN
                DO ibnd = 1, nbndsub
                   DO jbnd = 1, nbndsub
                      WRITE(iunit, '(I8,3F10.5,2I4,3ES16.8)') ik, xk_cryst, &
                           ibnd, jbnd, &
                           ABS(edmatf_b(ibnd, jbnd))**2, &
                           REAL(edmatf_b(ibnd, jbnd)), AIMAG(edmatf_b(ibnd, jbnd))
                   ENDDO
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    IF (ionode) THEN
       CLOSE(iunit)
       WRITE(stdout, '(5X,A,I8,A)') 'Interpolated M at ', nktotf, ' fine k-points'
       WRITE(stdout, '(5X,A,A)') 'Written to ', TRIM(fname)
       WRITE(stdout, '(5X,A)') REPEAT('=', 60)
       FLUSH(stdout)
    ENDIF

  END SUBROUTINE ed_fine_interp

  SUBROUTINE ed_fine_interp_offdiag(nbndsub, nrr, irvec, ndegen, chw, edmatw, &
                                     nk1f, nk2f, nk3f, prefix_in)
    !-----------------------------------------------------------------------
    ! Part B extension: Off-diagonal k-point interpolation.
    !
    ! For a fixed initial k_i, compute M(n,m, k_i, k_f) for ALL k_f
    ! on the fine grid. Uses momentum transfer q = k_f - k_i:
    !
    !   M_W(i,j,q) = Σ_R exp(iq·R)/ndegen(R) · M(i,j,R)
    !   M_B(n,m,k_i,k_f) = U†(k_i) · M_W(q) · U(k_f)
    !
    ! Output: prefix_edmat_scatter.dat
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode, stdout
    USE wan2bloch_edi, ONLY : get_cfac, hamwan2bloch_with_evec, edmatwan2bloch
    USE ep_constants, ONLY : czero, cone
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nbndsub, nrr
    INTEGER, INTENT(IN) :: irvec(3, nrr), ndegen(nrr)
    COMPLEX(dp), INTENT(IN) :: chw(nbndsub, nbndsub, nrr)
    COMPLEX(dp), INTENT(IN) :: edmatw(nbndsub, nbndsub, nrr)
    INTEGER, INTENT(IN) :: nk1f, nk2f, nk3f
    CHARACTER(LEN=*), INTENT(IN) :: prefix_in

    INTEGER :: ik, i1, i2, i3, nktotf, ibnd, jbnd, ir, iunit
    REAL(dp) :: xk_i(3), xk_f(3), xq(3)
    REAL(dp) :: eig_i(nbndsub), eig_f(nbndsub)
    COMPLEX(dp) :: evec_i(nbndsub, nbndsub), evec_f(nbndsub, nbndsub)
    COMPLEX(dp) :: cfac_i(nrr), cfac_f(nrr), cfac_q(nrr)
    COMPLEX(dp) :: edmatf_w(nbndsub, nbndsub)
    COMPLEX(dp) :: edmatf_b(nbndsub, nbndsub)
    COMPLEX(dp) :: tmp(nbndsub, nbndsub)
    CHARACTER(LEN=256) :: fname
    REAL(dp), PARAMETER :: twopi = 6.283185307179586_dp
    COMPLEX(dp), PARAMETER :: ci = (0.0_dp, 1.0_dp)

    ! K point in crystal coordinates (MoS2): (1/3, 1/3, 0)
    xk_i = (/1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 0.0_dp/)

    nktotf = nk1f * nk2f * nk3f

    IF (ionode) THEN
       WRITE(stdout, '(/,5X,A)') REPEAT('=', 60)
       WRITE(stdout, '(5X,A)')   'Part B: Off-diagonal scattering from K point'
       WRITE(stdout, '(5X,A)')   REPEAT('=', 60)
       WRITE(stdout, '(5X,A,3F10.5)') 'Initial k (crystal): ', xk_i
       WRITE(stdout, '(5X,A,3I4,A,I8)') 'Fine k-grid: ', nk1f, nk2f, nk3f, &
            '  total: ', nktotf
       FLUSH(stdout)
    ENDIF

    ! Diagonalize H(k_i) once → eigenvalues + eigenvectors at initial k
    CALL get_cfac(nrr, irvec, xk_i, cfac_i)
    CALL hamwan2bloch_with_evec(nbndsub, nrr, ndegen, eig_i, evec_i, chw, cfac_i)

    ! Find CBM index at K (band with smallest eigenvalue above gap)
    ! For MoS2 with 11 bands, CBM at K is typically band 8
    IF (ionode) THEN
       WRITE(stdout, '(5X,A)') 'Eigenvalues at K (Ry):'
       DO ibnd = 1, nbndsub
          WRITE(stdout, '(5X,A,I3,A,F12.6)') '  band ', ibnd, ': ', eig_i(ibnd)
       ENDDO
       FLUSH(stdout)
    ENDIF

    ! Open output file
    iunit = 85
    IF (ionode) THEN
       fname = TRIM(prefix_in) // '_edmat_scatter.dat'
       OPEN(iunit, FILE=TRIM(fname), FORM='formatted')
       WRITE(iunit, '(A)') '# Off-diagonal scattering: M(n, m, k_i=K, k_f)'
       WRITE(iunit, '(A,3F10.5)') '# k_initial (crystal) = ', xk_i
       WRITE(iunit, '(A)') '# ik  kfx  kfy  kfz  ibnd  jbnd  |M|^2  Re(M)  Im(M)'
    ENDIF

    ik = 0
    DO i3 = 0, nk3f - 1
       DO i2 = 0, nk2f - 1
          DO i1 = 0, nk1f - 1
             ik = ik + 1
             xk_f(1) = DBLE(i1) / DBLE(nk1f)
             xk_f(2) = DBLE(i2) / DBLE(nk2f)
             xk_f(3) = DBLE(i3) / DBLE(nk3f)

             ! Momentum transfer q = k_f - k_i
             xq = xk_f - xk_i

             ! Phase factors for q: exp(iq·R)
             DO ir = 1, nrr
                cfac_q(ir) = EXP(ci * twopi * DOT_PRODUCT(xq, DBLE(irvec(:, ir))))
             ENDDO

             ! Interpolate M(R) at q: M_W(q) = Σ_R exp(iq·R)/ndegen · M(R)
             CALL edmatwan2bloch(nbndsub, nrr, ndegen, edmatw, cfac_q, edmatf_w)

             ! Diagonalize H(k_f) → eigenvectors U(k_f)
             CALL get_cfac(nrr, irvec, xk_f, cfac_f)
             CALL hamwan2bloch_with_evec(nbndsub, nrr, ndegen, eig_f, evec_f, chw, cfac_f)

             ! Rotate to Bloch: M_B = U†(k_i) · M_W(q) · U(k_f)
             ! tmp = M_W · U(k_f)
             CALL ZGEMM('N', 'N', nbndsub, nbndsub, nbndsub, &
                          cone, edmatf_w, nbndsub, evec_f, nbndsub, czero, tmp, nbndsub)
             ! M_B = U†(k_i) · tmp
             CALL ZGEMM('C', 'N', nbndsub, nbndsub, nbndsub, &
                          cone, evec_i, nbndsub, tmp, nbndsub, czero, edmatf_b, nbndsub)

             ! Write all band pairs
             IF (ionode) THEN
                DO ibnd = 1, nbndsub
                   DO jbnd = 1, nbndsub
                      WRITE(iunit, '(I8,3F10.5,2I4,3ES16.8)') ik, xk_f, &
                           ibnd, jbnd, &
                           ABS(edmatf_b(ibnd, jbnd))**2, &
                           REAL(edmatf_b(ibnd, jbnd)), AIMAG(edmatf_b(ibnd, jbnd))
                   ENDDO
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    IF (ionode) THEN
       CLOSE(iunit)
       WRITE(stdout, '(5X,A,I8,A)') 'Computed M(K -> k_f) at ', nktotf, ' fine k-points'
       WRITE(stdout, '(5X,A,A)') 'Written to ', TRIM(fname)
       WRITE(stdout, '(5X,A)') REPEAT('=', 60)
       FLUSH(stdout)
    ENDIF

  END SUBROUTINE ed_fine_interp_offdiag


  SUBROUTINE ed_coarse_full_q(nbndsub, nrr_k, irvec_k, ndegen_k, &
                               edmatw_2d, prefix_in)
    !-----------------------------------------------------------------------
    ! Full double-FT: M(k_i, k_f) for ALL (k_i, k_f) pairs on the coarse grid.
    ! Computes M(R,Rp) via (paper Eq.5):
    !   M_W(k_i, k_f) = U†(k_i) · M_B(k_i, k_f) · U(k_f)
    !   M(R,Rp) = (1/Nk²) Σ_{ki,kf} exp(+iki·R) exp(-ikf·Rp) M_W
    !
    ! Pool-parallel: each pool handles its local k_i, broadcasts ψ(k_f).
    ! No clean_pw, no read_file_new — QE state preserved.
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode, ionode_id, stdout
    USE mp, ONLY : mp_sum
    USE mp_world, ONLY : world_comm
    USE mp_global, ONLY : inter_pool_comm
    USE mp_pools, ONLY : npool, my_pool_id, intra_pool_comm
    USE cell_base, ONLY : at
    USE klist, ONLY : xk, igk_k, ngk, nkstot, nks
    USE wvfct, ONLY : npwx, nbnd
    USE noncollin_module, ONLY : npol, lspinorb
    USE fft_base, ONLY : dffts
    USE fft_interfaces, ONLY : invfft
    USE io_files, ONLY : restart_dir
    USE pw_restart_new, ONLY : read_collected_wfc
    USE wann_common, ONLY : u_mat, n_wannier, num_bands, excluded_band
    USE edic_mod, ONLY : V_d, V_p, V_colin
    USE uspp_param, ONLY : nh
    USE uspp, ONLY : dvan, dvan_so
    USE becmod, ONLY : bec_type, calbec, allocate_bec_type, deallocate_bec_type
    USE ep_constants, ONLY : czero, cone
    USE constants, ONLY : tpi
    USE parallelism, ONLY : fkbounds
    USE wigner_seitz_edi, ONLY : wigner_seitz_k
    USE parallel_include
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nbndsub, nrr_k
    INTEGER, INTENT(IN) :: irvec_k(3, nrr_k), ndegen_k(nrr_k)
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: edmatw_2d(:,:,:,:)
    CHARACTER(LEN=*), INTENT(IN) :: prefix_in

    INTEGER :: ik, ib_i, ib_j, ir, ir_e, ir_p
    INTEGER :: irx, iry, irz, ir1mod, ir2mod, ir3mod, inr, irnmod
    INTEGER :: ibnd, nbnd_kept, lower_bnd, upper_bnd
    INTEGER :: ikf_p, ikf_l, ierr_mpi
    INTEGER :: my_rank, num_ranks, nkf, nrr_p
    INTEGER :: na, nt, ih, jh, ikb, jkb, ijkb0, nkb_d, nkb_p
    REAL(dp) :: xkf_cryst(3), rdotk, rdotkf, arg, d1, d2, d3
    REAL(dp), ALLOCATABLE :: xk_cryst_loc(:,:), xk_cryst_all(:,:)
    COMPLEX(dp) :: mlocal, mnl_d, mnl_p, phase
    COMPLEX(dp), ALLOCATABLE :: V_folded_kf(:)
    COMPLEX(dp), ALLOCATABLE :: psir_k(:,:,:), psir_kq(:,:,:), psic_tmp(:)
    COMPLEX(dp), ALLOCATABLE :: evc_tmp(:,:)
    COMPLEX(dp), ALLOCATABLE :: edmatkq(:,:), edmatkq_loc(:,:), edmatkq_nl(:,:)
    COMPLEX(dp), ALLOCATABLE :: edms(:,:), eptmp(:,:)
    COMPLEX(dp), ALLOCATABLE :: edmat_bloch(:,:,:,:)
    COMPLEX(dp), ALLOCATABLE :: edmat_bloch_loc(:,:,:,:)
    COMPLEX(dp), ALLOCATABLE :: edmat_bloch_nl(:,:,:,:)
    COMPLEX(dp), ALLOCATABLE :: vkb_d(:,:), vkb_p(:,:)
    COMPLEX(dp), ALLOCATABLE :: cu_k(:,:), cu_kq(:,:)
    COMPLEX(dp), ALLOCATABLE :: cu_all(:,:,:)
    INTEGER, ALLOCATABLE :: band_map(:)

    ! q-grid = k-grid (every k-point is a valid momentum transfer)
    nkf = nkstot
    nrr_p = nrr_k  ! R and R' on the same WS grid
    CALL MPI_Comm_rank(world_comm, my_rank, ierr_mpi)
    CALL MPI_Comm_size(world_comm, num_ranks, ierr_mpi)
    CALL fkbounds(nkstot, lower_bnd, upper_bnd)

    IF (ionode) THEN
       WRITE(stdout, '(/,5X,A)') REPEAT('=', 60)
       WRITE(stdout, '(5X,A)')   'Full double-FT: M(k, k'') for all (k,k'') pairs'
       WRITE(stdout, '(5X,A)')   REPEAT('=', 60)
       WRITE(stdout, '(5X,A,I6)') 'Nk_i = Nk_f = ', nkstot
       WRITE(stdout, '(5X,A,I4,A,I4)') 'nks=', nks, '  npool=', npool
       WRITE(stdout, '(5X,A,I6)') 'nrr (R = Rp)=', nrr_k
       FLUSH(stdout)
    ENDIF

    ! Build band mapping
    nbnd_kept = 0
    ALLOCATE(band_map(nbnd))
    DO ibnd = 1, nbnd
       IF (.NOT. excluded_band(ibnd)) THEN
          nbnd_kept = nbnd_kept + 1
          band_map(nbnd_kept) = ibnd
       ENDIF
    ENDDO

    ! Gather xk in crystal coords for ALL k-points (small array)
    ALLOCATE(xk_cryst_loc(3, nks), xk_cryst_all(3, nkstot))
    xk_cryst_loc(:, 1:nks) = xk(:, 1:nks)
    CALL cryst_to_cart(nks, xk_cryst_loc, at, -1)
    xk_cryst_all = 0.0_dp
    xk_cryst_all(:, lower_bnd:upper_bnd) = xk_cryst_loc(:, 1:nks)
    CALL mp_sum(xk_cryst_all, inter_pool_comm)
    CALL cryst_to_cart(nks, xk_cryst_loc, at, 1)  ! restore to Cartesian

    ! Gather cu (combined rotation) for ALL k-points (small array)
    ALLOCATE(cu_all(nbnd_kept, nbndsub, nkstot))
    cu_all = czero
    cu_all(:, :, lower_bnd:upper_bnd) = u_mat(:, :, 1:nks)
    CALL mp_sum(cu_all, inter_pool_comm)

    ! Cache pool-local wavefunctions in real space ONCE (avoid repeated disk I/O)
    ALLOCATE(psir_k(dffts%nnr, npol, nbnd_kept))
    ALLOCATE(psir_kq(dffts%nnr, npol, nbnd_kept))
    ALLOCATE(psic_tmp(dffts%nnr))
    ALLOCATE(evc_tmp(npwx * npol, nbnd))
    ALLOCATE(V_folded_kf(dffts%nnr))
    ALLOCATE(edmatkq(nbnd_kept, nbnd_kept))
    ALLOCATE(edmatkq_loc(nbnd_kept, nbnd_kept))
    ALLOCATE(edmatkq_nl(nbnd_kept, nbnd_kept))
    ALLOCATE(cu_k(nbnd_kept, nbndsub), cu_kq(nbnd_kept, nbndsub))
    ALLOCATE(edms(nbndsub, nbndsub), eptmp(nbnd_kept, nbndsub))
    ALLOCATE(edmatw_2d(nbndsub, nbndsub, nrr_k, nrr_p))
    edmatw_2d = czero
    ! Bloch-basis storage for output (before Wannier rotation)
    ALLOCATE(edmat_bloch(nbnd_kept, nbnd_kept, nkstot, nkstot))
    ALLOCATE(edmat_bloch_loc(nbnd_kept, nbnd_kept, nkstot, nkstot))
    ALLOCATE(edmat_bloch_nl(nbnd_kept, nbnd_kept, nkstot, nkstot))
    edmat_bloch = czero
    edmat_bloch_loc = czero
    edmat_bloch_nl = czero

    ! Count nonlocal projectors for defect and pristine supercells
    nkb_d = 0
    DO nt = 1, V_d%ntyp
       DO na = 1, V_d%nat
          IF (V_d%ityp(na) == nt) nkb_d = nkb_d + nh(nt)
       ENDDO
    ENDDO
    nkb_p = 0
    DO nt = 1, V_p%ntyp
       DO na = 1, V_p%nat
          IF (V_p%ityp(na) == nt) nkb_p = nkb_p + nh(nt)
       ENDDO
    ENDDO
    ALLOCATE(vkb_d(npwx, nkb_d))
    ALLOCATE(vkb_p(npwx, nkb_p))

    ! Panel broadcast: cache local psi(k), compute V_folded on-the-fly per (ki,kf)
    BLOCK
       COMPLEX(dp), ALLOCATABLE :: psir_cache(:,:,:,:), psir_recv(:,:,:,:)
       ! Cached becp arrays: becp_d_cache(nkb_d, npol_or_1, nbnd, nks) etc.
       ! For scalar-rel: becp%k(nkb, nbnd) — store as (nkb, 1, nbnd, nks)
       ! For SOC:        becp%nc(nkb, npol, nbnd) — store as (nkb, npol, nbnd, nks)
       COMPLEX(dp), ALLOCATABLE :: becd_cache(:,:,:,:), becd_recv(:,:,:,:)
       COMPLEX(dp), ALLOCATABLE :: becp_cache(:,:,:,:), becpr_recv(:,:,:,:)
       COMPLEX(dp), ALLOCATABLE :: phase_ki(:), phase_kf(:)
       INTEGER :: ik_cache, ib_cache, ig_cache, ik_global, ipol_cache
       INTEGER :: ip, src_lower, src_nks, nkbase, nkrest, nks_max
       INTEGER :: ikf_local, ikf_global, npairs_done
       REAL(dp) :: qcryst(3), xki_cryst(3)
       TYPE(bec_type) :: becp_tmp_d, becp_tmp_p

       ! Maximum k-points per pool (for receive buffer sizing)
       nkbase = nkstot / npool
       nkrest = MOD(nkstot, npool)
       nks_max = nkbase
       IF (nkrest > 0) nks_max = nkbase + 1

       ! Step 1: Cache pool-local wavefunctions in real space + beta projections
       ! For SOC (npol=2): cache both spinor components
       ALLOCATE(psir_cache(dffts%nnr, npol, nbnd_kept, nks))
       ALLOCATE(becd_cache(nkb_d, npol, nbnd, nks))
       ALLOCATE(becp_cache(nkb_p, npol, nbnd, nks))
       becd_cache = czero
       becp_cache = czero

       DO ik_cache = 1, nks
          CALL read_collected_wfc(restart_dir(), ik_cache, evc_tmp)
          ! FFT to real space
          DO ib_cache = 1, nbnd_kept
             DO ipol_cache = 1, npol
                psic_tmp = (0.0_dp, 0.0_dp)
                DO ig_cache = 1, ngk(ik_cache)
                   psic_tmp(dffts%nl(igk_k(ig_cache, ik_cache))) = &
                        evc_tmp(ig_cache + (ipol_cache-1)*npwx, band_map(ib_cache))
                ENDDO
                CALL invfft('Wave', psic_tmp, dffts)
                psir_cache(:, ipol_cache, ib_cache, ik_cache) = psic_tmp
             ENDDO
          ENDDO
          ! Compute and cache beta projections for nonlocal part
          CALL get_betavkb(ngk(ik_cache), igk_k(1,ik_cache), xk(1,ik_cache), &
                            vkb_d, V_d%nat, V_d%ityp, V_d%tau, nkb_d)
          CALL allocate_bec_type(nkb_d, nbnd, becp_tmp_d)
          CALL calbec(ngk(ik_cache), vkb_d, evc_tmp, becp_tmp_d)
          IF (lspinorb) THEN
             becd_cache(:, :, :, ik_cache) = becp_tmp_d%nc(:, :, :)
          ELSE
             becd_cache(:, 1, :, ik_cache) = becp_tmp_d%k(:, :)
          ENDIF
          CALL deallocate_bec_type(becp_tmp_d)

          CALL get_betavkb(ngk(ik_cache), igk_k(1,ik_cache), xk(1,ik_cache), &
                            vkb_p, V_p%nat, V_p%ityp, V_p%tau, nkb_p)
          CALL allocate_bec_type(nkb_p, nbnd, becp_tmp_p)
          CALL calbec(ngk(ik_cache), vkb_p, evc_tmp, becp_tmp_p)
          IF (lspinorb) THEN
             becp_cache(:, :, :, ik_cache) = becp_tmp_p%nc(:, :, :)
          ELSE
             becp_cache(:, 1, :, ik_cache) = becp_tmp_p%k(:, :)
          ENDIF
          CALL deallocate_bec_type(becp_tmp_p)
       ENDDO

       IF (ionode) THEN
          WRITE(stdout, '(5X,A,I4,A)') &
               'Cached ', nks, ' pool-local k-points (psir + becp)'
          FLUSH(stdout)
       ENDIF

       ! Step 2: Panel broadcast — iterate over source pools in lockstep
       ! All pools call MPI_Bcast with the SAME root, preventing deadlock.
       ! V_folded is computed ON-THE-FLY for each (ki,kf) pair using the
       ! EXACT q = kf - ki (no BZ wrapping) to avoid exp(iG*r) phase errors.
       ALLOCATE(psir_recv(dffts%nnr, npol, nbnd_kept, nks_max))
       ALLOCATE(becd_recv(nkb_d, npol, nbnd, nks_max))
       ALLOCATE(becpr_recv(nkb_p, npol, nbnd, nks_max))
       ALLOCATE(phase_ki(nrr_k), phase_kf(nrr_k))
       edmatw_2d = czero
       npairs_done = 0

       DO ip = 0, npool - 1
          ! Determine k-point range of source pool ip
          IF (ip < nkrest) THEN
             src_lower = ip * (nkbase + 1) + 1
             src_nks = nkbase + 1
          ELSE
             src_lower = nkrest * (nkbase + 1) + (ip - nkrest) * nkbase + 1
             src_nks = nkbase
          ENDIF

          ! Source pool broadcasts its cached wavefunctions + becp to all pools
          IF (my_pool_id == ip) THEN
             psir_recv(:, :, :, 1:src_nks) = psir_cache(:, :, :, 1:nks)
             becd_recv(:, :, :, 1:src_nks) = becd_cache(:, :, :, 1:nks)
             becpr_recv(:, :, :, 1:src_nks) = becp_cache(:, :, :, 1:nks)
          ENDIF
          CALL MPI_Bcast(psir_recv, dffts%nnr * npol * nbnd_kept * src_nks * 2, &
               MPI_DOUBLE_PRECISION, ip, inter_pool_comm, ierr_mpi)
          CALL MPI_Bcast(becd_recv, nkb_d * npol * nbnd * src_nks * 2, &
               MPI_DOUBLE_PRECISION, ip, inter_pool_comm, ierr_mpi)
          CALL MPI_Bcast(becpr_recv, nkb_p * npol * nbnd * src_nks * 2, &
               MPI_DOUBLE_PRECISION, ip, inter_pool_comm, ierr_mpi)

          ! Process all (ki_local, kf_from_src) pairs
          DO ikf_local = 1, src_nks
             ikf_global = src_lower + ikf_local - 1

             DO ik = 1, nks
                ik_global = ik + lower_bnd - 1

                ! Compute V_folded on-the-fly with EXACT q = kf - ki (no wrapping!)
                ! Using wrapped q would give V_folded(q+G) = exp(iG*r) * V_folded(q),
                ! introducing an unphysical phase that corrupts the double-FT.
                qcryst = xk_cryst_all(:, ikf_global) - xk_cryst_all(:, ik_global)
                V_folded_kf = (0.0_dp, 0.0_dp)
                IF (ABS(qcryst(1)) < 1.0d-8 .AND. ABS(qcryst(2)) < 1.0d-8 &
                     .AND. ABS(qcryst(3)) < 1.0d-8) THEN
                   inr = 0
                   DO irz = 0, V_d%nr3 - 1
                      ir3mod = irz - (irz / dffts%nr3) * dffts%nr3
                      DO iry = 0, V_d%nr2 - 1
                         ir2mod = iry - (iry / dffts%nr2) * dffts%nr2
                         DO irx = 0, V_d%nr1 - 1
                            ir1mod = irx - (irx / dffts%nr1) * dffts%nr1
                            inr = inr + 1
                            irnmod = ir3mod * dffts%nr1 * dffts%nr2 + &
                                     ir2mod * dffts%nr1 + ir1mod + 1
                            V_folded_kf(irnmod) = V_folded_kf(irnmod) + V_colin(inr)
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE
                   ! Phase per grid step: q·r = 2π × (q1·irx/nr1 + q2·iry/nr2 + q3·irz/nr3)
                   ! q is in crystal coords, grid indices are in crystal-lattice units
                   ! DO NOT mix with at (Cartesian lattice vectors) — that gives wrong
                   ! phases for non-orthogonal lattices (hexagonal, etc.)
                   d1 = tpi * qcryst(1) / dffts%nr1
                   d2 = tpi * qcryst(2) / dffts%nr2
                   d3 = tpi * qcryst(3) / dffts%nr3
                   inr = 0
                   DO irz = 0, V_d%nr3 - 1
                      ir3mod = irz - (irz / dffts%nr3) * dffts%nr3
                      DO iry = 0, V_d%nr2 - 1
                         ir2mod = iry - (iry / dffts%nr2) * dffts%nr2
                         DO irx = 0, V_d%nr1 - 1
                            ir1mod = irx - (irx / dffts%nr1) * dffts%nr1
                            inr = inr + 1
                            irnmod = ir3mod * dffts%nr1 * dffts%nr2 + &
                                     ir2mod * dffts%nr1 + ir1mod + 1
                            arg = irx * d1 + iry * d2 + irz * d3
                            phase = CMPLX(COS(arg), SIN(arg), dp)
                            V_folded_kf(irnmod) = V_folded_kf(irnmod) + V_colin(inr) * phase
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

                ! Matrix element: M_B = sum_{r,σ} u*_{σ}(ki,r) u_{σ}(kf,r) V_q(r) / nnr
                psir_k(:,:,:) = psir_cache(:,:,:, ik)
                psir_kq(:,:,:) = psir_recv(:,:,:, ikf_local)

                DO ib_i = 1, nbnd_kept
                   DO ib_j = 1, nbnd_kept
                      ! LOCAL part
                      mlocal = czero
                      DO ipol_cache = 1, npol
                         DO ir = 1, dffts%nnr
                            mlocal = mlocal + CONJG(psir_k(ir, ipol_cache, ib_i)) * &
                                 psir_kq(ir, ipol_cache, ib_j) * V_folded_kf(ir)
                         ENDDO
                      ENDDO
                      mlocal = mlocal / DBLE(dffts%nnr)

                      ! NONLOCAL defect: becp_ki from cache, becp_kf from recv
                      mnl_d = czero
                      ijkb0 = 0
                      DO nt = 1, V_d%ntyp
                         DO na = 1, V_d%nat
                            IF (V_d%ityp(na) == nt) THEN
                               IF (lspinorb) THEN
                                  DO ih = 1, nh(nt)
                                     ikb = ijkb0 + ih
                                     mnl_d = mnl_d + &
                                          CONJG(becd_cache(ikb,1,band_map(ib_i),ik)) * becd_recv(ikb,1,band_map(ib_j),ikf_local) * dvan_so(ih,ih,1,nt) + &
                                          CONJG(becd_cache(ikb,1,band_map(ib_i),ik)) * becd_recv(ikb,2,band_map(ib_j),ikf_local) * dvan_so(ih,ih,2,nt) + &
                                          CONJG(becd_cache(ikb,2,band_map(ib_i),ik)) * becd_recv(ikb,1,band_map(ib_j),ikf_local) * dvan_so(ih,ih,3,nt) + &
                                          CONJG(becd_cache(ikb,2,band_map(ib_i),ik)) * becd_recv(ikb,2,band_map(ib_j),ikf_local) * dvan_so(ih,ih,4,nt)
                                     DO jh = ih + 1, nh(nt)
                                        jkb = ijkb0 + jh
                                        mnl_d = mnl_d + &
                                             CONJG(becd_cache(ikb,1,band_map(ib_i),ik)) * becd_recv(jkb,1,band_map(ib_j),ikf_local) * dvan_so(ih,jh,1,nt) + &
                                             CONJG(becd_cache(jkb,1,band_map(ib_i),ik)) * becd_recv(ikb,1,band_map(ib_j),ikf_local) * dvan_so(jh,ih,1,nt) + &
                                             CONJG(becd_cache(ikb,1,band_map(ib_i),ik)) * becd_recv(jkb,2,band_map(ib_j),ikf_local) * dvan_so(ih,jh,2,nt) + &
                                             CONJG(becd_cache(jkb,1,band_map(ib_i),ik)) * becd_recv(ikb,2,band_map(ib_j),ikf_local) * dvan_so(jh,ih,2,nt) + &
                                             CONJG(becd_cache(ikb,2,band_map(ib_i),ik)) * becd_recv(jkb,1,band_map(ib_j),ikf_local) * dvan_so(ih,jh,3,nt) + &
                                             CONJG(becd_cache(jkb,2,band_map(ib_i),ik)) * becd_recv(ikb,1,band_map(ib_j),ikf_local) * dvan_so(jh,ih,3,nt) + &
                                             CONJG(becd_cache(ikb,2,band_map(ib_i),ik)) * becd_recv(jkb,2,band_map(ib_j),ikf_local) * dvan_so(ih,jh,4,nt) + &
                                             CONJG(becd_cache(jkb,2,band_map(ib_i),ik)) * becd_recv(ikb,2,band_map(ib_j),ikf_local) * dvan_so(jh,ih,4,nt)
                                     ENDDO
                                  ENDDO
                               ELSE
                                  DO ih = 1, nh(nt)
                                     ikb = ijkb0 + ih
                                     mnl_d = mnl_d + CONJG(becd_cache(ikb,1,band_map(ib_i),ik)) * &
                                          becd_recv(ikb,1,band_map(ib_j),ikf_local) * dvan(ih, ih, nt)
                                     DO jh = ih + 1, nh(nt)
                                        jkb = ijkb0 + jh
                                        mnl_d = mnl_d + &
                                             (CONJG(becd_cache(ikb,1,band_map(ib_i),ik)) * becd_recv(jkb,1,band_map(ib_j),ikf_local) + &
                                              CONJG(becd_cache(jkb,1,band_map(ib_i),ik)) * becd_recv(ikb,1,band_map(ib_j),ikf_local)) * &
                                             dvan(ih, jh, nt)
                                     ENDDO
                                  ENDDO
                               ENDIF
                               ijkb0 = ijkb0 + nh(nt)
                            ENDIF
                         ENDDO
                      ENDDO

                      ! NONLOCAL pristine
                      mnl_p = czero
                      ijkb0 = 0
                      DO nt = 1, V_p%ntyp
                         DO na = 1, V_p%nat
                            IF (V_p%ityp(na) == nt) THEN
                               IF (lspinorb) THEN
                                  DO ih = 1, nh(nt)
                                     ikb = ijkb0 + ih
                                     mnl_p = mnl_p + &
                                          CONJG(becp_cache(ikb,1,band_map(ib_i),ik)) * becpr_recv(ikb,1,band_map(ib_j),ikf_local) * dvan_so(ih,ih,1,nt) + &
                                          CONJG(becp_cache(ikb,1,band_map(ib_i),ik)) * becpr_recv(ikb,2,band_map(ib_j),ikf_local) * dvan_so(ih,ih,2,nt) + &
                                          CONJG(becp_cache(ikb,2,band_map(ib_i),ik)) * becpr_recv(ikb,1,band_map(ib_j),ikf_local) * dvan_so(ih,ih,3,nt) + &
                                          CONJG(becp_cache(ikb,2,band_map(ib_i),ik)) * becpr_recv(ikb,2,band_map(ib_j),ikf_local) * dvan_so(ih,ih,4,nt)
                                     DO jh = ih + 1, nh(nt)
                                        jkb = ijkb0 + jh
                                        mnl_p = mnl_p + &
                                             CONJG(becp_cache(ikb,1,band_map(ib_i),ik)) * becpr_recv(jkb,1,band_map(ib_j),ikf_local) * dvan_so(ih,jh,1,nt) + &
                                             CONJG(becp_cache(jkb,1,band_map(ib_i),ik)) * becpr_recv(ikb,1,band_map(ib_j),ikf_local) * dvan_so(jh,ih,1,nt) + &
                                             CONJG(becp_cache(ikb,1,band_map(ib_i),ik)) * becpr_recv(jkb,2,band_map(ib_j),ikf_local) * dvan_so(ih,jh,2,nt) + &
                                             CONJG(becp_cache(jkb,1,band_map(ib_i),ik)) * becpr_recv(ikb,2,band_map(ib_j),ikf_local) * dvan_so(jh,ih,2,nt) + &
                                             CONJG(becp_cache(ikb,2,band_map(ib_i),ik)) * becpr_recv(jkb,1,band_map(ib_j),ikf_local) * dvan_so(ih,jh,3,nt) + &
                                             CONJG(becp_cache(jkb,2,band_map(ib_i),ik)) * becpr_recv(ikb,1,band_map(ib_j),ikf_local) * dvan_so(jh,ih,3,nt) + &
                                             CONJG(becp_cache(ikb,2,band_map(ib_i),ik)) * becpr_recv(jkb,2,band_map(ib_j),ikf_local) * dvan_so(ih,jh,4,nt) + &
                                             CONJG(becp_cache(jkb,2,band_map(ib_i),ik)) * becpr_recv(ikb,2,band_map(ib_j),ikf_local) * dvan_so(jh,ih,4,nt)
                                     ENDDO
                                  ENDDO
                               ELSE
                                  DO ih = 1, nh(nt)
                                     ikb = ijkb0 + ih
                                     mnl_p = mnl_p + CONJG(becp_cache(ikb,1,band_map(ib_i),ik)) * &
                                          becpr_recv(ikb,1,band_map(ib_j),ikf_local) * dvan(ih, ih, nt)
                                     DO jh = ih + 1, nh(nt)
                                        jkb = ijkb0 + jh
                                        mnl_p = mnl_p + &
                                             (CONJG(becp_cache(ikb,1,band_map(ib_i),ik)) * becpr_recv(jkb,1,band_map(ib_j),ikf_local) + &
                                              CONJG(becp_cache(jkb,1,band_map(ib_i),ik)) * becpr_recv(ikb,1,band_map(ib_j),ikf_local)) * &
                                             dvan(ih, jh, nt)
                                     ENDDO
                                  ENDDO
                               ENDIF
                               ijkb0 = ijkb0 + nh(nt)
                            ENDIF
                         ENDDO
                      ENDDO

                      edmatkq(ib_i, ib_j) = mlocal + mnl_d - mnl_p
                      edmatkq_loc(ib_i, ib_j) = mlocal
                      edmatkq_nl(ib_i, ib_j) = mnl_d - mnl_p
                   ENDDO
                ENDDO

                ! Store Bloch-basis M for output (before Wannier rotation)
                edmat_bloch(:, :, ik_global, ikf_global) = edmatkq(:, :)
                edmat_bloch_loc(:, :, ik_global, ikf_global) = edmatkq_loc(:, :)
                edmat_bloch_nl(:, :, ik_global, ikf_global) = edmatkq_nl(:, :)

                ! Wannier rotation: M_W = U_dag(ki) M_B U(kf)
                cu_k  = cu_all(:, :, ik_global)
                cu_kq = cu_all(:, :, ikf_global)
                CALL ZGEMM('N', 'N', nbnd_kept, nbndsub, nbnd_kept, &
                     cone, edmatkq, nbnd_kept, cu_kq, nbnd_kept, czero, eptmp, nbnd_kept)
                CALL ZGEMM('C', 'N', nbndsub, nbndsub, nbnd_kept, &
                     cone, cu_k, nbnd_kept, eptmp, nbnd_kept, czero, edms, nbndsub)

                ! Double FT: M(Re,Rp) += exp(+iki*Re) exp(-ikf*Rp) / Nk^2 * M_W
                xki_cryst = xk_cryst_all(:, ik_global)
                xkf_cryst = xk_cryst_all(:, ikf_global)
                DO ir = 1, nrr_k
                   rdotk  = tpi * DOT_PRODUCT(xki_cryst, DBLE(irvec_k(:, ir)))
                   phase_ki(ir) = EXP(CMPLX(0.0_dp, +rdotk, dp))
                   rdotkf = tpi * DOT_PRODUCT(xkf_cryst, DBLE(irvec_k(:, ir)))
                   phase_kf(ir) = EXP(CMPLX(0.0_dp, -rdotkf, dp))
                ENDDO
                DO ir_p = 1, nrr_k
                   DO ir_e = 1, nrr_k
                      edmatw_2d(:,:,ir_e,ir_p) = edmatw_2d(:,:,ir_e,ir_p) + &
                           phase_ki(ir_e) * phase_kf(ir_p) / DBLE(nkstot * nkstot) * edms(:,:)
                   ENDDO
                ENDDO

                npairs_done = npairs_done + 1
             ENDDO  ! ik (local ki)
          ENDDO  ! ikf_local (kf from source pool)

          IF (ionode) THEN
             WRITE(stdout, '(5X,A,I4,A,I4,A,I6,A,I6)') &
                  'kf pool ', ip, ' / ', npool - 1, &
                  '  pairs: ', npairs_done, ' / ', nks * nkstot
             FLUSH(stdout)
          ENDIF
       ENDDO  ! ip (source pool)

       DEALLOCATE(phase_ki, phase_kf, psir_recv, psir_cache)
       DEALLOCATE(becd_cache, becd_recv, becp_cache, becpr_recv)

       ! Gather FT contributions from all pools
       CALL mp_sum(edmatw_2d, inter_pool_comm)
    END BLOCK

    ! Gather Bloch-basis M from all pools and write to file
    CALL mp_sum(edmat_bloch, inter_pool_comm)
    CALL mp_sum(edmat_bloch_loc, inter_pool_comm)
    CALL mp_sum(edmat_bloch_nl, inter_pool_comm)

    IF (ionode) THEN
       BLOCK
          INTEGER :: iunit_b, iki_b, ikf_b, ib1, ib2
          CHARACTER(LEN=256) :: fname_b

          fname_b = TRIM(prefix_in) // '_edmat_bloch.dat'
          iunit_b = 89
          OPEN(iunit_b, FILE=TRIM(fname_b), FORM='formatted')
          WRITE(iunit_b, '(A)') '# Bloch-basis M(n,m,ki,kf) from ed_coarse_full_q (before Wannier rotation)'
          WRITE(iunit_b, '(A,I6,A,I4)') '# nkstot=', nkstot, '  nbnd_kept=', nbnd_kept
          WRITE(iunit_b, '(A)') '# iki ikf  kix kiy kiz  kfx kfy kfz  ibnd jbnd  |M|^2  Re(M)  Im(M)  |M_loc|^2  |M_nl|^2'

          DO iki_b = 1, nkstot
             DO ikf_b = 1, nkstot
                DO ib1 = 1, nbnd_kept
                   DO ib2 = 1, nbnd_kept
                      IF (ABS(edmat_bloch(ib1, ib2, iki_b, ikf_b)) < 1.0d-30) CYCLE
                      WRITE(iunit_b, '(2I6,6F10.5,2I4,5ES16.8)') &
                           iki_b, ikf_b, &
                           xk_cryst_all(:, iki_b), xk_cryst_all(:, ikf_b), &
                           ib1, ib2, &
                           ABS(edmat_bloch(ib1, ib2, iki_b, ikf_b))**2, &
                           REAL(edmat_bloch(ib1, ib2, iki_b, ikf_b)), &
                           AIMAG(edmat_bloch(ib1, ib2, iki_b, ikf_b)), &
                           ABS(edmat_bloch_loc(ib1, ib2, iki_b, ikf_b))**2, &
                           ABS(edmat_bloch_nl(ib1, ib2, iki_b, ikf_b))**2
                   ENDDO
                ENDDO
             ENDDO
          ENDDO

          CLOSE(iunit_b)
          WRITE(stdout, '(5X,A,A)') 'Bloch-basis M written to ', TRIM(fname_b)
       END BLOCK
    ENDIF

    DEALLOCATE(edmat_bloch, edmat_bloch_loc, edmat_bloch_nl)
    DEALLOCATE(edmatkq_loc, edmatkq_nl)



    DEALLOCATE(psir_k, psir_kq, psic_tmp, evc_tmp, V_folded_kf, vkb_d, vkb_p)
    DEALLOCATE(edmatkq, cu_k, cu_kq, edms, eptmp)
    DEALLOCATE(xk_cryst_loc, xk_cryst_all, cu_all, band_map)

    IF (ionode) THEN
       WRITE(stdout, '(5X,A)') 'All ki-kf pairs completed. Writing output...'
       CALL write_edmatw_2d_v2(nbndsub, nrr_k, nrr_p, irvec_k, irvec_k, &
                                ndegen_k, ndegen_k, edmatw_2d, prefix_in)
       ! Write full M(R,Rp) binary for edwread
       CALL write_edmatw_2d_bin(nbndsub, nrr_k, irvec_k, ndegen_k, &
                                 edmatw_2d, prefix_in)
       WRITE(stdout, '(5X,A)') REPEAT('=', 60)
    ENDIF

  CONTAINS

    FUNCTION find_q_index(qvec, nk, xk_c) RESULT(iq_out)
      ! Find k-grid index closest to qvec (modulo reciprocal lattice)
      INTEGER, INTENT(IN) :: nk
      REAL(dp), INTENT(IN) :: qvec(3), xk_c(3, nk)
      INTEGER :: iq_out, jk
      REAL(dp) :: diff(3), dist, min_dist
      min_dist = 1.0d10
      iq_out = 1
      DO jk = 1, nk
         diff = qvec - xk_c(:, jk)
         diff = diff - NINT(diff)
         dist = SUM(diff**2)
         IF (dist < min_dist) THEN
            min_dist = dist
            iq_out = jk
         ENDIF
      ENDDO
    END FUNCTION

    FUNCTION find_kf_cryst(iki, dkf, nk, xk_c) RESULT(ikf_out)
      INTEGER, INTENT(IN) :: iki, nk
      REAL(dp), INTENT(IN) :: dkf(3), xk_c(3, nk)
      INTEGER :: ikf_out, jk
      REAL(dp) :: xkf_target(3), diff(3), dist, min_dist
      xkf_target = xk_c(:, iki) + dkf
      xkf_target = xkf_target - FLOOR(xkf_target)
      min_dist = 1.0d10
      ikf_out = 1
      DO jk = 1, nk
         diff = xkf_target - xk_c(:, jk)
         diff = diff - NINT(diff)
         dist = SUM(diff**2)
         IF (dist < min_dist) THEN
            min_dist = dist
            ikf_out = jk
         ENDIF
      ENDDO
    END FUNCTION

    FUNCTION pool_of_k(ik_global, nktot, np) RESULT(ipool)
      INTEGER, INTENT(IN) :: ik_global, nktot, np
      INTEGER :: ipool, nk_per_pool, remainder
      nk_per_pool = nktot / np
      remainder = MOD(nktot, np)
      IF (ik_global <= (nk_per_pool + 1) * remainder) THEN
         ipool = (ik_global - 1) / (nk_per_pool + 1)
      ELSE
         ipool = remainder + (ik_global - (nk_per_pool + 1) * remainder - 1) / nk_per_pool
      ENDIF
    END FUNCTION

    FUNCTION local_k_index(ik_global, nktot, np, ipool) RESULT(ik_local)
      INTEGER, INTENT(IN) :: ik_global, nktot, np, ipool
      INTEGER :: ik_local, nk_per_pool, remainder, offset
      nk_per_pool = nktot / np
      remainder = MOD(nktot, np)
      IF (ipool < remainder) THEN
         offset = ipool * (nk_per_pool + 1)
      ELSE
         offset = remainder * (nk_per_pool + 1) + (ipool - remainder) * nk_per_pool
      ENDIF
      ik_local = ik_global - offset
    END FUNCTION

    SUBROUTINE fetch_wfc_from_pool(src_pool, ik_src, evc_out, npwx_npol, nb, &
                                    pool_comm, my_pid, np)
      INTEGER, INTENT(IN) :: src_pool, ik_src, npwx_npol, nb, pool_comm, my_pid, np
      COMPLEX(dp), INTENT(OUT) :: evc_out(npwx_npol, nb)
      INTEGER :: tag, ierr, status(MPI_STATUS_SIZE)

      tag = ik_src + 10000
      IF (my_pid == src_pool) THEN
         ! I own this k-point — read and send to all who need it
         CALL read_collected_wfc(restart_dir(), ik_src, evc_out)
      ENDIF
      ! Broadcast from src_pool to all pools
      CALL MPI_Bcast(evc_out, npwx_npol * nb * 2, MPI_DOUBLE_PRECISION, &
                      src_pool, pool_comm, ierr)
    END SUBROUTINE

  END SUBROUTINE ed_coarse_full_q


  SUBROUTINE write_edmatw_2d(nbndsub, nrr, irvec, ndegen, edmatw_2d, prefix_in)
    USE io_global, ONLY : ionode, stdout
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nbndsub, nrr
    INTEGER, INTENT(IN) :: irvec(3, nrr), ndegen(nrr)
    COMPLEX(dp), INTENT(IN) :: edmatw_2d(nbndsub, nbndsub, nrr, nrr)
    CHARACTER(LEN=*), INTENT(IN) :: prefix_in

    INTEGER :: ir_e, ir_p, ii, jj, iunit
    REAL(dp) :: max_abs
    CHARACTER(LEN=256) :: fname

    iunit = 81
    fname = TRIM(prefix_in) // '_edmatw_2d.dat'
    OPEN(iunit, FILE=TRIM(fname), FORM='formatted')
    WRITE(iunit, '(A)') '# M(R,Rp) from full double-FT'
    WRITE(iunit, '(I6, I6)') nbndsub, nrr

    DO ir_p = 1, nrr
       DO ir_e = 1, nrr
          max_abs = MAXVAL(ABS(edmatw_2d(:,:,ir_e,ir_p)))
          WRITE(iunit, '(6I5, ES16.8)') irvec(:,ir_e), irvec(:,ir_p), max_abs
       ENDDO
    ENDDO
    CLOSE(iunit)
    WRITE(stdout, '(5X,A,A)') 'M(R,Rp) written to ', TRIM(fname)

    ! Also write decay comparison
    fname = TRIM(prefix_in) // '_decay_2d.dat'
    OPEN(iunit, FILE=TRIM(fname), FORM='formatted')
    WRITE(iunit, '(A)') '# |R|  |Rp|  max|M(R,Rp)|  (for comparison with diagonal)'
    DO ir_p = 1, nrr
       DO ir_e = 1, nrr
          max_abs = MAXVAL(ABS(edmatw_2d(:,:,ir_e,ir_p)))
          IF (max_abs > 1.0d-10) THEN
             WRITE(iunit, '(2F14.6, ES16.8)') &
                  SQRT(DBLE(SUM(irvec(:,ir_e)**2))), &
                  SQRT(DBLE(SUM(irvec(:,ir_p)**2))), max_abs
          ENDIF
       ENDDO
    ENDDO
    CLOSE(iunit)
    WRITE(stdout, '(5X,A,A)') 'Decay data written to ', TRIM(fname)

  END SUBROUTINE write_edmatw_2d

  SUBROUTINE write_edmatw_2d_v2(nbndsub, nrr_e, nrr_p, irvec_e, irvec_p, &
                                 ndegen_e, ndegen_p, edmatw_2d, prefix_in)
    USE io_global, ONLY : ionode, stdout
    USE cell_base, ONLY : at, alat
    USE ep_constants, ONLY : bohr2ang
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nbndsub, nrr_e, nrr_p
    INTEGER, INTENT(IN) :: irvec_e(3, nrr_e), irvec_p(3, nrr_p)
    INTEGER, INTENT(IN) :: ndegen_e(nrr_e), ndegen_p(nrr_p)
    COMPLEX(dp), INTENT(IN) :: edmatw_2d(nbndsub, nbndsub, nrr_e, nrr_p)
    CHARACTER(LEN=*), INTENT(IN) :: prefix_in

    INTEGER :: ir_e, ir_p, iunit
    REAL(dp) :: max_abs, re_len, rp_len, rvec_e(3), rvec_p(3)
    CHARACTER(LEN=256) :: fname

    iunit = 81

    ! Write decay data (|R| in Angstrom, |Rp| in Angstrom, max|M|)
    fname = TRIM(prefix_in) // '_decay_2d.dat'
    OPEN(iunit, FILE=TRIM(fname), FORM='formatted')
    WRITE(iunit, '(A)') '# |R|(Ang)  |Rp|(Ang)  max|M(R,Rp)|'
    DO ir_p = 1, nrr_p
       rvec_p = MATMUL(at, DBLE(irvec_p(:, ir_p))) * alat * bohr2ang
       rp_len = SQRT(DOT_PRODUCT(rvec_p, rvec_p))
       DO ir_e = 1, nrr_e
          max_abs = MAXVAL(ABS(edmatw_2d(:,:,ir_e,ir_p)))
          IF (max_abs > 1.0d-10) THEN
             rvec_e = MATMUL(at, DBLE(irvec_e(:, ir_e))) * alat * bohr2ang
             re_len = SQRT(DOT_PRODUCT(rvec_e, rvec_e))
             WRITE(iunit, '(2F14.6, ES16.8)') re_len, rp_len, max_abs
          ENDIF
       ENDDO
    ENDDO
    CLOSE(iunit)
    WRITE(stdout, '(5X,A,A)') 'Decay data written to ', TRIM(fname)
    WRITE(stdout, '(5X,A,I6,A,I6)') 'nrr_e=', nrr_e, '  nrr_p=', nrr_p
  END SUBROUTINE write_edmatw_2d_v2


  SUBROUTINE write_edmatw_2d_bin(nbndsub, nrr, irvec, ndegen, edmatw_2d, prefix_in)
    !-----------------------------------------------------------------------
    ! Write full M(R,Rp) in binary format for edwread restart.
    !-----------------------------------------------------------------------
    USE io_global, ONLY : stdout
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nbndsub, nrr
    INTEGER, INTENT(IN) :: irvec(3, nrr), ndegen(nrr)
    COMPLEX(dp), INTENT(IN) :: edmatw_2d(nbndsub, nbndsub, nrr, nrr)
    CHARACTER(LEN=*), INTENT(IN) :: prefix_in

    INTEGER :: iunit, ir
    CHARACTER(LEN=256) :: fname

    iunit = 82
    fname = TRIM(prefix_in) // '_edmatw_2d.bin'
    OPEN(iunit, FILE=TRIM(fname), FORM='unformatted')
    WRITE(iunit) nbndsub, nrr
    WRITE(iunit) ndegen(1:nrr)
    DO ir = 1, nrr
       WRITE(iunit) irvec(:, ir)
    ENDDO
    WRITE(iunit) edmatw_2d
    CLOSE(iunit)
    WRITE(stdout, '(5X,A,A)') 'M(R,Rp) binary written to ', TRIM(fname)
  END SUBROUTINE write_edmatw_2d_bin


  SUBROUTINE read_edmatw_2d_file(prefix_in, nbndsub, nrr, ndegen, irvec, edmatw_2d)
    !-----------------------------------------------------------------------
    ! Read full M(R,Rp) from binary file.
    ! Ionode reads, then broadcasts to all ranks (avoids filesystem races).
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode, ionode_id, stdout
    USE mp, ONLY : mp_bcast
    USE mp_world, ONLY : world_comm
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: prefix_in
    INTEGER, INTENT(OUT) :: nbndsub, nrr
    INTEGER, ALLOCATABLE, INTENT(OUT) :: ndegen(:), irvec(:,:)
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: edmatw_2d(:,:,:,:)

    INTEGER :: iunit, ir
    CHARACTER(LEN=256) :: fname

    nbndsub = 0
    nrr = 0

    IF (ionode) THEN
       iunit = 82
       fname = TRIM(prefix_in) // '_edmatw_2d.bin'
       OPEN(iunit, FILE=TRIM(fname), FORM='unformatted', STATUS='old')
       READ(iunit) nbndsub, nrr
    ENDIF
    CALL mp_bcast(nbndsub, ionode_id, world_comm)
    CALL mp_bcast(nrr, ionode_id, world_comm)

    ALLOCATE(ndegen(nrr), irvec(3, nrr))

    IF (ionode) THEN
       READ(iunit) ndegen(1:nrr)
       DO ir = 1, nrr
          READ(iunit) irvec(:, ir)
       ENDDO
    ENDIF
    CALL mp_bcast(ndegen, ionode_id, world_comm)
    CALL mp_bcast(irvec, ionode_id, world_comm)

    ALLOCATE(edmatw_2d(nbndsub, nbndsub, nrr, nrr))

    IF (ionode) THEN
       READ(iunit) edmatw_2d
       CLOSE(iunit)
       WRITE(stdout, '(5X,A,A)') 'M(R,Rp) binary read from ', TRIM(fname)
       WRITE(stdout, '(5X,A,I4,A,I6)') '  nbndsub=', nbndsub, '  nrr=', nrr
    ENDIF
    CALL mp_bcast(edmatw_2d, ionode_id, world_comm)
  END SUBROUTINE read_edmatw_2d_file


  SUBROUTINE ed_fine_interp_2d(nbndsub, nrr, irvec, ndegen, chw, edmatw_2d, &
                                nk1f, nk2f, nk3f, prefix_in)
    !-----------------------------------------------------------------------
    ! Part B with full double-FT M(R,Rp).
    !
    ! For each (k_i, k_f) on the fine grid:
    !   M_W(k_i,k_f) = Σ_{R,R'} exp(ik_i·R) exp(-ik_f·Rp) M(R,Rp)
    !   Diagonalize H(k_i) → U(k_i), H(k_f) → U(k_f)
    !   M_B = U†(k_i) · M_W · U(k_f)
    !
    ! For diagonal (k_i = k_f): outputs prefix_edmatf_2d.dat
    ! For off-diagonal from K: outputs prefix_edmat_scatter_2d.dat
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode, stdout
    USE wan2bloch_edi, ONLY : get_cfac, hamwan2bloch_with_evec, edmatwan2bloch_2d
    USE ep_constants, ONLY : czero, cone
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nbndsub, nrr
    INTEGER, INTENT(IN) :: irvec(3, nrr), ndegen(nrr)
    COMPLEX(dp), INTENT(IN) :: chw(nbndsub, nbndsub, nrr)
    COMPLEX(dp), INTENT(IN) :: edmatw_2d(nbndsub, nbndsub, nrr, nrr)
    INTEGER, INTENT(IN) :: nk1f, nk2f, nk3f
    CHARACTER(LEN=*), INTENT(IN) :: prefix_in

    INTEGER :: ik, i1, i2, i3, nktotf, ibnd, jbnd, iunit_diag, iunit_scatter
    REAL(dp) :: xk_i(3), xk_f(3), eig_i(nbndsub), eig_f(nbndsub)
    COMPLEX(dp) :: evec_i(nbndsub, nbndsub), evec_f(nbndsub, nbndsub)
    COMPLEX(dp) :: cfac_ki(nrr), cfac_kf(nrr)
    COMPLEX(dp) :: edmatf_w(nbndsub, nbndsub), edmatf_b(nbndsub, nbndsub)
    COMPLEX(dp) :: tmp(nbndsub, nbndsub)
    CHARACTER(LEN=256) :: fname
    REAL(dp), PARAMETER :: twopi = 6.283185307179586_dp
    COMPLEX(dp), PARAMETER :: ci = (0.0_dp, 1.0_dp)

    ! K point for off-diagonal scattering
    REAL(dp) :: xk_K(3)
    COMPLEX(dp) :: cfac_K(nrr), evec_K(nbndsub, nbndsub)
    REAL(dp) :: eig_K(nbndsub)

    nktotf = nk1f * nk2f * nk3f
    xk_K = (/1.0_dp/3.0_dp, 1.0_dp/3.0_dp, 0.0_dp/)

    IF (ionode) THEN
       WRITE(stdout, '(/,5X,A)') REPEAT('=', 60)
       WRITE(stdout, '(5X,A)')   'Part B (2D): Interpolation using M(R,Rp)'
       WRITE(stdout, '(5X,A)')   REPEAT('=', 60)
       WRITE(stdout, '(5X,A,3I4,A,I8)') 'Fine k-grid: ', nk1f, nk2f, nk3f, &
            '  total: ', nktotf
       FLUSH(stdout)
    ENDIF

    ! Diagonalize at K point once
    CALL get_cfac(nrr, irvec, xk_K, cfac_K)
    CALL hamwan2bloch_with_evec(nbndsub, nrr, ndegen, eig_K, evec_K, chw, cfac_K)

    iunit_diag = 80
    iunit_scatter = 79
    IF (ionode) THEN
       fname = TRIM(prefix_in) // '_edmatf_2d.dat'
       OPEN(iunit_diag, FILE=TRIM(fname), FORM='formatted')
       WRITE(iunit_diag, '(A)') '# Interpolated M from full double-FT M(R,Rp) — diagonal k_i=k_f'
       WRITE(iunit_diag, '(A)') '# ik  kx  ky  kz  ibnd  jbnd  |M|^2  Re(M)  Im(M)'

       fname = TRIM(prefix_in) // '_edmat_scatter_2d.dat'
       OPEN(iunit_scatter, FILE=TRIM(fname), FORM='formatted')
       WRITE(iunit_scatter, '(A)') '# Interpolated M from full double-FT — K -> k_f'
       WRITE(iunit_scatter, '(A,3F10.5)') '# k_initial (crystal) = ', xk_K
       WRITE(iunit_scatter, '(A)') '# ik  kfx  kfy  kfz  ibnd  jbnd  |M|^2  Re(M)  Im(M)'
    ENDIF

    ik = 0
    DO i3 = 0, nk3f - 1
       DO i2 = 0, nk2f - 1
          DO i1 = 0, nk1f - 1
             ik = ik + 1
             xk_f(1) = DBLE(i1) / DBLE(nk1f)
             xk_f(2) = DBLE(i2) / DBLE(nk2f)
             xk_f(3) = DBLE(i3) / DBLE(nk3f)

             CALL get_cfac(nrr, irvec, xk_f, cfac_kf)
             CALL hamwan2bloch_with_evec(nbndsub, nrr, ndegen, eig_f, evec_f, chw, cfac_kf)

             ! === Diagonal: k_i = k_f ===
             cfac_ki = cfac_kf
             CALL edmatwan2bloch_2d(nbndsub, nrr, ndegen, edmatw_2d, cfac_ki, cfac_kf, edmatf_w)
             CALL ZGEMM('N', 'N', nbndsub, nbndsub, nbndsub, &
                         cone, edmatf_w, nbndsub, evec_f, nbndsub, czero, tmp, nbndsub)
             CALL ZGEMM('C', 'N', nbndsub, nbndsub, nbndsub, &
                         cone, evec_f, nbndsub, tmp, nbndsub, czero, edmatf_b, nbndsub)

             IF (ionode) THEN
                DO ibnd = 1, nbndsub
                   DO jbnd = 1, nbndsub
                      WRITE(iunit_diag, '(I8,3F10.5,2I4,3ES16.8)') ik, xk_f, &
                           ibnd, jbnd, ABS(edmatf_b(ibnd,jbnd))**2, &
                           REAL(edmatf_b(ibnd,jbnd)), AIMAG(edmatf_b(ibnd,jbnd))
                   ENDDO
                ENDDO
             ENDIF

             ! === Off-diagonal: k_i = K, k_f varies ===
             CALL edmatwan2bloch_2d(nbndsub, nrr, ndegen, edmatw_2d, cfac_K, cfac_kf, edmatf_w)
             CALL ZGEMM('N', 'N', nbndsub, nbndsub, nbndsub, &
                         cone, edmatf_w, nbndsub, evec_f, nbndsub, czero, tmp, nbndsub)
             CALL ZGEMM('C', 'N', nbndsub, nbndsub, nbndsub, &
                         cone, evec_K, nbndsub, tmp, nbndsub, czero, edmatf_b, nbndsub)

             IF (ionode) THEN
                DO ibnd = 1, nbndsub
                   DO jbnd = 1, nbndsub
                      WRITE(iunit_scatter, '(I8,3F10.5,2I4,3ES16.8)') ik, xk_f, &
                           ibnd, jbnd, ABS(edmatf_b(ibnd,jbnd))**2, &
                           REAL(edmatf_b(ibnd,jbnd)), AIMAG(edmatf_b(ibnd,jbnd))
                   ENDDO
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    IF (ionode) THEN
       CLOSE(iunit_diag)
       CLOSE(iunit_scatter)
       WRITE(stdout, '(5X,A,I8,A)') 'Interpolated M(2D) at ', nktotf, ' fine k-points'
       WRITE(stdout, '(5X,A)') 'Diagonal: ' // TRIM(prefix_in) // '_edmatf_2d.dat'
       WRITE(stdout, '(5X,A)') 'K->k_f:  ' // TRIM(prefix_in) // '_edmat_scatter_2d.dat'
       WRITE(stdout, '(5X,A)') REPEAT('=', 60)
       FLUSH(stdout)
    ENDIF

  END SUBROUTINE ed_fine_interp_2d

  SUBROUTINE write_edmatw(nbndsub, nrr, irvec, ndegen, edmatw, prefix_in)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nbndsub, nrr
    INTEGER, INTENT(IN) :: irvec(3, nrr), ndegen(nrr)
    COMPLEX(dp), INTENT(IN) :: edmatw(nbndsub, nbndsub, nrr)
    CHARACTER(LEN=*), INTENT(IN) :: prefix_in

    INTEGER :: ir, ii, jj, iunit
    CHARACTER(LEN=256) :: fname

    iunit = 87
    fname = TRIM(prefix_in) // '_edmatw.dat'
    OPEN(iunit, FILE=TRIM(fname), FORM='formatted')
    WRITE(iunit, '(A)') '# Electron-defect matrix elements in Wannier basis'
    WRITE(iunit, '(I6)') nbndsub
    WRITE(iunit, '(I6)') nrr

    ir = 0
    DO WHILE (ir < nrr)
       WRITE(iunit, '(15I5)') ndegen(ir+1:MIN(ir+15, nrr))
       ir = MIN(ir + 15, nrr)
    ENDDO

    DO ir = 1, nrr
       DO jj = 1, nbndsub
          DO ii = 1, nbndsub
             WRITE(iunit, '(3I5, 2I5, 2E20.12)') irvec(:, ir), ii, jj, &
                  REAL(edmatw(ii, jj, ir), dp), AIMAG(edmatw(ii, jj, ir))
          ENDDO
       ENDDO
    ENDDO
    CLOSE(iunit)
  END SUBROUTINE write_edmatw

  SUBROUTINE read_edmatw_file(prefix_in, nbndsub, nrr, ndegen, irvec, edmatw)
    !-----------------------------------------------------------------------
    ! Read M(R) from prefix_edmatw.dat (written by write_edmatw).
    ! Allows skipping Part A and going directly to Part B.
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode, ionode_id, stdout
    USE mp, ONLY : mp_bcast
    USE mp_world, ONLY : world_comm
    USE parallel_include
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: prefix_in
    INTEGER, INTENT(OUT) :: nbndsub, nrr
    INTEGER, ALLOCATABLE, INTENT(OUT) :: ndegen(:), irvec(:,:)
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: edmatw(:,:,:)

    INTEGER :: iunit, ir, ii, jj, r1, r2, r3, idx_i, idx_j, ios
    INTEGER :: ntot, ierr
    REAL(dp) :: re_part, im_part
    REAL(dp), ALLOCATABLE :: rbuf(:)
    CHARACTER(LEN=256) :: fname, line

    iunit = 88
    fname = TRIM(prefix_in) // '_edmatw.dat'
    nbndsub = 0; nrr = 0

    IF (ionode) THEN
       OPEN(iunit, FILE=TRIM(fname), STATUS='old', FORM='formatted', IOSTAT=ios)
       IF (ios /= 0) CALL errore('read_edmatw_file', 'Cannot open ' // TRIM(fname), 1)

       READ(iunit, '(A)') line  ! header
       READ(iunit, *) nbndsub
       READ(iunit, *) nrr

       ALLOCATE(ndegen(nrr))
       ir = 0
       DO WHILE (ir < nrr)
          READ(iunit, *, IOSTAT=ios) ndegen(ir+1:MIN(ir+15, nrr))
          IF (ios /= 0) EXIT
          ir = MIN(ir + 15, nrr)
       ENDDO

       ALLOCATE(irvec(3, nrr))
       ALLOCATE(edmatw(nbndsub, nbndsub, nrr))
       edmatw = (0.0_dp, 0.0_dp)

       DO ir = 1, nrr
          DO jj = 1, nbndsub
             DO ii = 1, nbndsub
                READ(iunit, *) r1, r2, r3, idx_i, idx_j, re_part, im_part
                IF (ii == 1 .AND. jj == 1) irvec(:, ir) = (/r1, r2, r3/)
                edmatw(idx_i, idx_j, ir) = CMPLX(re_part, im_part, dp)
             ENDDO
          ENDDO
       ENDDO

       CLOSE(iunit)
       WRITE(stdout, '(5X,A,A)') 'Read M(R) from ', TRIM(fname)
       WRITE(stdout, '(5X,A,I4,A,I6)') 'nbndsub = ', nbndsub, ', nrr = ', nrr
    ENDIF

    CALL mp_bcast(nbndsub, ionode_id, world_comm)
    CALL mp_bcast(nrr, ionode_id, world_comm)

    IF (.NOT. ionode) THEN
       ALLOCATE(ndegen(nrr))
       ALLOCATE(irvec(3, nrr))
       ALLOCATE(edmatw(nbndsub, nbndsub, nrr))
    ENDIF

    CALL mp_bcast(ndegen, ionode_id, world_comm)
    CALL mp_bcast(irvec, ionode_id, world_comm)

    ! Broadcast complex 3D array as two real 1D arrays
    ntot = nbndsub * nbndsub * nrr
    ALLOCATE(rbuf(ntot))
    IF (ionode) rbuf = RESHAPE(REAL(edmatw, dp), (/ntot/))
    CALL MPI_Bcast(rbuf, ntot, MPI_DOUBLE_PRECISION, ionode_id, world_comm, ierr)
    IF (.NOT. ionode) edmatw = RESHAPE(CMPLX(rbuf, 0.0_dp, dp), (/nbndsub, nbndsub, nrr/))
    IF (ionode) rbuf = RESHAPE(AIMAG(edmatw), (/ntot/))
    CALL MPI_Bcast(rbuf, ntot, MPI_DOUBLE_PRECISION, ionode_id, world_comm, ierr)
    IF (.NOT. ionode) &
         edmatw = CMPLX(REAL(edmatw, dp), RESHAPE(rbuf, (/nbndsub, nbndsub, nrr/)), dp)
    DEALLOCATE(rbuf)

  END SUBROUTINE read_edmatw_file

  SUBROUTINE load_supercell_pot(prefix_in, outdir_in, vf_out)
    !-----------------------------------------------------------------------
    ! Load supercell KS potential via get_vloc_onthefly on ionode.
    !
    ! Uses MPI_COMM_SELF so ionode runs as a standalone single-rank QE.
    ! Does NOT call clean_pw — the caller is responsible for cleanup.
    ! This allows consecutive calls (V_d then V_p) where the second call
    ! overwrites the first call's QE state in-place (same dimensions).
    ! Calling clean_pw between loads would corrupt QE's internal state.
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode, ionode_id, stdout
    USE mp_world, ONLY : world_comm
    USE mp_pools, ONLY : intra_pool_comm, inter_pool_comm, npool, my_pool_id
    USE mp_bands, ONLY : intra_bgrp_comm, inter_bgrp_comm, nbgrp
    USE mp_images, ONLY : intra_image_comm, nimage
    USE mp, ONLY : mp_bcast
    USE io_files, ONLY : prefix, tmp_dir
    USE edic_mod, ONLY : V_file
    USE parallel_include
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: prefix_in, outdir_in
    TYPE(V_file), INTENT(INOUT) :: vf_out

    INTERFACE
       SUBROUTINE get_vloc_onthefly(prefix_in, outdir_in, vf_out, is_noncolin, bxc_out)
         USE kinds, ONLY : dp
         USE edic_mod, ONLY : V_file
         CHARACTER(LEN=*), INTENT(IN) :: prefix_in, outdir_in
         TYPE(V_file), INTENT(INOUT) :: vf_out
         LOGICAL, INTENT(IN) :: is_noncolin
         REAL(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL :: bxc_out(:,:)
       END SUBROUTINE
    END INTERFACE

    INTEGER :: save_world, save_pool, save_inter_pool
    INTEGER :: save_bgrp, save_inter_bgrp, save_image
    INTEGER :: save_npool, save_pool_id, save_nbgrp, save_nimage
    INTEGER :: nrtot, ierr
    REAL(dp), ALLOCATABLE :: pot_tmp(:)
    CHARACTER(LEN=256) :: save_prefix, save_tmpdir

    ! Save communicator state
    save_world = world_comm
    save_pool = intra_pool_comm
    save_inter_pool = inter_pool_comm
    save_bgrp = intra_bgrp_comm
    save_inter_bgrp = inter_bgrp_comm
    save_image = intra_image_comm
    save_npool = npool
    save_pool_id = my_pool_id
    save_nbgrp = nbgrp
    save_nimage = nimage
    save_prefix = prefix
    save_tmpdir = tmp_dir

    ! ionode reads supercell as standalone single-rank QE
    ! NO clean_pw before — let read_file_new overwrite previous state in-place
    IF (ionode) THEN
       world_comm = MPI_COMM_SELF
       intra_pool_comm = MPI_COMM_SELF
       inter_pool_comm = MPI_COMM_SELF
       intra_bgrp_comm = MPI_COMM_SELF
       inter_bgrp_comm = MPI_COMM_SELF
       intra_image_comm = MPI_COMM_SELF
       npool = 1
       my_pool_id = 0
       nbgrp = 1
       nimage = 1

       CALL get_vloc_onthefly(prefix_in, outdir_in, vf_out, .FALSE.)
       ! NO clean_pw after — preserve state for potential next load
    ENDIF

    ! Restore communicator state
    world_comm = save_world
    intra_pool_comm = save_pool
    inter_pool_comm = save_inter_pool
    intra_bgrp_comm = save_bgrp
    inter_bgrp_comm = save_inter_bgrp
    intra_image_comm = save_image
    npool = save_npool
    my_pool_id = save_pool_id
    nbgrp = save_nbgrp
    nimage = save_nimage
    prefix = save_prefix
    tmp_dir = save_tmpdir

    ! Broadcast from ionode to all ranks
    CALL mp_bcast(vf_out%nr1, ionode_id, world_comm)
    CALL mp_bcast(vf_out%nr2, ionode_id, world_comm)
    CALL mp_bcast(vf_out%nr3, ionode_id, world_comm)
    CALL mp_bcast(vf_out%nr1x, ionode_id, world_comm)
    CALL mp_bcast(vf_out%nr2x, ionode_id, world_comm)
    CALL mp_bcast(vf_out%nr3x, ionode_id, world_comm)
    CALL mp_bcast(vf_out%nat, ionode_id, world_comm)
    CALL mp_bcast(vf_out%ntyp, ionode_id, world_comm)
    CALL mp_bcast(vf_out%alat, ionode_id, world_comm)
    CALL mp_bcast(vf_out%ibrav, ionode_id, world_comm)
    CALL mp_bcast(vf_out%celldm, ionode_id, world_comm)
    CALL mp_bcast(vf_out%at, ionode_id, world_comm)
    CALL mp_bcast(vf_out%omega, ionode_id, world_comm)

    nrtot = vf_out%nr1 * vf_out%nr2 * vf_out%nr3

    ! Trim ionode's pot to logical grid (remove FFT padding if any)
    IF (ionode) THEN
       ALLOCATE(pot_tmp(nrtot))
       pot_tmp(:) = vf_out%pot(1:nrtot)
       DEALLOCATE(vf_out%pot)
       ALLOCATE(vf_out%pot(nrtot))
       vf_out%pot(:) = pot_tmp(:)
       DEALLOCATE(pot_tmp)
    ELSE
       ALLOCATE(vf_out%pot(nrtot))
       ALLOCATE(vf_out%ityp(vf_out%nat))
       ALLOCATE(vf_out%tau(3, vf_out%nat))
    ENDIF
    CALL mp_bcast(vf_out%pot, ionode_id, world_comm)
    CALL mp_bcast(vf_out%ityp, ionode_id, world_comm)
    CALL mp_bcast(vf_out%tau, ionode_id, world_comm)

    IF (ionode) WRITE(stdout, '(5X,A,3I6,A,I6)') &
         '  grid: ', vf_out%nr1, vf_out%nr2, vf_out%nr3, &
         '  nat:', vf_out%nat

  END SUBROUTINE load_supercell_pot


  SUBROUTINE load_pot_from_file(potfile, vf_out)
    !-----------------------------------------------------------------------
    ! Read supercell KS potential from a Gaussian cube file written by
    ! extract_pot.x. ionode reads, broadcasts to all ranks.
    !
    ! The cube file contains V_total = V_Hxc + V_ionic in Ry on the
    ! supercell FFT grid. Also extracts atomic positions (nat, ityp, tau)
    ! needed for the nonlocal matrix element computation.
    !
    ! Cube format: x outermost loop, z innermost.
    ! Stored in vf_out%pot with x-fastest ordering (QE convention):
    !   pot(ix + iy*nr1 + iz*nr1*nr2 + 1)
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode, ionode_id, stdout
    USE mp_world, ONLY : world_comm
    USE mp, ONLY : mp_bcast
    USE cell_base, ONLY : alat_prim => alat
    USE edic_mod, ONLY : V_file
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: potfile
    TYPE(V_file), INTENT(INOUT) :: vf_out

    INTEGER :: iunit, nrtot, ia, irx, iry, irz, idx
    INTEGER :: nr1, nr2, nr3
    REAL(dp) :: voxel1(3), voxel2(3), voxel3(3), origin(3)
    REAL(dp) :: at_charge, val
    CHARACTER(LEN=256) :: line

    iunit = 98

    IF (ionode) THEN
       OPEN(iunit, FILE=TRIM(potfile), STATUS='old', ACTION='read')

       ! Header: 2 comment lines
       READ(iunit, '(A)') line
       READ(iunit, '(A)') line

       ! Number of atoms and origin
       READ(iunit, *) vf_out%nat, origin(1), origin(2), origin(3)

       ! Grid dimensions and voxel vectors
       READ(iunit, *) nr1, voxel1(1), voxel1(2), voxel1(3)
       READ(iunit, *) nr2, voxel2(1), voxel2(2), voxel2(3)
       READ(iunit, *) nr3, voxel3(1), voxel3(2), voxel3(3)
       vf_out%nr1 = nr1; vf_out%nr2 = nr2; vf_out%nr3 = nr3
       vf_out%nr1x = nr1; vf_out%nr2x = nr2; vf_out%nr3x = nr3
       nrtot = nr1 * nr2 * nr3

       ! Reconstruct lattice vectors and alat from voxels
       ! Cube stores positions in bohr. Convert to QE convention:
       !   at in units of alat, tau in units of alat
       ! Use the primitive cell's alat (from cell_base, still loaded)
       vf_out%alat = alat_prim
       vf_out%at(:,1) = voxel1(:) * DBLE(nr1) / alat_prim
       vf_out%at(:,2) = voxel2(:) * DBLE(nr2) / alat_prim
       vf_out%at(:,3) = voxel3(:) * DBLE(nr3) / alat_prim
       vf_out%omega = ABS( &
            vf_out%at(1,1)*(vf_out%at(2,2)*vf_out%at(3,3) - vf_out%at(3,2)*vf_out%at(2,3)) - &
            vf_out%at(2,1)*(vf_out%at(1,2)*vf_out%at(3,3) - vf_out%at(3,2)*vf_out%at(1,3)) + &
            vf_out%at(3,1)*(vf_out%at(1,2)*vf_out%at(2,3) - vf_out%at(2,2)*vf_out%at(1,3)) ) &
            * alat_prim**3
       vf_out%ibrav = 0
       vf_out%celldm = 0.0_dp
       vf_out%celldm(1) = alat_prim

       ! Atomic positions (cube format: ityp, charge, x, y, z in bohr)
       ! Convert from bohr to units of alat (QE convention for get_betavkb)
       vf_out%ntyp = 0
       ALLOCATE(vf_out%ityp(vf_out%nat))
       ALLOCATE(vf_out%tau(3, vf_out%nat))
       DO ia = 1, vf_out%nat
          READ(iunit, *) vf_out%ityp(ia), at_charge, &
               vf_out%tau(1,ia), vf_out%tau(2,ia), vf_out%tau(3,ia)
          ! Convert tau from bohr to alat units
          vf_out%tau(:,ia) = vf_out%tau(:,ia) / alat_prim
          IF (vf_out%ityp(ia) > vf_out%ntyp) vf_out%ntyp = vf_out%ityp(ia)
       ENDDO

       ! Read volumetric data: cube format packs up to 6 values per line.
       ! Cube order: x outermost, z innermost.
       ! Store in x-fastest order: pot(ix + iy*nr1 + iz*nr1*nr2 + 1)
       BLOCK
          REAL(dp), ALLOCATABLE :: raw(:)
          INTEGER :: iraw, nlines, nvals_line, ipos
          REAL(dp) :: buf(6)

          ALLOCATE(raw(nrtot))
          iraw = 0
          DO WHILE (iraw < nrtot)
             ! Read one line — may contain up to 6 values
             nvals_line = MIN(6, nrtot - iraw)
             READ(iunit, *, END=900) buf(1:nvals_line)
             DO ipos = 1, nvals_line
                iraw = iraw + 1
                raw(iraw) = buf(ipos)
             ENDDO
          ENDDO
900       CONTINUE

          ! Reorder: cube is (x outer, z inner) → QE is (x fastest)
          ALLOCATE(vf_out%pot(nrtot))
          iraw = 0
          DO irx = 0, nr1 - 1
             DO iry = 0, nr2 - 1
                DO irz = 0, nr3 - 1
                   iraw = iraw + 1
                   idx = irx + iry * nr1 + irz * nr1 * nr2 + 1
                   vf_out%pot(idx) = raw(iraw)
                ENDDO
             ENDDO
          ENDDO
          DEALLOCATE(raw)
       END BLOCK

       CLOSE(iunit)

       WRITE(stdout, '(5X,A,A)') 'Read potential from cube: ', TRIM(potfile)
       WRITE(stdout, '(5X,A,3I6,A,I6)') '  grid: ', nr1, nr2, nr3, '  nat:', vf_out%nat
       WRITE(stdout, '(5X,A,2ES14.6)') '  V min/max = ', MINVAL(vf_out%pot), MAXVAL(vf_out%pot)
    ENDIF

    ! Broadcast to all ranks
    CALL mp_bcast(vf_out%nr1, ionode_id, world_comm)
    CALL mp_bcast(vf_out%nr2, ionode_id, world_comm)
    CALL mp_bcast(vf_out%nr3, ionode_id, world_comm)
    CALL mp_bcast(vf_out%nr1x, ionode_id, world_comm)
    CALL mp_bcast(vf_out%nr2x, ionode_id, world_comm)
    CALL mp_bcast(vf_out%nr3x, ionode_id, world_comm)
    CALL mp_bcast(vf_out%nat, ionode_id, world_comm)
    CALL mp_bcast(vf_out%ntyp, ionode_id, world_comm)
    CALL mp_bcast(vf_out%alat, ionode_id, world_comm)
    CALL mp_bcast(vf_out%ibrav, ionode_id, world_comm)
    CALL mp_bcast(vf_out%omega, ionode_id, world_comm)
    CALL mp_bcast(vf_out%celldm, ionode_id, world_comm)
    CALL mp_bcast(vf_out%at, ionode_id, world_comm)

    nrtot = vf_out%nr1 * vf_out%nr2 * vf_out%nr3
    IF (.NOT. ionode) THEN
       ALLOCATE(vf_out%pot(nrtot))
       ALLOCATE(vf_out%ityp(vf_out%nat))
       ALLOCATE(vf_out%tau(3, vf_out%nat))
    ENDIF

    BLOCK
       USE parallel_include
       INTEGER :: ierr
       CALL MPI_Bcast(vf_out%pot, nrtot, MPI_DOUBLE_PRECISION, &
            ionode_id, world_comm, ierr)
       CALL MPI_Bcast(vf_out%ityp, vf_out%nat, MPI_INTEGER, &
            ionode_id, world_comm, ierr)
       CALL MPI_Bcast(vf_out%tau, 3*vf_out%nat, MPI_DOUBLE_PRECISION, &
            ionode_id, world_comm, ierr)
    END BLOCK

  END SUBROUTINE load_pot_from_file




  SUBROUTINE read_filukk_edi(fname, nbnd, nkstot, nks, nbndsub)
    !-----------------------------------------------------------------------
    ! Read Wannier rotation matrices from filukk file.
    ! Sets u_mat, u_mat_opt, lwindow, excluded_band, num_bands, n_wannier
    ! in wann_common module, so ed_coarse_calc can use them.
    ! The filukk stores the COMBINED rotation u_kc = u_opt * u.
    ! We store u_kc in u_mat and set num_bands = n_wannier so
    ! ed_coarse_calc uses cu = u_mat directly (no u_mat_opt needed).
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode, ionode_id, stdout
    USE mp, ONLY : mp_bcast
    USE mp_world, ONLY : world_comm
    USE wann_common, ONLY : u_mat, u_mat_opt, lwindow, excluded_band, &
                             n_wannier, num_bands, iknum, wann_centers
    USE global_var, ONLY : nbndep, nbndskip, ibndkept
    USE parallelism, ONLY : fkbounds
    USE parallel_include
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: fname
    INTEGER, INTENT(IN) :: nbnd, nkstot, nks, nbndsub

    INTEGER :: ik, ibnd, iw, iunit, nbndep_f, nbndskip_f
    INTEGER :: lower_bnd, upper_bnd, ntot, ierr
    COMPLEX(dp), ALLOCATABLE :: u_kc_all(:,:,:)
    REAL(dp), ALLOCATABLE :: rbuf(:)
    LOGICAL, ALLOCATABLE :: lwindow_all(:,:)
    INTEGER, ALLOCATABLE :: lwin_int(:)

    iunit = 89
    iknum = nkstot

    ! nbndsub may not be broadcast to non-ionode ranks yet.
    ! Use ionode's value and broadcast.
    IF (ionode) n_wannier = nbndsub
    CALL MPI_Bcast(n_wannier, 1, MPI_INTEGER, ionode_id, world_comm, ierr)

    IF (ionode) THEN
       OPEN(UNIT=iunit, FILE=TRIM(fname), FORM='formatted', STATUS='old')

       ! Read nbndep, nbndskip
       READ(iunit, *) nbndep_f, nbndskip_f
    ENDIF
    CALL MPI_Barrier(world_comm, ierr)
    CALL MPI_Bcast(nbndep_f, 1, MPI_INTEGER, ionode_id, world_comm, ierr)
    CALL MPI_Bcast(nbndskip_f, 1, MPI_INTEGER, ionode_id, world_comm, ierr)

    nbndep = nbndep_f
    nbndskip = nbndskip_f
    num_bands = nbndep

    ! Read ibndkept
    IF (.NOT. ALLOCATED(ibndkept)) ALLOCATE(ibndkept(nbndep))
    IF (ionode) THEN
       DO ibnd = 1, nbndep
          READ(iunit, *) ibndkept(ibnd)
       ENDDO
    ENDIF
    CALL MPI_Barrier(world_comm, ierr)
    CALL MPI_Bcast(ibndkept, nbndep, MPI_INTEGER, ionode_id, world_comm, ierr)

    ! Read u_kc (combined rotation) for all k-points
    ALLOCATE(u_kc_all(nbndep, n_wannier, nkstot))
    u_kc_all = (0.0_dp, 0.0_dp)
    IF (ionode) THEN
       DO ik = 1, nkstot
          DO ibnd = 1, nbndep
             DO iw = 1, n_wannier
                READ(iunit, *) u_kc_all(ibnd, iw, ik)
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    ! Broadcast u_kc using direct MPI_Bcast (QE wrapper has issues with large arrays)
    ntot = nbndep * n_wannier * nkstot
    ALLOCATE(rbuf(2*ntot))
    IF (ionode) THEN
       DO iw = 1, ntot
          rbuf(2*iw-1) = REAL(u_kc_all(MOD(iw-1,nbndep)+1, &
               MOD((iw-1)/nbndep,n_wannier)+1, (iw-1)/(nbndep*n_wannier)+1), dp)
          rbuf(2*iw) = AIMAG(u_kc_all(MOD(iw-1,nbndep)+1, &
               MOD((iw-1)/nbndep,n_wannier)+1, (iw-1)/(nbndep*n_wannier)+1))
       ENDDO
    ENDIF
    CALL MPI_Barrier(world_comm, ierr)
    CALL MPI_Bcast(rbuf, 2*ntot, MPI_DOUBLE_PRECISION, ionode_id, world_comm, ierr)
    IF (.NOT. ionode) THEN
       DO iw = 1, ntot
          u_kc_all(MOD(iw-1,nbndep)+1, MOD((iw-1)/nbndep,n_wannier)+1, &
               (iw-1)/(nbndep*n_wannier)+1) = CMPLX(rbuf(2*iw-1), rbuf(2*iw), dp)
       ENDDO
    ENDIF
    DEALLOCATE(rbuf)

    ! Read lwindow for all k-points (broadcast as integer via MPI_Bcast)
    ALLOCATE(lwindow_all(nbndep, nkstot))
    ntot = nbndep * nkstot
    ALLOCATE(lwin_int(ntot))
    lwin_int = 0
    IF (ionode) THEN
       DO ik = 1, nkstot
          DO ibnd = 1, nbndep
             READ(iunit, *) lwindow_all(ibnd, ik)
             IF (lwindow_all(ibnd, ik)) lwin_int((ik-1)*nbndep + ibnd) = 1
          ENDDO
       ENDDO
    ENDIF
    CALL MPI_Bcast(lwin_int, ntot, MPI_INTEGER, ionode_id, world_comm, ierr)
    DO ik = 1, nkstot
       DO ibnd = 1, nbndep
          lwindow_all(ibnd, ik) = lwin_int((ik-1)*nbndep + ibnd) > 0
       ENDDO
    ENDDO
    DEALLOCATE(lwin_int)

    ! Read excluded_band (broadcast as integer via MPI_Bcast)
    IF (.NOT. ALLOCATED(excluded_band)) ALLOCATE(excluded_band(nbnd))
    ALLOCATE(lwin_int(nbnd))
    lwin_int = 0
    IF (ionode) THEN
       DO ibnd = 1, nbnd
          READ(iunit, *) excluded_band(ibnd)
          IF (excluded_band(ibnd)) lwin_int(ibnd) = 1
       ENDDO
    ENDIF
    CALL MPI_Bcast(lwin_int, nbnd, MPI_INTEGER, ionode_id, world_comm, ierr)
    excluded_band = lwin_int > 0
    DEALLOCATE(lwin_int)

    ! Read wann_centers
    IF (.NOT. ALLOCATED(wann_centers)) ALLOCATE(wann_centers(3, n_wannier))
    IF (ionode) THEN
       DO iw = 1, n_wannier
          READ(iunit, *) wann_centers(:, iw)
       ENDDO
       CLOSE(iunit)
    ENDIF
    CALL MPI_Bcast(wann_centers, 3*n_wannier, MPI_DOUBLE_PRECISION, ionode_id, world_comm, ierr)

    ! Extract pool-local k-points
    CALL fkbounds(nkstot, lower_bnd, upper_bnd)

    ! Store u_kc in u_mat (pool-local), set num_bands = n_wannier
    ! so ed_coarse_calc uses cu = u_mat directly
    IF (ALLOCATED(u_mat)) DEALLOCATE(u_mat)
    ALLOCATE(u_mat(nbndep, n_wannier, nks))
    u_mat(:,:,1:nks) = u_kc_all(:,:,lower_bnd:upper_bnd)

    ! u_mat_opt not needed (num_bands = n_wannier path), but allocate for safety
    IF (ALLOCATED(u_mat_opt)) DEALLOCATE(u_mat_opt)
    ALLOCATE(u_mat_opt(nbndep, nbndep, nks))
    u_mat_opt = (0.0_dp, 0.0_dp)

    ! lwindow (pool-local)
    IF (ALLOCATED(lwindow)) DEALLOCATE(lwindow)
    ALLOCATE(lwindow(nbndep, nks))
    lwindow(:,1:nks) = lwindow_all(:,lower_bnd:upper_bnd)

    DEALLOCATE(u_kc_all, lwindow_all)

    IF (ionode) THEN
       WRITE(stdout, '(5X,A,A)') 'Read Wannier data from ', TRIM(fname)
       WRITE(stdout, '(5X,A,I4,A,I4,A,I4)') &
            'nbndep=', nbndep, ' n_wannier=', n_wannier, ' nkstot=', nkstot
    ENDIF

  END SUBROUTINE read_filukk_edi


  SUBROUTINE generate_nscf_input(filki, filkf, edi_outdir_in)
    !-----------------------------------------------------------------------
    ! Auto-generate nscf_custom.in for k-points not in the current NSCF data.
    ! Reads parameters from QE modules (populated by read_file).
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode, stdout
    USE io_files, ONLY : prefix, pseudo_dir
    USE cell_base, ONLY : at, alat, tpiba2
    USE ions_base, ONLY : nat, ntyp => nsp, ityp, tau, atm
    USE wvfct, ONLY : nbnd
    USE gvect, ONLY : gcutm
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: filki, filkf, edi_outdir_in
    REAL(dp) :: ecutwfc_val, ecutrho_val

    INTEGER :: nki, nkf, nk_merged, ik, jk, iunit
    REAL(dp), ALLOCATABLE :: xki(:,:), xkf(:,:), xk_merged(:,:)
    REAL(dp) :: diff(3), dist
    LOGICAL :: is_dup
    INTEGER :: ia

    CALL read_kpoints_file(filki, nki, xki)
    CALL read_kpoints_file(filkf, nkf, xkf)

    IF (.NOT. ionode) THEN
       IF (ALLOCATED(xki)) DEALLOCATE(xki)
       IF (ALLOCATED(xkf)) DEALLOCATE(xkf)
       RETURN
    ENDIF

    ! Merge ki and kf lists, removing duplicates
    ALLOCATE(xk_merged(3, nki + nkf))
    nk_merged = 0
    DO ik = 1, nki
       is_dup = .FALSE.
       DO jk = 1, nk_merged
          diff = xki(:,ik) - xk_merged(:,jk)
          diff = diff - NINT(diff)
          IF (SUM(diff**2) < 1.0d-10) THEN
             is_dup = .TRUE.
             EXIT
          ENDIF
       ENDDO
       IF (.NOT. is_dup) THEN
          nk_merged = nk_merged + 1
          xk_merged(:, nk_merged) = xki(:, ik)
       ENDIF
    ENDDO
    DO ik = 1, nkf
       is_dup = .FALSE.
       DO jk = 1, nk_merged
          diff = xkf(:,ik) - xk_merged(:,jk)
          diff = diff - NINT(diff)
          IF (SUM(diff**2) < 1.0d-10) THEN
             is_dup = .TRUE.
             EXIT
          ENDIF
       ENDDO
       IF (.NOT. is_dup) THEN
          nk_merged = nk_merged + 1
          xk_merged(:, nk_merged) = xkf(:, ik)
       ENDIF
    ENDDO

    ! Compute cutoffs from QE internal variables (Ry units)
    ! In QE: ecutwfc = ecfixed or from input; gcutm = ecutrho/tpiba2
    ! We recover ecutrho = gcutm * tpiba2
    ecutrho_val = gcutm * tpiba2
    ! ecutwfc: use the ratio ecutrho/ecutwfc which is typically 4 (NC) or 8-12 (US/PAW)
    ! We can get ecutwfc from gvecw module
    BLOCK
       USE gvecw, ONLY : ecutwfc_in => ecutwfc
       ecutwfc_val = ecutwfc_in
    END BLOCK

    ! Write nscf_custom.in
    iunit = 90
    OPEN(iunit, FILE='nscf_custom.in', FORM='formatted')
    WRITE(iunit, '(A)')       '&CONTROL'
    WRITE(iunit, '(A)')       "  calculation = 'bands'"
    WRITE(iunit, '(A,A,A)')   "  prefix = '", TRIM(prefix), "'"
    WRITE(iunit, '(A,A,A)')   "  outdir = '", TRIM(edi_outdir_in), "'"
    WRITE(iunit, '(A,A,A)')   "  pseudo_dir = '", TRIM(pseudo_dir), "'"
    WRITE(iunit, '(A)')       '/'
    WRITE(iunit, '(A)')       '&SYSTEM'
    WRITE(iunit, '(A)')       '  ibrav = 0'
    WRITE(iunit, '(A,F12.6)') '  ecutwfc = ', ecutwfc_val
    WRITE(iunit, '(A,F12.6)') '  ecutrho = ', ecutrho_val
    WRITE(iunit, '(A,I4)')    '  nat = ', nat
    WRITE(iunit, '(A,I4)')    '  ntyp = ', ntyp
    WRITE(iunit, '(A,I4)')    '  nbnd = ', nbnd
    WRITE(iunit, '(A)')       '/'
    WRITE(iunit, '(A)')       '&ELECTRONS'
    WRITE(iunit, '(A)')       '/'
    WRITE(iunit, '(A)')       'ATOMIC_SPECIES'
    DO ia = 1, ntyp
       WRITE(iunit, '(A,A,A)') TRIM(atm(ia)), '  0.0  ', TRIM(atm(ia)) // '.upf'
    ENDDO
    WRITE(iunit, '(A)')       'CELL_PARAMETERS alat'
    WRITE(iunit, '(3F16.10)') at(:,1)
    WRITE(iunit, '(3F16.10)') at(:,2)
    WRITE(iunit, '(3F16.10)') at(:,3)
    WRITE(iunit, '(A)')       'ATOMIC_POSITIONS crystal'
    DO ia = 1, nat
       WRITE(iunit, '(A,3F16.10)') TRIM(atm(ityp(ia))), tau(:,ia)
    ENDDO
    WRITE(iunit, '(A)')       'K_POINTS crystal'
    WRITE(iunit, '(I6)')      nk_merged
    DO ik = 1, nk_merged
       WRITE(iunit, '(3F16.10,A)') xk_merged(:,ik), '  1.0'
    ENDDO
    CLOSE(iunit)

    WRITE(stdout, '(/,5X,A)') REPEAT('=', 60)
    WRITE(stdout, '(5X,A)')   'NSCF wavefunctions not found for required k-points.'
    WRITE(stdout, '(5X,A,I6,A)') 'Auto-generated: nscf_custom.in (', nk_merged, ' k-points)'
    WRITE(stdout, '(5X,A)')   ''
    WRITE(stdout, '(5X,A)')   'Please run:  pw.x < nscf_custom.in > nscf_custom.out'
    WRITE(stdout, '(5X,A)')   'Then re-run EDI with the same input.'
    WRITE(stdout, '(5X,A)')   REPEAT('=', 60)

    DEALLOCATE(xki, xkf, xk_merged)
  END SUBROUTINE generate_nscf_input


  SUBROUTINE ed_direct_from_files(filki, filkf, prefix_in, band_ed_str)
    !-----------------------------------------------------------------------
    ! Direct calculation mode: compute M(k_i, k_f) from wavefunctions
    ! for arbitrary k-points read from files.
    !
    ! band_ed_str: band selection string, e.g. '13-17' for bands 13 to 17.
    !   If empty, all bands are used.
    !
    ! Prerequisites: primitive cell state restored (read_file_new done),
    !   supercell potentials loaded in V_d, V_p, V_colin.
    ! Requires NSCF wavefunctions at ALL requested k-points.
    ! Supercell potentials are loaded internally.
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode, ionode_id, stdout
    USE mp, ONLY : mp_bcast, mp_sum
    USE mp_world, ONLY : world_comm
    USE mp_global, ONLY : inter_pool_comm
    USE mp_pools, ONLY : npool, my_pool_id
    USE cell_base, ONLY : at
    USE klist, ONLY : xk, igk_k, ngk, nkstot, nks
    USE wvfct, ONLY : npwx, nbnd
    USE noncollin_module, ONLY : npol, lspinorb
    USE fft_base, ONLY : dffts
    USE fft_interfaces, ONLY : invfft
    USE io_files, ONLY : restart_dir
    USE pw_restart_new, ONLY : read_collected_wfc
    USE edic_mod, ONLY : V_d, V_p, V_colin
    USE uspp_param, ONLY : nh
    USE uspp, ONLY : dvan, dvan_so
    USE becmod, ONLY : bec_type, calbec, allocate_bec_type, deallocate_bec_type
    USE ep_constants, ONLY : czero, cone
    USE constants, ONLY : tpi
    USE parallelism, ONLY : fkbounds
    USE parallel_include
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: filki, filkf, prefix_in, band_ed_str

    INTEGER :: nki, nkf, iki, ikf, ik_nscf, ibnd, jbnd, ierr_mpi
    INTEGER :: ig, ir, inr, irx, iry, irz, ir1mod, ir2mod, ir3mod, irnmod
    INTEGER :: iunit, lower_bnd, upper_bnd, nbnd_kept
    INTEGER :: ibnd_min, ibnd_max, ipos, ios_parse
    INTEGER :: na, nt, ih, jh, ikb, jkb, ijkb0, nkb_d, nkb_p
    REAL(dp), ALLOCATABLE :: xki(:,:), xkf(:,:), xk_cryst_nscf(:,:)
    REAL(dp) :: qcryst(3), d1, d2, d3, arg, diff(3)
    COMPLEX(dp) :: phase, mlocal, mnl_d, mnl_p
    COMPLEX(dp), ALLOCATABLE :: psir_ki(:,:,:), psir_kf(:,:,:)
    COMPLEX(dp), ALLOCATABLE :: evc_tmp(:,:), psic_tmp(:)
    COMPLEX(dp), ALLOCATABLE :: V_folded(:), edmatkq(:,:)
    COMPLEX(dp), ALLOCATABLE :: edmatkq_loc(:,:), edmatkq_nl(:,:)
    COMPLEX(dp), ALLOCATABLE :: vkb_d(:,:), vkb_p(:,:)
    INTEGER :: ipol_d
    INTEGER, ALLOCATABLE :: iki_map(:), ikf_map(:)
    CHARACTER(LEN=256) :: fname

    ! Read k-point files
    CALL read_kpoints_file(filki, nki, xki)
    CALL read_kpoints_file(filkf, nkf, xkf)

    IF (ionode) THEN
       WRITE(stdout, '(/,5X,A)') REPEAT('=', 60)
       WRITE(stdout, '(5X,A)')   'Direct calculation: M(k_i, k_f) from wavefunctions'
       WRITE(stdout, '(5X,A)')   REPEAT('=', 60)
       WRITE(stdout, '(5X,A,I6,A,I6,A,I10)') &
            'nki=', nki, '  nkf=', nkf, '  pairs=', nki * nkf
       FLUSH(stdout)
    ENDIF

    ! Build crystal-coord map of NSCF k-points
    CALL fkbounds(nkstot, lower_bnd, upper_bnd)
    ALLOCATE(xk_cryst_nscf(3, nkstot))
    xk_cryst_nscf = 0.0_dp
    xk_cryst_nscf(:, lower_bnd:upper_bnd) = xk(:, 1:nks)
    CALL cryst_to_cart(nks, xk_cryst_nscf(:, lower_bnd:upper_bnd), at, -1)
    CALL mp_sum(xk_cryst_nscf, inter_pool_comm)

    ! Map user k-points to NSCF indices
    ALLOCATE(iki_map(nki), ikf_map(nkf))
    iki_map = 0
    ikf_map = 0

    DO iki = 1, nki
       DO ik_nscf = 1, nkstot
          diff = xki(:,iki) - xk_cryst_nscf(:,ik_nscf)
          diff = diff - NINT(diff)
          IF (SUM(diff**2) < 1.0d-8) THEN
             iki_map(iki) = ik_nscf
             EXIT
          ENDIF
       ENDDO
    ENDDO
    DO ikf = 1, nkf
       DO ik_nscf = 1, nkstot
          diff = xkf(:,ikf) - xk_cryst_nscf(:,ik_nscf)
          diff = diff - NINT(diff)
          IF (SUM(diff**2) < 1.0d-8) THEN
             ikf_map(ikf) = ik_nscf
             EXIT
          ENDIF
       ENDDO
    ENDDO

    ! Check if all k-points are found
    IF (ANY(iki_map == 0) .OR. ANY(ikf_map == 0)) THEN
       IF (ionode) THEN
          WRITE(stdout, '(5X,A)') 'Some k-points not found in NSCF data.'
          DO iki = 1, nki
             IF (iki_map(iki) == 0) WRITE(stdout, '(5X,A,3F10.5)') &
                  '  Missing ki: ', xki(:,iki)
          ENDDO
          DO ikf = 1, nkf
             IF (ikf_map(ikf) == 0) WRITE(stdout, '(5X,A,3F10.5)') &
                  '  Missing kf: ', xkf(:,ikf)
          ENDDO
       ENDIF
       ! Auto-generate NSCF input and stop
       BLOCK
          USE io_files, ONLY : tmp_dir_loc => tmp_dir
          CALL generate_nscf_input(filki, filkf, tmp_dir_loc)
       END BLOCK
       DEALLOCATE(xki, xkf, xk_cryst_nscf, iki_map, ikf_map)
       CALL errore('ed_direct_from_files', &
            'k-points not in NSCF data. Run pw.x with nscf_custom.in first.', 1)
    ENDIF

    IF (ionode) THEN
       WRITE(stdout, '(5X,A)') 'All k-points found in NSCF data.'
       FLUSH(stdout)
    ENDIF

    ! Parse band_ed string (e.g. '13-17') to get band range
    IF (LEN_TRIM(band_ed_str) > 0) THEN
       ipos = INDEX(band_ed_str, '-')
       IF (ipos > 0) THEN
          READ(band_ed_str(1:ipos-1), *, IOSTAT=ios_parse) ibnd_min
          READ(band_ed_str(ipos+1:), *, IOSTAT=ios_parse) ibnd_max
       ELSE
          READ(band_ed_str, *, IOSTAT=ios_parse) ibnd_min
          ibnd_max = ibnd_min
       ENDIF
       IF (ibnd_min < 1) ibnd_min = 1
       IF (ibnd_max > nbnd) ibnd_max = nbnd
       IF (ionode) WRITE(stdout, '(5X,A,I4,A,I4,A,I4,A)') &
            'band_ed = ', ibnd_min, ' - ', ibnd_max, &
            '  (', ibnd_max - ibnd_min + 1, ' bands selected)'
    ELSE
       ibnd_min = 1
       ibnd_max = nbnd
       IF (ionode) WRITE(stdout, '(5X,A,I4,A)') &
            'band_ed not set, using all ', nbnd, ' bands'
    ENDIF
    nbnd_kept = ibnd_max - ibnd_min + 1

    ! Potentials (V_colin) and primitive cell state already set up by caller
    ! Panel broadcast: cache local wfcs, broadcast per pool (same as double-FT)
    ALLOCATE(psic_tmp(dffts%nnr))
    ALLOCATE(evc_tmp(npwx * npol, nbnd))
    ! Count nonlocal projectors
    nkb_d = 0
    DO nt = 1, V_d%ntyp
       DO na = 1, V_d%nat
          IF (V_d%ityp(na) == nt) nkb_d = nkb_d + nh(nt)
       ENDDO
    ENDDO
    nkb_p = 0
    DO nt = 1, V_p%ntyp
       DO na = 1, V_p%nat
          IF (V_p%ityp(na) == nt) nkb_p = nkb_p + nh(nt)
       ENDDO
    ENDDO
    ALLOCATE(vkb_d(npwx, nkb_d))
    ALLOCATE(vkb_p(npwx, nkb_p))
    ALLOCATE(V_folded(dffts%nnr))
    ALLOCATE(edmatkq(nbnd_kept, nbnd_kept))
    ALLOCATE(edmatkq_loc(nbnd_kept, nbnd_kept))
    ALLOCATE(edmatkq_nl(nbnd_kept, nbnd_kept))
    ALLOCATE(psir_ki(dffts%nnr, npol, nbnd_kept))
    ALLOCATE(psir_kf(dffts%nnr, npol, nbnd_kept))

    BLOCK
       COMPLEX(dp), ALLOCATABLE :: psir_cache(:,:,:,:), psir_recv(:,:,:,:)
       COMPLEX(dp), ALLOCATABLE :: becd_cache(:,:,:,:), becd_recv(:,:,:,:)
       COMPLEX(dp), ALLOCATABLE :: becpc_cache(:,:,:,:), becpc_recv(:,:,:,:)
       COMPLEX(dp), ALLOCATABLE :: becd_ki(:,:,:), becpc_ki(:,:,:)
       TYPE(bec_type) :: becp_tmp_d, becp_tmp_p
       INTEGER :: ik_cache, ib_cache, ig_cache, npairs_done, ipol_c
       INTEGER :: ip, src_lower, src_nks, nkbase, nkrest, nks_max
       INTEGER :: ikf_local, ik_pool, ik_l

       nkbase = nkstot / npool
       nkrest = MOD(nkstot, npool)
       nks_max = nkbase
       IF (nkrest > 0) nks_max = nkbase + 1

       ! Step 1: Cache pool-local wavefunctions + becp in real space (read disk ONCE)
       ALLOCATE(psir_cache(dffts%nnr, npol, nbnd_kept, nks))
       ALLOCATE(becd_cache(nkb_d, npol, nbnd, nks))
       ALLOCATE(becpc_cache(nkb_p, npol, nbnd, nks))
       becd_cache = czero
       becpc_cache = czero

       DO ik_cache = 1, nks
          CALL read_collected_wfc(restart_dir(), ik_cache, evc_tmp)
          ! FFT to real space
          DO ib_cache = 1, nbnd_kept
             DO ipol_c = 1, npol
                psic_tmp = (0.0_dp, 0.0_dp)
                DO ig_cache = 1, ngk(ik_cache)
                   psic_tmp(dffts%nl(igk_k(ig_cache, ik_cache))) = &
                        evc_tmp(ig_cache + (ipol_c-1)*npwx, ibnd_min + ib_cache - 1)
                ENDDO
                CALL invfft('Wave', psic_tmp, dffts)
                psir_cache(:, ipol_c, ib_cache, ik_cache) = psic_tmp
             ENDDO
          ENDDO
          ! Compute and cache beta projections
          CALL get_betavkb(ngk(ik_cache), igk_k(1,ik_cache), xk(1,ik_cache), &
                            vkb_d, V_d%nat, V_d%ityp, V_d%tau, nkb_d)
          CALL allocate_bec_type(nkb_d, nbnd, becp_tmp_d)
          CALL calbec(ngk(ik_cache), vkb_d, evc_tmp, becp_tmp_d)
          IF (lspinorb) THEN
             becd_cache(:, :, :, ik_cache) = becp_tmp_d%nc(:, :, :)
          ELSE
             becd_cache(:, 1, :, ik_cache) = becp_tmp_d%k(:, :)
          ENDIF
          CALL deallocate_bec_type(becp_tmp_d)

          CALL get_betavkb(ngk(ik_cache), igk_k(1,ik_cache), xk(1,ik_cache), &
                            vkb_p, V_p%nat, V_p%ityp, V_p%tau, nkb_p)
          CALL allocate_bec_type(nkb_p, nbnd, becp_tmp_p)
          CALL calbec(ngk(ik_cache), vkb_p, evc_tmp, becp_tmp_p)
          IF (lspinorb) THEN
             becpc_cache(:, :, :, ik_cache) = becp_tmp_p%nc(:, :, :)
          ELSE
             becpc_cache(:, 1, :, ik_cache) = becp_tmp_p%k(:, :)
          ENDIF
          CALL deallocate_bec_type(becp_tmp_p)
       ENDDO

       IF (ionode) THEN
          WRITE(stdout, '(5X,A,I4,A)') &
               'Cached ', nks, ' pool-local k-points (psir + becp)'
          FLUSH(stdout)
       ENDIF

       ! Step 2: Get psi(ki) + becp(ki) — broadcast from owning pool
       ALLOCATE(becd_ki(nkb_d, npol, nbnd))
       ALLOCATE(becpc_ki(nkb_p, npol, nbnd))
       DO iki = 1, nki
          CALL get_pool_and_local(iki_map(iki), nkstot, npool, ik_pool, ik_l)
          psir_ki = (0.0_dp, 0.0_dp)
          becd_ki = czero
          becpc_ki = czero
          IF (my_pool_id == ik_pool) THEN
             psir_ki(:,:,:) = psir_cache(:,:,:,ik_l)
             becd_ki(:,:,:) = becd_cache(:,:,:,ik_l)
             becpc_ki(:,:,:) = becpc_cache(:,:,:,ik_l)
          ENDIF
          CALL MPI_Bcast(psir_ki, dffts%nnr * npol * nbnd_kept * 2, &
               MPI_DOUBLE_PRECISION, ik_pool, inter_pool_comm, ierr_mpi)
          CALL MPI_Bcast(becd_ki, nkb_d * npol * nbnd * 2, &
               MPI_DOUBLE_PRECISION, ik_pool, inter_pool_comm, ierr_mpi)
          CALL MPI_Bcast(becpc_ki, nkb_p * npol * nbnd * 2, &
               MPI_DOUBLE_PRECISION, ik_pool, inter_pool_comm, ierr_mpi)
       ENDDO

       ! Step 3: Panel broadcast for kf — iterate over source pools
       ALLOCATE(psir_recv(dffts%nnr, npol, nbnd_kept, nks_max))
       ALLOCATE(becd_recv(nkb_d, npol, nbnd, nks_max))
       ALLOCATE(becpc_recv(nkb_p, npol, nbnd, nks_max))

       iunit = 87
       IF (ionode) THEN
          fname = TRIM(prefix_in) // '_edmat_direct.dat'
          OPEN(iunit, FILE=TRIM(fname), FORM='formatted')
          WRITE(iunit, '(A)') '# Direct e-d matrix elements M(n,m, k_i, k_f)'
          WRITE(iunit, '(A,I6,A,I6)') '# nki=', nki, ' nkf=', nkf
          WRITE(iunit, '(A,I4,A,I4)') '# band range: ', ibnd_min, ' - ', ibnd_max
          WRITE(iunit, '(A)') '# iki ikf  kix kiy kiz  kfx kfy kfz  ibnd jbnd  |M|^2  Re(M)  Im(M)  |M_loc|^2  |M_nl|^2'
       ENDIF

       npairs_done = 0

       ! Distribute (iki,ikf) pairs across pools for parallel computation
       BLOCK
          INTEGER :: npairs_total, ipair, ipair_global, pair_iki, pair_ikf
          INTEGER :: my_pair_lo, my_pair_hi, npairs_local
          INTEGER :: pairs_per_pool, pairs_extra
          COMPLEX(dp), ALLOCATABLE :: edmatkq_all(:,:,:,:)
          COMPLEX(dp), ALLOCATABLE :: edmatkq_loc_all(:,:,:,:)
          COMPLEX(dp), ALLOCATABLE :: edmatkq_nl_all(:,:,:,:)

          ! Allocate full result arrays (zero-initialized)
          ALLOCATE(edmatkq_all(nbnd_kept, nbnd_kept, nki, nkf))
          ALLOCATE(edmatkq_loc_all(nbnd_kept, nbnd_kept, nki, nkf))
          ALLOCATE(edmatkq_nl_all(nbnd_kept, nbnd_kept, nki, nkf))
          edmatkq_all = czero
          edmatkq_loc_all = czero
          edmatkq_nl_all = czero

          npairs_total = nki * nkf
          pairs_per_pool = npairs_total / npool
          pairs_extra = MOD(npairs_total, npool)
          IF (my_pool_id < pairs_extra) THEN
             my_pair_lo = my_pool_id * (pairs_per_pool + 1) + 1
             npairs_local = pairs_per_pool + 1
          ELSE
             my_pair_lo = pairs_extra * (pairs_per_pool + 1) + &
                          (my_pool_id - pairs_extra) * pairs_per_pool + 1
             npairs_local = pairs_per_pool
          ENDIF
          my_pair_hi = my_pair_lo + npairs_local - 1

          IF (ionode) THEN
             WRITE(stdout, '(5X,A,I8,A,I4)') &
                  'Total pairs: ', npairs_total, '  npool: ', npool
             FLUSH(stdout)
          ENDIF

       DO ip = 0, npool - 1
          ! Determine k-point range of source pool ip
          IF (ip < nkrest) THEN
             src_lower = ip * (nkbase + 1) + 1
             src_nks = nkbase + 1
          ELSE
             src_lower = nkrest * (nkbase + 1) + (ip - nkrest) * nkbase + 1
             src_nks = nkbase
          ENDIF

          ! Source pool broadcasts its cached wavefunctions + becp
          IF (my_pool_id == ip) THEN
             psir_recv(:, :, :, 1:src_nks) = psir_cache(:, :, :, 1:nks)
             becd_recv(:, :, :, 1:src_nks) = becd_cache(:, :, :, 1:nks)
             becpc_recv(:, :, :, 1:src_nks) = becpc_cache(:, :, :, 1:nks)
          ENDIF
          CALL MPI_Bcast(psir_recv, dffts%nnr * npol * nbnd_kept * src_nks * 2, &
               MPI_DOUBLE_PRECISION, ip, inter_pool_comm, ierr_mpi)
          CALL MPI_Bcast(becd_recv, nkb_d * npol * nbnd * src_nks * 2, &
               MPI_DOUBLE_PRECISION, ip, inter_pool_comm, ierr_mpi)
          CALL MPI_Bcast(becpc_recv, nkb_p * npol * nbnd * src_nks * 2, &
               MPI_DOUBLE_PRECISION, ip, inter_pool_comm, ierr_mpi)

          ! Process pairs assigned to this pool whose kf is in source pool ip
          DO iki = 1, nki
             DO ikf = 1, nkf
                ! Check if this kf's NSCF index is in source pool ip's range
                IF (ikf_map(ikf) < src_lower .OR. &
                    ikf_map(ikf) >= src_lower + src_nks) CYCLE

                ! Check if this pair is assigned to my pool
                ipair_global = (iki - 1) * nkf + ikf
                IF (ipair_global < my_pair_lo .OR. ipair_global > my_pair_hi) CYCLE

                ikf_local = ikf_map(ikf) - src_lower + 1
                psir_kf(:,:,:) = psir_recv(:,:,:, ikf_local)

                ! Compute V_folded and matrix element
                qcryst = xkf(:,ikf) - xki(:,iki)
                V_folded = (0.0_dp, 0.0_dp)
                IF (ABS(qcryst(1)) < 1.0d-8 .AND. ABS(qcryst(2)) < 1.0d-8 &
                     .AND. ABS(qcryst(3)) < 1.0d-8) THEN
                   inr = 0
                   DO irz = 0, V_d%nr3 - 1
                      ir3mod = irz - (irz / dffts%nr3) * dffts%nr3
                      DO iry = 0, V_d%nr2 - 1
                         ir2mod = iry - (iry / dffts%nr2) * dffts%nr2
                         DO irx = 0, V_d%nr1 - 1
                            ir1mod = irx - (irx / dffts%nr1) * dffts%nr1
                            inr = inr + 1
                            irnmod = ir3mod * dffts%nr1 * dffts%nr2 + ir2mod * dffts%nr1 + ir1mod + 1
                            V_folded(irnmod) = V_folded(irnmod) + V_colin(inr)
                         ENDDO
                      ENDDO
                   ENDDO
                ELSE
                   d1 = tpi * qcryst(1) / dffts%nr1
                   d2 = tpi * qcryst(2) / dffts%nr2
                   d3 = tpi * qcryst(3) / dffts%nr3
                   inr = 0
                   DO irz = 0, V_d%nr3 - 1
                      ir3mod = irz - (irz / dffts%nr3) * dffts%nr3
                      DO iry = 0, V_d%nr2 - 1
                         ir2mod = iry - (iry / dffts%nr2) * dffts%nr2
                         DO irx = 0, V_d%nr1 - 1
                            ir1mod = irx - (irx / dffts%nr1) * dffts%nr1
                            inr = inr + 1
                            irnmod = ir3mod * dffts%nr1 * dffts%nr2 + ir2mod * dffts%nr1 + ir1mod + 1
                            arg = irx * d1 + iry * d2 + irz * d3
                            phase = CMPLX(COS(arg), SIN(arg), dp)
                            V_folded(irnmod) = V_folded(irnmod) + V_colin(inr) * phase
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF

                ! Matrix element: local + nonlocal
                DO ibnd = 1, nbnd_kept
                   DO jbnd = 1, nbnd_kept
                      ! LOCAL part
                      mlocal = czero
                      DO ipol_d = 1, npol
                         DO ir = 1, dffts%nnr
                            mlocal = mlocal + CONJG(psir_ki(ir, ipol_d, ibnd)) * &
                                 psir_kf(ir, ipol_d, jbnd) * V_folded(ir)
                         ENDDO
                      ENDDO
                      mlocal = mlocal / DBLE(dffts%nnr)

                      ! NONLOCAL defect
                      mnl_d = czero
                      ijkb0 = 0
                      DO nt = 1, V_d%ntyp
                         DO na = 1, V_d%nat
                            IF (V_d%ityp(na) == nt) THEN
                               IF (lspinorb) THEN
                                  DO ih = 1, nh(nt)
                                     ikb = ijkb0 + ih
                                     mnl_d = mnl_d + &
                                          CONJG(becd_ki(ikb,1,ibnd_min+ibnd-1)) * becd_recv(ikb,1,ibnd_min+jbnd-1,ikf_local) * dvan_so(ih,ih,1,nt) + &
                                          CONJG(becd_ki(ikb,1,ibnd_min+ibnd-1)) * becd_recv(ikb,2,ibnd_min+jbnd-1,ikf_local) * dvan_so(ih,ih,2,nt) + &
                                          CONJG(becd_ki(ikb,2,ibnd_min+ibnd-1)) * becd_recv(ikb,1,ibnd_min+jbnd-1,ikf_local) * dvan_so(ih,ih,3,nt) + &
                                          CONJG(becd_ki(ikb,2,ibnd_min+ibnd-1)) * becd_recv(ikb,2,ibnd_min+jbnd-1,ikf_local) * dvan_so(ih,ih,4,nt)
                                     DO jh = ih + 1, nh(nt)
                                        jkb = ijkb0 + jh
                                        mnl_d = mnl_d + &
                                             CONJG(becd_ki(ikb,1,ibnd_min+ibnd-1)) * becd_recv(jkb,1,ibnd_min+jbnd-1,ikf_local) * dvan_so(ih,jh,1,nt) + &
                                             CONJG(becd_ki(jkb,1,ibnd_min+ibnd-1)) * becd_recv(ikb,1,ibnd_min+jbnd-1,ikf_local) * dvan_so(jh,ih,1,nt) + &
                                             CONJG(becd_ki(ikb,1,ibnd_min+ibnd-1)) * becd_recv(jkb,2,ibnd_min+jbnd-1,ikf_local) * dvan_so(ih,jh,2,nt) + &
                                             CONJG(becd_ki(jkb,1,ibnd_min+ibnd-1)) * becd_recv(ikb,2,ibnd_min+jbnd-1,ikf_local) * dvan_so(jh,ih,2,nt) + &
                                             CONJG(becd_ki(ikb,2,ibnd_min+ibnd-1)) * becd_recv(jkb,1,ibnd_min+jbnd-1,ikf_local) * dvan_so(ih,jh,3,nt) + &
                                             CONJG(becd_ki(jkb,2,ibnd_min+ibnd-1)) * becd_recv(ikb,1,ibnd_min+jbnd-1,ikf_local) * dvan_so(jh,ih,3,nt) + &
                                             CONJG(becd_ki(ikb,2,ibnd_min+ibnd-1)) * becd_recv(jkb,2,ibnd_min+jbnd-1,ikf_local) * dvan_so(ih,jh,4,nt) + &
                                             CONJG(becd_ki(jkb,2,ibnd_min+ibnd-1)) * becd_recv(ikb,2,ibnd_min+jbnd-1,ikf_local) * dvan_so(jh,ih,4,nt)
                                     ENDDO
                                  ENDDO
                               ELSE
                                  DO ih = 1, nh(nt)
                                     ikb = ijkb0 + ih
                                     mnl_d = mnl_d + CONJG(becd_ki(ikb,1,ibnd_min+ibnd-1)) * &
                                          becd_recv(ikb,1,ibnd_min+jbnd-1,ikf_local) * dvan(ih, ih, nt)
                                     DO jh = ih + 1, nh(nt)
                                        jkb = ijkb0 + jh
                                        mnl_d = mnl_d + &
                                             (CONJG(becd_ki(ikb,1,ibnd_min+ibnd-1)) * becd_recv(jkb,1,ibnd_min+jbnd-1,ikf_local) + &
                                              CONJG(becd_ki(jkb,1,ibnd_min+ibnd-1)) * becd_recv(ikb,1,ibnd_min+jbnd-1,ikf_local)) * &
                                             dvan(ih, jh, nt)
                                     ENDDO
                                  ENDDO
                               ENDIF
                               ijkb0 = ijkb0 + nh(nt)
                            ENDIF
                         ENDDO
                      ENDDO

                      ! NONLOCAL pristine
                      mnl_p = czero
                      ijkb0 = 0
                      DO nt = 1, V_p%ntyp
                         DO na = 1, V_p%nat
                            IF (V_p%ityp(na) == nt) THEN
                               IF (lspinorb) THEN
                                  DO ih = 1, nh(nt)
                                     ikb = ijkb0 + ih
                                     mnl_p = mnl_p + &
                                          CONJG(becpc_ki(ikb,1,ibnd_min+ibnd-1)) * becpc_recv(ikb,1,ibnd_min+jbnd-1,ikf_local) * dvan_so(ih,ih,1,nt) + &
                                          CONJG(becpc_ki(ikb,1,ibnd_min+ibnd-1)) * becpc_recv(ikb,2,ibnd_min+jbnd-1,ikf_local) * dvan_so(ih,ih,2,nt) + &
                                          CONJG(becpc_ki(ikb,2,ibnd_min+ibnd-1)) * becpc_recv(ikb,1,ibnd_min+jbnd-1,ikf_local) * dvan_so(ih,ih,3,nt) + &
                                          CONJG(becpc_ki(ikb,2,ibnd_min+ibnd-1)) * becpc_recv(ikb,2,ibnd_min+jbnd-1,ikf_local) * dvan_so(ih,ih,4,nt)
                                     DO jh = ih + 1, nh(nt)
                                        jkb = ijkb0 + jh
                                        mnl_p = mnl_p + &
                                             CONJG(becpc_ki(ikb,1,ibnd_min+ibnd-1)) * becpc_recv(jkb,1,ibnd_min+jbnd-1,ikf_local) * dvan_so(ih,jh,1,nt) + &
                                             CONJG(becpc_ki(jkb,1,ibnd_min+ibnd-1)) * becpc_recv(ikb,1,ibnd_min+jbnd-1,ikf_local) * dvan_so(jh,ih,1,nt) + &
                                             CONJG(becpc_ki(ikb,1,ibnd_min+ibnd-1)) * becpc_recv(jkb,2,ibnd_min+jbnd-1,ikf_local) * dvan_so(ih,jh,2,nt) + &
                                             CONJG(becpc_ki(jkb,1,ibnd_min+ibnd-1)) * becpc_recv(ikb,2,ibnd_min+jbnd-1,ikf_local) * dvan_so(jh,ih,2,nt) + &
                                             CONJG(becpc_ki(ikb,2,ibnd_min+ibnd-1)) * becpc_recv(jkb,1,ibnd_min+jbnd-1,ikf_local) * dvan_so(ih,jh,3,nt) + &
                                             CONJG(becpc_ki(jkb,2,ibnd_min+ibnd-1)) * becpc_recv(ikb,1,ibnd_min+jbnd-1,ikf_local) * dvan_so(jh,ih,3,nt) + &
                                             CONJG(becpc_ki(ikb,2,ibnd_min+ibnd-1)) * becpc_recv(jkb,2,ibnd_min+jbnd-1,ikf_local) * dvan_so(ih,jh,4,nt) + &
                                             CONJG(becpc_ki(jkb,2,ibnd_min+ibnd-1)) * becpc_recv(ikb,2,ibnd_min+jbnd-1,ikf_local) * dvan_so(jh,ih,4,nt)
                                     ENDDO
                                  ENDDO
                               ELSE
                                  DO ih = 1, nh(nt)
                                     ikb = ijkb0 + ih
                                     mnl_p = mnl_p + CONJG(becpc_ki(ikb,1,ibnd_min+ibnd-1)) * &
                                          becpc_recv(ikb,1,ibnd_min+jbnd-1,ikf_local) * dvan(ih, ih, nt)
                                     DO jh = ih + 1, nh(nt)
                                        jkb = ijkb0 + jh
                                        mnl_p = mnl_p + &
                                             (CONJG(becpc_ki(ikb,1,ibnd_min+ibnd-1)) * becpc_recv(jkb,1,ibnd_min+jbnd-1,ikf_local) + &
                                              CONJG(becpc_ki(jkb,1,ibnd_min+ibnd-1)) * becpc_recv(ikb,1,ibnd_min+jbnd-1,ikf_local)) * &
                                             dvan(ih, jh, nt)
                                     ENDDO
                                  ENDDO
                               ENDIF
                               ijkb0 = ijkb0 + nh(nt)
                            ENDIF
                         ENDDO
                      ENDDO

                      edmatkq(ibnd, jbnd) = mlocal + mnl_d - mnl_p
                      edmatkq_loc(ibnd, jbnd) = mlocal
                      edmatkq_nl(ibnd, jbnd) = mnl_d - mnl_p
                      edmatkq_all(ibnd, jbnd, iki, ikf) = mlocal + mnl_d - mnl_p
                      edmatkq_loc_all(ibnd, jbnd, iki, ikf) = mlocal
                      edmatkq_nl_all(ibnd, jbnd, iki, ikf) = mnl_d - mnl_p
                   ENDDO
                ENDDO

                npairs_done = npairs_done + 1
             ENDDO  ! ikf
          ENDDO  ! iki

          IF (ionode) THEN
             WRITE(stdout, '(5X,A,I4,A,I4,A,I6)') &
                  'kf pool ', ip, ' / ', npool - 1, &
                  '  my pairs so far: ', npairs_done
             FLUSH(stdout)
          ENDIF
       ENDDO  ! ip (source pool)

       ! Gather results: use mp_sum on a full (nki,nkf,nbnd,nbnd) array
       ! Each pool computed a non-overlapping subset of pairs, others are zero
       CALL mp_sum(edmatkq_all, inter_pool_comm)
       CALL mp_sum(edmatkq_loc_all, inter_pool_comm)
       CALL mp_sum(edmatkq_nl_all, inter_pool_comm)

       ! ionode writes ordered output
       IF (ionode) THEN
          DO iki = 1, nki
             DO ikf = 1, nkf
                DO ibnd = 1, nbnd_kept
                   DO jbnd = 1, nbnd_kept
                      WRITE(iunit, '(2I6,6F10.5,2I4,5ES16.8)') &
                           iki, ikf, xki(:,iki), xkf(:,ikf), &
                           ibnd_min + ibnd - 1, ibnd_min + jbnd - 1, &
                           ABS(edmatkq_all(ibnd, jbnd, iki, ikf))**2, &
                           REAL(edmatkq_all(ibnd, jbnd, iki, ikf)), &
                           AIMAG(edmatkq_all(ibnd, jbnd, iki, ikf)), &
                           ABS(edmatkq_loc_all(ibnd, jbnd, iki, ikf))**2, &
                           ABS(edmatkq_nl_all(ibnd, jbnd, iki, ikf))**2
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
          CLOSE(iunit)
          fname = TRIM(prefix_in) // '_edmat_direct.dat'
          WRITE(stdout, '(5X,A,I6,A,A)') 'Total pairs: ', nki * nkf, &
               '  Written to ', TRIM(fname)
          WRITE(stdout, '(5X,A)') REPEAT('=', 60)
       ENDIF
       DEALLOCATE(edmatkq_all, edmatkq_loc_all, edmatkq_nl_all)
       END BLOCK  ! pair distribution block

       DEALLOCATE(psir_cache, psir_recv)
       DEALLOCATE(becd_cache, becd_recv, becpc_cache, becpc_recv)
       DEALLOCATE(becd_ki, becpc_ki)
    END BLOCK

    DEALLOCATE(xki, xkf, xk_cryst_nscf, iki_map, ikf_map)
    DEALLOCATE(evc_tmp, psic_tmp, psir_ki, psir_kf, V_folded, edmatkq, edmatkq_loc, edmatkq_nl, vkb_d, vkb_p)

  CONTAINS
    SUBROUTINE get_pool_and_local(ik_g, nktot, np, ipool, ik_l)
      ! Map global k-index to owning pool and local index
      INTEGER, INTENT(IN) :: ik_g, nktot, np
      INTEGER, INTENT(OUT) :: ipool, ik_l
      INTEGER :: nkb, nkr, off
      nkb = nktot / np
      nkr = MOD(nktot, np)
      IF (ik_g <= (nkb + 1) * nkr) THEN
         ipool = (ik_g - 1) / (nkb + 1)
         ik_l = ik_g - ipool * (nkb + 1)
      ELSE
         ipool = nkr + (ik_g - (nkb + 1) * nkr - 1) / nkb
         off = nkr * (nkb + 1) + (ipool - nkr) * nkb
         ik_l = ik_g - off
      ENDIF
    END SUBROUTINE
  END SUBROUTINE ed_direct_from_files


  SUBROUTINE read_kpoints_file(fname, nk, xk_cryst)
    USE io_global, ONLY : ionode, ionode_id, stdout
    USE mp, ONLY : mp_bcast
    USE mp_world, ONLY : world_comm
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: fname
    INTEGER, INTENT(OUT) :: nk
    REAL(dp), ALLOCATABLE, INTENT(OUT) :: xk_cryst(:,:)
    INTEGER :: ik, iunit, ios

    iunit = 84
    nk = 0
    IF (ionode) THEN
       OPEN(iunit, FILE=TRIM(fname), STATUS='old', FORM='formatted', IOSTAT=ios)
       IF (ios /= 0) CALL errore('read_kpoints_file', 'Cannot open ' // TRIM(fname), 1)
       READ(iunit, *) nk
    ENDIF
    CALL mp_bcast(nk, ionode_id, world_comm)
    ALLOCATE(xk_cryst(3, nk))
    IF (ionode) THEN
       DO ik = 1, nk
          READ(iunit, *) xk_cryst(:, ik)
       ENDDO
       CLOSE(iunit)
       WRITE(stdout, '(5X,A,A,A,I6,A)') 'Read ', TRIM(fname), ': ', nk, ' k-points'
    ENDIF
    CALL mp_bcast(xk_cryst, ionode_id, world_comm)
  END SUBROUTINE read_kpoints_file


  SUBROUTINE ed_interp_from_file(nbndsub, nrr, irvec, ndegen, chw, edmatw_2d, &
                                  filki, filkf, prefix_in)
    !-----------------------------------------------------------------------
    ! Wannier interpolation from k-point files using full double-FT
    ! M(R,Rp) via edmatwan2bloch_2d.
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode, stdout
    USE wan2bloch_edi, ONLY : get_cfac, hamwan2bloch_with_evec, edmatwan2bloch_2d
    USE ep_constants, ONLY : czero, cone
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nbndsub, nrr
    INTEGER, INTENT(IN) :: irvec(3, nrr), ndegen(nrr)
    COMPLEX(dp), INTENT(IN) :: chw(nbndsub, nbndsub, nrr)
    COMPLEX(dp), INTENT(IN) :: edmatw_2d(nbndsub, nbndsub, nrr, nrr)
    CHARACTER(LEN=*), INTENT(IN) :: filki, filkf, prefix_in

    INTEGER :: iki, ikf, nki, nkf, ibnd, jbnd, iunit
    REAL(dp), ALLOCATABLE :: xki(:,:), xkf(:,:)
    REAL(dp) :: eig_i(nbndsub), eig_f(nbndsub)
    COMPLEX(dp) :: evec_i(nbndsub, nbndsub), evec_f(nbndsub, nbndsub)
    COMPLEX(dp) :: cfac_ki(nrr), cfac_kf(nrr)
    COMPLEX(dp) :: edmatf_w(nbndsub, nbndsub)
    COMPLEX(dp) :: edmatf_b(nbndsub, nbndsub), tmp(nbndsub, nbndsub)
    CHARACTER(LEN=256) :: fname

    CALL read_kpoints_file(filki, nki, xki)
    CALL read_kpoints_file(filkf, nkf, xkf)

    IF (ionode) THEN
       WRITE(stdout, '(/,5X,A)') REPEAT('=', 60)
       WRITE(stdout, '(5X,A)') 'Wannier interpolation from k-point files [double-FT]'
       WRITE(stdout, '(5X,A)') '  Using M(R,Rp) with edmatwan2bloch_2d'
       WRITE(stdout, '(5X,A)') REPEAT('=', 60)
       WRITE(stdout, '(5X,A,I6,A,I6,A,I10)') &
            'nki=', nki, ' nkf=', nkf, ' pairs=', nki*nkf

       iunit = 83
       fname = TRIM(prefix_in) // '_edmat_interp.dat'
       OPEN(iunit, FILE=TRIM(fname), FORM='formatted')
       WRITE(iunit, '(A)') '# Wannier-interpolated M(n,m, k_i, k_f) [double-FT]'
       WRITE(iunit, '(A,I6,A,I6)') '# nki=', nki, ' nkf=', nkf
       WRITE(iunit, '(A)') '# iki ikf  kix kiy kiz  kfx kfy kfz  ibnd jbnd  |M|^2  Re(M)  Im(M)'
    ENDIF

    DO iki = 1, nki
       ! Diag H(ki) -> U(ki)
       CALL get_cfac(nrr, irvec, xki(:,iki), cfac_ki)
       CALL hamwan2bloch_with_evec(nbndsub, nrr, ndegen, eig_i, evec_i, chw, cfac_ki)

       DO ikf = 1, nkf
          ! Diag H(kf) -> U(kf)
          CALL get_cfac(nrr, irvec, xkf(:,ikf), cfac_kf)
          CALL hamwan2bloch_with_evec(nbndsub, nrr, ndegen, eig_f, evec_f, chw, cfac_kf)

          ! Double-FT inverse: M_W(ki,kf) from M(Re, Rp)
          CALL edmatwan2bloch_2d(nbndsub, nrr, ndegen, edmatw_2d, cfac_ki, cfac_kf, edmatf_w)

          ! Rotate to Bloch: M_B = U_dag(ki) M_W U(kf)
          CALL ZGEMM('N', 'N', nbndsub, nbndsub, nbndsub, &
                      cone, edmatf_w, nbndsub, evec_f, nbndsub, czero, tmp, nbndsub)
          CALL ZGEMM('C', 'N', nbndsub, nbndsub, nbndsub, &
                      cone, evec_i, nbndsub, tmp, nbndsub, czero, edmatf_b, nbndsub)

          IF (ionode) THEN
             DO ibnd = 1, nbndsub
                DO jbnd = 1, nbndsub
                   WRITE(iunit, '(2I6,6F10.5,2I4,3ES16.8)') &
                        iki, ikf, xki(:,iki), xkf(:,ikf), ibnd, jbnd, &
                        ABS(edmatf_b(ibnd, jbnd))**2, &
                        REAL(edmatf_b(ibnd, jbnd)), AIMAG(edmatf_b(ibnd, jbnd))
                ENDDO
             ENDDO
          ENDIF
       ENDDO

       IF (ionode .AND. (MOD(iki, MAX(1,nki/10)) == 0 .OR. iki == nki)) THEN
          WRITE(stdout, '(5X,A,I6,A,I6)') '  k_i ', iki, ' / ', nki
          FLUSH(stdout)
       ENDIF
    ENDDO

    IF (ionode) THEN
       CLOSE(iunit)
       WRITE(stdout, '(5X,A,A)') 'Written to ', TRIM(fname)
       WRITE(stdout, '(5X,A)') REPEAT('=', 60)
    ENDIF

    DEALLOCATE(xki, xkf)
  END SUBROUTINE ed_interp_from_file


  SUBROUTINE build_vcolin_aligned(vf_d, vf_p, vcolin, nrtot, shift_out)
    !-----------------------------------------------------------------------
    ! Build V_colin = V_d - V_p with vacuum alignment for 2D systems.
    !
    ! For a slab geometry with vacuum along z:
    !   1. Compute planar-averaged V_d(z) and V_p(z)
    !   2. Find the z-plane furthest from any atom (vacuum)
    !   3. shift = V_p(z_vac) - V_d(z_vac)
    !   4. V_colin = V_d - V_p + shift
    !
    ! Uses pristine cell atoms to locate the vacuum (more atoms → cleaner).
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode, stdout
    USE edic_mod, ONLY : V_file
    IMPLICIT NONE

    TYPE(V_file), INTENT(IN) :: vf_d, vf_p
    INTEGER, INTENT(IN) :: nrtot
    REAL(dp), INTENT(OUT) :: vcolin(nrtot)
    REAL(dp), INTENT(OUT) :: shift_out

    INTEGER :: irz, ia, nr1, nr2, nr3, nxy, iz_vac
    REAL(dp) :: z_frac_grid, z_frac_atom, d, dmin, max_dmin
    REAL(dp), ALLOCATABLE :: vavg_d(:), vavg_p(:)

    nr1 = vf_d%nr1
    nr2 = vf_d%nr2
    nr3 = vf_d%nr3
    nxy = nr1 * nr2

    ! Planar-averaged potential along z for both supercells
    ALLOCATE(vavg_d(0:nr3-1), vavg_p(0:nr3-1))
    DO irz = 0, nr3 - 1
       vavg_d(irz) = SUM(vf_d%pot(irz*nxy+1 : (irz+1)*nxy)) / DBLE(nxy)
       vavg_p(irz) = SUM(vf_p%pot(irz*nxy+1 : (irz+1)*nxy)) / DBLE(nxy)
    ENDDO

    ! Find z-plane furthest from any atom (pristine cell)
    iz_vac = 0
    max_dmin = -1.0_dp
    DO irz = 0, nr3 - 1
       z_frac_grid = DBLE(irz) / DBLE(nr3)
       dmin = 1.0_dp
       DO ia = 1, vf_p%nat
          z_frac_atom = vf_p%tau(3, ia) / vf_p%at(3, 3)
          d = ABS(z_frac_grid - z_frac_atom)
          d = MIN(d, 1.0_dp - d)
          IF (d < dmin) dmin = d
       ENDDO
       IF (dmin > max_dmin) THEN
          max_dmin = dmin
          iz_vac = irz
       ENDIF
    ENDDO

    shift_out = vavg_p(iz_vac) - vavg_d(iz_vac)
    vcolin(:) = vf_d%pot(1:nrtot) - vf_p%pot(1:nrtot) + shift_out

    IF (ionode) THEN
       WRITE(stdout, '(5X,A)')          '--- Vacuum alignment (2D) ---'
       WRITE(stdout, '(5X,A,I6,A,F8.4)') &
            '  Vacuum z-index: ', iz_vac, '  z_frac = ', DBLE(iz_vac)/DBLE(nr3)
       WRITE(stdout, '(5X,A,F14.8,A)')  '  V_d(vac) planar avg = ', vavg_d(iz_vac), ' Ry'
       WRITE(stdout, '(5X,A,F14.8,A)')  '  V_p(vac) planar avg = ', vavg_p(iz_vac), ' Ry'
       WRITE(stdout, '(5X,A,F14.8,A)')  '  Alignment shift     = ', shift_out, ' Ry'
       WRITE(stdout, '(5X,A)')          '-----------------------------'
       FLUSH(stdout)
    ENDIF

    DEALLOCATE(vavg_d, vavg_p)

  END SUBROUTINE build_vcolin_aligned


  SUBROUTINE build_vcolin_corealign(vf_d, vf_p, vcolin, nrtot, shift_out, &
                                     defect_center, r_avg)
    !-----------------------------------------------------------------------
    ! Build V_colin = V_d - V_p with core-site alignment (2D & 3D).
    !
    ! 1. Find the atom in the pristine cell furthest from defect_center
    !    (minimum-image distance in the periodic supercell).
    ! 2. Average V_d and V_p over grid points within a sphere of radius
    !    r_avg (Bohr) around that atom.
    ! 3. shift = V_p_avg - V_d_avg
    ! 4. V_colin = V_d - V_p + shift
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode, stdout
    USE edic_mod, ONLY : V_file
    IMPLICIT NONE

    TYPE(V_file), INTENT(IN) :: vf_d, vf_p
    INTEGER, INTENT(IN) :: nrtot
    REAL(dp), INTENT(OUT) :: vcolin(nrtot)
    REAL(dp), INTENT(OUT) :: shift_out
    REAL(dp), INTENT(IN) :: defect_center(3)
    REAL(dp), INTENT(IN) :: r_avg

    INTEGER :: ia, ia_far, irx, iry, irz, inr, nsph
    INTEGER :: nr1, nr2, nr3
    REAL(dp) :: r_def_cart(3), dr(3), s(3), dist, max_dist
    REAL(dp) :: r_far(3), r_grid(3), dr_grid(3), s_grid(3), dist_grid
    REAL(dp) :: r_avg_alat, vd_sum, vp_sum
    REAL(dp) :: at(3,3), bg(3,3), cross(3), det

    at = vf_p%at
    nr1 = vf_d%nr1
    nr2 = vf_d%nr2
    nr3 = vf_d%nr3

    ! Compute reciprocal lattice vectors bg (without 2pi, in 1/alat)
    ! bg(:,i) satisfies sum_j at(j,i)*bg(j,k) = delta(i,k)
    cross(1) = at(2,2)*at(3,3) - at(3,2)*at(2,3)
    cross(2) = at(3,2)*at(1,3) - at(1,2)*at(3,3)
    cross(3) = at(1,2)*at(2,3) - at(2,2)*at(1,3)
    det = at(1,1)*cross(1) + at(2,1)*cross(2) + at(3,1)*cross(3)

    bg(1,1) = cross(1) / det
    bg(2,1) = cross(2) / det
    bg(3,1) = cross(3) / det

    cross(1) = at(2,3)*at(3,1) - at(3,3)*at(2,1)
    cross(2) = at(3,3)*at(1,1) - at(1,3)*at(3,1)
    cross(3) = at(1,3)*at(2,1) - at(2,3)*at(1,1)
    bg(1,2) = cross(1) / det
    bg(2,2) = cross(2) / det
    bg(3,2) = cross(3) / det

    cross(1) = at(2,1)*at(3,2) - at(3,1)*at(2,2)
    cross(2) = at(3,1)*at(1,2) - at(1,1)*at(3,2)
    cross(3) = at(1,1)*at(2,2) - at(2,1)*at(1,2)
    bg(1,3) = cross(1) / det
    bg(2,3) = cross(2) / det
    bg(3,3) = cross(3) / det

    ! Convert defect center: fractional -> Cartesian (alat)
    r_def_cart(:) = defect_center(1) * at(:,1) &
                  + defect_center(2) * at(:,2) &
                  + defect_center(3) * at(:,3)

    ! Find pristine atom furthest from defect center (minimum image)
    ia_far = 1
    max_dist = -1.0_dp
    DO ia = 1, vf_p%nat
       dr(:) = vf_p%tau(:, ia) - r_def_cart(:)
       ! Minimum image: convert to fractional, wrap, convert back
       s(1) = bg(1,1)*dr(1) + bg(2,1)*dr(2) + bg(3,1)*dr(3)
       s(2) = bg(1,2)*dr(1) + bg(2,2)*dr(2) + bg(3,2)*dr(3)
       s(3) = bg(1,3)*dr(1) + bg(2,3)*dr(2) + bg(3,3)*dr(3)
       s(:) = s(:) - NINT(s(:))
       dr(:) = s(1)*at(:,1) + s(2)*at(:,2) + s(3)*at(:,3)
       dist = SQRT(SUM(dr**2)) * vf_p%alat
       IF (dist > max_dist) THEN
          max_dist = dist
          ia_far = ia
          r_far(:) = vf_p%tau(:, ia)
       ENDIF
    ENDDO

    ! Average V_d and V_p over grid points within r_avg of r_far
    r_avg_alat = r_avg / vf_p%alat
    vd_sum = 0.0_dp
    vp_sum = 0.0_dp
    nsph = 0

    DO irz = 0, nr3 - 1
       DO iry = 0, nr2 - 1
          DO irx = 0, nr1 - 1
             inr = irz * nr1 * nr2 + iry * nr1 + irx + 1
             r_grid(:) = (DBLE(irx)/DBLE(nr1)) * at(:,1) &
                        + (DBLE(iry)/DBLE(nr2)) * at(:,2) &
                        + (DBLE(irz)/DBLE(nr3)) * at(:,3)
             dr_grid(:) = r_grid(:) - r_far(:)
             ! Minimum image
             s_grid(1) = bg(1,1)*dr_grid(1) + bg(2,1)*dr_grid(2) + bg(3,1)*dr_grid(3)
             s_grid(2) = bg(1,2)*dr_grid(1) + bg(2,2)*dr_grid(2) + bg(3,2)*dr_grid(3)
             s_grid(3) = bg(1,3)*dr_grid(1) + bg(2,3)*dr_grid(2) + bg(3,3)*dr_grid(3)
             s_grid(:) = s_grid(:) - NINT(s_grid(:))
             dr_grid(:) = s_grid(1)*at(:,1) + s_grid(2)*at(:,2) + s_grid(3)*at(:,3)
             dist_grid = SQRT(SUM(dr_grid**2))
             IF (dist_grid <= r_avg_alat) THEN
                vd_sum = vd_sum + vf_d%pot(inr)
                vp_sum = vp_sum + vf_p%pot(inr)
                nsph = nsph + 1
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    IF (nsph == 0) THEN
       IF (ionode) WRITE(stdout, '(5X,A)') &
            'WARNING: core_align_radius too small, no grid points found. Using shift=0.'
       shift_out = 0.0_dp
    ELSE
       shift_out = (vp_sum - vd_sum) / DBLE(nsph)
    ENDIF

    vcolin(:) = vf_d%pot(1:nrtot) - vf_p%pot(1:nrtot) + shift_out

    IF (ionode) THEN
       WRITE(stdout, '(5X,A)')          '--- Core-site alignment ---'
       WRITE(stdout, '(5X,A,I6,A,F10.4,A)') &
            '  Reference atom: ', ia_far, &
            '  dist from defect = ', max_dist, ' Bohr'
       WRITE(stdout, '(5X,A,3F10.5)')   '  Atom position (alat): ', r_far
       WRITE(stdout, '(5X,A,F8.3,A,I8)') &
            '  Averaging radius = ', r_avg, ' Bohr,  grid points = ', nsph
       WRITE(stdout, '(5X,A,F14.8,A)')  '  V_d(core) avg = ', vd_sum/DBLE(MAX(nsph,1)), ' Ry'
       WRITE(stdout, '(5X,A,F14.8,A)')  '  V_p(core) avg = ', vp_sum/DBLE(MAX(nsph,1)), ' Ry'
       WRITE(stdout, '(5X,A,F14.8,A)')  '  Alignment shift = ', shift_out, ' Ry'
       WRITE(stdout, '(5X,A)')          '---------------------------'
       FLUSH(stdout)
    ENDIF

  END SUBROUTINE build_vcolin_corealign


  SUBROUTINE write_vcolin_cube(vf_ref, vcolin, nrtot, prefix_in)
    !-----------------------------------------------------------------------
    ! Write V_colin to a Gaussian cube file for visualization (VESTA, VMD).
    ! Uses vf_ref (defect cell) for lattice vectors and atomic positions.
    ! Coordinates are in Bohr; potential values are in Ry.
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode, stdout
    USE edic_mod, ONLY : V_file
    IMPLICIT NONE

    TYPE(V_file), INTENT(IN) :: vf_ref
    INTEGER, INTENT(IN) :: nrtot
    REAL(dp), INTENT(IN) :: vcolin(nrtot)
    CHARACTER(LEN=*), INTENT(IN) :: prefix_in

    INTEGER :: iunit, irx, iry, irz, ia, inr, nval
    INTEGER :: nr1, nr2, nr3
    REAL(dp) :: voxel(3)
    CHARACTER(LEN=256) :: fname

    IF (.NOT. ionode) RETURN

    nr1 = vf_ref%nr1
    nr2 = vf_ref%nr2
    nr3 = vf_ref%nr3

    fname = TRIM(prefix_in) // '_vcolin.cube'
    iunit = 98
    OPEN(iunit, FILE=TRIM(fname), STATUS='REPLACE', ACTION='WRITE')

    WRITE(iunit, '(A)') 'V_colin = V_defect - V_pristine (vacuum-aligned) [Ry]'
    WRITE(iunit, '(A)') 'Generated by EDI code'

    WRITE(iunit, '(I5, 3F12.6)') vf_ref%nat, 0.0_dp, 0.0_dp, 0.0_dp

    voxel(:) = vf_ref%at(:, 1) * vf_ref%alat / DBLE(nr1)
    WRITE(iunit, '(I5, 3F12.6)') nr1, voxel
    voxel(:) = vf_ref%at(:, 2) * vf_ref%alat / DBLE(nr2)
    WRITE(iunit, '(I5, 3F12.6)') nr2, voxel
    voxel(:) = vf_ref%at(:, 3) * vf_ref%alat / DBLE(nr3)
    WRITE(iunit, '(I5, 3F12.6)') nr3, voxel

    DO ia = 1, vf_ref%nat
       WRITE(iunit, '(I5, 4F12.6)') vf_ref%ityp(ia), 0.0_dp, &
            vf_ref%tau(1, ia) * vf_ref%alat, &
            vf_ref%tau(2, ia) * vf_ref%alat, &
            vf_ref%tau(3, ia) * vf_ref%alat
    ENDDO

    ! Cube format: outermost = axis 1, innermost = axis 3
    DO irx = 0, nr1 - 1
       DO iry = 0, nr2 - 1
          nval = 0
          DO irz = 0, nr3 - 1
             inr = irz * nr1 * nr2 + iry * nr1 + irx + 1
             WRITE(iunit, '(ES13.5)', ADVANCE='NO') vcolin(inr)
             nval = nval + 1
             IF (MOD(nval, 6) == 0) WRITE(iunit, *)
          ENDDO
          IF (MOD(nval, 6) /= 0) WRITE(iunit, *)
       ENDDO
    ENDDO

    CLOSE(iunit)
    WRITE(stdout, '(5X,A,A)') 'V_colin cube file written: ', TRIM(fname)
    FLUSH(stdout)

  END SUBROUTINE write_vcolin_cube


END MODULE ed_coarse
