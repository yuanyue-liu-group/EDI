MODULE transport_edi
  !-----------------------------------------------------------------------
  ! SERTA/MRTA transport for electron-defect scattering.
  ! Full IBZ approach with symmetry (following EPW's compute_sigma_sym).
  !
  ! Strategy:
  !   1. Build IBZ from full BZ using crystal point group
  !   2. Compute scattering rates at IBZ ki only (nsym× speedup)
  !   3. Unfold inv_tau to full BZ (scalar, invariant under symmetry)
  !   4. Compute conductivity with symmetry-rotated velocities
  !   5. Enforces σ_xx = σ_yy for hexagonal systems automatically
  !-----------------------------------------------------------------------
  USE kinds, ONLY : dp
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: compute_transport, compute_dos_validation

  REAL(dp), PARAMETER :: twopi = 6.283185307179586_dp
  REAL(dp), PARAMETER :: kb_eV = 8.617333262e-5_dp
  REAL(dp), PARAMETER :: ryd2ev = 13.6056980659_dp
  REAL(dp), PARAMETER :: bohr2ang = 0.52917720859_dp
  REAL(dp), PARAMETER :: ang2cm = 1.0e-8_dp
  REAL(dp), PARAMETER :: electron_SI = 1.602176634e-19_dp
  REAL(dp), PARAMETER :: hbarJ = 1.054571800e-34_dp

CONTAINS

  SUBROUTINE compute_transport(nbndsub, nrr, irvec, ndegen, chw, edmatw_2d, &
                                nk1f, nk2f, nk3f, nstemp, temps, &
                                defect_conc, prefix_in)
    USE io_global, ONLY : ionode, stdout
    USE mp, ONLY : mp_sum
    USE mp_global, ONLY : inter_pool_comm
    USE mp_pools, ONLY : npool, my_pool_id
    USE mp_world, ONLY : nproc
    USE cell_base, ONLY : at, alat, bg, tpiba, omega
    USE symm_base, ONLY : s, nsym
    USE noncollin_module, ONLY : noncolin
    USE wan2bloch_edi, ONLY : get_cfac, hamwan2bloch_with_evec, edmatwan2bloch_2d
    USE ep_constants, ONLY : czero, cone
    USE delta_weights, ONLY : get_delta_weight, kindex_sub_2d
    USE edi_input, ONLY : delta_method, delta_sigma, &
                           transport_win_min, transport_win_max, carrier_conc
    USE parallel_include
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nbndsub, nrr, nstemp
    INTEGER, INTENT(IN) :: irvec(3, nrr), ndegen(nrr)
    INTEGER, INTENT(IN) :: nk1f, nk2f, nk3f
    COMPLEX(dp), INTENT(IN) :: chw(nbndsub, nbndsub, nrr)
    COMPLEX(dp), INTENT(IN) :: edmatw_2d(nbndsub, nbndsub, nrr, nrr)
    REAL(dp), INTENT(IN) :: temps(nstemp), defect_conc
    CHARACTER(LEN=*), INTENT(IN) :: prefix_in
    LOGICAL :: use_gaussian, use_adaptive
    REAL(dp) :: sfac  ! spin degeneracy: 2 for collinear, 1 for noncolin/SOC

    INTEGER :: nktotf, iki, ikf, ibnd, jbnd, itemp, i1, i2, i3, ndim
    INTEGER :: iunit, npairs_computed
    INTEGER :: ki_lower, ki_upper, nki_local
    REAL(dp) :: temp, ef_eV, ef_Ry, beta_eV
    REAL(dp) :: win_min_Ry, win_max_Ry
    COMPLEX(dp) :: evec_i(nbndsub, nbndsub), evec_f(nbndsub, nbndsub)
    COMPLEX(dp) :: edmatf_w(nbndsub, nbndsub), edmatf_b(nbndsub, nbndsub)
    COMPLEX(dp) :: tmp_mat(nbndsub, nbndsub)

    ! Pre-cached arrays
    REAL(dp), ALLOCATABLE :: eig_all_eV(:,:), eig_all_Ry(:,:)
    REAL(dp), ALLOCATABLE :: vel(:,:,:)
    COMPLEX(dp), ALLOCATABLE :: cfac_all(:,:)
    COMPLEX(dp), ALLOCATABLE :: evec_all(:,:,:)
    LOGICAL, ALLOCATABLE :: kf_active(:)

    ! Scattering rates
    REAL(dp), ALLOCATABLE :: inv_tau_serta(:,:), inv_tau_mrta(:,:)

    ! Adaptive smearing: pre-computed sigma per (ibnd, ik) in Ry
    REAL(dp), ALLOCATABLE :: sigma_adapt(:,:)

    ! Active ki/kf lists and pair distribution
    INTEGER, ALLOCATABLE :: active_ibz_list(:), active_kf_list(:)
    INTEGER :: n_active_ki, n_active_kf
    INTEGER :: my_pair_lower, my_pair_upper, my_npairs_local

    ! Delta function workspace
    REAL(dp), ALLOCATABLE :: bande_Ry(:,:), freq_zero(:,:)

    ! IBZ mapping
    INTEGER, ALLOCATABLE :: bztoibz(:)      ! (nktotf) full BZ → IBZ index
    INTEGER, ALLOCATABLE :: s_bztoibz(:)    ! (nktotf) symmetry op index for each BZ point
    INTEGER :: nibz                          ! number of IBZ points
    INTEGER, ALLOCATABLE :: ibz_idx(:)      ! (nibz) global indices of IBZ points
    INTEGER, ALLOCATABLE :: ibz_weight(:)   ! (nibz) number of BZ equivalents

    ! Symmetry rotation matrices in Cartesian
    REAL(dp) :: sr(3,3,48)  ! Cartesian rotation matrices

    REAL(dp) :: w_delta, costheta, vi(3), vf(3), vi_norm, vf_norm
    REAL(dp) :: wqf, n_d_cell, cell_area_bohr2, cell_area_cm2, inv_cell
    REAL(dp) :: dfdE, sigma_serta(3,3), sigma_mrta(3,3)
    REAL(dp) :: mobility_serta(3,3), mobility_mrta(3,3)
    REAL(dp) :: carrier_conc_cell, v_rot(3)
    CHARACTER(LEN=256) :: fname

    nktotf = nk1f * nk2f * nk3f
    ndim = 2
    IF (nk3f > 1) ndim = 3
    ! Spin degeneracy: 2 for collinear (each band is spin-degenerate), 1 for SOC
    IF (noncolin) THEN
       sfac = 1.0_dp
    ELSE
       sfac = 2.0_dp
    ENDIF
    win_min_Ry = transport_win_min / ryd2ev
    win_max_Ry = transport_win_max / ryd2ev
    wqf = 1.0_dp / DBLE(nktotf)

    ! Cell area/volume
    IF (ndim == 2) THEN
       cell_area_bohr2 = ABS(at(1,1)*at(2,2) - at(2,1)*at(1,2)) * alat**2
       cell_area_cm2 = cell_area_bohr2 * (bohr2ang * ang2cm)**2
       n_d_cell = defect_conc * cell_area_cm2
       inv_cell = 1.0_dp / cell_area_bohr2
    ELSE
       cell_area_cm2 = omega * (bohr2ang * ang2cm)**3
       n_d_cell = defect_conc * cell_area_cm2
       inv_cell = 1.0_dp / omega
    ENDIF

    ! Check nproc divisible by npool
    IF (MOD(nproc, npool) /= 0) THEN
       CALL errore('compute_transport', 'nproc must be divisible by npool', nproc)
    ENDIF

    IF (ionode) THEN
       WRITE(stdout, '(/,5X,A)') REPEAT('=', 60)
       WRITE(stdout, '(5X,A)')   'Phase 4: Transport (SERTA/MRTA) with IBZ symmetry'
       WRITE(stdout, '(5X,A)')   REPEAT('=', 60)
       WRITE(stdout, '(5X,A,3I4,A,I8)') 'Fine k-grid: ', nk1f, nk2f, nk3f, &
            '  total: ', nktotf
       WRITE(stdout, '(5X,A,F10.4,A,F10.4,A)') 'Energy window: ', transport_win_min, ' to ', transport_win_max, ' eV'
       WRITE(stdout, '(5X,A,ES10.3,A)') 'Carrier conc:  ', carrier_conc, ' cm^-2'
       WRITE(stdout, '(5X,A,ES10.4)') 'n_d per cell: ', n_d_cell
       WRITE(stdout, '(5X,A,I4,A,I4)') 'nsym = ', nsym, '  npool = ', npool
       WRITE(stdout, '(5X,A,F3.1)') 'Spin degeneracy (sfac): ', sfac
       WRITE(stdout, '(5X,A,A)') 'Delta method: ', TRIM(delta_method)
       IF (use_gaussian) &
            WRITE(stdout, '(5X,A,F8.4,A)') 'Gaussian sigma: ', delta_sigma, ' eV (fixed)'
       IF (use_adaptive) &
            WRITE(stdout, '(5X,A)') 'Adaptive sigma: velocity-dependent (EPW-style, no user parameter)'
       FLUSH(stdout)
    ENDIF

    ! =========================================================
    ! Determine delta function method
    use_gaussian = (TRIM(delta_method) == 'gaussian')
    use_adaptive = (TRIM(delta_method) == 'adaptive')

    ! Step 1: Build symmetry rotation matrices (crystal → Cartesian)
    ! Following EPW: sa = s, sb = bg*sa, sr = at*sb^T, sr = sr^T
    ! =========================================================
    CALL build_symm_cart(nsym, s, at, bg, sr)

    ! =========================================================
    ! Step 2: Pre-compute eigenvalues, eigenvectors, cfac, velocities
    ! =========================================================
    ALLOCATE(eig_all_eV(nbndsub, nktotf), eig_all_Ry(nbndsub, nktotf))
    ALLOCATE(evec_all(nbndsub, nbndsub, nktotf))
    ALLOCATE(cfac_all(nrr, nktotf))
    ALLOCATE(vel(3, nbndsub, nktotf))

    BLOCK
       REAL(dp) :: xk_cryst(3)
       INTEGER :: ik_count
       ik_count = 0
       DO i3 = 0, nk3f - 1
          DO i2 = 0, nk2f - 1
             DO i1 = 0, nk1f - 1
                ik_count = ik_count + 1
                xk_cryst(1) = DBLE(i1) / DBLE(nk1f)
                xk_cryst(2) = DBLE(i2) / DBLE(nk2f)
                xk_cryst(3) = DBLE(i3) / DBLE(nk3f)
                CALL get_cfac(nrr, irvec, xk_cryst, cfac_all(:, ik_count))
                CALL hamwan2bloch_with_evec(nbndsub, nrr, ndegen, &
                     eig_all_eV(:, ik_count), evec_all(:,:, ik_count), chw, cfac_all(:, ik_count))
             ENDDO
          ENDDO
       ENDDO
    END BLOCK

    eig_all_Ry = eig_all_eV / ryd2ev
    CALL compute_velocities_fast(nbndsub, nk1f, nk2f, nk3f, nktotf, eig_all_Ry, vel)

    ! Pre-compute adaptive smearing widths if requested
    ! σ(n,k) = 0.5 × sqrt(Σ_α [(2π/alat) × |v_nk · bg_α| / nkf_α]²) / sqrt(12)
    ! Depends only on (ibnd, ik), reused for all kf
    IF (use_adaptive) THEN
       ALLOCATE(sigma_adapt(nbndsub, nktotf))
       BLOCK
          INTEGER :: ik_a, ib_a
          REAL(dp) :: eta_tmp(3), vi_loc(3)
          REAL(dp), PARAMETER :: sigma_floor = 1.0d-3 / ryd2ev  ! 1 meV floor in Ry
          DO ik_a = 1, nktotf
             DO ib_a = 1, nbndsub
                vi_loc = vel(:, ib_a, ik_a)  ! velocity in Ry*bohr
                ! Project velocity onto reciprocal lattice directions, scale by grid spacing
                eta_tmp(1) = (twopi / alat) * ABS(DOT_PRODUCT(vi_loc, bg(:,1))) / DBLE(nk1f)
                eta_tmp(2) = (twopi / alat) * ABS(DOT_PRODUCT(vi_loc, bg(:,2))) / DBLE(nk2f)
                IF (nk3f > 1) THEN
                   eta_tmp(3) = (twopi / alat) * ABS(DOT_PRODUCT(vi_loc, bg(:,3))) / DBLE(nk3f)
                ELSE
                   eta_tmp(3) = 0.0_dp
                ENDIF
                sigma_adapt(ib_a, ik_a) = 0.5_dp * &
                     SQRT(eta_tmp(1)**2 + eta_tmp(2)**2 + eta_tmp(3)**2) / SQRT(12.0_dp)
                sigma_adapt(ib_a, ik_a) = MAX(sigma_adapt(ib_a, ik_a), sigma_floor)
             ENDDO
          ENDDO
       END BLOCK
       IF (ionode) THEN
          WRITE(stdout, '(5X,A,F8.4,A,F8.4,A)') 'Adaptive sigma range: ', &
               MINVAL(sigma_adapt) * ryd2ev * 1000.0_dp, ' - ', &
               MAXVAL(sigma_adapt) * ryd2ev * 1000.0_dp, ' meV'
          FLUSH(stdout)
       ENDIF
    ENDIF

    ! Delta function arrays (Ry)
    ALLOCATE(bande_Ry(0:nktotf-1, nbndsub), freq_zero(0:nktotf-1, 1))
    freq_zero = 0.0_dp
    DO iki = 1, nktotf
       bande_Ry(iki - 1, :) = eig_all_Ry(:, iki)
    ENDDO

    ! =========================================================
    ! Step 3: Build IBZ mapping
    ! =========================================================
    ALLOCATE(bztoibz(nktotf), s_bztoibz(nktotf))
    CALL build_ibz_mapping(nk1f, nk2f, nk3f, nktotf, nsym, s, &
                            bztoibz, s_bztoibz, nibz, ibz_idx, ibz_weight)

    IF (ionode) THEN
       WRITE(stdout, '(5X,A,I8,A,I8)') 'IBZ k-points: ', nibz, ' / ', nktotf
       WRITE(stdout, '(5X,A,F6.1,A)') 'Symmetry speedup: ', DBLE(nktotf)/DBLE(nibz), 'x'
       WRITE(stdout, '(5X,A)') 'Band energy ranges (eV):'
       DO ibnd = 1, nbndsub
          WRITE(stdout, '(5X,A,I2,A,F10.4,A,F10.4)') &
               '  Band ', ibnd, ':  min=', MINVAL(eig_all_eV(ibnd,:)), &
               '  max=', MAXVAL(eig_all_eV(ibnd,:))
       ENDDO
       FLUSH(stdout)
    ENDIF

    ! Fermi window mask and count active IBZ points
    ALLOCATE(kf_active(nktotf))
    kf_active = .FALSE.
    DO ikf = 1, nktotf
       IF (ANY(eig_all_Ry(:,ikf) >= win_min_Ry .AND. eig_all_Ry(:,ikf) <= win_max_Ry)) kf_active(ikf) = .TRUE.
    ENDDO

    ! Build list of active (ki_IBZ, kf_fullBZ) pairs after Fermi filtering
    ! ki: IBZ only (inv_tau is symmetry-invariant for both SERTA and MRTA)
    ! kf: full BZ (the sum over kf must cover the entire BZ)
    BLOCK
       INTEGER :: iki_tmp, ikf_tmp, nkf_active, npairs_total
       INTEGER :: pair_lower, pair_upper, npairs_local

       ! Build active ki list (IBZ, Fermi-filtered)
       n_active_ki = 0
       DO iki_tmp = 1, nibz
          IF (ANY(eig_all_Ry(:, ibz_idx(iki_tmp)) >= win_min_Ry .AND. eig_all_Ry(:, ibz_idx(iki_tmp)) <= win_max_Ry)) &
               n_active_ki = n_active_ki + 1
       ENDDO
       ALLOCATE(active_ibz_list(n_active_ki))
       n_active_ki = 0
       DO iki_tmp = 1, nibz
          IF (ANY(eig_all_Ry(:, ibz_idx(iki_tmp)) >= win_min_Ry .AND. eig_all_Ry(:, ibz_idx(iki_tmp)) <= win_max_Ry)) THEN
             n_active_ki = n_active_ki + 1
             active_ibz_list(n_active_ki) = iki_tmp  ! IBZ index
          ENDIF
       ENDDO

       ! Build active kf list
       nkf_active = COUNT(kf_active)
       ALLOCATE(active_kf_list(nkf_active))
       n_active_kf = 0
       DO ikf_tmp = 1, nktotf
          IF (kf_active(ikf_tmp)) THEN
             n_active_kf = n_active_kf + 1
             active_kf_list(n_active_kf) = ikf_tmp
          ENDIF
       ENDDO

       npairs_total = n_active_ki * n_active_kf

       ! Distribute pairs across pools
       pair_lower = my_pool_id * (npairs_total / npool) + MIN(my_pool_id, MOD(npairs_total, npool)) + 1
       IF (my_pool_id < MOD(npairs_total, npool)) THEN
          npairs_local = npairs_total / npool + 1
       ELSE
          npairs_local = npairs_total / npool
       ENDIF
       pair_upper = pair_lower + npairs_local - 1
       my_pair_lower = pair_lower
       my_pair_upper = pair_upper
       my_npairs_local = npairs_local

       IF (ionode) THEN
          WRITE(stdout, '(5X,A,F10.4,A,F10.4,A)') &
               'Transport window: [', transport_win_min, ', ', transport_win_max, '] eV'
          WRITE(stdout, '(5X,A,I8,A,I8)') 'Active kf (full BZ): ', n_active_kf, ' / ', nktotf
          WRITE(stdout, '(5X,A,I8,A,I8)') 'Active ki (IBZ):     ', n_active_ki, ' / ', nibz
          WRITE(stdout, '(5X,A,I10)') 'Total (ki,kf) pairs: ', npairs_total
          WRITE(stdout, '(5X,A,F8.1)') 'Pairs per pool:      ', &
               DBLE(npairs_total) / DBLE(MAX(1, npool))
          FLUSH(stdout)
       ENDIF
    END BLOCK

    ALLOCATE(inv_tau_serta(nbndsub, nktotf))
    ALLOCATE(inv_tau_mrta(nbndsub, nktotf))

    iunit = 88
    IF (ionode) THEN
       fname = TRIM(prefix_in) // '_transport.dat'
       OPEN(iunit, FILE=TRIM(fname), FORM='formatted')
       WRITE(iunit, '(A)') '# EDI Transport with IBZ symmetry'
       WRITE(iunit, '(A,3I4,A,I8,A,I8)') '# Grid: ', nk1f, nk2f, nk3f, &
            '  full=', nktotf, '  IBZ=', nibz
       WRITE(iunit, '(A,F10.4,A,F10.4,A)') '# Window=', transport_win_min, ' to ', transport_win_max, ' eV'
       WRITE(iunit, '(A)') '# T(K)  mu_SERTA_xx  mu_MRTA_xx  mu_SERTA_yy  mu_MRTA_yy  (cm^2/Vs)'
    ENDIF

    ! =========================================================
    ! Step 4: Temperature loop
    ! =========================================================
    DO itemp = 1, nstemp
       temp = temps(itemp)
       beta_eV = 1.0_dp / (kb_eV * temp)

       ! Compute Fermi level from carrier_conc via bisection
       ! n_target = carrier_conc × Ω_cell_cm
       ! Solve: (1/Nk) Σ_{n,k∈window} f(ε_nk, E_F, T) = n_target
       BLOCK
          REAL(dp) :: n_target, n_lo, n_hi, n_mid, ef_lo, ef_hi, ef_mid
          INTEGER :: ibis, ik_bis, ib_bis
          n_target = carrier_conc * cell_area_cm2  ! per cell (dimensionless)
          ! E_F can be outside the energy window (e.g. in the band gap for
          ! dilute carrier concentrations), so search a wider range
          ef_lo = transport_win_min - 2.0_dp  ! 2 eV below window
          ef_hi = transport_win_max + 2.0_dp  ! 2 eV above window
          DO ibis = 1, 100  ! max bisection iterations
             ef_mid = 0.5_dp * (ef_lo + ef_hi)
             n_mid = 0.0_dp
             DO ik_bis = 1, nktotf
                DO ib_bis = 1, nbndsub
                   IF (eig_all_eV(ib_bis, ik_bis) < transport_win_min .OR. &
                       eig_all_eV(ib_bis, ik_bis) > transport_win_max) CYCLE
                   n_mid = n_mid + sfac * wqf * fermi_func(eig_all_eV(ib_bis, ik_bis), ef_mid, beta_eV)
                ENDDO
             ENDDO
             IF (n_mid < n_target) THEN
                ef_lo = ef_mid  ! need more carriers → raise E_F
             ELSE
                ef_hi = ef_mid  ! too many carriers → lower E_F
             ENDIF
             IF (ABS(ef_hi - ef_lo) < 1.0d-8) EXIT
          ENDDO
          ef_eV = 0.5_dp * (ef_lo + ef_hi)
          ef_Ry = ef_eV / ryd2ev
       END BLOCK

       IF (ionode) THEN
          WRITE(stdout, '(/,5X,A,F8.1,A,F10.4,A)') &
               'T = ', temp, ' K  E_F = ', ef_eV, ' eV'
          FLUSH(stdout)
       ENDIF

       inv_tau_serta = 0.0_dp
       inv_tau_mrta = 0.0_dp
       npairs_computed = 0

       ! =========================================================
       ! Step 5: Scattering rate loop — (ki,kf) PAIRS distributed across pools
       ! Pairs are flattened: pair_id = (iki_idx-1)*n_active_kf + ikf_idx
       ! This gives optimal load balance — every pool processes equal work.
       ! =========================================================
       BLOCK
          INTEGER :: ipair, iki_idx, ikf_idx, iki_bz

          DO ipair = my_pair_lower, my_pair_upper
             ! Decode pair index -> (iki_idx, ikf_idx) in active lists
             iki_idx = (ipair - 1) / n_active_kf + 1
             ikf_idx = MOD(ipair - 1, n_active_kf) + 1

             iki_bz = ibz_idx(active_ibz_list(iki_idx))  ! IBZ index → full BZ index
             ikf = active_kf_list(ikf_idx)

             ! Progress print from ionode every 20%
             IF (ionode .AND. MOD(ipair - my_pair_lower + 1, MAX(1, my_npairs_local/5)) == 0) THEN
                WRITE(stdout, '(5X,A,I8,A,I8,A,I8,A)') &
                     '  pairs: ', ipair - my_pair_lower + 1, ' / ', my_npairs_local, &
                     '  (total: ', n_active_ki * n_active_kf, ')'
                FLUSH(stdout)
             ENDIF

             CALL edmatwan2bloch_2d(nbndsub, nrr, ndegen, edmatw_2d, &
                  cfac_all(:,iki_bz), cfac_all(:,ikf), edmatf_w)
             CALL ZGEMM('N', 'N', nbndsub, nbndsub, nbndsub, &
                         cone, edmatf_w, nbndsub, evec_all(:,:,ikf), nbndsub, &
                         czero, tmp_mat, nbndsub)
             CALL ZGEMM('C', 'N', nbndsub, nbndsub, nbndsub, &
                         cone, evec_all(:,:,iki_bz), nbndsub, tmp_mat, nbndsub, &
                         czero, edmatf_b, nbndsub)

             DO ibnd = 1, nbndsub
                IF (eig_all_Ry(ibnd, iki_bz) < win_min_Ry .OR. eig_all_Ry(ibnd, iki_bz) > win_max_Ry) CYCLE
                DO jbnd = 1, nbndsub
                   IF (eig_all_Ry(jbnd, ikf) < win_min_Ry .OR. eig_all_Ry(jbnd, ikf) > win_max_Ry) CYCLE

                   IF (use_gaussian .OR. use_adaptive) THEN
                      ! Gaussian broadening: δ(E) ≈ (1/σ√2π) exp(-E²/2σ²)
                      ! PURE delta function in 1/Ry — no BZ normalization (wqf provides 1/Nk)
                      BLOCK
                         REAL(dp) :: dE_gauss, sigma_Ry, xg
                         IF (use_adaptive) THEN
                            sigma_Ry = sigma_adapt(ibnd, iki_bz)
                         ELSE
                            sigma_Ry = delta_sigma / ryd2ev
                         ENDIF
                         dE_gauss = eig_all_Ry(ibnd, iki_bz) - eig_all_Ry(jbnd, ikf)
                         xg = dE_gauss / sigma_Ry
                         IF (ABS(xg) < 6.0_dp) THEN
                            w_delta = (1.0_dp / (sigma_Ry * SQRT(twopi))) * &
                                 EXP(-0.5_dp * xg**2)
                         ELSE
                            w_delta = 0.0_dp
                         ENDIF
                      END BLOCK
                   ELSE
                      ! Triangular (2D) or tetrahedral (3D)
                      ! get_delta_weight returns π/Nk × area_frac
                      ! Divide by π AND multiply by Nk to get PURE δ(E) in 1/Ry
                      ! (wqf provides the 1/Nk BZ normalization separately)
                      BLOCK
                         INTEGER :: iq_transfer
                         REAL(dp), PARAMETER :: pi_norm = 3.141592653589793_dp
                         iq_transfer = kindex_sub_2d(ikf-1, iki_bz-1, nk1f, nk2f)
                         w_delta = get_delta_weight(iq_transfer, iki_bz-1, ibnd, jbnd, 1, &
                                    freq_zero, bande_Ry, nk1f, nk2f, nk3f, &
                                    1.0_dp, 0.0_dp, nktotf, nbndsub, 1, ndim) &
                                    / pi_norm * DBLE(nktotf)
                      END BLOCK
                   ENDIF
                   IF (ABS(w_delta) < 1.0d-8) CYCLE

                   ! inv_tau = 2π × n_d_cell × (1/Nk) × Σ_kf |M|² × δ(ε)
                   ! w_delta is PURE δ(ε) in 1/Ry (no BZ normalization), wqf provides 1/Nk
                   inv_tau_serta(ibnd, iki_bz) = inv_tau_serta(ibnd, iki_bz) + &
                        twopi * n_d_cell * wqf * ABS(edmatf_b(ibnd, jbnd))**2 * w_delta

                   vi = vel(:, ibnd, iki_bz)
                   vf = vel(:, jbnd, ikf)
                   vi_norm = SQRT(SUM(vi**2))
                   vf_norm = SQRT(SUM(vf**2))
                   IF (vi_norm > 1.0d-10 .AND. vf_norm > 1.0d-10) THEN
                      costheta = DOT_PRODUCT(vi, vf) / (vi_norm * vf_norm)
                   ELSE
                      costheta = 0.0_dp
                   ENDIF
                   inv_tau_mrta(ibnd, iki_bz) = inv_tau_mrta(ibnd, iki_bz) + &
                        twopi * n_d_cell * wqf * ABS(edmatf_b(ibnd, jbnd))**2 * &
                        w_delta * (1.0_dp - costheta)
                ENDDO
             ENDDO
             npairs_computed = npairs_computed + 1
          ENDDO  ! ipair
       END BLOCK


       ! Gather pair count and inv_tau across pools
       CALL mp_sum(inv_tau_serta, inter_pool_comm)
       CALL mp_sum(inv_tau_mrta, inter_pool_comm)

       ! Unfold inv_tau from IBZ to full BZ (both SERTA and MRTA are
       ! symmetry-invariant scalars: the kf sum remaps under rotation)
       DO iki = 1, nktotf
          IF (bztoibz(iki) /= iki) THEN
             inv_tau_serta(:, iki) = inv_tau_serta(:, bztoibz(iki))
             inv_tau_mrta(:, iki) = inv_tau_mrta(:, bztoibz(iki))
          ENDIF
       ENDDO

       ! Diagnostics: inv_tau and velocity statistics
       IF (ionode) THEN
          BLOCK
             INTEGER :: n_act, ik_d, ib_d
             REAL(dp) :: sum_s, sum_m, max_s, min_s, avg_v2
             n_act = 0; sum_s = 0; sum_m = 0; max_s = 0; min_s = 1d30; avg_v2 = 0
             DO ik_d = 1, nktotf
                DO ib_d = 1, nbndsub
                   IF (inv_tau_serta(ib_d, ik_d) > 1.0d-10) THEN
                      n_act = n_act + 1
                      sum_s = sum_s + inv_tau_serta(ib_d, ik_d)
                      sum_m = sum_m + inv_tau_mrta(ib_d, ik_d)
                      max_s = MAX(max_s, inv_tau_serta(ib_d, ik_d))
                      min_s = MIN(min_s, inv_tau_serta(ib_d, ik_d))
                      avg_v2 = avg_v2 + SUM(vel(:,ib_d,ik_d)**2)
                   ENDIF
                ENDDO
             ENDDO
             IF (n_act > 0) THEN
                WRITE(stdout, '(5X,A,I6)') '  States with nonzero inv_tau: ', n_act
                WRITE(stdout, '(5X,A,ES12.4,A,ES12.4)') '  avg inv_tau: SERTA=', sum_s/n_act, '  MRTA=', sum_m/n_act
                WRITE(stdout, '(5X,A,ES12.4,A,ES12.4)') '  inv_tau range: ', min_s, ' - ', max_s
                WRITE(stdout, '(5X,A,ES12.4,A)') '  avg |v|^2 = ', avg_v2/n_act, ' (Ry*bohr)^2'
                WRITE(stdout, '(5X,A,ES12.4,A)') '  avg tau SERTA = ', 0.04838_dp/(sum_s/n_act), ' fs'
             ENDIF
             FLUSH(stdout)
          END BLOCK
       ENDIF

       ! =========================================================
       ! Step 6: Compute conductivity
       ! SERTA: use IBZ + symmetry-rotated velocities (inv_tau_SERTA is scalar-invariant)
       ! MRTA: use FULL BZ directly (inv_tau_MRTA depends on velocity direction)
       ! Both enforce σ_xx = σ_yy for hexagonal systems
       ! =========================================================
       sigma_serta = 0.0_dp
       sigma_mrta = 0.0_dp
       carrier_conc_cell = 0.0_dp

       ! --- SERTA σ: IBZ with symmetry rotation ---
       DO iki = 1, nibz
          BLOCK
             INTEGER :: iki_bz, isym, nb

             iki_bz = ibz_idx(iki)
             DO ibnd = 1, nbndsub
                IF (eig_all_eV(ibnd, iki_bz) < transport_win_min .OR. eig_all_eV(ibnd, iki_bz) > transport_win_max) CYCLE

                ! -df/dE in 1/eV → -df/dε in 1/Ry: multiply by ryd2ev (chain rule)
                dfdE = fermi_deriv(eig_all_eV(ibnd, iki_bz), ef_eV, beta_eV) * ryd2ev
                vi = vel(:, ibnd, iki_bz)

                carrier_conc_cell = carrier_conc_cell + &
                     sfac * (DBLE(ibz_weight(iki)) / DBLE(nktotf)) * &
                     fermi_func(eig_all_eV(ibnd, iki_bz), ef_eV, beta_eV)

                IF (inv_tau_serta(ibnd, iki_bz) > 1.0d-10) THEN
                   DO nb = 1, nktotf
                      IF (bztoibz(nb) /= iki_bz) CYCLE
                      isym = s_bztoibz(nb)
                      v_rot(1) = sr(1,1,isym)*vi(1) + sr(1,2,isym)*vi(2) + sr(1,3,isym)*vi(3)
                      v_rot(2) = sr(2,1,isym)*vi(1) + sr(2,2,isym)*vi(2) + sr(2,3,isym)*vi(3)
                      v_rot(3) = sr(3,1,isym)*vi(1) + sr(3,2,isym)*vi(2) + sr(3,3,isym)*vi(3)
                      DO i1 = 1, 3
                         DO i2 = 1, 3
                            sigma_serta(i1, i2) = sigma_serta(i1, i2) + &
                                 sfac * wqf * v_rot(i1) * v_rot(i2) / inv_tau_serta(ibnd, iki_bz) * dfdE
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
             ENDDO
          END BLOCK
       ENDDO

       ! --- MRTA σ: IBZ with symmetry-rotated velocities (same as SERTA) ---
       ! inv_tau_MRTA is also symmetry-invariant (kf sum remaps under rotation)
       DO iki = 1, nibz
          BLOCK
             INTEGER :: iki_bz2, isym2, nb2
             iki_bz2 = ibz_idx(iki)
             DO ibnd = 1, nbndsub
                IF (eig_all_eV(ibnd, iki_bz2) < transport_win_min .OR. eig_all_eV(ibnd, iki_bz2) > transport_win_max) CYCLE
                IF (inv_tau_mrta(ibnd, iki_bz2) > 1.0d-10) THEN
                   dfdE = fermi_deriv(eig_all_eV(ibnd, iki_bz2), ef_eV, beta_eV) * ryd2ev
                   vi = vel(:, ibnd, iki_bz2)
                   DO nb2 = 1, nktotf
                      IF (bztoibz(nb2) /= iki_bz2) CYCLE
                      isym2 = s_bztoibz(nb2)
                      v_rot(1) = sr(1,1,isym2)*vi(1) + sr(1,2,isym2)*vi(2) + sr(1,3,isym2)*vi(3)
                      v_rot(2) = sr(2,1,isym2)*vi(1) + sr(2,2,isym2)*vi(2) + sr(2,3,isym2)*vi(3)
                      v_rot(3) = sr(3,1,isym2)*vi(1) + sr(3,2,isym2)*vi(2) + sr(3,3,isym2)*vi(3)
                      DO i1 = 1, 3
                         DO i2 = 1, 3
                            sigma_mrta(i1, i2) = sigma_mrta(i1, i2) + &
                                 sfac * wqf * v_rot(i1) * v_rot(i2) / inv_tau_mrta(ibnd, iki_bz2) * dfdE
                         ENDDO
                      ENDDO
                   ENDDO
                ENDIF
             ENDDO
          END BLOCK
       ENDDO

       ! Convert σ → μ
       mobility_serta = 0.0_dp
       mobility_mrta = 0.0_dp
       IF (carrier_conc_cell > 1.0d-30) THEN
          BLOCK
             REAL(dp) :: conv
             conv = electron_SI * (bohr2ang * ang2cm)**2 / hbarJ
             mobility_serta = sigma_serta * conv / carrier_conc_cell
             mobility_mrta = sigma_mrta * conv / carrier_conc_cell
          END BLOCK
       ENDIF

       ! Gather total pair count for reporting
       CALL mp_sum(npairs_computed, inter_pool_comm)

       IF (ionode) THEN
          WRITE(stdout, '(5X,A,F8.1,A,I10,A)') 'T = ', temp, ' K  (', npairs_computed, ' ki-kf pairs computed)'
          WRITE(stdout, '(5X,A,ES12.4,A)') '  mu_SERTA_xx = ', mobility_serta(1,1), ' cm^2/Vs'
          WRITE(stdout, '(5X,A,ES12.4,A)') '  mu_SERTA_yy = ', mobility_serta(2,2), ' cm^2/Vs'
          WRITE(stdout, '(5X,A,ES12.4,A)') '  mu_MRTA_xx  = ', mobility_mrta(1,1), ' cm^2/Vs'
          WRITE(stdout, '(5X,A,ES12.4,A)') '  mu_MRTA_yy  = ', mobility_mrta(2,2), ' cm^2/Vs'
          FLUSH(stdout)
          WRITE(iunit, '(F10.2, 4ES16.6)') temp, &
               mobility_serta(1,1), mobility_mrta(1,1), &
               mobility_serta(2,2), mobility_mrta(2,2)
       ENDIF
    ENDDO  ! itemp

    ! Output
    IF (ionode) THEN
       CLOSE(iunit)
       fname = TRIM(prefix_in) // '_transport.dat'
       WRITE(stdout, '(/,5X,A,A)') 'Mobility written to ', TRIM(fname)

       fname = TRIM(prefix_in) // '_inv_tau.dat'
       OPEN(iunit, FILE=TRIM(fname), FORM='formatted')
       WRITE(iunit, '(A)') '# ik  ibnd  E(eV)  inv_tau_SERTA(Ry)  inv_tau_MRTA(Ry)  tau_SERTA(fs)  tau_MRTA(fs)'
       DO iki = 1, nktotf
          DO ibnd = 1, nbndsub
             IF (inv_tau_serta(ibnd, iki) > 1.0d-10 .OR. inv_tau_mrta(ibnd, iki) > 1.0d-10) THEN
                WRITE(iunit, '(I8, I4, F12.6, 4ES16.6)') iki, ibnd, &
                     eig_all_eV(ibnd, iki), &
                     inv_tau_serta(ibnd, iki), inv_tau_mrta(ibnd, iki), &
                     0.04838_dp / MAX(inv_tau_serta(ibnd, iki), 1.0d-30), &
                     0.04838_dp / MAX(inv_tau_mrta(ibnd, iki), 1.0d-30)
             ENDIF
          ENDDO
       ENDDO
       CLOSE(iunit)
       WRITE(stdout, '(5X,A,A)') 'Scattering rates written to ', TRIM(fname)
       WRITE(stdout, '(5X,A)') REPEAT('=', 60)
    ENDIF

    DEALLOCATE(eig_all_eV, eig_all_Ry, evec_all, cfac_all, vel, kf_active)
    DEALLOCATE(inv_tau_serta, inv_tau_mrta, bande_Ry, freq_zero)
    DEALLOCATE(bztoibz, s_bztoibz, ibz_idx, ibz_weight, active_ibz_list, active_kf_list)
    IF (ALLOCATED(sigma_adapt)) DEALLOCATE(sigma_adapt)

  CONTAINS

    REAL(dp) FUNCTION fermi_func(E, Ef, beta)
      REAL(dp), INTENT(IN) :: E, Ef, beta
      REAL(dp) :: x
      x = (E - Ef) * beta
      IF (x > 40.0_dp) THEN; fermi_func = 0.0_dp
      ELSEIF (x < -40.0_dp) THEN; fermi_func = 1.0_dp
      ELSE; fermi_func = 1.0_dp / (EXP(x) + 1.0_dp)
      ENDIF
    END FUNCTION

    REAL(dp) FUNCTION fermi_deriv(E, Ef, beta)
      REAL(dp), INTENT(IN) :: E, Ef, beta
      REAL(dp) :: f
      f = fermi_func(E, Ef, beta)
      fermi_deriv = beta * f * (1.0_dp - f)
    END FUNCTION

  END SUBROUTINE compute_transport


  SUBROUTINE build_symm_cart(nsym, s, at, bg, sr)
    !-----------------------------------------------------------------------
    ! Convert symmetry matrices from crystal to Cartesian (EPW convention)
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nsym, s(3,3,48)
    REAL(dp), INTENT(IN) :: at(3,3), bg(3,3)
    REAL(dp), INTENT(OUT) :: sr(3,3,48)
    INTEGER :: isym
    REAL(dp) :: sa(3,3), sb(3,3)

    DO isym = 1, nsym
       sa = DBLE(s(:,:,isym))
       sb = MATMUL(bg, sa)
       sr(:,:,isym) = TRANSPOSE(MATMUL(at, TRANSPOSE(sb)))
    ENDDO
  END SUBROUTINE build_symm_cart


  SUBROUTINE build_ibz_mapping(nk1, nk2, nk3, nktotf, nsym, s, &
                                bztoibz, s_bztoibz, nibz, ibz_idx, ibz_weight)
    !-----------------------------------------------------------------------
    ! Build full BZ → IBZ mapping using crystal symmetry operations.
    ! For each k-point, apply all symmetry ops and map to the smallest index.
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nk1, nk2, nk3, nktotf, nsym, s(3,3,48)
    INTEGER, INTENT(OUT) :: bztoibz(nktotf), s_bztoibz(nktotf), nibz
    INTEGER, ALLOCATABLE, INTENT(OUT) :: ibz_idx(:), ibz_weight(:)

    INTEGER :: ik, isym, i1, i2, i3, j1, j2, j3, jk
    INTEGER :: sk(3), kv(3)
    INTEGER, ALLOCATABLE :: equiv(:)

    ALLOCATE(equiv(nktotf))

    ! Initialize: each point maps to itself
    DO ik = 1, nktotf
       equiv(ik) = ik
       s_bztoibz(ik) = 1  ! identity
    ENDDO

    ! For each k-point, find the smallest equivalent under all symmetries
    ik = 0
    DO i3 = 0, nk3 - 1
       DO i2 = 0, nk2 - 1
          DO i1 = 0, nk1 - 1
             ik = ik + 1
             kv = (/i1, i2, i3/)

             DO isym = 1, nsym
                ! Apply symmetry: sk = s × kv (in grid units)
                sk(1) = s(1,1,isym)*kv(1) + s(1,2,isym)*kv(2) + s(1,3,isym)*kv(3)
                sk(2) = s(2,1,isym)*kv(1) + s(2,2,isym)*kv(2) + s(2,3,isym)*kv(3)
                sk(3) = s(3,1,isym)*kv(1) + s(3,2,isym)*kv(2) + s(3,3,isym)*kv(3)

                ! Wrap to [0, nki)
                j1 = MOD(sk(1), nk1); IF (j1 < 0) j1 = j1 + nk1
                j2 = MOD(sk(2), nk2); IF (j2 < 0) j2 = j2 + nk2
                j3 = MOD(sk(3), nk3); IF (j3 < 0) j3 = j3 + nk3

                jk = j3 * nk2 * nk1 + j2 * nk1 + j1 + 1

                ! Map both to the smaller index
                IF (jk < equiv(ik)) THEN
                   equiv(ik) = jk
                   s_bztoibz(ik) = isym
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    ! Propagate: ensure equiv(ik) points to the true representative
    ! (handle chains: if equiv(a) = b and equiv(b) = c, then equiv(a) should be c)
    DO ik = 1, nktotf
       DO WHILE (equiv(equiv(ik)) /= equiv(ik))
          equiv(ik) = equiv(equiv(ik))
       ENDDO
    ENDDO

    ! Build bztoibz
    bztoibz = equiv

    ! Find the symmetry op for each BZ point → its IBZ representative
    ! (re-scan to get the correct s for the final representative)
    ik = 0
    DO i3 = 0, nk3 - 1
       DO i2 = 0, nk2 - 1
          DO i1 = 0, nk1 - 1
             ik = ik + 1
             IF (bztoibz(ik) == ik) THEN
                s_bztoibz(ik) = 1  ! identity for IBZ points
             ELSE
                ! Find which symmetry maps ik → bztoibz(ik)
                kv = (/i1, i2, i3/)
                DO isym = 1, nsym
                   sk(1) = s(1,1,isym)*kv(1) + s(1,2,isym)*kv(2) + s(1,3,isym)*kv(3)
                   sk(2) = s(2,1,isym)*kv(1) + s(2,2,isym)*kv(2) + s(2,3,isym)*kv(3)
                   sk(3) = s(3,1,isym)*kv(1) + s(3,2,isym)*kv(2) + s(3,3,isym)*kv(3)
                   j1 = MOD(sk(1), nk1); IF (j1 < 0) j1 = j1 + nk1
                   j2 = MOD(sk(2), nk2); IF (j2 < 0) j2 = j2 + nk2
                   j3 = MOD(sk(3), nk3); IF (j3 < 0) j3 = j3 + nk3
                   jk = j3 * nk2 * nk1 + j2 * nk1 + j1 + 1
                   IF (jk == bztoibz(ik)) THEN
                      s_bztoibz(ik) = isym
                      EXIT
                   ENDIF
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    ! Count IBZ points and build index/weight arrays
    nibz = 0
    DO ik = 1, nktotf
       IF (bztoibz(ik) == ik) nibz = nibz + 1
    ENDDO

    ALLOCATE(ibz_idx(nibz), ibz_weight(nibz))
    nibz = 0
    DO ik = 1, nktotf
       IF (bztoibz(ik) == ik) THEN
          nibz = nibz + 1
          ibz_idx(nibz) = ik
          ibz_weight(nibz) = 0
       ENDIF
    ENDDO

    ! Count weight (number of BZ equivalents) for each IBZ point
    DO ik = 1, nktotf
       DO i1 = 1, nibz
          IF (bztoibz(ik) == ibz_idx(i1)) THEN
             ibz_weight(i1) = ibz_weight(i1) + 1
             EXIT
          ENDIF
       ENDDO
    ENDDO

    DEALLOCATE(equiv)
  END SUBROUTINE build_ibz_mapping


  SUBROUTINE compute_velocities_fast(nbndsub, nk1f, nk2f, nk3f, nktotf, eig_all, vel)
    !-----------------------------------------------------------------------
    ! Group velocities via finite differences on the eigenvalue grid.
    ! v_α = dE/dk_cart_α in Ry·bohr (same units as EPW's vmef)
    ! Uses: v = (J^{-T}) × dE/dk_cryst, where J = bg × tpiba
    ! Equivalently: v_α = Σ_i at(α,i) × alat × dEdk_i / (2π)
    !             = Σ_i at(α,i) × dEdk_i / tpiba  (since alat/(2π) = 1/tpiba)
    !-----------------------------------------------------------------------
    USE cell_base, ONLY : at, alat, bg, tpiba
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nbndsub, nk1f, nk2f, nk3f, nktotf
    REAL(dp), INTENT(IN) :: eig_all(nbndsub, nktotf)
    REAL(dp), INTENT(OUT) :: vel(3, nbndsub, nktotf)

    INTEGER :: ik, i1, i2, i3, ibnd, ik_p, ik_m
    REAL(dp) :: dk1, dk2, dk3, dEdk(3)

    dk1 = 1.0_dp / DBLE(nk1f)
    dk2 = 1.0_dp / DBLE(nk2f)
    dk3 = 1.0_dp / DBLE(nk3f)

    ik = 0
    DO i3 = 0, nk3f - 1
       DO i2 = 0, nk2f - 1
          DO i1 = 0, nk1f - 1
             ik = ik + 1
             DO ibnd = 1, nbndsub
                ik_p = i3*nk2f*nk1f + i2*nk1f + MOD(i1+1,nk1f) + 1
                ik_m = i3*nk2f*nk1f + i2*nk1f + MOD(i1-1+nk1f,nk1f) + 1
                dEdk(1) = (eig_all(ibnd, ik_p) - eig_all(ibnd, ik_m)) / (2.0_dp * dk1)

                ik_p = i3*nk2f*nk1f + MOD(i2+1,nk2f)*nk1f + i1 + 1
                ik_m = i3*nk2f*nk1f + MOD(i2-1+nk2f,nk2f)*nk1f + i1 + 1
                dEdk(2) = (eig_all(ibnd, ik_p) - eig_all(ibnd, ik_m)) / (2.0_dp * dk2)

                IF (nk3f > 1) THEN
                   ik_p = MOD(i3+1,nk3f)*nk2f*nk1f + i2*nk1f + i1 + 1
                   ik_m = MOD(i3-1+nk3f,nk3f)*nk2f*nk1f + i2*nk1f + i1 + 1
                   dEdk(3) = (eig_all(ibnd, ik_p) - eig_all(ibnd, ik_m)) / (2.0_dp * dk3)
                ELSE
                   dEdk(3) = 0.0_dp
                ENDIF

                ! v_α = Σ_i at(α,i) × dEdk(i) × alat / (2π)  [Ry·bohr]
                ! = Σ_i at(α,i) × dEdk(i) / tpiba  (since tpiba = 2π/alat)
                vel(1,ibnd,ik) = (at(1,1)*dEdk(1)+at(1,2)*dEdk(2)+at(1,3)*dEdk(3))/tpiba
                vel(2,ibnd,ik) = (at(2,1)*dEdk(1)+at(2,2)*dEdk(2)+at(2,3)*dEdk(3))/tpiba
                vel(3,ibnd,ik) = (at(3,1)*dEdk(1)+at(3,2)*dEdk(2)+at(3,3)*dEdk(3))/tpiba
             ENDDO
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE compute_velocities_fast

  SUBROUTINE compute_dos_validation(nbndsub, nrr, irvec, ndegen, chw, &
                                      nk1f, nk2f, nk3f, prefix_in)
    !-----------------------------------------------------------------------
    ! For each ki in the Fermi window, compute the local DOS:
    !   dos(ki) = Σ_{m,kf} δ(ε_nki - ε_mkf)
    ! using ALL three methods: triangular, fixed Gaussian, adaptive Gaussian.
    !-----------------------------------------------------------------------
    USE io_global, ONLY : ionode, stdout
    USE cell_base, ONLY : at, alat, bg, tpiba
    USE wan2bloch_edi, ONLY : get_cfac, hamwan2bloch_with_evec
    USE delta_weights, ONLY : get_delta_weight, kindex_sub_2d
    USE edi_input, ONLY : transport_win_min, transport_win_max, delta_sigma
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nbndsub, nrr, nk1f, nk2f, nk3f
    INTEGER, INTENT(IN) :: irvec(3, nrr), ndegen(nrr)
    COMPLEX(dp), INTENT(IN) :: chw(nbndsub, nbndsub, nrr)
    CHARACTER(LEN=*), INTENT(IN) :: prefix_in

    INTEGER :: nktotf, ik, ikf, ibnd, jbnd, ndim, iunit
    INTEGER :: i1, i2, i3
    REAL(dp) :: win_min_Ry_dos, win_max_Ry_dos, w_tri, w_gauss, w_adapt
    REAL(dp) :: dos_tri, dos_gauss, dos_adapt
    REAL(dp) :: sigma_fix_Ry, sigma_ad, dE_ry, xarg
    REAL(dp), ALLOCATABLE :: eig_all_eV(:,:), eig_all_Ry(:,:)
    REAL(dp), ALLOCATABLE :: vel_dos(:,:,:)
    REAL(dp), ALLOCATABLE :: sigma_adapt_dos(:,:)
    REAL(dp), ALLOCATABLE :: bande_Ry(:,:), freq_zero(:,:)
    COMPLEX(dp), ALLOCATABLE :: cfac(:)
    COMPLEX(dp) :: evec_tmp(nbndsub, nbndsub)
    CHARACTER(LEN=256) :: fname
    REAL(dp), PARAMETER :: pi_val = 3.141592653589793_dp

    nktotf = nk1f * nk2f * nk3f
    ndim = 2
    IF (nk3f > 1) ndim = 3
    win_min_Ry_dos = transport_win_min / ryd2ev
    win_max_Ry_dos = transport_win_max / ryd2ev
    sigma_fix_Ry = delta_sigma / ryd2ev

    IF (ionode) THEN
       WRITE(stdout, '(/,5X,A)') REPEAT('=', 60)
       WRITE(stdout, '(5X,A)')   'DOS validation: triangular vs Gaussian vs adaptive'
       WRITE(stdout, '(5X,A,F8.4,A)') 'Fixed Gaussian sigma: ', delta_sigma, ' eV'
       WRITE(stdout, '(5X,A)')   REPEAT('=', 60)
       FLUSH(stdout)
    ENDIF

    ! Interpolate eigenvalues
    ALLOCATE(eig_all_eV(nbndsub, nktotf), eig_all_Ry(nbndsub, nktotf), cfac(nrr))
    ik = 0
    DO i3 = 0, nk3f - 1
       DO i2 = 0, nk2f - 1
          DO i1 = 0, nk1f - 1
             ik = ik + 1
             CALL get_cfac(nrr, irvec, &
                  (/DBLE(i1)/DBLE(nk1f), DBLE(i2)/DBLE(nk2f), DBLE(i3)/DBLE(nk3f)/), cfac)
             CALL hamwan2bloch_with_evec(nbndsub, nrr, ndegen, eig_all_eV(:,ik), &
                  evec_tmp, chw, cfac)
          ENDDO
       ENDDO
    ENDDO
    DEALLOCATE(cfac)
    eig_all_Ry = eig_all_eV / ryd2ev

    ! Compute velocities for adaptive sigma
    ALLOCATE(vel_dos(3, nbndsub, nktotf))
    CALL compute_velocities_fast(nbndsub, nk1f, nk2f, nk3f, nktotf, eig_all_Ry, vel_dos)

    ! Pre-compute adaptive sigma
    ALLOCATE(sigma_adapt_dos(nbndsub, nktotf))
    BLOCK
       INTEGER :: ik_a, ib_a
       REAL(dp) :: eta_tmp(3), vi_loc(3)
       REAL(dp), PARAMETER :: sigma_floor = 1.0d-3 / 13.6056980659_dp
       DO ik_a = 1, nktotf
          DO ib_a = 1, nbndsub
             vi_loc = vel_dos(:, ib_a, ik_a)
             eta_tmp(1) = (2.0_dp*pi_val / alat) * ABS(DOT_PRODUCT(vi_loc, bg(:,1))) / DBLE(nk1f)
             eta_tmp(2) = (2.0_dp*pi_val / alat) * ABS(DOT_PRODUCT(vi_loc, bg(:,2))) / DBLE(nk2f)
             IF (nk3f > 1) THEN
                eta_tmp(3) = (2.0_dp*pi_val / alat) * ABS(DOT_PRODUCT(vi_loc, bg(:,3))) / DBLE(nk3f)
             ELSE
                eta_tmp(3) = 0.0_dp
             ENDIF
             sigma_adapt_dos(ib_a, ik_a) = 0.5_dp * &
                  SQRT(eta_tmp(1)**2 + eta_tmp(2)**2 + eta_tmp(3)**2) / SQRT(12.0_dp)
             sigma_adapt_dos(ib_a, ik_a) = MAX(sigma_adapt_dos(ib_a, ik_a), sigma_floor)
          ENDDO
       ENDDO
    END BLOCK
    DEALLOCATE(vel_dos)

    IF (ionode) THEN
       WRITE(stdout, '(5X,A,F8.2,A,F8.2,A)') 'Adaptive sigma range: ', &
            MINVAL(sigma_adapt_dos) * ryd2ev * 1000.0_dp, ' - ', &
            MAXVAL(sigma_adapt_dos) * ryd2ev * 1000.0_dp, ' meV'
       FLUSH(stdout)
    ENDIF

    ! Prepare delta function arrays
    ALLOCATE(bande_Ry(0:nktotf-1, nbndsub), freq_zero(0:nktotf-1, 1))
    freq_zero = 0.0_dp
    DO ik = 1, nktotf
       bande_Ry(ik-1, :) = eig_all_Ry(:, ik)
    ENDDO

    ! Open output
    iunit = 89
    IF (ionode) THEN
       fname = TRIM(prefix_in) // '_dos_validation.dat'
       OPEN(iunit, FILE=TRIM(fname), FORM='formatted')
       WRITE(iunit, '(A)') '# DOS validation: for each (ki,ibnd) in Fermi window'
       WRITE(iunit, '(A)') '# dos = SUM_{m,kf} delta(eps_nki - eps_mkf)'
       WRITE(iunit, '(A)') '# ik  ibnd  E(eV)  dos_triangular  dos_gaussian_fixed  dos_gaussian_adaptive'
    ENDIF

    ! Compute DOS per ki with all three methods
    DO ik = 1, nktotf
       DO ibnd = 1, nbndsub
          IF (eig_all_Ry(ibnd, ik) < win_min_Ry_dos .OR. eig_all_Ry(ibnd, ik) > win_max_Ry_dos) CYCLE

          dos_tri = 0.0_dp
          dos_gauss = 0.0_dp
          dos_adapt = 0.0_dp
          sigma_ad = sigma_adapt_dos(ibnd, ik)

          DO ikf = 1, nktotf
             DO jbnd = 1, nbndsub
                IF (eig_all_Ry(jbnd, ikf) < win_min_Ry_dos .OR. eig_all_Ry(jbnd, ikf) > win_max_Ry_dos) CYCLE

                dE_ry = eig_all_Ry(ibnd, ik) - eig_all_Ry(jbnd, ikf)

                ! Triangular (iq = kf - ki); divide by π to match Gaussian normalization
                BLOCK
                   INTEGER :: iq_t
                   iq_t = kindex_sub_2d(ikf-1, ik-1, nk1f, nk2f)
                   w_tri = get_delta_weight(iq_t, ik-1, ibnd, jbnd, 1, &
                            freq_zero, bande_Ry, nk1f, nk2f, nk3f, &
                            1.0_dp, 0.0_dp, nktotf, nbndsub, 1, ndim) / pi_val
                END BLOCK
                dos_tri = dos_tri + w_tri

                ! Fixed Gaussian
                xarg = dE_ry / sigma_fix_Ry
                IF (ABS(xarg) < 6.0_dp) THEN
                   dos_gauss = dos_gauss + &
                        (1.0_dp / (sigma_fix_Ry * SQRT(2.0_dp * pi_val))) * &
                        EXP(-0.5_dp * xarg**2) / DBLE(nktotf)
                ENDIF

                ! Adaptive Gaussian
                xarg = dE_ry / sigma_ad
                IF (ABS(xarg) < 6.0_dp) THEN
                   dos_adapt = dos_adapt + &
                        (1.0_dp / (sigma_ad * SQRT(2.0_dp * pi_val))) * &
                        EXP(-0.5_dp * xarg**2) / DBLE(nktotf)
                ENDIF
             ENDDO
          ENDDO

          IF (ionode) THEN
             WRITE(iunit, '(I8, I4, F12.6, 3ES16.6)') &
                  ik, ibnd, eig_all_eV(ibnd, ik), dos_tri, dos_gauss, dos_adapt
          ENDIF
       ENDDO

       IF (ionode .AND. MOD(ik, MAX(1, nktotf/10)) == 0) THEN
          WRITE(stdout, '(5X,A,I8,A,I8)') '  ki: ', ik, ' / ', nktotf
          FLUSH(stdout)
       ENDIF
    ENDDO

    IF (ionode) THEN
       CLOSE(iunit)
       WRITE(stdout, '(5X,A,A)') 'DOS validation written to ', TRIM(fname)
       WRITE(stdout, '(5X,A)') REPEAT('=', 60)
    ENDIF

    DEALLOCATE(eig_all_eV, eig_all_Ry, sigma_adapt_dos, bande_Ry, freq_zero)
  END SUBROUTINE compute_dos_validation




END MODULE transport_edi
