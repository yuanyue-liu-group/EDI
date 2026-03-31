MODULE edi_pw2wan
  USE kinds, ONLY : dp
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC :: edi_run_wannier90
  PUBLIC :: edi_build_ham_r
  PUBLIC :: edi_interp_bands

CONTAINS

  SUBROUTINE edi_run_wannier90(seedname_in, nk1, nk2, nk3)
    USE io_global, ONLY : ionode, ionode_id, stdout
    USE mp, ONLY : mp_bcast
    USE mp_world, ONLY : world_comm
    USE cell_base, ONLY : at, bg
    USE klist, ONLY : nkstot, nks, xk
    USE io_files, ONLY : prefix
    USE noncollin_module, ONLY : noncolin
    USE wann_common
    USE input, ONLY : xk_all, xk_cryst, scdm_proj
    USE edi_input, ONLY : write_winfil_edi
    USE pw2wan, ONLY : pw2wan90epw, setup_nnkp, ylm_expansion, &
                        compute_amn_para, compute_mmn_para, write_band, run_wannier
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: seedname_in
    INTEGER, INTENT(IN) :: nk1, nk2, nk3
    INTEGER :: ierr

    WRITE(stdout, '(/,5X,A)') REPEAT('=', 60)
    WRITE(stdout, '(5X,A)') 'EDI: Running Wannier90 internally'
    WRITE(stdout, '(5X,A)') REPEAT('=', 60)

    CALL write_winfil_edi(TRIM(seedname_in))
    WRITE(stdout, '(5X,A)') 'Wrote ' // TRIM(seedname_in) // '.win'

    seedname = TRIM(seedname_in)
    seedname2 = TRIM(seedname_in)
    wan_mode = 'library'
    write_amn = .TRUE.
    write_mmn = .TRUE.
    write_eig = .TRUE.
    write_unk = .FALSE.
    write_spn = .FALSE.
    write_uhu = .FALSE.
    write_uIu = .FALSE.
    write_sHu = .FALSE.
    write_sIu = .FALSE.
    write_dmn = .FALSE.
    read_sym  = .FALSE.
    irr_bz    = .FALSE.
    wvfn_formatted = .FALSE.
    reduce_unk = .FALSE.
    reduce_unk_factor = 1
    logwann = .TRUE.
    scdm_proj = .FALSE.

    ispinw  = 0
    ikstart = 1
    ikstop  = nkstot
    iknum   = nkstot

    mp_grid(1) = nk1
    mp_grid(2) = nk2
    mp_grid(3) = nk3

    IF (nk1 * nk2 * nk3 /= iknum) &
       CALL errore('edi_run_wannier90', 'k-grid does not match nkstot', iknum)

    IF (.NOT. ALLOCATED(xk_all)) THEN
       ALLOCATE(xk_all(3, nkstot), STAT=ierr)
       IF (ierr /= 0) CALL errore('edi_run_wannier90', 'Error allocating xk_all', 1)
       CALL poolcollect(3, nks, xk, nkstot, xk_all)
    ENDIF

    WRITE(stdout, '(5X,A,I3,A,I3,A,I3,A,I6,A)') &
         'k-grid: ', nk1, ' x ', nk2, ' x ', nk3, &
         ' (', iknum, ' k-points)'

    ALLOCATE(kpt_latt(3, iknum), STAT=ierr)
    IF (ierr /= 0) CALL errore('edi_run_wannier90', 'Error allocating kpt_latt', 1)
    kpt_latt(:, 1:iknum) = xk_cryst(:, 1:iknum)
    CALL mp_bcast(kpt_latt, ionode_id, world_comm)

    WRITE(stdout, '(5X,A)') 'Setting up nearest-neighbor k-points...'
    CALL setup_nnkp()

    WRITE(stdout, '(5X,A)') 'Computing Ylm expansion...'
    CALL ylm_expansion()

    WRITE(stdout, '(5X,A)') 'Computing A (projection) matrix...'
    CALL compute_amn_para()

    WRITE(stdout, '(5X,A)') 'Computing M (overlap) matrix...'
    CALL compute_mmn_para()

    WRITE(stdout, '(5X,A)') 'Writing eigenvalues...'
    CALL write_band()

    WRITE(stdout, '(5X,A)') 'Running Wannier90 minimization...'
    CALL run_wannier()

    WRITE(stdout, '(5X,A)') 'Wannier90 completed.'
    WRITE(stdout, '(5X,A,I4)') '  Number of Wannier functions: ', n_wannier
    WRITE(stdout, '(5X,A,I4)') '  Number of bands: ', num_bands
    WRITE(stdout, '(5X,A,F12.4,A)') '  Total spread: ', spreads(1), ' Ang^2'
    WRITE(stdout, '(5X,A)') REPEAT('=', 60)
  END SUBROUTINE edi_run_wannier90

  SUBROUTINE edi_build_ham_r(nk1, nk2, nk3, nbndsub_out, nrr_out, &
                              irvec_out, ndegen_out, chw_out)
    USE io_global, ONLY : ionode, ionode_id, stdout
    USE kinds, ONLY : dp
    USE constants, ONLY : rytoev
    USE cell_base, ONLY : at
    USE klist, ONLY : nkstot, nks, xk
    USE wvfct, ONLY : nbnd, et
    USE mp_world, ONLY : world_comm
    USE wann_common, ONLY : u_mat, u_mat_opt, lwindow, n_wannier, &
                         num_bands, excluded_band, wann_centers, iknum
    USE input, ONLY : xk_cryst_input => xk_cryst
    USE wigner_seitz_edi, ONLY : wigner_seitz_k
    USE wan2bloch_edi, ONLY : hambloch2wan_edi
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nk1, nk2, nk3
    INTEGER, INTENT(OUT) :: nbndsub_out, nrr_out
    INTEGER, ALLOCATABLE, INTENT(OUT) :: irvec_out(:,:), ndegen_out(:)
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: chw_out(:,:,:)

    INTEGER :: ik, ib, ib_kept, ib_win, nbndsub
    REAL(dp), ALLOCATABLE :: et_sub(:,:), xk_cryst_loc(:,:), wslen(:)
    COMPLEX(dp), ALLOCATABLE :: cu(:,:,:)
    LOGICAL :: use_disentangle

    nbndsub = n_wannier
    nbndsub_out = nbndsub
    use_disentangle = (num_bands > n_wannier)

    ALLOCATE(et_sub(nbndsub, nkstot))
    ALLOCATE(cu(nbndsub, nbndsub, nkstot))
    ALLOCATE(xk_cryst_loc(3, nkstot))

    xk_cryst_loc(:,:) = xk_cryst_input(:,:)

    IF (use_disentangle) THEN
       DO ik = 1, nkstot
          ib_kept = 0
          ib_win = 0
          DO ib = 1, nbnd
             IF (excluded_band(ib)) CYCLE
             ib_kept = ib_kept + 1
             IF (lwindow(ib_kept, ik)) THEN
                ib_win = ib_win + 1
                et_sub(ib_win, ik) = et(ib, ik)
             ENDIF
          ENDDO
       ENDDO
       DO ik = 1, nkstot
          cu(:,:,ik) = MATMUL(u_mat_opt(:, :, ik), u_mat(:, :, ik))
       ENDDO
    ELSE
       DO ik = 1, nkstot
          ib_kept = 0
          DO ib = 1, nbnd
             IF (excluded_band(ib)) CYCLE
             ib_kept = ib_kept + 1
             IF (ib_kept <= nbndsub) et_sub(ib_kept, ik) = et(ib, ik)
          ENDDO
          cu(:,:,ik) = u_mat(:,:,ik)
       ENDDO
    ENDIF

    CALL wigner_seitz_k(nk1, nk2, nk3, at, irvec_out, ndegen_out, wslen, nrr_out)
    IF (ionode) WRITE(stdout, '(5X,A,I6)') 'Wigner-Seitz R-points: ', nrr_out

    ALLOCATE(chw_out(nbndsub, nbndsub, nrr_out))
    chw_out = (0.0_dp, 0.0_dp)

    BLOCK
      INTEGER :: ir_l, ik_l
      REAL(dp) :: rdotk_l
      COMPLEX(dp) :: cfac_l, ctmp_l
      COMPLEX(dp), ALLOCATABLE :: chs_l(:,:,:)
      REAL(dp), PARAMETER :: twopi_l = 6.283185307179586_dp
      COMPLEX(dp), PARAMETER :: ci_l = (0.0_dp, 1.0_dp)

      ALLOCATE(chs_l(nbndsub, nbndsub, nkstot))
      chs_l = (0.0_dp, 0.0_dp)
      DO ik_l = 1, nkstot
         DO ib_win = 1, nbndsub
            DO ib_kept = 1, ib_win
               ctmp_l = (0.0_dp, 0.0_dp)
               DO ib = 1, nbndsub
                  ctmp_l = ctmp_l + CONJG(cu(ib, ib_kept, ik_l)) * &
                           et_sub(ib, ik_l) * cu(ib, ib_win, ik_l)
               ENDDO
               chs_l(ib_kept, ib_win, ik_l) = ctmp_l
               chs_l(ib_win, ib_kept, ik_l) = CONJG(ctmp_l)
            ENDDO
         ENDDO
      ENDDO

      DO ir_l = 1, nrr_out
         DO ik_l = 1, nkstot
            rdotk_l = twopi_l * DOT_PRODUCT(xk_cryst_loc(:, ik_l), &
                                             DBLE(irvec_out(:, ir_l)))
            cfac_l = EXP(-ci_l * rdotk_l) / DBLE(nkstot)
            chw_out(:, :, ir_l) = chw_out(:, :, ir_l) + cfac_l * chs_l(:, :, ik_l)
         ENDDO
      ENDDO
      DEALLOCATE(chs_l)
    END BLOCK

    IF (ionode) WRITE(stdout, '(5X,A)') 'H(R) in Wannier basis constructed.'

    DEALLOCATE(et_sub, cu, xk_cryst_loc, wslen)
  END SUBROUTINE edi_build_ham_r

  SUBROUTINE edi_interp_bands(nbndsub, nrr, irvec, ndegen, chw, &
                               nk_fine, xk_fine, eig_fine)
    USE io_global, ONLY : stdout
    USE wan2bloch_edi, ONLY : hamwan2bloch_edi, get_cfac
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nbndsub, nrr, nk_fine
    INTEGER, INTENT(IN) :: irvec(3, nrr), ndegen(nrr)
    COMPLEX(dp), INTENT(IN) :: chw(nbndsub, nbndsub, nrr)
    REAL(dp), INTENT(IN) :: xk_fine(3, nk_fine)
    REAL(dp), INTENT(OUT) :: eig_fine(nbndsub, nk_fine)

    INTEGER :: ik
    COMPLEX(dp), ALLOCATABLE :: cfac(:)

    ALLOCATE(cfac(nrr))
    DO ik = 1, nk_fine
       CALL get_cfac(nrr, irvec, xk_fine(:, ik), cfac)
       CALL hamwan2bloch_edi(nbndsub, nrr, ndegen, eig_fine(:, ik), chw, cfac)
    ENDDO
    DEALLOCATE(cfac)
  END SUBROUTINE edi_interp_bands

END MODULE edi_pw2wan
