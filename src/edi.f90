Program edi
  USE kinds, ONLY : dp
  USE io_files, ONLY : prefix, tmp_dir, restart_dir
  USE wvfct, ONLY : nbnd, et
  USE klist, ONLY : nkstot, nks, xk
  USE lsda_mod, ONLY : isk
  USE cell_base, ONLY : at
  USE environment, ONLY : environment_start, environment_end
  USE mp_global, ONLY : mp_startup, mp_global_end, inter_pool_comm
  USE mp_pools, ONLY : npool
  USE mp_bands, ONLY : nbgrp
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp, ONLY : mp_bcast, mp_sum
  USE mp_world, ONLY : world_comm
  USE constants, ONLY : rytoev
  USE noncollin_module, ONLY : noncolin, lspinorb, npol

  USE wann_common
  USE global_var, ONLY : nbndep, nbndskip, ibndkept
  USE input, ONLY : xk_all, xk_cryst, xk_loc, nkc1, nkc2, nkc3
  USE edi_input, ONLY : read_edinput, write_winfil_edi, &
                         sync_input_module, &
                         coarse_nk1, coarse_nk2, coarse_nk3, &
                         fine_nk1, fine_nk2, fine_nk3, &
                         edi_prefix, edi_outdir, edinput_nml, &
                         nbndsub, proj, wdata, bands_skipped, &
                         dis_win_min, dis_win_max, dis_froz_min, dis_froz_max, &
                         num_iter, iprint_w90, auto_projections, scdm_proj, filukk, &
                         wannierize, do_edmat, edwwrite, edwread, &
                         edmat_direct_from_file, edmat_interp_from_file, &
                         filki_direct, filkf_direct, filki_interp, filkf_interp, &
                         sc_nk1, sc_nk2, sc_nk3, band_ed, &
                         do_transport, nstemp, temps, &
                         transport_win_min, transport_win_max, &
                         carrier_conc, defect_conc, delta_method, delta_sigma, &
                         pristine_prefix, pristine_outdir, &
                         defect_prefix, defect_outdir, &
                         potfile_d, potfile_p, &
                         pot_align, defect_center, core_align_radius
  USE edi_pw2wan, ONLY : edi_run_wannier90, edi_interp_bands
  USE edi_read_hr, ONLY : read_hr_file
  USE ed_coarse, ONLY : load_supercell_pot, load_pot_from_file, read_filukk_edi, &
                         ed_interp_from_file, ed_coarse_full_q, ed_fine_interp_2d, &
                         read_edmatw_2d_file, &
                         ed_direct_from_files, generate_nscf_input, &
                         build_vcolin_aligned, build_vcolin_corealign, &
                         write_vcolin_cube
  USE edic_mod, ONLY : V_d, V_p, V_colin

  IMPLICIT NONE

  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  INTEGER :: ios, ik, ibnd, ib, ib_kept, ierr
  INTEGER :: nbndsub_loc, nrr
  INTEGER, ALLOCATABLE :: irvec(:,:), ndegen(:)
  COMPLEX(dp), ALLOCATABLE :: chw(:,:,:)
  REAL(dp), ALLOCATABLE :: eig_interp(:,:), xk_fine(:,:), et_kept(:,:)
  REAL(dp) :: max_err, err, vac_shift
  LOGICAL :: need_wf
  CHARACTER(LEN=256) :: outdir

  CALL mp_startup(start_images=.TRUE.)
  CALL environment_start('EDI')

  IF (nbgrp > 1) CALL errore('edi', 'band groups not supported', nbgrp)

  ios = 0

  IF (ionode) THEN
     CALL input_from_file()
     READ(5, NML=edinput_nml, IOSTAT=ios)
     IF (ios < 0) ios = 0
  ENDIF

  CALL mp_bcast(ios, ionode_id, world_comm)
  IF (ios > 0) CALL errore('edi', 'error reading edinput_nml namelist', ios)

  CALL mp_bcast(edi_prefix, ionode_id, world_comm)
  CALL mp_bcast(edi_outdir, ionode_id, world_comm)
  CALL mp_bcast(wannierize, ionode_id, world_comm)
  CALL mp_bcast(nbndsub, ionode_id, world_comm)
  CALL mp_bcast(coarse_nk1, ionode_id, world_comm)
  CALL mp_bcast(coarse_nk2, ionode_id, world_comm)
  CALL mp_bcast(coarse_nk3, ionode_id, world_comm)
  CALL mp_bcast(dis_win_min, ionode_id, world_comm)
  CALL mp_bcast(dis_win_max, ionode_id, world_comm)
  CALL mp_bcast(dis_froz_min, ionode_id, world_comm)
  CALL mp_bcast(dis_froz_max, ionode_id, world_comm)
  CALL mp_bcast(num_iter, ionode_id, world_comm)
  CALL mp_bcast(auto_projections, ionode_id, world_comm)
  CALL mp_bcast(scdm_proj, ionode_id, world_comm)
  CALL mp_bcast(bands_skipped, ionode_id, world_comm)
  CALL mp_bcast(filukk, ionode_id, world_comm)
  CALL mp_bcast(do_edmat, ionode_id, world_comm)
  CALL mp_bcast(pristine_prefix, ionode_id, world_comm)
  CALL mp_bcast(pristine_outdir, ionode_id, world_comm)
  CALL mp_bcast(defect_prefix, ionode_id, world_comm)
  CALL mp_bcast(defect_outdir, ionode_id, world_comm)
  CALL mp_bcast(edwwrite, ionode_id, world_comm)
  CALL mp_bcast(edwread, ionode_id, world_comm)
  CALL mp_bcast(edmat_direct_from_file, ionode_id, world_comm)
  CALL mp_bcast(edmat_interp_from_file, ionode_id, world_comm)
  CALL mp_bcast(filki_direct, ionode_id, world_comm)
  CALL mp_bcast(filkf_direct, ionode_id, world_comm)
  CALL mp_bcast(filki_interp, ionode_id, world_comm)
  CALL mp_bcast(filkf_interp, ionode_id, world_comm)
  CALL mp_bcast(sc_nk1, ionode_id, world_comm)
  CALL mp_bcast(sc_nk2, ionode_id, world_comm)
  CALL mp_bcast(sc_nk3, ionode_id, world_comm)
  CALL mp_bcast(band_ed, ionode_id, world_comm)
  CALL mp_bcast(fine_nk1, ionode_id, world_comm)
  CALL mp_bcast(fine_nk2, ionode_id, world_comm)
  CALL mp_bcast(fine_nk3, ionode_id, world_comm)
  CALL mp_bcast(do_transport, ionode_id, world_comm)
  CALL mp_bcast(nstemp, ionode_id, world_comm)
  CALL mp_bcast(temps, ionode_id, world_comm)
  CALL mp_bcast(transport_win_min, ionode_id, world_comm)
  CALL mp_bcast(transport_win_max, ionode_id, world_comm)
  CALL mp_bcast(carrier_conc, ionode_id, world_comm)
  CALL mp_bcast(defect_conc, ionode_id, world_comm)
  CALL mp_bcast(delta_method, ionode_id, world_comm)
  CALL mp_bcast(delta_sigma, ionode_id, world_comm)
  CALL mp_bcast(potfile_d, ionode_id, world_comm)
  CALL mp_bcast(potfile_p, ionode_id, world_comm)
  CALL mp_bcast(pot_align, ionode_id, world_comm)
  CALL mp_bcast(defect_center, ionode_id, world_comm)
  CALL mp_bcast(core_align_radius, ionode_id, world_comm)

  prefix = TRIM(edi_prefix)
  tmp_dir = trimcheck(TRIM(edi_outdir))
  CALL mp_bcast(prefix, ionode_id, world_comm)
  CALL mp_bcast(tmp_dir, ionode_id, world_comm)

  CALL sync_input_module()

  !============================================================
  ! Direct calculation mode: EXCLUSIVE — skips all other modes
  !============================================================
  IF (edmat_direct_from_file) THEN
     wannierize = .FALSE.
     do_edmat = .FALSE.
     edmat_interp_from_file = .FALSE.
     edwread = .FALSE.

     IF (ionode) THEN
        WRITE(stdout, '(/,5X,A)') REPEAT('=', 60)
        WRITE(stdout, '(5X,A)')   'EDI — Direct Calculation Mode'
        WRITE(stdout, '(5X,A)')   REPEAT('=', 60)
        WRITE(stdout, '(5X,A)')   'WARNING: Wannierization, Part A (Bloch->Wannier),'
        WRITE(stdout, '(5X,A)')   '  Part B (Wannier interpolation) are DISABLED.'
        WRITE(stdout, '(5X,A)')   '  Computing M(k_i,k_f) directly from wavefunctions.'
        WRITE(stdout, '(5X,A)')   REPEAT('=', 60)
        WRITE(stdout, '(5X,A,A)') 'prefix  = ', TRIM(prefix)
        WRITE(stdout, '(5X,A,A)') 'outdir  = ', TRIM(tmp_dir)
        WRITE(stdout, '(5X,A,A)') 'filki   = ', TRIM(filki_direct)
        WRITE(stdout, '(5X,A,A)') 'filkf   = ', TRIM(filkf_direct)
        FLUSH(stdout)
     ENDIF

     IF (ionode) WRITE(stdout, '(5X,A)') 'Reading primitive cell data...'
     CALL read_file()
     IF (ionode) THEN
        WRITE(stdout, '(5X,A,I6)') 'nkstot = ', nkstot
        WRITE(stdout, '(5X,A,I6)') 'nbnd   = ', nbnd
        WRITE(stdout, '(5X,A,I6)') 'npol   = ', npol
        IF (lspinorb) THEN
           WRITE(stdout, '(5X,A)') 'SOC: enabled (noncolin=T, lspinorb=T)'
        ELSE IF (noncolin) THEN
           WRITE(stdout, '(5X,A)') 'Noncollinear: enabled (noncolin=T, lspinorb=F)'
        ELSE
           WRITE(stdout, '(5X,A)') 'SOC: disabled (collinear calculation)'
        ENDIF
     ENDIF

     ! Load supercell potentials
     IF (LEN_TRIM(potfile_d) > 0 .AND. LEN_TRIM(potfile_p) > 0) THEN
        ! Use pre-extracted cube files (from extract_pot.x)
        IF (ionode) WRITE(stdout, '(5X,A)') 'Loading potentials from cube files...'
        CALL load_pot_from_file(TRIM(potfile_d), V_d)
        CALL load_pot_from_file(TRIM(potfile_p), V_p)
     ELSE
        ! On-the-fly extraction (MPI_COMM_SELF on ionode)
        IF (ionode) WRITE(stdout, '(5X,A)') 'Loading defect supercell potential...'
        CALL load_supercell_pot(TRIM(defect_prefix), TRIM(defect_outdir), V_d)
        IF (ionode) WRITE(stdout, '(5X,A)') 'Loading pristine supercell potential...'
        CALL load_supercell_pot(TRIM(pristine_prefix), TRIM(pristine_outdir), V_p)
     ENDIF

     ALLOCATE(V_colin(V_d%nr1 * V_d%nr2 * V_d%nr3))
     SELECT CASE (TRIM(pot_align))
     CASE ('vacuum')
        CALL build_vcolin_aligned(V_d, V_p, V_colin, SIZE(V_colin), vac_shift)
     CASE ('core')
        CALL build_vcolin_corealign(V_d, V_p, V_colin, SIZE(V_colin), vac_shift, &
                                     defect_center, core_align_radius)
     CASE ('none')
        V_colin(:) = V_d%pot(:) - V_p%pot(:)
        vac_shift = 0.0_dp
        IF (ionode) WRITE(stdout, '(5X,A)') 'pot_align = none: no alignment applied'
     CASE DEFAULT
        CALL errore('edi', 'Unknown pot_align: '//TRIM(pot_align), 1)
     END SELECT
     IF (ionode) THEN
        WRITE(stdout, '(5X,A)') '--- V_colin DIAGNOSTIC ---'
        WRITE(stdout, '(5X,A,2ES14.6)') 'V_d%pot  min/max = ', MINVAL(V_d%pot), MAXVAL(V_d%pot)
        WRITE(stdout, '(5X,A,2ES14.6)') 'V_p%pot  min/max = ', MINVAL(V_p%pot), MAXVAL(V_p%pot)
        WRITE(stdout, '(5X,A,2ES14.6)') 'V_d-V_p  min/max = ', &
             MINVAL(V_d%pot - V_p%pot), MAXVAL(V_d%pot - V_p%pot)
        WRITE(stdout, '(5X,A,2ES14.6)') 'V_colin  min/max = ', MINVAL(V_colin), MAXVAL(V_colin)
        WRITE(stdout, '(5X,A,I12)') 'SIZE(V_d%pot) = ', SIZE(V_d%pot)
        WRITE(stdout, '(5X,A,I12)') 'SIZE(V_colin) = ', SIZE(V_colin)
        WRITE(stdout, '(5X,A)') '--------------------------'
        ! Check grid ordering: sample specific (ix,iy,iz) points
        ! Vacancy site should be near (120,120,170) in the 240x240x300 grid
        ! Point in vacuum should be near (0,0,0)
        BLOCK
           INTEGER :: nr1_d, nr2_d, nr3_d, idx
           nr1_d = V_d%nr1; nr2_d = V_d%nr2; nr3_d = V_d%nr3
           ! Point at (0,0,0) — vacuum
           idx = 0 + 0*nr1_d + 0*nr1_d*nr2_d + 1
           WRITE(stdout, '(5X,A,I10,A,2F12.6)') 'idx=', idx, &
                '  (0,0,0) V_d/V_p = ', V_d%pot(idx), V_p%pot(idx)
           ! Point at (120,120,0) — vacuum xy-center
           idx = 120 + 120*nr1_d + 0*nr1_d*nr2_d + 1
           WRITE(stdout, '(5X,A,I10,A,2F12.6)') 'idx=', idx, &
                '  (120,120,0) V_d/V_p = ', V_d%pot(idx), V_p%pot(idx)
           ! Point at (120,120,170) — near vacancy
           idx = 120 + 120*nr1_d + 170*nr1_d*nr2_d + 1
           WRITE(stdout, '(5X,A,I10,A,2F12.6)') 'idx=', idx, &
                '  (120,120,170) V_d/V_p = ', V_d%pot(idx), V_p%pot(idx)
           ! Point at (0,0,170) — same z-plane, different xy
           idx = 0 + 0*nr1_d + 170*nr1_d*nr2_d + 1
           WRITE(stdout, '(5X,A,I10,A,2F12.6)') 'idx=', idx, &
                '  (0,0,170) V_d/V_p = ', V_d%pot(idx), V_p%pot(idx)
           ! Point at (120,120,150) — Mo layer
           idx = 120 + 120*nr1_d + 150*nr1_d*nr2_d + 1
           WRITE(stdout, '(5X,A,I10,A,2F12.6)') 'idx=', idx, &
                '  (120,120,150) V_d/V_p = ', V_d%pot(idx), V_p%pot(idx)
           ! Check where V_d-V_p is largest
           idx = MAXLOC(V_d%pot - V_p%pot, DIM=1)
           WRITE(stdout, '(5X,A,I10,A,3I6)') 'maxdiff at idx=', idx, &
                '  (ix,iy,iz) = ', MOD(idx-1, nr1_d), MOD((idx-1)/nr1_d, nr2_d), &
                (idx-1)/(nr1_d*nr2_d)
           WRITE(stdout, '(5X,A,3F14.6)') '  V_d, V_p, diff = ', &
                V_d%pot(idx), V_p%pot(idx), V_d%pot(idx) - V_p%pot(idx)
        END BLOCK
        WRITE(stdout, '(5X,A)') '--------------------------'
        FLUSH(stdout)
     ENDIF
     CALL write_vcolin_cube(V_d, V_colin, SIZE(V_colin), TRIM(edi_prefix))

     ! Restore primitive cell state (same as Part A)
     IF (ionode) WRITE(stdout, '(5X,A)') 'Restoring primitive cell state...'
     CALL clean_pw(.TRUE.)
     prefix = TRIM(edi_prefix)
     tmp_dir = trimcheck(TRIM(edi_outdir))
     need_wf = .TRUE.
     CALL read_file_new(need_wf)
     IF (ionode) WRITE(stdout, '(5X,A)') 'Primitive cell restored.'

     ! Compute M directly
     CALL ed_direct_from_files(TRIM(filki_direct), TRIM(filkf_direct), &
          TRIM(edi_prefix), TRIM(band_ed))

     CALL environment_end('EDI')
     CALL mp_global_end()
     STOP
  ENDIF

  IF (ionode) THEN
     WRITE(stdout, '(/,5X,A)') REPEAT('=', 60)
     WRITE(stdout, '(5X,A)')   'EDI — Electron-Defect Interaction'
     WRITE(stdout, '(5X,A)')   'Wannierization mode'
     WRITE(stdout, '(5X,A)')   REPEAT('=', 60)
     WRITE(stdout, '(5X,A,A)') 'prefix  = ', TRIM(prefix)
     WRITE(stdout, '(5X,A,A)') 'outdir  = ', TRIM(tmp_dir)
     WRITE(stdout, '(5X,A,I6)') 'nbndsub = ', nbndsub
     WRITE(stdout, '(5X,A,3I4)') 'k-grid  = ', coarse_nk1, coarse_nk2, coarse_nk3
     WRITE(stdout, '(5X,A,I4)') 'npool   = ', npool
  ENDIF

  IF (ionode) WRITE(stdout, '(5X,A)') 'Reading NSCF save data...'
  CALL read_file()
  IF (ionode) THEN
     WRITE(stdout, '(5X,A,I6)') 'nkstot = ', nkstot
     WRITE(stdout, '(5X,A,I6)') 'nks    = ', nks
     WRITE(stdout, '(5X,A,I6)') 'nbnd   = ', nbnd
  ENDIF

  IF (coarse_nk1 * coarse_nk2 * coarse_nk3 /= nkstot) THEN
     CALL errore('edi', 'k-grid dimensions do not match nkstot', nkstot)
  ENDIF

  nbndep = nbnd
  nbndskip = 0

  ALLOCATE(xk_all(3, nkstot), STAT=ierr)
  IF (ierr /= 0) CALL errore('edi', 'Error allocating xk_all', 1)
  CALL poolcollect(3, nks, xk, nkstot, xk_all)

  ALLOCATE(xk_cryst(3, nkstot), STAT=ierr)
  IF (ierr /= 0) CALL errore('edi', 'Error allocating xk_cryst', 1)
  xk_cryst(:, :) = xk_all(:, :)
  CALL cryst_to_cart(nkstot, xk_cryst, at, -1)
  CALL mp_bcast(xk_cryst, ionode_id, world_comm)

  ALLOCATE(xk_loc(3, nks), STAT=ierr)
  IF (ierr /= 0) CALL errore('edi', 'Error allocating xk_loc', 1)
  xk_loc(:, 1:nks) = xk(:, 1:nks)

  CALL openfilepw()

  IF (wannierize) THEN
     CALL edi_run_wannier90(TRIM(prefix), coarse_nk1, coarse_nk2, coarse_nk3)
     ! After Wannier90, read back the combined rotation u_kc from filukk.
     ! This ensures cu = u_kc (already combined) with correct dimensions,
     ! avoiding basis mismatch between u_mat_opt/lwindow and edmatkq indices.
     IF (ionode) WRITE(stdout, '(5X,A)') 'Reading combined rotation from filukk...'
     CALL read_filukk_edi(filukk, nbnd, nkstot, nks, nbndsub)
  ELSE
     IF (ionode) WRITE(stdout, '(5X,A)') 'wannierize = .FALSE., reading Wannier data from filukk'
     CALL read_filukk_edi(filukk, nbnd, nkstot, nks, nbndsub)
  ENDIF

  CALL read_hr_file(TRIM(prefix), nbndsub_loc, nrr, ndegen, irvec, chw)

  IF (ionode) THEN
     WRITE(stdout, '(/,5X,A)') 'Validating: Wannier interpolation vs NSCF eigenvalues'
     WRITE(stdout, '(5X,A)') REPEAT('-', 56)
  ENDIF

  ALLOCATE(xk_fine(3, nkstot))
  xk_fine(:,:) = xk_cryst(:,:)

  ALLOCATE(eig_interp(nbndsub_loc, nkstot))
  CALL edi_interp_bands(nbndsub_loc, nrr, irvec, ndegen, chw, &
                         nkstot, xk_fine, eig_interp)

  IF (ionode) THEN
     ALLOCATE(et_kept(nbndsub_loc, nkstot))
     DO ik = 1, nkstot
        ib_kept = 0
        DO ib = 1, nbnd
           IF (excluded_band(ib)) CYCLE
           ib_kept = ib_kept + 1
           IF (ib_kept <= nbndsub_loc) et_kept(ib_kept, ik) = et(ib, ik) * rytoev
        ENDDO
     ENDDO

     WRITE(stdout, '(5X,A4,A5,A16,A16,A16)') &
          'ik', 'ib', 'E_nscf(eV)', 'E_w90(eV)', 'Error(eV)'

     max_err = 0.0_dp
     DO ik = 1, nkstot
        DO ibnd = 1, nbndsub_loc
           err = ABS(eig_interp(ibnd, ik) - et_kept(ibnd, ik))
           IF (err > max_err) max_err = err
           IF (ik <= 3 .OR. ik == nkstot) THEN
              WRITE(stdout, '(5X,I4,I5,3F16.6)') ik, ibnd, &
                   et_kept(ibnd, ik), eig_interp(ibnd, ik), err
           ENDIF
        ENDDO
     ENDDO

     WRITE(stdout, '(/,5X,A,ES12.4,A)') 'Max error: ', max_err, ' eV'
     IF (max_err < 1.0d-4) THEN
        WRITE(stdout, '(5X,A)') 'PASSED: Wannier interpolation < 1e-4 eV'
     ELSE
        WRITE(stdout, '(5X,A)') 'Note: error exceeds 1e-4 eV'
     ENDIF
     WRITE(stdout, '(5X,A,/)') REPEAT('=', 60)

     DEALLOCATE(et_kept)
  ENDIF

  IF (ALLOCATED(eig_interp)) DEALLOCATE(eig_interp)
  IF (ALLOCATED(xk_fine)) DEALLOCATE(xk_fine)

  !============================================================
  ! Electron-defect matrix element interpolation
  !============================================================
  IF (do_edmat) THEN

     IF (edwread) THEN
        !=========================================================
        ! edwread = .true.: Read M(R,Rp) from binary, skip Part A
        !=========================================================
        IF (ionode) THEN
           WRITE(stdout, '(/,5X,A)') REPEAT('-', 60)
           WRITE(stdout, '(5X,A)') 'edwread = .TRUE. : reading M(R,Rp) from binary file'
           WRITE(stdout, '(5X,A)') REPEAT('-', 60)
           FLUSH(stdout)
        ENDIF

     ELSE
        !=========================================================
        ! Part A: Load potentials and compute M(R,Rp)
        !=========================================================
        IF (LEN_TRIM(potfile_d) > 0 .AND. LEN_TRIM(potfile_p) > 0) THEN
           IF (ionode) WRITE(stdout, '(5X,A)') 'Loading potentials from cube files...'
           CALL load_pot_from_file(TRIM(potfile_d), V_d)
           CALL load_pot_from_file(TRIM(potfile_p), V_p)
        ELSE
           IF (ionode) WRITE(stdout, '(5X,A)') 'Loading defect supercell potential...'
           CALL load_supercell_pot(TRIM(defect_prefix), TRIM(defect_outdir), V_d)
           IF (ionode) WRITE(stdout, '(5X,A)') 'Loading pristine supercell potential...'
           CALL load_supercell_pot(TRIM(pristine_prefix), TRIM(pristine_outdir), V_p)
        ENDIF

        ALLOCATE(V_colin(V_d%nr1 * V_d%nr2 * V_d%nr3))
        SELECT CASE (TRIM(pot_align))
        CASE ('vacuum')
           CALL build_vcolin_aligned(V_d, V_p, V_colin, SIZE(V_colin), vac_shift)
        CASE ('core')
           CALL build_vcolin_corealign(V_d, V_p, V_colin, SIZE(V_colin), vac_shift, &
                                        defect_center, core_align_radius)
        CASE ('none')
           V_colin(:) = V_d%pot(:) - V_p%pot(:)
           vac_shift = 0.0_dp
           IF (ionode) WRITE(stdout, '(5X,A)') 'pot_align = none: no alignment applied'
        CASE DEFAULT
           CALL errore('edi', 'Unknown pot_align: '//TRIM(pot_align), 1)
        END SELECT
        CALL write_vcolin_cube(V_d, V_colin, SIZE(V_colin), TRIM(edi_prefix))

        IF (ionode) WRITE(stdout, '(5X,A)') 'Restoring primitive cell state...'
        CALL clean_pw(.TRUE.)
        prefix = TRIM(edi_prefix)
        tmp_dir = trimcheck(TRIM(edi_outdir))
        need_wf = .TRUE.
        CALL read_file_new(need_wf)

        BLOCK
           COMPLEX(dp), ALLOCATABLE :: edmatw_2d(:,:,:,:)
           CALL ed_coarse_full_q(nbndsub_loc, nrr, irvec, ndegen, &
                                  edmatw_2d, TRIM(edi_prefix))
           IF (ALLOCATED(edmatw_2d)) DEALLOCATE(edmatw_2d)
        END BLOCK
        IF (ALLOCATED(V_colin)) DEALLOCATE(V_colin)
     ENDIF

     !=========================================================
     ! Part B / Transport: read _edmatw_2d.bin and either
     !   - do_transport: run transport (Fermi-filtered, no file output)
     !   - else: run Part B interpolation (full grid, file output)
     !=========================================================
     BLOCK
        USE transport_edi, ONLY : compute_transport, compute_dos_validation
        INTEGER :: nbndsub_2d, nrr_2d, nk1_use, nk2_use, nk3_use
        INTEGER, ALLOCATABLE :: ndegen_2d(:), irvec_2d(:,:)
        COMPLEX(dp), ALLOCATABLE :: edmatw_2d_r(:,:,:,:)
        LOGICAL :: file_exists

        IF (fine_nk1 > 0 .AND. fine_nk2 > 0 .AND. fine_nk3 > 0) THEN
           nk1_use = fine_nk1; nk2_use = fine_nk2; nk3_use = fine_nk3
        ELSE
           nk1_use = coarse_nk1; nk2_use = coarse_nk2; nk3_use = coarse_nk3
        ENDIF

        IF (ionode) INQUIRE(FILE=TRIM(edi_prefix)//'_edmatw_2d.bin', EXIST=file_exists)
        CALL mp_bcast(file_exists, ionode_id, world_comm)
        IF (file_exists) THEN
           CALL read_edmatw_2d_file(TRIM(edi_prefix), nbndsub_2d, nrr_2d, &
                                     ndegen_2d, irvec_2d, edmatw_2d_r)

           IF (do_transport) THEN
              ! DOS validation first
              CALL compute_dos_validation(nbndsub_2d, nrr_2d, irvec_2d, ndegen_2d, chw, &
                   nk1_use, nk2_use, nk3_use, TRIM(edi_prefix))

              ! Transport: Fermi-filtered interpolation → scattering rates → mobility
              IF (ionode) THEN
                 WRITE(stdout, '(/,5X,A)') REPEAT('-', 60)
                 WRITE(stdout, '(5X,A)') 'Transport mode: skipping Part B file output'
                 WRITE(stdout, '(5X,A)') '  Computing scattering rates with Fermi filtering'
                 WRITE(stdout, '(5X,A)') REPEAT('-', 60)
                 FLUSH(stdout)
              ENDIF
              CALL compute_transport(nbndsub_2d, nrr_2d, irvec_2d, ndegen_2d, chw, &
                   edmatw_2d_r, nk1_use, nk2_use, nk3_use, nstemp, temps, &
                   defect_conc, TRIM(edi_prefix))
           ELSE
              ! Part B: full-grid interpolation → output files
              IF (ionode) THEN
                 WRITE(stdout, '(/,5X,A)') REPEAT('-', 60)
                 WRITE(stdout, '(5X,A)') 'Part B: interpolating M(R,Rp) on fine grid'
                 WRITE(stdout, '(5X,A)') REPEAT('-', 60)
                 FLUSH(stdout)
              ENDIF
              CALL ed_fine_interp_2d(nbndsub_2d, nrr_2d, irvec_2d, ndegen_2d, chw, &
                   edmatw_2d_r, nk1_use, nk2_use, nk3_use, TRIM(edi_prefix))
           ENDIF

           IF (ALLOCATED(edmatw_2d_r)) DEALLOCATE(edmatw_2d_r)
           IF (ALLOCATED(ndegen_2d)) DEALLOCATE(ndegen_2d)
           IF (ALLOCATED(irvec_2d)) DEALLOCATE(irvec_2d)
        ELSE
           IF (ionode) THEN
              WRITE(stdout, '(5X,A)') 'ERROR: _edmatw_2d.bin not found'
              WRITE(stdout, '(5X,A)') '  Run first with edwread=.false. to generate it'
           ENDIF
        ENDIF
     END BLOCK
  ENDIF

  !============================================================
  ! Section 3.D: Custom k-point interpolation from files
  !============================================================
  IF (edmat_interp_from_file .AND. filki_interp /= ' ' .AND. filkf_interp /= ' ') THEN
     BLOCK
        INTEGER :: nbndsub_2d, nrr_2d
        INTEGER, ALLOCATABLE :: ndegen_2d(:), irvec_2d(:,:)
        COMPLEX(dp), ALLOCATABLE :: edmatw_2d_f(:,:,:,:)
        LOGICAL :: file_2d_exists

        IF (.NOT. ALLOCATED(chw)) THEN
           CALL read_hr_file(TRIM(edi_prefix), nbndsub_2d, nrr_2d, ndegen_2d, irvec_2d, chw)
        ENDIF

        IF (ionode) INQUIRE(FILE=TRIM(edi_prefix)//'_edmatw_2d.bin', EXIST=file_2d_exists)
        CALL mp_bcast(file_2d_exists, ionode_id, world_comm)
        IF (file_2d_exists) THEN
           CALL read_edmatw_2d_file(TRIM(edi_prefix), nbndsub_2d, nrr_2d, &
                                     ndegen_2d, irvec_2d, edmatw_2d_f)
           CALL ed_interp_from_file(nbndsub_2d, nrr_2d, irvec_2d, ndegen_2d, chw, &
                                     edmatw_2d_f, &
                                     TRIM(filki_interp), TRIM(filkf_interp), TRIM(edi_prefix))
           DEALLOCATE(edmatw_2d_f)
           IF (ALLOCATED(ndegen_2d)) DEALLOCATE(ndegen_2d)
           IF (ALLOCATED(irvec_2d)) DEALLOCATE(irvec_2d)
        ELSE
           IF (ionode) THEN
              WRITE(stdout, '(5X,A)') 'ERROR: edmat_interp_from_file requires _edmatw_2d.bin'
              WRITE(stdout, '(5X,A)') '  Run first with edwread=.false. to generate it'
           ENDIF
        ENDIF
     END BLOCK
  ENDIF

  IF (ALLOCATED(irvec)) DEALLOCATE(irvec)
  IF (ALLOCATED(ndegen)) DEALLOCATE(ndegen)
  IF (ALLOCATED(chw)) DEALLOCATE(chw)
  IF (ALLOCATED(xk_all)) DEALLOCATE(xk_all)
  IF (ALLOCATED(xk_cryst)) DEALLOCATE(xk_cryst)
  IF (ALLOCATED(xk_loc)) DEALLOCATE(xk_loc)
  IF (ALLOCATED(ibndkept)) DEALLOCATE(ibndkept)

  CALL environment_end('EDI')
  CALL mp_global_end()

End Program edi
