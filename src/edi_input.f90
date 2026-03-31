MODULE edi_input
  USE kinds, ONLY : dp
  IMPLICIT NONE
  SAVE

  INTEGER, PARAMETER :: nwanxx = 200

  ! Wannier90 control
  LOGICAL :: wannierize = .FALSE.
  INTEGER :: nbndsub = 0
  REAL(dp) :: dis_win_min = -9999.0_dp
  REAL(dp) :: dis_win_max =  9999.0_dp
  REAL(dp) :: dis_froz_min = -9999.0_dp
  REAL(dp) :: dis_froz_max =  9999.0_dp
  INTEGER :: num_iter = 200
  INTEGER :: iprint_w90 = 2
  CHARACTER(LEN=256) :: w90_seedname = 'wannier90'
  CHARACTER(LEN=256) :: filukk = 'filukk'
  LOGICAL :: auto_projections = .TRUE.
  LOGICAL :: scdm_proj = .FALSE.
  CHARACTER(LEN=256) :: proj(nwanxx) = ' '
  CHARACTER(LEN=256) :: wdata(nwanxx) = ' '
  CHARACTER(LEN=256) :: bands_skipped = ' '

  ! Transport control
  LOGICAL :: do_transport = .FALSE.
  INTEGER :: nstemp = 1
  REAL(dp) :: temps(100) = 300.0_dp
  REAL(dp) :: transport_win_min = -9999.0_dp  ! energy window lower bound (eV)
  REAL(dp) :: transport_win_max =  9999.0_dp  ! energy window upper bound (eV)
  REAL(dp) :: carrier_conc = 1.0d13   ! carrier concentration (cm^-2 for 2D, cm^-3 for 3D)
  REAL(dp) :: defect_conc = 1.0d12    ! defect concentration (cm^-2 for 2D, cm^-3 for 3D)
  CHARACTER(LEN=32) :: delta_method = 'triangular'
  REAL(dp) :: delta_sigma = 0.1_dp    ! Gaussian sigma (eV), for delta_method='gaussian'
  LOGICAL :: iterative_bte = .FALSE.
  INTEGER :: bte_max_iter = 100
  REAL(dp) :: bte_conv_thr = 1.0d-4

  ! Wannier interpolation of M
  LOGICAL :: wannier_interp_m = .FALSE.
  REAL(dp) :: m_interp_rcut = -1.0_dp

  ! k-grid dimensions
  INTEGER :: coarse_nk1 = 0
  INTEGER :: coarse_nk2 = 0
  INTEGER :: coarse_nk3 = 0
  INTEGER :: fine_nk1 = 0
  INTEGER :: fine_nk2 = 0
  INTEGER :: fine_nk3 = 0

  CHARACTER(LEN=256) :: edi_prefix = 'pwscf'
  CHARACTER(LEN=256) :: edi_outdir = './'

  ! Control for reading/writing Wannier-basis matrix elements
  LOGICAL :: edwwrite = .TRUE.   ! write M(R) to prefix_edmatw.dat after Part A
  LOGICAL :: edwread  = .FALSE.  ! read M(R) from prefix_edmatw.dat, skip Part A

  INTEGER :: sc_nk1 = 0, sc_nk2 = 0, sc_nk3 = 0  ! supercell dimensions (legacy, unused)

  ! Band selection for direct calculation mode (e.g. '13-17')
  CHARACTER(LEN=256) :: band_ed = ' '

  ! Custom k-point mode: direct calculation or interpolation from file
  LOGICAL :: edmat_direct_from_file = .FALSE.
  LOGICAL :: edmat_interp_from_file = .FALSE.
  CHARACTER(LEN=256) :: filki_direct = ' '
  CHARACTER(LEN=256) :: filkf_direct = ' '
  CHARACTER(LEN=256) :: filki_interp = ' '
  CHARACTER(LEN=256) :: filkf_interp = ' '

  ! Supercell parameters for electron-defect matrix elements
  LOGICAL :: do_edmat = .FALSE.
  CHARACTER(LEN=256) :: pristine_prefix = ' '
  CHARACTER(LEN=256) :: pristine_outdir = ' '
  CHARACTER(LEN=256) :: defect_prefix = ' '
  CHARACTER(LEN=256) :: defect_outdir = ' '

  ! Pre-extracted potential cube files (from extract_pot.x)
  ! If set, these are used instead of on-the-fly loading
  CHARACTER(LEN=256) :: potfile_d = ' '
  CHARACTER(LEN=256) :: potfile_p = ' '

  ! Potential alignment: 'vacuum' (2D), 'core' (2D/3D), 'none'
  CHARACTER(LEN=16) :: pot_align = 'vacuum'
  ! Defect center in fractional coordinates of the supercell (required for 'core')
  REAL(dp) :: defect_center(3) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
  ! Averaging sphere radius in Bohr for core alignment
  REAL(dp) :: core_align_radius = 2.0_dp

  NAMELIST / edinput / &
       wannierize, nbndsub, dis_win_min, dis_win_max, &
       dis_froz_min, dis_froz_max, num_iter, iprint_w90, &
       w90_seedname, filukk, auto_projections, scdm_proj, &
       proj, wdata, bands_skipped, &
       do_transport, nstemp, temps, transport_win_min, transport_win_max, &
       carrier_conc, defect_conc, delta_method, delta_sigma, &
       iterative_bte, bte_max_iter, bte_conv_thr, &
       wannier_interp_m, m_interp_rcut, &
       coarse_nk1, coarse_nk2, coarse_nk3, &
       fine_nk1, fine_nk2, fine_nk3, &
       do_edmat, pristine_prefix, pristine_outdir, &
       defect_prefix, defect_outdir, &
       potfile_d, potfile_p, &
       pot_align, defect_center, core_align_radius, &
       edwwrite, edwread, &
       edmat_direct_from_file, edmat_interp_from_file, &
       filki_direct, filkf_direct, filki_interp, filkf_interp, &
       sc_nk1, sc_nk2, sc_nk3, &
       band_ed

  NAMELIST / edinput_nml / &
       edi_prefix, edi_outdir, &
       wannierize, nbndsub, dis_win_min, dis_win_max, &
       dis_froz_min, dis_froz_max, num_iter, iprint_w90, &
       w90_seedname, filukk, auto_projections, scdm_proj, &
       proj, wdata, bands_skipped, &
       do_transport, nstemp, temps, transport_win_min, transport_win_max, &
       carrier_conc, defect_conc, delta_method, delta_sigma, &
       iterative_bte, bte_max_iter, bte_conv_thr, &
       wannier_interp_m, m_interp_rcut, &
       coarse_nk1, coarse_nk2, coarse_nk3, &
       fine_nk1, fine_nk2, fine_nk3, &
       do_edmat, pristine_prefix, pristine_outdir, &
       defect_prefix, defect_outdir, &
       potfile_d, potfile_p, &
       pot_align, defect_center, core_align_radius, &
       edwwrite, edwread, &
       edmat_direct_from_file, edmat_interp_from_file, &
       filki_direct, filkf_direct, filki_interp, filkf_interp, &
       sc_nk1, sc_nk2, sc_nk3, &
       band_ed

CONTAINS

  SUBROUTINE read_edinput(ios)
    USE io_global, ONLY : ionode
    USE mp, ONLY : mp_bcast
    USE mp_images, ONLY : intra_image_comm
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: ios

    ios = 0

    IF (ionode) THEN
       READ(5, edinput, IOSTAT=ios)
       IF (ios < 0) ios = 0
    ENDIF

    CALL mp_bcast(ios, 0, intra_image_comm)
    CALL mp_bcast(wannierize, 0, intra_image_comm)
    CALL mp_bcast(nbndsub, 0, intra_image_comm)
    CALL mp_bcast(dis_win_min, 0, intra_image_comm)
    CALL mp_bcast(dis_win_max, 0, intra_image_comm)
    CALL mp_bcast(dis_froz_min, 0, intra_image_comm)
    CALL mp_bcast(dis_froz_max, 0, intra_image_comm)
    CALL mp_bcast(num_iter, 0, intra_image_comm)
    CALL mp_bcast(iprint_w90, 0, intra_image_comm)
    CALL mp_bcast(w90_seedname, 0, intra_image_comm)
    CALL mp_bcast(filukk, 0, intra_image_comm)
    CALL mp_bcast(auto_projections, 0, intra_image_comm)
    CALL mp_bcast(scdm_proj, 0, intra_image_comm)
    CALL mp_bcast(bands_skipped, 0, intra_image_comm)
    CALL mp_bcast(do_transport, 0, intra_image_comm)
    CALL mp_bcast(nstemp, 0, intra_image_comm)
    CALL mp_bcast(temps, 0, intra_image_comm)
    CALL mp_bcast(transport_win_min, 0, intra_image_comm)
    CALL mp_bcast(transport_win_max, 0, intra_image_comm)
    CALL mp_bcast(carrier_conc, 0, intra_image_comm)
    CALL mp_bcast(defect_conc, 0, intra_image_comm)
    CALL mp_bcast(iterative_bte, 0, intra_image_comm)
    CALL mp_bcast(bte_max_iter, 0, intra_image_comm)
    CALL mp_bcast(bte_conv_thr, 0, intra_image_comm)
    CALL mp_bcast(wannier_interp_m, 0, intra_image_comm)
    CALL mp_bcast(m_interp_rcut, 0, intra_image_comm)
    CALL mp_bcast(coarse_nk1, 0, intra_image_comm)
    CALL mp_bcast(coarse_nk2, 0, intra_image_comm)
    CALL mp_bcast(coarse_nk3, 0, intra_image_comm)
    CALL mp_bcast(fine_nk1, 0, intra_image_comm)
    CALL mp_bcast(fine_nk2, 0, intra_image_comm)
    CALL mp_bcast(fine_nk3, 0, intra_image_comm)
    CALL mp_bcast(do_edmat, 0, intra_image_comm)
    CALL mp_bcast(pristine_prefix, 0, intra_image_comm)
    CALL mp_bcast(pristine_outdir, 0, intra_image_comm)
    CALL mp_bcast(defect_prefix, 0, intra_image_comm)
    CALL mp_bcast(defect_outdir, 0, intra_image_comm)
    CALL mp_bcast(pot_align, 0, intra_image_comm)
    CALL mp_bcast(potfile_d, 0, intra_image_comm)
    CALL mp_bcast(potfile_p, 0, intra_image_comm)
    CALL mp_bcast(defect_center, 0, intra_image_comm)
    CALL mp_bcast(core_align_radius, 0, intra_image_comm)
  END SUBROUTINE read_edinput

  SUBROUTINE sync_input_module()
    USE input, ONLY : nbndsub_inp => nbndsub, proj_inp => proj, &
                       wdata_inp => wdata, bands_skipped_inp => bands_skipped, &
                       dis_win_min_inp => dis_win_min, dis_win_max_inp => dis_win_max, &
                       dis_froz_min_inp => dis_froz_min, dis_froz_max_inp => dis_froz_max, &
                       num_iter_inp => num_iter, iprint_inp => iprint, &
                       auto_projections_inp => auto_projections, &
                       scdm_proj_inp => scdm_proj, filukk_inp => filukk, &
                       nkc1, nkc2, nkc3
    IMPLICIT NONE
    nbndsub_inp = nbndsub
    proj_inp(1:nwanxx) = proj(1:nwanxx)
    wdata_inp(1:nwanxx) = wdata(1:nwanxx)
    bands_skipped_inp = bands_skipped
    dis_win_min_inp = dis_win_min
    dis_win_max_inp = dis_win_max
    dis_froz_min_inp = dis_froz_min
    dis_froz_max_inp = dis_froz_max
    num_iter_inp = num_iter
    iprint_inp = iprint_w90
    auto_projections_inp = auto_projections
    scdm_proj_inp = scdm_proj
    filukk_inp = filukk
    nkc1 = coarse_nk1
    nkc2 = coarse_nk2
    nkc3 = coarse_nk3
  END SUBROUTINE sync_input_module

  SUBROUTINE write_winfil_edi(prefix_in)
    USE io_global, ONLY : ionode
    USE wvfct, ONLY : nbnd, et
    USE klist, ONLY : nks, nkstot
    USE constants, ONLY : rytoev
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: prefix_in

    INTEGER :: i, iunit
    LOGICAL :: random, notfound
    CHARACTER(LEN=255) :: dummy, lvalue
    INTEGER :: pos, pos1, pos2

    IF (.NOT. ionode) RETURN

    iunit = 97
    OPEN(UNIT=iunit, FILE=TRIM(prefix_in)//'.win', FORM='formatted')

    IF (auto_projections) THEN
       WRITE(iunit, '(a)') 'auto_projections = .true.'
    ELSE
       WRITE(iunit, '(a)') 'begin projections'
       random = .TRUE.
       DO i = 1, nbndsub
          IF (proj(i) /= ' ') THEN
             WRITE(iunit, *) TRIM(proj(i))
             random = .FALSE.
          ENDIF
       ENDDO
       IF (random) WRITE(iunit, '(a)') 'random'
       WRITE(iunit, '(a)') 'end projections'
    ENDIF

    IF (bands_skipped /= ' ') WRITE(iunit, *) bands_skipped

    WRITE(iunit, '("num_wann = ", i3)') nbndsub
    WRITE(iunit, '("iprint = ", i3)') iprint_w90

    IF (dis_win_min > -9000.0_dp) WRITE(iunit, '("dis_win_min = ", f18.12)') dis_win_min
    IF (dis_win_max <  9000.0_dp) WRITE(iunit, '("dis_win_max = ", f18.12)') dis_win_max
    IF (dis_froz_min > -9000.0_dp) WRITE(iunit, '("dis_froz_min = ", f18.12)') dis_froz_min
    IF (dis_froz_max <  9000.0_dp) WRITE(iunit, '("dis_froz_max = ", f18.12)') dis_froz_max
    WRITE(iunit, '("num_iter = ", i7)') num_iter
    WRITE(iunit, '(a)') 'write_bvec = .true.'

    notfound = .TRUE.
    DO i = 1, nwanxx
       IF (wdata(i) /= ' ') THEN
          pos = INDEX(TRIM(ADJUSTL(wdata(i))), 'write_hr')
          IF (pos == 1) notfound = .FALSE.
       ENDIF
    ENDDO

    DO i = 1, nwanxx
       IF (wdata(i) /= ' ') WRITE(iunit, *) TRIM(wdata(i))
    ENDDO

    IF (notfound) WRITE(iunit, *) 'write_hr = .true.'

    CLOSE(iunit)
  END SUBROUTINE write_winfil_edi

END MODULE edi_input
