Subroutine initialization()
  !-----------------------------------------------------------------------
  !
  ! See files INPUT_BANDS.* in Doc/ directory for usage
  !
  USE control_flags,    ONLY : io_level
  USE io_files,  ONLY : prefix, tmp_dir, nwordwfc, iunwfc, restart_dir
  USE mp_global, ONLY : mp_startup
  USE mp_pools,  ONLY : npool
  USE control_flags, ONLY : gamma_only
  USE environment,   ONLY : environment_start, environment_end
  USE wvfct,     ONLY : nbnd, npwx
  USE klist,     ONLY : nkstot, two_fermi_energies, nks
  USE noncollin_module, ONLY : noncolin, i_cons, npol
  USE lsda_mod,  ONLY : nspin
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp,        ONLY : mp_bcast
  USE mp_images, ONLY : intra_image_comm
  USE parameters,ONLY : npk
  USE wavefunctions,    ONLY : evc

  USE pw_restart_new, ONLY : read_xml_file

  USE pw_restart_new,   ONLY : read_collected_wfc
  !
  USE wavefunctions_gpum, ONLY : using_evc
  USE buffers,        ONLY : open_buffer, close_buffer, save_buffer
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  CHARACTER (len=256) :: filband, filp, outdir
  LOGICAL :: lsigma(4), lsym, lp, no_overlap, plot_2d, wfc_is_collected, exst
  INTEGER :: spin_component, firstk, lastk
  INTEGER :: ios
  integer :: ik
  !
  NAMELIST / bands / outdir, prefix, filband, filp, spin_component, lsigma,&
                       lsym, lp, filp, firstk, lastk, no_overlap, plot_2d
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'BANDS' )
  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  filband = 'bands.out'
  lsym=.true.
  no_overlap=.true.
  plot_2d=.false.
  lsigma=.false.
  lp=.false.
  filp='p_avg.dat'
  firstk=0
  lastk=npk
  spin_component = 1
  !
  ios = 0
  !
  IF ( ionode )  THEN
     !
     CALL input_from_file ( )
     !
     READ (5, bands, iostat = ios)
     !
     lsigma(4)=.false.
     tmp_dir = trimcheck (outdir)
     !
  ENDIF

  wfc_is_collected=.true.
  CALL read_file_new(wfc_is_collected)
  write(*,1003) 'Number of Bands', nbnd
  write(*,1003) 'Number of Plane Waves', npwx
  write(*,1003) 'Number of K-points', nks
  1003 format(A24," = ",I6)
  nwordwfc = nbnd*npwx*npol
  io_level = 1
  
 

  
END Subroutine initialization

