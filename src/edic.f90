Program edic
  USE environment,   ONLY : environment_start, environment_end

      Use kinds,    only: dp
      USE io_files,  ONLY : prefix, tmp_dir, nwordwfc, iunwfc, restart_dir
      Use edic_mod,   only: V_file, V_loc, V_0, Bxc_1, Bxc_2, Bxc_3, V_p
      Use edic_mod
      Use wvfct, ONLY: npwx, nbnd, wg, et, g2kin
      Use fft_base,  ONLY: dfftp, dffts
      Use fft_interfaces, ONLY : fwfft, invfft
      Use edic_mod, Only : evc1,evc2,evc3,evc4,&
                               psic1, psic2, psic3, psic4
      USE pw_restart_new,   ONLY : read_collected_wfc
  !    Use klist,  Only : nks
      Use edic_mod, Only: m_loc, m_nloc

  USE control_flags,    ONLY : io_level
  USE mp_global, ONLY : mp_startup
  USE mp_pools,  ONLY : npool
  USE control_flags, ONLY : gamma_only
  USE environment,   ONLY : environment_start, environment_end
  USE klist,     ONLY : nkstot, two_fermi_energies, nks
  USE noncollin_module, ONLY : noncolin, i_cons, npol
  USE lsda_mod,  ONLY : nspin
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp,        ONLY : mp_bcast
  USE mp_images, ONLY : intra_image_comm
  USE parameters,ONLY : npk
  USE wavefunctions,    ONLY : evc

  USE pw_restart_new, ONLY : read_xml_file

  !
  USE wavefunctions_gpum, ONLY : using_evc
  USE buffers,        ONLY : open_buffer, close_buffer, save_buffer
 
      Implicit none
      
      
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
      integer :: tmp_unit
      integer, external :: find_free_unit
      integer :: ios
      integer :: ikk, ikk0, ibnd, ibnd0, ik, ik0, nk
  CHARACTER (len=256) :: filband, filp, outdir
  LOGICAL :: lsigma(4), lsym, lp, no_overlap, plot_2d, wfc_is_collected, exst
  INTEGER :: spin_component, firstk, lastk
  NAMELIST / bands / outdir, prefix, filband, filp, spin_component, lsigma,&
                       lsym, lp, filp, firstk, lastk, no_overlap, plot_2d

  CALL mp_startup ( )
write(*,*) '0start '
  CALL environment_start ( 'BANDS' )
write(*,*) '1start '
  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
!  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
!  IF ( trim( outdir ) == ' ' ) outdir = './'
  !filband = 'bands.out'
  !lsym=.true.
  !no_overlap=.true.
  !plot_2d=.false.
  !lsigma=.false.
  !lp=.false.
  !filp='p_avg.dat'
  !firstk=0
  !lastk=npk
  !spin_component = 1
  !
  ios = 0
  !
!  IF ( ionode )  THEN
!     !
write(*,*) '0'
     CALL input_from_file ( )
write(*,*) '1'
!     !
     READ (5, bands, iostat = ios)
write(*,*) '2',prefix,outdir
!     !
!     lsigma(4)=.false.
     tmp_dir = trimcheck (outdir)
!     !
!  ENDIF


  wfc_is_collected=.true.
  CALL read_file_new(wfc_is_collected)
  write(*,1003) 'Number of Bands', nbnd
  write(*,1003) 'Number of Plane Waves', npwx
  write(*,1003) 'Number of K-points', nks
  1000 format(A24," = ",I6)
  nwordwfc = nbnd*npwx*npol
  io_level = 1
  

      tmp_unit = find_free_unit()

      OPEN(unit=tmp_unit,file = 'calcmdefect.dat',status='old',err=20)
      20 continue
         READ(tmp_unit,calcmcontrol,iostat=ios)
      CLOSE(tmp_unit)
call getvloc()

call getepsdat()

!       V_0%filename = V_0_filename
!      Bxc_1%filename = Bxc_1_filename
!      Bxc_2%filename = Bxc_2_filename
!      Bxc_3%filename = Bxc_3_filename
!      V_p%filename = V_p_filename
!      
!
!      call read_perturb_file(V_0)
!      
!      ! call read_perturb_file(Bxc_1)
!      ! call read_perturb_file(Bxc_2)
!       call read_perturb_file(Bxc_3)
!       
!       call read_perturb_file(V_p)
      
      

!       allocate(V_loc( V_0%nr1 * V_0%nr2 * V_0%nr3, 2))
!       call get_vloc_colin()
       
      allocate(evc1(2*npwx,nbnd))
      allocate(evc2(2*npwx,nbnd))
      !allocate(evc3(2*npwx,nbnd))
      !allocate(evc4(2*npwx,nbnd))
      allocate(psic1(dfftp%nnr))
      allocate(psic2(dfftp%nnr))
      allocate(psic3(dfftp%nnr))
      allocate(psic4(dfftp%nnr))
      write(*,*)'ML2',size(psic2)

      write(*,"(///A56)")'----------------------------'
      write (*,"(/A55/)") 'Start M calculation k loop'
      write(*,"(A56//)")'----------------------------'
      nk = nks
      ibnd0=bnd_initial
      ibnd=bnd_final

      do ik0 = kpoint_initial,kpoint_final
            do ik = 1, nk
      write(*,*)'ik',ik
                  ikk = ik 
                  ikk0 = ik0
           
                  CALL read_collected_wfc ( restart_dir(), ikk, evc2 )
      write(*,*)'evc2',evc2(1,2),size(evc2)
                  CALL read_collected_wfc ( restart_dir(), ikk0, evc1 )
                  call calcmdefect_ml_rs(ibnd0,ibnd,ik0,ik)
      write(*,*)'evc2',evc2(1,1)
                  call calcmdefect_mnl_ks(ibnd0,ibnd,ik0,ik)
                  
                  !call calcmdefect_ml_rs_noncolin(ibnd0,ibnd,ikk0,ikk)
                  !call calcmdefect_mnl_ks_soc(ibnd0,ibnd,ikk0,ikk)

1003 format(A24,I6,I6,A6,I6,I6 " ( ",e17.9," , ",e17.9," ) ",e17.9//)
            write (*,1003) 'M_tot ni ki --> nf kf ', ibnd0,ikk0, '-->', ibnd,ikk, &
            m_loc+m_nloc, abs(m_loc+m_nloc)
            end do
      end do
  call environment_end('BANDS')
End Program edic






