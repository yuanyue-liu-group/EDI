Program edic
!  USE environment,   ONLY : environment_start, environment_end

      Use kinds,    only: dp

      USE io_files,  ONLY : prefix, tmp_dir, nwordwfc, iunwfc, restart_dir
!      Use edic_mod,   only: V_file, V_nc, V_colin, V_d, Bxc_1, Bxc_2, Bxc_3, V_p
!      Use edic_mod,   only: gw_epsq1_data,gw_epsq0_data
!      Use edic_mod, only: v_p_shift,v_d_shift
      Use edic_mod
!      use edic_mod, only: bndkp_pair
!      Use edic_mod, Only : evc1,evc2,evc3,evc4,psic1, psic2, psic3, psic4
!      Use edic_mod, Only: m_loc, m_nloc

      Use wvfct, ONLY: npwx, nbnd!, wg, et, g2kin
      Use fft_base,  ONLY: dfftp!, dffts
!      Use fft_interfaces, ONLY : fwfft, invfft
      USE pw_restart_new,   ONLY : read_collected_wfc
      Use klist,  Only : nks
!  USE parameters,ONLY : npk
!  USE pw_restart_new, ONLY : read_xml_file
  USE wavefunctions,    ONLY : evc
!  USE wavefunctions_gpum, ONLY : using_evc
!  USE buffers,        ONLY : open_buffer, close_buffer, save_buffer


!  USE control_flags,    ONLY : io_level
!  USE control_flags, ONLY : gamma_only
  USE environment,   ONLY : environment_start, environment_end
!  USE klist,     ONLY : nkstot, two_fermi_energies, nks
  !USE noncollin_module, ONLY : noncolin, i_cons, npol
  !USE noncollin_module, ONLY :  i_cons, npol
  USE noncollin_module, ONLY :  npol
  !USE lsda_mod,  ONLY : nspin

  USE mp_global, ONLY : mp_startup,mp_global_end
  USE io_global, ONLY : ionode
  !USE io_global, ONLY : ionode, ionode_id, stdout
  !USE mp,        ONLY : mp_bcast
  !USE mp_images, ONLY : intra_image_comm
  USE mp_pools, ONLY : npool, my_pool_id
  USE mp_images, ONLY : nimage, my_image_id

use hdf5
  USE parallel_include
 
      Implicit none
      
      

  CHARACTER(LEN=256), EXTERNAL :: trimcheck
      integer :: tmp_unit
      integer, external :: find_free_unit
      integer :: ios
      integer :: ig,ikk, ikk0, ibnd, ibnd0, ik, ik0, nk
integer:: p_rank,p_size
!  CHARACTER (len=256) :: filband, filp, outdir
!  LOGICAL :: lsigma(4), lsym, lp, no_overlap, plot_2d, wfc_is_collected, exst
!  INTEGER :: spin_component, firstk, lastk
!integer :: nr1x_perturb, nr2x_perturb, nr3x_perturb, nr1_perturb, nr2_perturb, nr3_perturb, &
!nat_perturb, ntyp_perturb, ibrav_perturb, plot_num_perturb,  i_perturb,nkb_perturb
!  NAMELIST / bands / outdir, prefix, filband, filp, spin_component, lsigma,&
!                       lsym, lp, filp, firstk, lastk, no_overlap, plot_2d


!!!!!!!!!!!! hdf5 debug
!INTEGER(HID_T)                               :: loc_id, attr_id, data_type, mem_type
!integer :: ierr
!CALL H5Tcopy_f( H5T_NATIVE_INTEGER, mem_type, ierr )      
!!!!!!!!!!!! hdf5 debug

  CALL mp_startup( start_images=.TRUE. )
  !CALL mp_startup ( )
  IF ( ionode )  THEN
      write(*,"(///A56)")'----------------------------'
      write (*,"(/A55/)") 'Start EDIC program '
      write(*,"(A56//)")'----------------------------'
  ENDIF
!      write(*,*)' pool_id',my_pool_id
!      write(*,*)' pool_id',my_pool_id
      write(*,*)' 1image_id',my_image_id,nimage
call  mpi_comm_rank(mpi_comm_world,ig,ik)
call  mpi_barrier(mpi_comm_world)
      write(*,*)' 2image_id',mpi_comm_world,my_image_id,ig,ik
call  mpi_barrier(mpi_comm_world)



  CALL environment_start ( 'EDIC' )

!#if defined(__INTEL_COMPILER)
!    CALL remove_stack_limit ( )
!#endif
!    CALL init_clocks(.TRUE.) 
!    CALL start_clock( TRIM(code) )



call  mpi_comm_rank(mpi_comm_world,p_rank,ik)
call  mpi_comm_size(mpi_comm_world,p_size,ik)

      write(*,*)' 2image_id',my_image_id,p_rank,p_size,ik
write(*,*) '1start '
  !
  !   set default values for variables in namelist
  !
!  prefix = 'pwscf'
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
  IF ( ionode )  THEN
!     !
write(*,*) '0'
     CALL input_from_file ( )
write(*,*) '1'
!     !
!     READ (5, bands, iostat = ios)
!write(*,*) '2',prefix,outdir
!     !
!     lsigma(4)=.false.
!     !
  ENDIF

      tmp_unit = find_free_unit()
      OPEN(unit=tmp_unit,file = 'calcmdefect.dat',status='old',err=20)
      20 continue
         READ(tmp_unit,calcmcontrol,iostat=ios)
      CLOSE(tmp_unit)

     tmp_dir = trimcheck (outdir)

  wfc_is_collected=.true.
  CALL read_file_new(wfc_is_collected)
  write(*,1003) 'Number of Bands', nbnd
  write(*,1003) 'Number of Plane Waves', npwx
  write(*,1003) 'Number of K-points', nks
  1000 format(A24," = ",I6)
  nwordwfc = nbnd*npwx*npol
  !io_level = 1
  
!write(*,*)'nr1,nat_perturb)',-dfftp%nr1,dfftp%nr1,nat_perturb
!call getvrsc()
!write(*,*)'nr1,nat_perturb)',-dfftp%nr1,dfftp%nr1,nat_perturb

call getwtdata()

if (calcmcharge) then

      write(*,"(///A56)")'----------------------------'
      write (*,"(/A55/)") 'Read QEH dielectric data '
call getepsdata()
 
      write(*,"(///A56)")'----------------------------'
      write (*,"(/A55/)") 'Read BGW dielectric data '
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 call gw_eps_init(gw_epsmat_filename,gw_epsq1_data)
 call gw_eps_init(gw_eps0mat_filename,gw_epsq0_data)
call get_gind_rhoandpsi_gw(gw_epsq1_data)
call get_gind_rhoandpsi_gw(gw_epsq0_data)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
endif



    V_d%filename = V_d_filename
    V_p%filename = V_p_filename
    call read_perturb_file(V_d)
    call read_perturb_file(V_p)
    if (lvacalign) then
        v_d_shift = sum(V_d%pot(vac_idx*V_d%nr1*V_d%nr2:(vac_idx+1)*V_d%nr1*V_d%nr2))/(V_d%nr1*V_d%nr2)
        v_p_shift = sum(v_p%pot(vac_idx*v_p%nr1*v_p%nr2:(vac_idx+1)*v_p%nr1*v_p%nr2))/(v_p%nr1*v_p%nr2)
    elseif (lcorealign) then
        v_d_shift = core_v_d
        v_p_shift = core_v_p
    endif
 
if (noncolin .or. lspinorb)then
    Bxc_1%filename = Bxc_1_filename
    Bxc_2%filename = Bxc_2_filename
    Bxc_3%filename = Bxc_3_filename
    !call read_perturb_file(Bxc_1)
    !call read_perturb_file(Bxc_2)
    call read_perturb_file(Bxc_3)

    allocate(V_nc( V_d%nr1 * V_d%nr2 * V_d%nr3, 2))
       V_nc(:, 1) = V_d%pot(:) -V_p%pot(:) + Bxc_3%pot(:) - V_d_shift + V_p_shift
        V_nc(:, 2) = V_d%pot(:) -V_p%pot(:) - Bxc_3%pot(:) - V_d_shift + V_p_shift
else
    allocate(V_colin( V_d%nr1 * V_d%nr2 * V_d%nr3))
    V_colin(:) = V_d%pot(:) -V_p%pot(:) - V_d_shift + V_p_shift
endif
write(*,*) 'vcolin' ,shape(v_colin)

!       allocate(V_loc( V_d%nr1 * V_d%nr2 * V_d%nr3, 2))
!       call get_vloc_colin()
       

if (noncolin .or. lspinorb)then
      allocate(evc1(2*npwx,nbnd))
      allocate(evc2(2*npwx,nbnd))
else
      allocate(evc1(1*npwx,nbnd))
      allocate(evc2(1*npwx,nbnd))
endif
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

! k pair parralellization


!      do ik0 = kpoint_initial,kpoint_final
!            do ik = 1, nk

do ig = 1,bndkp_pair%npairs
      write(*,*)'ig image_id_loop image',ig,p_rank,p_size 
if ( (ig<=bndkp_pair%npairs/p_size*(p_rank+1) .and. ig>bndkp_pair%npairs/p_size*(p_rank)) &
                        .or.(p_rank==p_size-1 .and. ig>bndkp_pair%npairs/p_size*(p_rank)) ) then

      write(*,*)'in loop, image_id',ig,my_image_id
                  ikk = bndkp_pair%kp_idx(ig,1) 
                  ikk0 = bndkp_pair%kp_idx(ig,2) 
                  ibnd = bndkp_pair%bnd_idx(ig,1) 
                  ibnd0 = bndkp_pair%bnd_idx(ig,2) 
           
  ! allocate(evc(1*npwx,nbnd))
write(*,*)'evc size',shape(evc),shape(evc1)
write(*,*)' restart_dir()', restart_dir()
write(*,*)'ikk,ikk0',ikk,ikk0
   !CALL read_collected_wfc ( restart_dir(), ikk, evc )

                  CALL read_collected_wfc ( restart_dir(), ikk, evc2 )
      write(*,*)'evc2',evc2(1,1),size(evc2)
                  CALL read_collected_wfc ( restart_dir(), ikk0, evc1 )
      write(*,*)'evc1',evc1(1,1),shape(evc1)

if (noncolin )then
!                  call calcmdefect_ml_rs_noncolin(ibnd0,ibnd,ikk0,ikk)
!                  call calcmdefect_mnl_ks_noncolin(ibnd0,ibnd,ikk0,ikk,v_d)
endif

if ( lspinorb)then
!                  call calcmdefect_ml_rs_noncolin(ibnd0,ibnd,ikk0,ikk)
!                  call calcmdefect_mnl_ks_soc(ibnd0,ibnd,ikk0,ikk,v_d)
endif
if ( (.not. lspinorb ).and. (.not. noncolin ))then
                  call calcmdefect_ml_rs(ibnd0,ibnd,ik0,ik,V_colin)
                  call calcmdefect_mnl_ks(ibnd0,ibnd,ik0,ik,v_d)
                  call calcmdefect_mnl_ks(ibnd0,ibnd,ik0,ik,v_p)
endif
!      write(*,*)'evc2',evc2(1,1)
!      write(*,*)'evc2',evc2(1,1)
                  
if (calcmcharge) then


if (mcharge_dolfa) then
                  call calcmdefect_charge_lfa(ibnd0,ibnd,ikk0,ikk)
else
                  call calcmdefect_charge_nolfa(ibnd0,ibnd,ikk0,ikk)
endif
endif

1003 format(A24,I6,I6,A6,I6,I6 " ( ",e17.9," , ",e17.9," ) ",e17.9//)
            write (*,1003) 'M_tot ni ki --> nf kf ', ibnd0,ikk0, '-->', ibnd,ikk, &
            m_loc+m_nloc, abs(m_loc+m_nloc)
endif
!            end do
!      end do
enddo

  call environment_end('EDIC')

  CALL mp_global_end()

End Program edic






