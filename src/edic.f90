Program edic
  !USE environment,   ONLY : environment_start, environment_end
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
  !USE klist,     ONLY : nkstot, two_fermi_energies, nks
  !USE noncollin_module, ONLY : noncolin, i_cons, npol
  !USE noncollin_module, ONLY :  i_cons, npol
  USE noncollin_module, ONLY :  npol
  !USE lsda_mod,  ONLY : nspin
  USE mp_global, ONLY : mp_startup,mp_global_end
  USE io_global, ONLY : ionode
  !USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp,        ONLY : mp_bcast
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
  integer :: kp_idx_i,kp_idx_f, bnd_idx_i,bnd_idx_f
  integer:: p_rank,p_size
  CHARACTER(LEN=6), EXTERNAL :: int_to_char  
  !  CHARACTER (len=256) :: filband, filp, outdir
  !  LOGICAL :: lsigma(4), lsym, lp, no_overlap, plot_2d, wfc_is_collected, exst
  !  INTEGER :: spin_component, firstk, lastk
  !integer :: nr1x_perturb, nr2x_perturb, nr3x_perturb, nr1_perturb, nr2_perturb, nr3_perturb, &
  !nat_perturb, ntyp_perturb, ibrav_perturb, plot_num_perturb,  i_perturb,nkb_perturb
  !  NAMELIST / bands / outdir, prefix, filband, filp, spin_component, lsigma,&
  !                       lsym, lp, filp, firstk, lastk, no_overlap, plot_2d
  
  
  !!!!!!!!!!!!! hdf5 debug
  !INTEGER(HID_T)                               :: loc_id, attr_id, data_type, mem_type
  integer :: ierr
  complex :: mlocal0,mlocal1,mlocal,mnonlocal0,mnonlocal1,mnonlocal,mcharge
  !CALL H5Tcopy_f( H5T_NATIVE_INTEGER, mem_type, ierr )      
  !write(*,*) 'ierr        ', ierr
  !!!!!!!!!!!!! hdf5 debug

  CALL mp_startup( start_images=.TRUE. )
  !CALL mp_startup ( )
  IF ( ionode )  THEN
      write(*,"(///A56)")'----------------------------'
      write (*,"(/A55/)") 'Start EDIC program '
      write(*,"(A56//)")'----------------------------'
  ENDIF
  !      write(*,*)' pool_id',my_pool_id
  !      write(*,*)' pool_id',my_pool_id
  !write(*,*)' 1image_id',my_image_id,nimage
  !call  mpi_comm_rank(mpi_comm_world,ig,ik)
  !call  mpi_barrier(mpi_comm_world)
  !write(*,*)' 2image_id',mpi_comm_world,my_image_id,ig,ik
  !call  mpi_barrier(mpi_comm_world)



  CALL environment_start ( 'EDIC' )

  !#if defined(__INTEL_COMPILER)
  !    CALL remove_stack_limit ( )
  !#endif
  !    CALL init_clocks(.TRUE.) 
  !    CALL start_clock( TRIM(code) )



  call  mpi_comm_rank(mpi_comm_world,p_rank,ik)
  call  mpi_comm_size(mpi_comm_world,p_size,ik)

  !write(*,*)' 2image_id',my_image_id,p_rank,p_size,ik
  !write(*,*) '1start '
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
      write(*,*) '0'
      CALL input_from_file ( )
      write(*,*) '1'
      !READ (5, bands, iostat = ios)
      !write(*,*) '2',prefix,outdir
      !lsigma(4)=.false.
  ENDIF

  tmp_unit = find_free_unit()

  OPEN(unit=tmp_unit,file = 'calcmdefect.dat',status='old',err=20)
  20 continue
  READ(tmp_unit,calcmcontrol,iostat=ios)
  CLOSE(tmp_unit)

  tmp_dir = trimcheck (outdir)
  !call  mpi_barrier(mpi_comm_world)

  wfc_is_collected=.true.
  CALL read_file_new(wfc_is_collected)
  write(*,1003) 'Number of Bands', nbnd
  write(*,1003) 'Number of Plane Waves', npwx
  write(*,1003) 'Number of K-points', nks
  1000 format(A24," = ",I6)
  nwordwfc = nbnd*npwx*npol
  !io_level = 1
  
  !call  mpi_barrier(mpi_comm_world)
  !write(*,*)'nr1,nat_perturb)',-dfftp%nr1,dfftp%nr1,nat_perturb
  !call getvrsc()
  !write(*,*)'nr1,nat_perturb)',-dfftp%nr1,dfftp%nr1,nat_perturb

  call getwtdata()

  !call  mpi_barrier(mpi_comm_world)
  if (calcmcharge) then

    !if(eps_type=='qeh') then
    if(doqeh) then
      write(*,"(///A56)")'----------------------------'
      write (*,"(/A55/)") 'Read QEH dielectric data '
      call get_qeh_eps_data()
    endif
   

    !elseif(eps_type=='gw') then
    if(dogwfull .or. dogwdiag) then
      write(*,"(///A56)")'----------------------------'
      write (*,"(/A55/)") 'Read BGW dielectric data '
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
      !!!!!!!!!!!!! test 1
      ! IF ( ionode )  THEN
      !do ig=0,p_size-1
      !if ( p_rank==ig) then
      ! call gw_eps_read(gw_epsmat_filename,gw_epsq1_data)
      ! call gw_eps_read(gw_eps0mat_filename,gw_epsq0_data)
      !endif
      !call  mpi_barrier(mpi_comm_world)
      !enddo
      
      !!!!!!!!!!!!! test 2
      !if ( p_rank==1) then
      ! call gw_eps_read(gw_epsmat_filename,gw_epsq1_data)
      ! call gw_eps_read(gw_eps0mat_filename,gw_epsq0_data)
      !      write (*,*) 'rank0',p_rank
      !endif
      !      write (*,*) 'rank1',p_rank
      !call  mpi_barrier(mpi_comm_world)
      !      write (*,*) 'rank2',p_rank
      !
      !
      !if ( p_rank==0) then
      !      write (*,*) 'rank3',p_rank
      ! call gw_eps_read(gw_epsmat_filename,gw_epsq1_data)
      ! call gw_eps_read(gw_eps0mat_filename,gw_epsq0_data)
      !      write (*,*) 'rank4',p_rank
      !endif
      !      write (*,*) 'rank5',p_rank
      !call  mpi_barrier(mpi_comm_world)
      !
      !
      
      !!!!!!!!!!!!!! test 3 work
      !if ( p_rank==0) then
      ! call gw_eps_read('epsmat.h5.0',gw_epsq1_data)
      !endif
      !if ( p_rank==1) then
      ! call gw_eps_read('epsmat.h5.1',gw_epsq1_data)
      !endif
      
      !!!!!!!!!!!!! test 4
      !if ( p_rank==0) then
      !      write (*,*) 'rank3',p_rank
      ! call gw_eps_read(gw_epsmat_filename,gw_epsq1_data)
      ! call gw_eps_read(gw_eps0mat_filename,gw_epsq0_data)
      !call get_gind_rhoandpsi_gw(gw_epsq1_data)
      !call get_gind_rhoandpsi_gw(gw_epsq0_data)
      !      write (*,*) 'rank4',p_rank
      !
      !endif
      !call  mpi_barrier(mpi_comm_world)
      !
      !      write (*,*) 'rank5',p_rank
      !mpi_bcast(gw_epsq1_data%q_g_commonsuset_size,1,mpi_integer,0,mpi_comm_world)
      !      write (*,*) 'rank6',p_rank
      !
      
      !!!!!!!!!!!!!!!!!!! test 3.2
      !do ig=0,p_size-1
      ! if ( p_rank==ig) then
      !!
      ! call gw_eps_read(trim(gw_epsmat_filename)//'.'//trim(int_to_char(p_rank)),gw_epsq1_data)
      ! call gw_eps_read(trim(gw_eps0mat_filename)//'.'//trim(int_to_char(p_rank)),gw_epsq0_data)
      !endif
      !enddo
      
      !!!!!!!!!!!!!!!!!! test 3.3
      !call flush(6)
      if ( p_rank==0) then
          write(*,*)'gw read 1', p_rank
          call gw_eps_read(gw_epsmat_filename,gw_epsq1_data)
          call gw_eps_read(gw_eps0mat_filename,gw_epsq0_data)
          !call gw_eps_read(trim(gw_epsmat_filename)//'.'//trim(int_to_char(p_rank)),gw_epsq1_data)
          !call gw_eps_read(trim(gw_eps0mat_filename)//'.'//trim(int_to_char(p_rank)),gw_epsq0_data)
      endif
  
  
      !write(*,*)'gw read 2', p_rank
      !call  mpi_barrier(mpi_comm_world)
      !write(*,*)'gw read 3', p_rank
      !call  mpi_barrier(mpi_comm_world)
      !write(*,*)'gw read 4', p_rank
      ! call gw_eps_bcast(gw_epsq0_data)
  
  
      ! CALL MPI_BCAST( array, n, MPI_DOUBLE_PRECISION, root, gid, ierr )
      ! if ( p_rank/=0) then
      !allocate(gw_epsq0_data%ng_data(1))
      !endif 
      !write(*,*)'gw read 5 rank1', gw_epsq0_data%ng_data, 1, MPI_integer, 0, mpi_comm_world, ierr
      !           CALL MPI_BCAST( gw_epsq0_data%ng_data, 1, MPI_integer, 0, mpi_comm_world, ierr )
      write(*,*)'gw read 5 rank', p_rank,gw_epsq0_data%ng_data,MPI_comm_world 
      call gw_eps_bcast(p_rank,0,gw_epsq1_data,MPI_comm_world,mpi_integer,MPI_DOUBLE_PRECISION)
      call gw_eps_bcast(p_rank,0,gw_epsq0_data,MPI_comm_world,mpi_integer,MPI_DOUBLE_PRECISION)
      write(*,*)'gw read 5 rank', p_rank,gw_epsq0_data%ng_data 
  
      ! call gw_eps_read(gw_epsmat_filename,gw_epsq1_data)
      ! call gw_eps_read(gw_eps0mat_filename,gw_epsq0_data)
      
      !  write(*,*) 'rank, eps', p_rank,gw_epsq1_data%gind_rho2psi(:)
      ! ENDIF
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      call gw_eps_init(gw_epsq1_data)
      call gw_eps_init(gw_epsq0_data)

      call get_gind_rhoandpsi_gw(gw_epsq1_data)
      call get_gind_rhoandpsi_gw(gw_epsq0_data)
    endif
    !elseif(eps_type=='tf')then
    !  if (k0screen_read<0) stop ('screening wavelength not correct')
    !  write(*,*) 'TF screening wavelength',k0screen_read
    !else
    !  stop ('eps_type incorrect')
    !endif
  
  endif


  if (calcmlocal .or. calcmnonlocal) then
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
          Bxc_1_p%filename = Bxc_1_p_filename
          Bxc_2_p%filename = Bxc_2_p_filename
          Bxc_3_p%filename = Bxc_3_p_filename
          call read_perturb_file(Bxc_1_p)
          call read_perturb_file(Bxc_2_p)
          call read_perturb_file(Bxc_3_p)
          Bxc_1_d%filename = Bxc_1_d_filename
          Bxc_2_d%filename = Bxc_2_d_filename
          Bxc_3_d%filename = Bxc_3_d_filename
          call read_perturb_file(Bxc_1_d)
          call read_perturb_file(Bxc_2_d)
          call read_perturb_file(Bxc_3_d)
      
          allocate(V_nc( V_d%nr1 * V_d%nr2 * V_d%nr3, 4))
          V_nc(:, 1) =     V_d%pot(:) -    V_p%pot(:) - V_d_shift + V_p_shift
          V_nc(:, 2) = Bxc_1_d%pot(:) -Bxc_1_p%pot(:) 
          V_nc(:, 3) = Bxc_2_d%pot(:) -Bxc_2_p%pot(:) 
          V_nc(:, 4) = Bxc_3_d%pot(:) -Bxc_3_p%pot(:) 
      else
          allocate(V_colin( V_d%nr1 * V_d%nr2 * V_d%nr3))
          V_colin(:) = V_d%pot(:) -V_p%pot(:) - V_d_shift + V_p_shift
      endif
      write(*,*) 'vcolin' ,shape(v_colin)
    
      !       allocate(V_loc( V_d%nr1 * V_d%nr2 * V_d%nr3, 2))
      !       call get_vloc_colin()
  endif
       
  call  mpi_barrier(mpi_comm_world)

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
  !      write(*,*)'ML2',size(psic2)

  write(*,"(///A56)")'----------------------------'
  write (*,"(/A55/)") 'Start M calculation k loop'
  write(*,"(A56//)")'----------------------------'
  !nk = nks
  !ibnd0=bnd_initial
  !ibnd=bnd_final

  ! k pair parralellization
  
  
  !      do ik0 = kpoint_initial,kpoint_final
  !            do ik = 1, nk
  
  call  mpi_barrier(mpi_comm_world)
  do ig = 1,bndkp_pair%npairs
    write(*,*)'ig image_id_loop image',ig,p_rank,p_size 
    if ( (ig<=bndkp_pair%npairs/p_size*(p_rank+1) .and. ig>bndkp_pair%npairs/p_size*(p_rank)) &
                          .or.(p_rank==p_size-1 .and. ig>bndkp_pair%npairs/p_size*(p_rank)) ) then
      write(*,*)'in loop, image_id',ig,my_image_id
      kp_idx_i = bndkp_pair%kp_idx(ig,1) 
      kp_idx_f = bndkp_pair%kp_idx(ig,2) 
      bnd_idx_i = bndkp_pair%bnd_idx(ig,1) 
      bnd_idx_f = bndkp_pair%bnd_idx(ig,2) 
      ! allocate(evc(1*npwx,nbnd))
      write(*,*)'evc1 size',shape(evc1)
      write(*,*)'restart_dir()', restart_dir()
      write(*,*)'kp_idx_i,kp_idx_f',kp_idx_i,kp_idx_f
      write(*,*)'bnd_idx_i,bnd_idx_f',bnd_idx_i,bnd_idx_f
      !CALL read_collected_wfc ( restart_dir(), kp_idx_i, evc )
      CALL read_collected_wfc ( restart_dir(), kp_idx_i, evc1 )
      write(*,*)'evc1',evc1(1,1),size(evc1)
      CALL read_collected_wfc ( restart_dir(), kp_idx_f, evc2 )
      write(*,*)'evc2',evc2(1,1),shape(evc2)

      if (calcmcharge .and. mcharge_dolfa) then
          call calcmdefect_charge_lfa(bnd_idx_f,bnd_idx_i,kp_idx_f,kp_idx_i,noncolin,mcharge)
          bndkp_pair%m(ig)=mcharge
      endif
      if (calcmcharge .and. .not. mcharge_dolfa) then
          !write(*,*) k0screen_read
          write(*,*) 'k0sc1',k0screen_read
          call calcmdefect_charge_nolfa(bnd_idx_f,bnd_idx_i,kp_idx_f,kp_idx_i,noncolin,mcharge)
          bndkp_pair%m(ig)=mcharge
      endif

      if (noncolin .and. .not. lspinorb .and. calcmlocal)then
          call calcmdefect_ml_rs_noncolin(bnd_idx_f,bnd_idx_i,kp_idx_f,kp_idx_i,mlocal)
      endif
      if (noncolin .and. .not. lspinorb .and. calcmnonlocal)then
          call calcmdefect_mnl_ks_noncolin(bnd_idx_f,bnd_idx_i,kp_idx_f,kp_idx_i,v_d,mnonlocal0)
          call calcmdefect_mnl_ks_noncolin(bnd_idx_f,bnd_idx_i,kp_idx_f,kp_idx_i,v_p,mnonlocal1)
      endif
      
      if (noncolin .and. lspinorb .and. calcmlocal)then
          call calcmdefect_ml_rs_noncolin(bnd_idx_f,bnd_idx_i,kp_idx_f,kp_idx_i,mlocal)
      endif
      if (noncolin .and. lspinorb .and. calcmnonlocal)then
          call calcmdefect_mnl_ks_soc(bnd_idx_f,bnd_idx_i,kp_idx_f,kp_idx_i,v_d,mnonlocal0)
          call calcmdefect_mnl_ks_soc(bnd_idx_f,bnd_idx_i,kp_idx_f,kp_idx_i,v_p,mnonlocal1)
      endif
      if ( .not. noncolin .and. calcmlocal)then
          call calcmdefect_ml_rs(bnd_idx_f,bnd_idx_i,kp_idx_f,kp_idx_i,V_colin,mlocal)
      endif
      if ( .not. noncolin .and. calcmnonlocal)then
          call calcmdefect_mnl_ks(bnd_idx_f,bnd_idx_i,kp_idx_f,kp_idx_i,v_d,mnonlocal0)
          call calcmdefect_mnl_ks(bnd_idx_f,bnd_idx_i,kp_idx_f,kp_idx_i,v_p,mnonlocal1)
          mnonlocal=mnonlocal0-mnonlocal1
      endif
      bndkp_pair%m(ig)=mlocal+mnonlocal0-mnonlocal1
      !      write(*,*)'evc2',evc2(1,1)
      !      write(*,*)'evc2',evc2(1,1)
   
  
      write (*,1003) 'M_tot ni ki --> nf kf ', bnd_idx_f,kp_idx_f, '-->', bnd_idx_i,kp_idx_i, &
      m_loc+m_nloc, abs(m_loc+m_nloc)
      
    endif
  !            end do
  !      end do
  enddo
  call environment_end('EDIC')
  CALL mp_global_end()
1003 format(A24,I6,I6,A6,I6,I6 " ( ",e17.9," , ",e17.9," ) ",e17.9//)
End Program edic






