Program edic
  Use kinds,    only: dp
  USE io_files,  ONLY : prefix, tmp_dir, nwordwfc, iunwfc, restart_dir
  Use edic_mod
  Use wvfct, ONLY: npwx, nbnd
  Use fft_base,  ONLY: dfftp
  USE pw_restart_new,   ONLY : read_collected_wfc
  Use klist,  Only : nks
  USE wavefunctions,    ONLY : evc
  USE environment,   ONLY : environment_start, environment_end
  USE noncollin_module, ONLY :  npol
  USE mp_global, ONLY : mp_startup,mp_global_end
  USE io_global, ONLY : ionode
  USE mp,        ONLY : mp_bcast
  USE mp_pools, ONLY : npool, my_pool_id
  USE mp_images, ONLY : nimage, my_image_id
  use hdf5
  USE parallel_include
 
  Implicit none

  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  integer :: tmp_unit
  integer, external :: find_free_unit
  integer :: ios
  integer :: nchunk
  integer :: ig,ikk, ikk0, ibnd, ibnd0, ik, ik0, nk
  integer :: kp_idx_i,kp_idx_f, bnd_idx_i,bnd_idx_f
  integer:: p_rank,p_size
  CHARACTER(LEN=6), EXTERNAL :: int_to_char  
  integer :: ierr
  COMPLEX(DP)  :: mlocal0,mlocal1,mlocal,mnonlocal0,mnonlocal1,mnonlocal,mcharge
    
  !!!!!!!!!!!!! hdf5 debug
  !INTEGER(HID_T)                               :: loc_id, attr_id, data_type, mem_type
  !CALL H5Tcopy_f( H5T_NATIVE_INTEGER, mem_type, ierr )      
  !write(*,*) 'ierr        ', ierr
  !!!!!!!!!!!!! hdf5 debug

  integer::p_source
  real(dp):: m_tmp1,m_tmp2,m_tmp3,m_tmp4

  CALL mp_startup( start_images=.TRUE. )
  !CALL mp_startup ( )
  IF ( ionode )  THEN
      write(*,"(///A56)")'----------------------------'
      write (*,"(/A55/)") 'Start EDIC program '
      write(*,"(A56//)")'----------------------------'
  ENDIF
  

  CALL environment_start ( 'EDIC' )

  call  mpi_comm_rank(mpi_comm_world,p_rank,ik)
  call  mpi_comm_size(mpi_comm_world,p_size,ik)

  ios = 0
  IF ( ionode )  THEN
      CALL input_from_file ( )
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
  write(*,*) 'dogw', dogwfull
  write(*,*) 'dogw', dogwdiag
  write(*,*) 'dogw', do2d,do3d
  1000 format(A24," = ",I6)
  nwordwfc = nbnd*npwx*npol

  call getwtdata()
  call  mpi_barrier(mpi_comm_world,ierr)

  if (calcmcharge) then

    if(doqeh) then
      write(*,"(///A56)")'----------------------------'
      write (*,"(/A55/)") 'Read QEH dielectric data '
      call get_qeh_eps_data()
    endif
   

    if(dogwfull .or. dogwdiag) then
      write(*,"(///A56)")'----------------------------'
      write (*,"(/A55/)") 'Read BGW dielectric data '
      if ( p_rank==0) then
          write(*,*)'gw read 1', p_rank
          call gw_eps_read(gw_epsmat_filename,gw_epsq1_data)
          call gw_eps_read(gw_eps0mat_filename,gw_epsq0_data)
      endif
 
      write(*,*)'gw read 5 rank', p_rank,gw_epsq0_data%ng_data,MPI_comm_world 
      call gw_eps_bcast(p_rank,0,gw_epsq1_data,MPI_comm_world,mpi_integer,MPI_DOUBLE_PRECISION)
      call gw_eps_bcast(p_rank,0,gw_epsq0_data,MPI_comm_world,mpi_integer,MPI_DOUBLE_PRECISION)
      write(*,*)'gw read 5 rank', p_rank,gw_epsq0_data%ng_data 
   
      call gw_eps_init(gw_epsq1_data)
      call gw_eps_init(gw_epsq0_data)

      call get_gind_rhoandpsi_gw(gw_epsq1_data)
      call get_gind_rhoandpsi_gw(gw_epsq0_data)
    endif
  
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
          allocate(V_nc1( V_d%nr1 * V_d%nr2 * V_d%nr3, 4))
          allocate(V_nc2( V_d%nr1 * V_d%nr2 * V_d%nr3, 4))
          allocate(V_nc3( V_d%nr1 * V_d%nr2 * V_d%nr3, 4))
          allocate(V_nc4( V_d%nr1 * V_d%nr2 * V_d%nr3, 4))
          V_nc(:, 1) =     V_d%pot(:) -    V_p%pot(:) - V_d_shift + V_p_shift
      else
          allocate(V_colin( V_d%nr1 * V_d%nr2 * V_d%nr3))
          V_colin(:) = V_d%pot(:) -V_p%pot(:) - V_d_shift + V_p_shift
      endif
      write(*,*) 'vcolin' ,shape(v_colin)
    
  endif
       
  call  mpi_barrier(mpi_comm_world,ierr)

  if (noncolin .or. lspinorb)then
      allocate(evc1(2*npwx,nbnd))
      allocate(evc2(2*npwx,nbnd))
  else
      allocate(evc1(1*npwx,nbnd))
      allocate(evc2(1*npwx,nbnd))
  endif
  allocate(psic1(dfftp%nnr))
  allocate(psic2(dfftp%nnr))
  allocate(psic3(dfftp%nnr))
  allocate(psic4(dfftp%nnr))

  write(*,"(///A56)")'----------------------------'
  write (*,"(/A55/)") 'Start M calculation k loop'
  write(*,"(A56//)")'----------------------------'

  call  mpi_barrier(mpi_comm_world,ierr)
  ! k pair parralellization
  
  nchunk=bndkp_pair%npairs/(p_size)+1
  do ig = 1,bndkp_pair%npairs
    write(*,*)'ig image_id_loop image',ig,p_rank,p_size 
    if ( (ig<=nchunk*(p_rank+1) .and. ig>nchunk*(p_rank)) )then
      write(*,*)'in loop, image_id',ig,my_image_id
      kp_idx_i = bndkp_pair%kp_idx(ig,1) 
      kp_idx_f = bndkp_pair%kp_idx(ig,2) 
      bnd_idx_i = bndkp_pair%bnd_idx(ig,1) 
      bnd_idx_f = bndkp_pair%bnd_idx(ig,2) 
      write(*,*)'evc1 size',shape(evc1)
      write(*,*)'restart_dir()', restart_dir()
      write(*,*)'kp_idx_i,kp_idx_f',kp_idx_i,kp_idx_f
      write(*,*)'bnd_idx_i,bnd_idx_f',bnd_idx_i,bnd_idx_f
      CALL read_collected_wfc ( restart_dir(), kp_idx_i, evc1 )
      write(*,*)'evc1',evc1(1,1),size(evc1)
      CALL read_collected_wfc ( restart_dir(), kp_idx_f, evc2 )
      write(*,*)'evc2',evc2(1,1),shape(evc2)

      if (calcmcharge .and. mcharge_dolfa) then
          call calcmdefect_charge_lfa(bnd_idx_f,bnd_idx_i,kp_idx_f,kp_idx_i,noncolin,mcharge)
          bndkp_pair%mc(ig)=mcharge
      endif
      if (calcmcharge .and. .not. mcharge_dolfa) then
          write(*,*) 'k0sc1',k0screen_read
          call calcmdefect_charge_nolfa(bnd_idx_f,bnd_idx_i,kp_idx_f,kp_idx_i,noncolin,mcharge)
          bndkp_pair%mc(ig)=mcharge
      endif

      if (noncolin .and. .not. lspinorb .and. calcmlocal)then
          call calcmdefect_ml_rs_noncolin(bnd_idx_f,bnd_idx_i,kp_idx_f,kp_idx_i,mlocal)
      endif
      if (noncolin .and. .not. lspinorb .and. calcmnonlocal)then
          call calcmdefect_mnl_ks_noncolin(bnd_idx_f,bnd_idx_i,kp_idx_f,kp_idx_i,v_d,mnonlocal0)
          call calcmdefect_mnl_ks_noncolin(bnd_idx_f,bnd_idx_i,kp_idx_f,kp_idx_i,v_p,mnonlocal1)
      endif
      
      if (noncolin .and. lspinorb .and. calcmlocal)then
          V_nc(:, 1) =     V_d%pot(:) - V_d_shift 
          call calcmdefect_ml_rs_noncolin(bnd_idx_f,bnd_idx_i,kp_idx_f,kp_idx_i,mlocal)
          V_nc(:, 1) =     V_p%pot(:) - V_p_shift
          call calcmdefect_ml_rs_noncolin(bnd_idx_f,bnd_idx_i,kp_idx_f,kp_idx_i,mlocal)
          V_nc(:, 1) =     V_d%pot(:) 
          call calcmdefect_ml_rs_noncolin(bnd_idx_f,bnd_idx_i,kp_idx_f,kp_idx_i,mlocal)
          V_nc(:, 1) =     V_p%pot(:)
          call calcmdefect_ml_rs_noncolin(bnd_idx_f,bnd_idx_i,kp_idx_f,kp_idx_i,mlocal)
          call calcmdefect_ml_rs_noncolin(bnd_idx_f,bnd_idx_i,kp_idx_f,kp_idx_i,mlocal)
          write (*,*)  'Mifl',mlocal,mnonlocal0,mnonlocal1
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
   

      write (*,*)    'ibnd, ik,ki(xyz),ei,vi(xyz);fbnd, fk,kf(xyz),ef,vf(xyz);  wt, m,mc,ml,mnl0,mnl1'
      write (*,*)  'Mif',  bndkp_pair%bnd_idx(ig,1),bndkp_pair%kp_idx(ig,1),&
                   bndkp_pair%k_coord(ig,1,1),bndkp_pair%k_coord(ig,2,1),bndkp_pair%k_coord(ig,3,1),&
                   bndkp_pair%e_pair(ig,1),&
                   bndkp_pair%v_pair(ig,1,1),bndkp_pair%v_pair(ig,2,1),bndkp_pair%v_pair(ig,3,1),&
                   bndkp_pair%bnd_idx(ig,2),bndkp_pair%kp_idx(ig,2),&
                   bndkp_pair%k_coord(ig,1,2),bndkp_pair%k_coord(ig,2,2),bndkp_pair%k_coord(ig,3,2),&
                   bndkp_pair%e_pair(ig,2),&
                   bndkp_pair%v_pair(ig,1,2),bndkp_pair%v_pair(ig,2,2),bndkp_pair%v_pair(ig,3,2), &
                   bndkp_pair%wt(ig),bndkp_pair%m(ig),bndkp_pair%mc(ig),mlocal,mnonlocal0,mnonlocal1
 
  
      write (*,1003) 'M_tot ni ki --> nf kf ', bnd_idx_f,kp_idx_f, '-->', bnd_idx_i,kp_idx_i, &
      m_loc+m_nloc, abs(m_loc+m_nloc)
      
    endif
  enddo


  call  mpi_barrier(mpi_comm_world,ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! bcast bndkp_pair%m and mc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nchunk=bndkp_pair%npairs/(p_size)+1
  do ig = 1,bndkp_pair%npairs
    write(*,*)'ig image_id_loop image',ig,p_rank,p_size 
    if ( (ig<=nchunk*(p_rank+1) .and. ig>nchunk*(p_rank)) )then
      write(*,*)'in loop, image_id',ig,my_image_id
      p_source=p_rank
      m_tmp1= real(bndkp_pair%m(ig))
      m_tmp2= aimag(bndkp_pair%m(ig))
      m_tmp3= real(bndkp_pair%mc(ig))
      m_tmp4= aimag(bndkp_pair%mc(ig))
    endif
    p_source=(ig-1)/nchunk
    write(*,*)'ps',p_rank,p_source
    write(*,*)'ig, image_id',ig,p_rank,p_source,my_image_id
    CALL MPI_BCAST(   m_tmp1, 1, MPI_DOUBLE_PRECISION,p_source, mpi_comm_world, ierr )
    CALL MPI_BCAST(   m_tmp2, 1, MPI_DOUBLE_PRECISION,p_source, mpi_comm_world, ierr )
    CALL MPI_BCAST(   m_tmp3, 1, MPI_DOUBLE_PRECISION,p_source, mpi_comm_world, ierr )
    CALL MPI_BCAST(   m_tmp4, 1, MPI_DOUBLE_PRECISION,p_source, mpi_comm_world, ierr )
    bndkp_pair%m(ig)=complex(m_tmp1,m_tmp2)
    bndkp_pair%mc(ig)=complex(m_tmp3,m_tmp4)
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! gamma mu calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if ( p_rank==0) call postprocess()

  call environment_end('EDIC')
  CALL mp_global_end()
1003 format(A24,I6,I6,A6,I6,I6 " ( ",e17.9," , ",e17.9," ) ",e17.9//)
End Program edic






