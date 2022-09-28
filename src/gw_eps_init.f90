!!!!!!!!!!!!!!!!1
! eps0mat.h5 should be from 0 to the smallest |q| of epsmat.h5 grids
! epsmat.h5 should be a fine grid of q points
!!!!!!!!!!!!!!!!

subroutine gw_eps_init(gw_)
  USE kinds, ONLY: DP,sgl
  USE HDF5
  use edic_mod, only: gw_eps_data
  
  CHARACTER(LEN=256) :: eps_filename_
  type(gw_eps_data),intent (inout) ,target:: gw_
  
  !!!!!!!!!!hdf5
  CHARACTER(LEN=256) :: h5filename      ! Dataset name
  CHARACTER(LEN=256) :: h5datasetname = "matrix-diagonal"     ! Dataset name
  INTEGER     ::   h5rank,h5error ! Error flag
  !INTEGER     ::  i, j
  !real(dp), DIMENSION(3,1465,2) :: h5dataset_data, data_out ! Data buffers
  !real(dp), DIMENSION(3,1465,2) :: h5dataset_data, data_out ! Data buffers
  real(dp), allocatable :: h5dataset_data_double(:), data_out(:)
  integer, allocatable :: h5dataset_data_integer(:)
  INTEGER(HSIZE_T), allocatable :: h5dims(:),h5maxdims(:)
  
  integer :: h5dims1(1),h5dims2(2),h5dims3(3),h5dims4(4),h5dims5(5),h5dims6(6)
  integer:: p_rank,p_size,ik



  !!!!!!!!!!!!!!!!!!!
  ! prep read gw h5 data
  
  ! qabs
  write(*,*) gw_%bvec_data(:,1)
  write(*,*) gw_%bvec_data(:,2)
  write(*,*) gw_%bvec_data(:,3)
  write(*,*) gw_%qpts_data(:,1)
  write(*,*) gw_%qpts_data(:,2)
  write(*,*) gw_%qpts_data(:,3)

  if (  allocated(gw_%qabs)) then
     deallocate(gw_%qabs)
  endif
  allocate(gw_%qabs(gw_%nq_data(1)))
  do ig1 = 1, gw_%nq_data(1)
      gw_%qabs(ig1)=norm2(&
              gw_%qpts_data(1,ig1)*gw_%bvec_data(:,1)+ &
              gw_%qpts_data(2,ig1)*gw_%bvec_data(:,2)+ &
              gw_%qpts_data(3,ig1)*gw_%bvec_data(:,3))
             !*gw_%blat_data(1)

      !write(*,*)              gw_%qpts_data(1,ig1)*gw_%bvec_data(1,:)
      !write(*,*)              gw_%qpts_data(2,ig1)*gw_%bvec_data(2,:)
      !write(*,*)              gw_%qpts_data(3,ig1)*gw_%bvec_data(3,:)

      !debug
      write(*,*)'gw_qabs debug ', gw_%qabs(ig1),gw_%epsmat_diag_data(:,1,ig1)
      !debug
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  select only common g vectors for eps with different q
  !!!!!!!!!!!!!!


      ! convert eps(q) g index to common gw-rho based g index
      !gw_%q_g_commonsubset_size
      !gw_%q_g_commonsubset2rho(:,:)
      !do ig = 1, gw_%ng_data(1)
      !  do iq=1,gw_%nq_data(1)
      !    gind_gw_%eps=gw_%gind_rho2eps_data(iq,ig)
      !       if      (gindgw_%_eps<gw_%nmtx(iq))  then
      !  enddo
      !enddo
      !eps(gw_%gind_rho2eps_data(iq,1:gw_%nmtx_data(iq)))


  !!!! initialize q_g_commonsubset_indinrho to gw's ind eps2rho
  !!!! q_g_comm
  if (  allocated(gw_%q_g_commonsubset_indinrho)) then
      deallocate(gw_%q_g_commonsubset_indinrho)
  endif
  allocate(gw_%q_g_commonsubset_indinrho(gw_%nmtx_max_data(1)))
  gw_%q_g_commonsubset_indinrho(:)=0
  gw_%q_g_commonsubset_indinrho(:)=gw_%gind_eps2rho_data(1:gw_%nmtx_data(1),1)
  write(*,*)  'gw_%q_g_commonsubset_indinrho',gw_%q_g_commonsubset_indinrho(1:40),'shape',shape(gw_%q_g_commonsubset_indinrho)

  !!!!!!!!!!!!!!
  !!!!  throw away index from q=1, where resulted gindinrho is not in other q points
  do iq=1,gw_%nq_data(1)
    do ig=1,gw_%nmtx_data(iq)
      if(gw_%q_g_commonsubset_indinrho(ig)>0) then
        if (gw_%gind_rho2eps_data(gw_%q_g_commonsubset_indinrho(ig),iq)>gw_%nmtx_data(iq) ) then
           gw_%q_g_commonsubset_indinrho(ig)=0
        endif
      endif
    enddo
  enddo
  write(*,*)  'gw_%q_g_commonsubset_indinrho',gw_%q_g_commonsubset_indinrho(1:40),'shape',shape(gw_%q_g_commonsubset_indinrho)


  ig=0
  do ig1=1,gw_%nmtx_max_data(1)
    if(gw_%q_g_commonsubset_indinrho(ig1)>0) ig=ig+1
  enddo

  !write(*,*)  'gw_%q_g_commonsubset_indinrho',gw_%q_g_commonsubset_indinrho(:)
  if (  allocated(gw_%q_g_commonsubset_indinrhotmp1)) then
     deallocate(gw_%q_g_commonsubset_indinrhotmp1)
  endif
  allocate(gw_%q_g_commonsubset_indinrhotmp1(ig))
  gw_%q_g_commonsubset_indinrhotmp1(:)=0


  ig1=1
  do ig=1,gw_%nmtx_max_data(1)
    if(gw_%q_g_commonsubset_indinrho(ig)>0) then 
       !write(*,*) gw_q_g_commonsubset_indinrhotmp1(ig1),gw_%q_g_commonsubset_indinrho(ig)
       gw_%q_g_commonsubset_indinrhotmp1(ig1)=gw_%q_g_commonsubset_indinrho(ig) 
       ig1=ig1+1
    endif
  enddo
  deallocate(gw_%q_g_commonsubset_indinrho)
  allocate(gw_%q_g_commonsubset_indinrho(size(gw_%q_g_commonsubset_indinrhotmp1)))
  gw_%q_g_commonsubset_indinrho(:)=gw_%q_g_commonsubset_indinrhotmp1(:) 
  
  write(*,*)  'gw_%q_g_commonsubset_indinrho finalized',&
         gw_%q_g_commonsubset_indinrho(1:40),'shape',shape(gw_%q_g_commonsubset_indinrho)
  gw_%q_g_commonsubset_size=size(gw_%q_g_commonsubset_indinrho)
  do ig=1,gw_%q_g_commonsubset_size
    if(gw_%q_g_commonsubset_indinrho(ig)==0) stop ('commonsubset indinrho has 0') 
  enddo
  


  !!!!!!!!!!!!!
  !  convert eps(q) g index to common gw-rho based g index
  !!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  ! prep read gw h5 data
  !!!!!!!!!!!!!!!!!!!!
  
  
  
  
  !select case( h5rank)
  !  case (1)
  !h5dims3=h5dims
  !allocate(gw_%eps0mat_diag_data(h5dims3(1),h5dims3(2),h5dims3(3)))
  !gw_%eps0mat_diag_data=reshape(h5dataset_data,h5dims3)
  !  case default
  !  write(*,*) 'h5 read error'
  !end select 
  
  
  !!!!!!!!!!!!!!!!!!!!!!!
  !!md5sum not working for non-text files
  !    CALL md5_from_file('t.tgz',epsmatf_md5_cksum)
  !    write (*,*) 'GW epsmat files:',trim(eps0mat_filename),'  MD5 sum:',epsmatf_md5_cksum
  !    CALL md5_from_file('t1.tgz',epsmatf_md5_cksum)
  !    write (*,*) 'GW epsmat files:',trim(eps0mat_filename),'  MD5 sum:',epsmatf_md5_cksum
  !
  !    CALL md5_from_file('eps0mat.10-epsilon_subsampling-cutoff10.h5',epsmatf_md5_cksum)
  !    write (*,*) 'GW epsmat files:',trim(eps0mat_filename),'  MD5 sum:',epsmatf_md5_cksum
  !    CALL md5_from_file(eps0mat_filename, epsmatf_md5_cksum)
  !    write (*,*) 'GW epsmat files:',trim(eps0mat_filename),'  MD5 sum:',epsmatf_md5_cksum
  !    CALL md5_from_file(epsmat_filename, epsmatf_md5_cksum)
  !    write (*,*) 'GW epsmat files:',trim(epsmat_filename),'  MD5 sum:',epsmatf_md5_cksum
  !!!!!!!!!!!!!!!!!!!!!!!

  !
  !!!!!! gweps read 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    

end subroutine gw_eps_init

