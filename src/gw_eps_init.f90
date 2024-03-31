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
  real(dp), allocatable :: h5dataset_data_double(:), data_out(:)
  integer, allocatable :: h5dataset_data_integer(:)
  INTEGER(HSIZE_T), allocatable :: h5dims(:),h5maxdims(:)
  
  integer :: h5dims1(1),h5dims2(2),h5dims3(3),h5dims4(4),h5dims5(5),h5dims6(6)
  integer:: p_rank,p_size,ik


  if (  allocated(gw_%qabs)) then
     deallocate(gw_%qabs)
  endif
  allocate(gw_%qabs(gw_%nq_data(1)))
  do ig1 = 1, gw_%nq_data(1)
      gw_%qabs(ig1)=norm2(&
              gw_%qpts_data(1,ig1)*gw_%bvec_data(:,1)+ &
              gw_%qpts_data(2,ig1)*gw_%bvec_data(:,2)+ &
              gw_%qpts_data(3,ig1)*gw_%bvec_data(:,3))
  enddo
  if (  allocated(gw_%q_g_commonsubset_indinrho)) then
      deallocate(gw_%q_g_commonsubset_indinrho)
  endif
  allocate(gw_%q_g_commonsubset_indinrho(gw_%nmtx_max_data(1)))
  gw_%q_g_commonsubset_indinrho(:)=0
  gw_%q_g_commonsubset_indinrho(:)=gw_%gind_eps2rho_data(1:gw_%nmtx_data(1),1)

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


  ig=0
  do ig1=1,gw_%nmtx_max_data(1)
    if(gw_%q_g_commonsubset_indinrho(ig1)>0) ig=ig+1
  enddo

  if (  allocated(gw_%q_g_commonsubset_indinrhotmp1)) then
     deallocate(gw_%q_g_commonsubset_indinrhotmp1)
  endif
  allocate(gw_%q_g_commonsubset_indinrhotmp1(ig))
  gw_%q_g_commonsubset_indinrhotmp1(:)=0


  ig1=1
  do ig=1,gw_%nmtx_max_data(1)
    if(gw_%q_g_commonsubset_indinrho(ig)>0) then 
       gw_%q_g_commonsubset_indinrhotmp1(ig1)=gw_%q_g_commonsubset_indinrho(ig) 
       ig1=ig1+1
    endif
  enddo
  deallocate(gw_%q_g_commonsubset_indinrho)
  allocate(gw_%q_g_commonsubset_indinrho(size(gw_%q_g_commonsubset_indinrhotmp1)))
  gw_%q_g_commonsubset_indinrho(:)=gw_%q_g_commonsubset_indinrhotmp1(:) 
  
  gw_%q_g_commonsubset_size=size(gw_%q_g_commonsubset_indinrho)
  do ig=1,gw_%q_g_commonsubset_size
    if(gw_%q_g_commonsubset_indinrho(ig)==0) stop ('commonsubset indinrho has 0') 
  enddo
  
end subroutine gw_eps_init

