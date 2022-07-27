 subroutine gw_eps_init(h5filename,gw_ng_data ,gw_nmtx_max_data ,gw_nmtx_data ,gw_gind_eps2rho_data ,gw_gind_rho2eps_data ,&
                                 gw_g_components_data ,gw_bvec_data ,gw_blat_data ,gw_qpts_data ,gw_nq_data ,&
              gw_epsmat_diag_data ,gw_epsmat_full_data ,gw_q_g_commonsubset_indinrho ,gw_q_g_commonsubset_size,gw_qabs)
    
USE kinds, ONLY: DP,sgl
USE HDF5


!!!!!!!!!!hdf5
  CHARACTER(LEN=256) :: h5filename      ! Dataset name
  CHARACTER(LEN=256) :: h5datasetname = "matrix-diagonal"     ! Dataset name
 INTEGER     ::   h5rank,h5error ! Error flag
  !INTEGER     ::  i, j
!  real(dp), DIMENSION(3,1465,2) :: h5dataset_data, data_out ! Data buffers
!  real(dp), DIMENSION(3,1465,2) :: h5dataset_data, data_out ! Data buffers
  real(dp), allocatable :: h5dataset_data_double(:), data_out(:)
  integer, allocatable :: h5dataset_data_integer(:)
  INTEGER(HSIZE_T), allocatable :: h5dims(:),h5maxdims(:)


!  real(dp), allocatable :: gw_epsmat_diag_data(:,:,2),  gw_eps0mat_diag_data(:,:,2)
  real(dp), allocatable :: gw_epsmat_diag_data_q1(:,:,:),  gw_epsmat_diag_data_q0(:,:,:)
  !complex(dp), allocatable :: gw_epsmat_diag_data(:,:,:),  gw_eps0mat_diag_data(:,:,:)
!  real(dp), allocatable :: gw_epsmat_full_data(:,1,1,:,:,2),  gw_eps0mat_full_data(:,1,1,:,:,2)
  real(dp), allocatable :: gw_epsmat_full_data_q1(:,:,:,:,:,:),  gw_epsmat_full_data_q0(:,:,:,:,:,:)
!  real(dp), allocatable :: gw_epsallmat_full_data(:,1,1,:,:,2)
  real(dp), allocatable :: gw_epsmat_full_data_qall(:,:,:,:,:,:)

  real(dp), allocatable :: gw_vcoul_data_q1(:,:),gw_qpts_Data_q1(:,:)
  real(dp), allocatable :: gw_blat_data_q1(:),gw_bvec_Data_q1(:,:)
  integer, allocatable :: gw_gind_eps2rho_data_q1(:,:), gw_gind_rho2eps_data_q1(:,:),gw_nmtx_data_q1(:)

!q0
  real(dp), allocatable :: gw_vcoul_data_q0(:,:),gw_qpts_Data_q0(:,:)
  real(dp), allocatable :: gw_blat_data_q0(:),gw_bvec_Data_q0(:,:)
  integer, allocatable :: gw_gind_eps2rho_data_q0(:,:), gw_gind_rho2eps_data_q0(:,:),gw_nmtx_data_q0(:)
!q0

  integer :: h5dims1(1),h5dims2(2),h5dims3(3),h5dims4(4),h5dims5(5),h5dims6(6)

   integer, allocatable :: gw_grho_data_q1(:),  gw_geps_data_q1(:),gw_g_components_data_q1(:,:)
  integer, allocatable :: gw_nq_data_q1(:),gw_nmtx_max_data_q1(:),gw_fftgrid_data_q1(:),gw_qgrid_data_q1(:),gw_ng_data_q1(:)

!q0
   integer, allocatable :: gw_grho_data_q0(:),  gw_geps_data_q0(:),gw_g_components_data_q0(:,:)
  integer, allocatable :: gw_nq_data_q0(:),gw_nmtx_max_data_q0(:),gw_fftgrid_data_q0(:),gw_qgrid_data_q0(:),gw_ng_data_q0(:)
!q0


!  integer(i8b), allocatable :: gw_nqi8(:)

    real(DP),allocatable ::gw_qabs_q1(:)
    INTEGER :: gw_q_g_commonsubset_size_q1
    integer(DP),allocatable ::gw_q_g_commonsubset_indinrho_q1(:)

!q0
    real(DP),allocatable ::gw_qabs_q0(:)
    INTEGER :: gw_q_g_commonsubset_size_q0
    integer(DP),allocatable ::gw_q_g_commonsubset_indinrho_q0(:)
!q0


!!!!!!!!!!!!!!!!!!
    integer(DP),allocatable ::gind_rho2psi_gw(:)
    real(DP) ::gvec_gw(3)
    integer(DP),allocatable ::gind_psi2rho_gw(:)

    integer(DP),allocatable ::gind_rho2psi_gw_q0(:)
    real(DP) ::gvec_gw_q0(3)
    integer(DP),allocatable ::gind_psi2rho_gw_q0(:)

    integer(DP),allocatable ::gind_rho2psi_gw_q1(:)
    real(DP) ::gvec_gw_q1(3)
    integer(DP),allocatable ::gind_psi2rho_gw_q1(:)
!!!!!!!!!!!!!!!!!!




!gw_ng_data 
!gw_nmtx_max_data 
!gw_nmtx_data 
!gw_gind_eps2rho_data 
!gw_gind_rho2eps_data 
!gw_g_components_data 
!gw_bvec_data 
!gw_blat_data 
!gw_qpts_data 
!gw_nq_data 
!gw_epsmat_diag_data 
!gw_epsmat_full_data 
!gw_q_g_commonsubset_indinrho 
!gw_q_g_commonsubset_size



  !real(dp) ,dimension(:,:), intent (inout) :: gw_vcoul_data,gw_qpts_data
  !real(dp) ,allocatable, intent (inout) :: gw_vcoul_data(:,:),gw_qpts_data(:,:)
  real(dp) ,allocatable, intent (inout) :: gw_qpts_data(:,:)
!  real(dp) ,allocatable, intent (inout) :: gw_vcoul_data(:,:)
  real(dp), allocatable,intent (inout)  :: gw_blat_data(:),gw_bvec_data(:,:)
  integer, allocatable,intent (inout)  :: gw_gind_eps2rho_data(:,:), gw_gind_rho2eps_data(:,:),gw_nmtx_data(:)
   integer, allocatable,intent (inout)  :: gw_g_components_data(:,:)
   integer, allocatable  :: gw_grho_data(:),  gw_geps_data(:)
   !integer, allocatable,intent (inout)  :: gw_grho_data(:),  gw_geps_data(:),gw_g_components_data(:,:)
  integer, allocatable  :: gw_qgrid_data(:),gw_fftgrid_data(:)
  integer, allocatable ,intent (inout) :: gw_nq_data(:),gw_nmtx_max_data(:),gw_ng_data(:)
    
    integer(DP),allocatable ::gw_q_g_commonsubset_indinrhotmp1(:)
    real(DP),allocatable ,intent (inout) ::gw_qabs(:)
    INTEGER ,intent (inout) :: gw_q_g_commonsubset_size
    integer(DP),allocatable ,intent (inout) ::gw_q_g_commonsubset_indinrho(:)
  real(dp), allocatable ,intent (inout) :: gw_epsmat_full_data(:,:,:,:,:,:)

  real(dp), allocatable,intent (inout)  :: gw_epsmat_diag_data(:,:,:)



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALL h5gn_members_f(file_id, "/mats", nmembers, error)
! write(*,*) "Number of root group member is " , nmembers
! do i = 0, nmembers - 1
!    CALL h5gget_obj_info_idx_f(file_id, "/mats", i, name_buffer, dtype, error)
! write(*,*) trim(name_buffer), dtype
! end do



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!! gweps read 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!
! inverse alll dimensions for description
  h5datasetname='/mf_header/crystal/blat'              !f8 
  h5datasetname='/mf_header/crystal/bvec'              !f8 (3,3)

  h5datasetname='/mf_header/gspace/components'         !I4 (ng,3) G pts within cutoff
  h5datasetname='/mf_header/gspace/ng'                 !
  h5datasetname='/mf_header/gspace/FFTgrid'            !i4 (3)
  h5datasetname='/mf_header/gspace/ecutrho'            !
  h5datasetname='/eps_header/gspace/gind_eps2rho'      !i4 (nq,ng)
  h5datasetname='/eps_header/gspace/gind_rho2eps'      !i4 (nq,ng)
  h5datasetname='/eps_header/gspace/nmtx_max'          !i4 
  h5datasetname='/eps_header/gspace/nmtx'              !i4 (nq)  G pts for eps
                                                        
  h5datasetname='/eps_header/gspace/vcoul'             !f8 (nq,nmtx_max)
  h5datasetname='/eps_header/qpoints/nq'               !
  h5datasetname='/eps_header/qpoints/qpts'             !f8 (nq,3)
  h5datasetname='/eps_header/qpoints/qgrid'            !i4 (3)
                                                        
                                                        
  h5datasetname='/mats/matrix'                         !f8 (nq, 1,1, nmtx_max,nmtx_max,2)
  h5datasetname='/mats/matrix-diagonal'                         !f8 (nq, 1,1, nmtx_max,nmtx_max,2)
  !hdf5



h5datasetname='/mf_header/gspace/ng'      !i4 (nq,ng)
call h5gw_read(h5filename,h5datasetname,h5dataset_data_double,h5dataset_Data_integer,h5dims,h5rank,h5error)
if (h5error<0)  write(*,*)  'h5error',h5error
if (h5rank/=1) then
 write(*,*)  'h5rank error(should be 1)',h5rank 
else
 h5dims1=h5dims
 allocate(gw_ng_data(h5dims1(1)))
 gw_ng_data=reshape(h5dataset_data_integer,h5dims1)
 write(*,*)  'shape h5dataset',shape(h5dataset_data_integer)
 write(*,*)  'ng()',gw_ng_data(:)
 deallocate(h5dims)
 deallocate(h5dataset_Data_integer)
 deallocate(h5dataset_Data_double)
endif


h5datasetname='/eps_header/gspace/nmtx_max'      !i4 (nq,ng)
call h5gw_read(h5filename,h5datasetname,h5dataset_data_double,h5dataset_Data_integer,h5dims,h5rank,h5error)
if (h5error<0)  write(*,*)  'h5error',h5error
if (h5rank/=1) then
 write(*,*)  'h5rank error(should be 1)',h5rank 
else
 h5dims1=h5dims
 allocate(gw_nmtx_max_data(h5dims1(1)))
 gw_nmtx_max_data=reshape(h5dataset_data_integer,h5dims1)
 write(*,*)  'shape h5dataset',shape(h5dataset_data_integer)
 write(*,*)  'nmtx_max()',gw_nmtx_max_data(:)
 deallocate(h5dims)
 deallocate(h5dataset_Data_integer)
 deallocate(h5dataset_Data_double)
endif




h5datasetname='/eps_header/gspace/nmtx'      !i4 (nq,ng)
call h5gw_read(h5filename,h5datasetname,h5dataset_data_double,h5dataset_Data_integer,h5dims,h5rank,h5error)
if (h5error<0)  write(*,*)  'h5error',h5error
if (h5rank/=1) then
 write(*,*)  'h5rank error(should be 1)',h5rank 
else
 h5dims1=h5dims
 allocate(gw_nmtx_data(h5dims1(1)))
 gw_nmtx_data=reshape(h5dataset_data_integer,h5dims1)
 write(*,*)  'shape h5dataset',shape(h5dataset_data_integer)
 write(*,*)  'nmtx()',gw_nmtx_data(:)
 deallocate(h5dims)
 deallocate(h5dataset_Data_integer)
 deallocate(h5dataset_Data_double)
endif



h5datasetname='/eps_header/gspace/gind_eps2rho'      !i4 (nq,ng)
call h5gw_read(h5filename,h5datasetname,h5dataset_data_double,h5dataset_Data_integer,h5dims,h5rank,h5error)
if (h5error<0)  write(*,*)  'h5error',h5error
if (h5rank/=2) then
 write(*,*)  'h5rank error(should be 2)',h5rank 
else
 h5dims2=h5dims
 allocate(gw_gind_eps2rho_data(h5dims2(1),h5dims2(2)))
 gw_gind_eps2rho_data=reshape(h5dataset_data_integer,h5dims2)
 write(*,*)  'shape h5dataset',shape(gw_gind_eps2rho_data)
 write(*,*)  'gw_gind_eps2rho_data()',gw_gind_eps2rho_data(1:100,1)
 deallocate(h5dims)
 deallocate(h5dataset_Data_integer)
 deallocate(h5dataset_Data_double)
endif


h5datasetname='/eps_header/gspace/gind_rho2eps'      !i4 (nq,ng)
call h5gw_read(h5filename,h5datasetname,h5dataset_data_double,h5dataset_Data_integer,h5dims,h5rank,h5error)
if (h5error<0)  write(*,*)  'h5error',h5error
if (h5rank/=2) then
 write(*,*)  'h5rank error(should be 2)',h5rank 
else
 h5dims2=h5dims
 allocate(gw_gind_rho2eps_data(h5dims2(1),h5dims2(2)))
 gw_gind_rho2eps_data=reshape(h5dataset_data_integer,h5dims2)
 write(*,*)  'shape h5dataset',shape(h5dataset_data_integer)
 write(*,*)  'gw_gind_rho2eps_data()',gw_gind_rho2eps_data(1:100,1)
 deallocate(h5dims)
 deallocate(h5dataset_Data_integer)
 deallocate(h5dataset_Data_double)
endif




h5datasetname='/mf_header/gspace/components'               !
call h5gw_read(h5filename,h5datasetname,h5dataset_data_double,h5dataset_Data_integer,h5dims,h5rank,h5error)
if (h5error<0)  write(*,*)  'h5error',h5error
if (h5rank/=2) then
 write(*,*)  'h5rank error(should be 2)',h5rank 
else
 h5dims2=h5dims
 allocate(gw_g_components_data(h5dims2(1),h5dims2(2)))
 gw_g_components_data=reshape(h5dataset_data_integer,h5dims2)
 write(*,*)  'shape h5dataset',shape(h5dataset_data_integer)
 write(*,*)  'gw_g_components_data()',gw_g_components_data(:,1:7)
 deallocate(h5dims)
 deallocate(h5dataset_Data_integer)
 deallocate(h5dataset_Data_double)
endif



h5datasetname='/mf_header/crystal/bvec'               !
call h5gw_read(h5filename,h5datasetname,h5dataset_data_double,h5dataset_Data_integer,h5dims,h5rank,h5error)
if (h5error<0)  write(*,*)  'h5error',h5error
if (h5rank/=2) then
 write(*,*)  'h5rank error(should be 2)',h5rank 
else
 h5dims2=h5dims
 allocate(gw_bvec_data(h5dims2(1),h5dims2(2)))
 gw_bvec_data=reshape(h5dataset_data_double,h5dims2)
 write(*,*)  'shape h5dataset',shape(h5dataset_data_double)
 write(*,*)  'gw_bvec_data()',gw_bvec_data(:,:)
 deallocate(h5dims)
 deallocate(h5dataset_Data_integer)
 deallocate(h5dataset_Data_double)
endif



h5datasetname='/mf_header/crystal/blat'               !
call h5gw_read(h5filename,h5datasetname,h5dataset_data_double,h5dataset_Data_integer,h5dims,h5rank,h5error)
if (h5error<0)  write(*,*)  'h5error',h5error
if (h5rank/=1) then
 write(*,*)  'h5rank error(should be 1)',h5rank 
else
 h5dims1=h5dims
 allocate(gw_blat_data(h5dims1(1)))
 gw_blat_data=reshape(h5dataset_data_double,h5dims1)
 write(*,*)  'shape h5dataset',shape(h5dataset_data_double)
 write(*,*)  'gw_blat_data()',gw_blat_data(:)
 deallocate(h5dims)
 deallocate(h5dataset_Data_integer)
 deallocate(h5dataset_Data_double)
endif




h5datasetname='/eps_header/qpoints/qpts'               !
call h5gw_read(h5filename,h5datasetname,h5dataset_data_double,h5dataset_Data_integer,h5dims,h5rank,h5error)
if (h5error<0)  write(*,*)  'h5error',h5error
if (h5rank/=2) then
 write(*,*)  'h5rank error(should be 2)',h5rank 
else
 h5dims2=h5dims
 allocate(gw_qpts_data(h5dims2(1),h5dims2(2)))
 gw_qpts_data=reshape(h5dataset_data_double,h5dims2)
 write(*,*)  'shape h5dataset',shape(h5dataset_data_double)
 write(*,*)  'gw_qpts_data()',gw_qpts_data(:,:)
 deallocate(h5dims)
 deallocate(h5dataset_Data_integer)
 deallocate(h5dataset_Data_double)
endif



h5datasetname='/eps_header/qpoints/nq'               !
call h5gw_read(h5filename,h5datasetname,h5dataset_data_double,h5dataset_Data_integer,h5dims,h5rank,h5error)
if (h5error<0)  write(*,*)  'h5error',h5error
if (h5rank/=1) then
 write(*,*)  'h5rank error(should be 3)',h5rank 
else
          !      write(*,*) 'sizeof(int(i4b)):',sizeof(gw_nq)
!                write(*,*) 'sizeof(int(i8b)):',sizeof(gw_nqi8)
                write(*,*) 'sizeof(int):',sizeof(h5rank)
 h5dims1=h5dims
 allocate(gw_nq_data(h5dims1(1)))
 gw_nq_data=reshape(h5dataset_data_integer,h5dims1)
 write(*,*)  'shape h5dataset',shape(h5dataset_data_integer)
 write(*,*)  'gw_nq_data()',gw_nq_data(:)
 deallocate(h5dims)
 deallocate(h5dataset_Data_integer)
 deallocate(h5dataset_Data_double)
endif


h5datasetname='/mats/matrix-diagonal'                         !f8 (nq, 1,1, nmtx_max,nmtx_max,2)
call h5gw_read(h5filename,h5datasetname,h5dataset_data_double,h5dataset_Data_integer,h5dims,h5rank,h5error)
if (h5error<0)  write(*,*)  'h5error',h5error
if (h5rank/=3) then
 write(*,*)  'h5rank error(should be 3)',h5rank 
else
 h5dims3=h5dims
 allocate(gw_epsmat_diag_data(h5dims3(1),h5dims3(2),h5dims3(3)))
 gw_epsmat_diag_data=reshape(h5dataset_data_double,h5dims3)
 write(*,*)  'shape h5dataset',shape(h5dataset_data_double)
 write(*,*)  'gw_epsmat_diag_data(:,1,1)',gw_epsmat_diag_data(:,1,:)
 deallocate(h5dims)
 deallocate(h5dataset_Data_integer)
 deallocate(h5dataset_Data_double)
endif

h5datasetname='/mats/matrix'                         !f8 (nq, 1,1, nmtx_max,nmtx_max,2)
call h5gw_read(h5filename,h5datasetname,h5dataset_data_double,h5dataset_Data_integer,h5dims,h5rank,h5error)
if (h5error<0)  write(*,*)  'h5error',h5error
if (h5rank/=6) then
 write(*,*)  'h5rank error(should be 6)',h5rank 
else
 h5dims6=h5dims
 write(*,*)  'hdims',h5dims 
 allocate(gw_epsmat_full_data(h5dims6(1),h5dims6(2),h5dims6(3),h5dims6(4),h5dims6(5),h5dims6(6)))
 gw_epsmat_full_data=reshape(h5dataset_data_double,h5dims6)
 write(*,*)  'shape h5dataset',shape(h5dataset_data_double)
 write(*,*)  'gw_epsmat_full_data(:,1,1)diag',gw_epsmat_full_data(:,1,1,1,1,1)
 write(*,*)  'gw_epsmat_full_data(:,1,1)diag',gw_epsmat_full_data(:,1,1,1,1,2)
 write(*,*)  'gw_epsmat_full_data(:,1,1)diag',gw_epsmat_full_data(:,1,1,1,1,3)
 write(*,*)  'gw_epsmat_full_data(:,1,1)wing',gw_epsmat_full_data(:,1,1,1,1,1)
 write(*,*)  'gw_epsmat_full_data(:,1,1)wing',gw_epsmat_full_data(:,1,2,1,1,1)
 write(*,*)  'gw_epsmat_full_data(:,1,1)wing',gw_epsmat_full_data(:,1,3,1,1,1)
 write(*,*)  'gw_epsmat_full_data(:,1,1)wing',gw_epsmat_full_data(:,1,4,1,1,1)
 write(*,*)  'gw_epsmat_full_data(:,1,1)wing',gw_epsmat_full_data(:,1,5,1,1,1)
 write(*,*)  'gw_epsmat_full_data(:,1,1)wing',gw_epsmat_full_data(:,1,6,1,1,1)
endif


!!!!!!!!!!!!!!!!!!!
! prep read gw h5 data

! qabs
write(*,*) gw_bvec_data(:,1)
write(*,*) gw_bvec_data(:,2)
write(*,*) gw_bvec_data(:,3)
write(*,*) gw_qpts_data(:,1)
write(*,*) gw_qpts_data(:,2)
write(*,*) gw_qpts_data(:,3)

    allocate(gw_qabs(gw_nq_data(1)))
    do ig1 = 1, gw_nq_data(1)
      gw_qabs(ig1)=norm2(&
              gw_qpts_data(1,ig1)*gw_bvec_data(:,1)+ &
              gw_qpts_data(2,ig1)*gw_bvec_data(:,2)+ &
              gw_qpts_data(3,ig1)*gw_bvec_data(:,3))
!*gw_blat_data(1)

!write(*,*)              gw_qpts_data(1,ig1)*gw_bvec_data(1,:)
!write(*,*)              gw_qpts_data(2,ig1)*gw_bvec_data(2,:)
!write(*,*)              gw_qpts_data(3,ig1)*gw_bvec_data(3,:)


!debug
write(*,*)'gw_qabs debug ', gw_qabs(ig1),gw_epsmat_diag_data(:,1,ig1)
!debug
    enddo

!!!!!!!!!!!!!!
!  convert eps(q) g index to common gw-rho based g index
!     gw_q_g_commonsubset_size
!    gw_q_g_commonsubset2rho(:,:)
!    do ig = 1, gw_ng_data(1)
!      do iq=1,gw_nq_data(1)
!        gind_gw_eps=gw_gind_rho2eps_data(iq,ig)
!           if      (gind_gw_eps<gw_nmtx(iq))  then
!      enddo
!    enddo
!eps(gw_gind_rho2eps_data(iq,1:gw_nmtx_data(iq)))

allocate(gw_q_g_commonsubset_indinrho(gw_nmtx_max_data(1)))
gw_q_g_commonsubset_indinrho(:)=0
gw_q_g_commonsubset_indinrho(:)=gw_gind_eps2rho_data(1:gw_nmtx_data(1),1)

!write(*,*)  'gw_q_g_commonsubset_indinrho',gw_q_g_commonsubset_indinrho(1:10),shape(gw_q_g_commonsubset_indinrho)

do iq=1,gw_nq_data(1)
  do ig=1,gw_nmtx_max_data(1)
    if(gw_q_g_commonsubset_indinrho(ig)>0) then
      if (gw_gind_rho2eps_data(gw_q_g_commonsubset_indinrho(ig),iq)>gw_nmtx_data(iq) ) then
         gw_q_g_commonsubset_indinrho(ig)=0
       endif
    endif
  enddo
enddo
!write(*,*)  'gw_q_g_commonsubset_indinrho',gw_q_g_commonsubset_indinrho(:)
ig=0
  do ig1=1,gw_nmtx_max_data(1)
    if(gw_q_g_commonsubset_indinrho(ig1)>0) ig=ig+1
  enddo

write(*,*)  'gw_q_g_commonsubset_indinrho',gw_q_g_commonsubset_indinrho(:)
allocate(gw_q_g_commonsubset_indinrhotmp1(ig))
ig1=1
do ig=1,gw_nmtx_max_data(1)
  if(gw_q_g_commonsubset_indinrho(ig)>0) then 
!     write(*,*) gw_q_g_commonsubset_indinrhotmp1(ig1),gw_q_g_commonsubset_indinrho(ig)
     gw_q_g_commonsubset_indinrhotmp1(ig1)=gw_q_g_commonsubset_indinrho(ig) 
     ig1=ig1+1
  endif
enddo
deallocate(gw_q_g_commonsubset_indinrho)
allocate(gw_q_g_commonsubset_indinrho(size(gw_q_g_commonsubset_indinrhotmp1)))
gw_q_g_commonsubset_indinrho(:)=gw_q_g_commonsubset_indinrhotmp1(:) 

write(*,*)  'gw_q_g_commonsubset_indinrho',gw_q_g_commonsubset_indinrho(:),shape(gw_q_g_commonsubset_indinrho)
gw_q_g_commonsubset_size=size(gw_q_g_commonsubset_indinrho)
!  convert eps(q) g index to common gw-rho based g index
!!!!!!!!!!!!!


! prep read gw h5 data
!!!!!!!!!!!!!!!!!!!




!select case( h5rank)
!  case (1)
!h5dims3=h5dims
!allocate(gw_eps0mat_diag_data(h5dims3(1),h5dims3(2),h5dims3(3)))
!gw_eps0mat_diag_data=reshape(h5dataset_data,h5dims3)
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
    
