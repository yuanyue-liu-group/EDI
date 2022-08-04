 !subroutine gw_eps_init(h5filename,gw_ng_data ,gw_nmtx_max_data ,gw_nmtx_data ,gw_gind_eps2rho_data ,gw_gind_rho2eps_data ,&
 !                                gw_g_components_data ,gw_bvec_data ,gw_blat_data ,gw_qpts_data ,gw_nq_data ,&
 !             gw_epsmat_diag_data ,gw_epsmat_full_data ,gw_q_g_commonsubset_indinrho ,gw_q_g_commonsubset_size,gw_qabs)
    
 subroutine gw_eps_init(eps_filename_,gw_)
USE kinds, ONLY: DP,sgl
USE HDF5
use edic_mod, only: gw_eps_data
   !   Use edic_mod,   only: gw_epsq1_data,gw_epsq0_data

CHARACTER(LEN=256) :: eps_filename_
type(gw_eps_data),intent (inout) ,target:: gw_

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

  integer :: h5dims1(1),h5dims2(2),h5dims3(3),h5dims4(4),h5dims5(5),h5dims6(6)
!
!!!  real(dp), allocatable :: gw_epsmat_diag_data(:,:,2),  gw_eps0mat_diag_data(:,:,2)
!!  real(dp), allocatable :: gw_epsmat_diag_data_q1(:,:,:),  gw_epsmat_diag_data_q0(:,:,:)
!!  !complex(dp), allocatable :: gw_epsmat_diag_data(:,:,:),  gw_eps0mat_diag_data(:,:,:)
!!!  real(dp), allocatable :: gw_epsmat_full_data(:,1,1,:,:,2),  gw_eps0mat_full_data(:,1,1,:,:,2)
!!  real(dp), allocatable :: gw_epsmat_full_data_q1(:,:,:,:,:,:),  gw_epsmat_full_data_q0(:,:,:,:,:,:)
!!!  real(dp), allocatable :: gw_epsallmat_full_data(:,1,1,:,:,2)
!!  real(dp), allocatable :: gw_epsmat_full_data_qall(:,:,:,:,:,:)
!!
!!  real(dp), allocatable :: gw_vcoul_data_q1(:,:),gw_qpts_Data_q1(:,:)
!!  real(dp), allocatable :: gw_blat_data_q1(:),gw_bvec_Data_q1(:,:)
!!  integer, allocatable :: gw_gind_eps2rho_data_q1(:,:), gw_gind_rho2eps_data_q1(:,:),gw_nmtx_data_q1(:)
!!
!!!q0
!!  real(dp), allocatable :: gw_vcoul_data_q0(:,:),gw_qpts_Data_q0(:,:)
!!  real(dp), allocatable :: gw_blat_data_q0(:),gw_bvec_Data_q0(:,:)
!!  integer, allocatable :: gw_gind_eps2rho_data_q0(:,:), gw_gind_rho2eps_data_q0(:,:),gw_nmtx_data_q0(:)
!!!q0
!
!
!!   integer, allocatable :: gw_grho_data_q1(:),  gw_geps_data_q1(:),gw_g_components_data_q1(:,:)
!!  integer, allocatable :: gw_nq_data_q1(:),gw_nmtx_max_data_q1(:),gw_fftgrid_data_q1(:),gw_qgrid_data_q1(:),gw_ng_data_q1(:)
!!
!!!q0
!!   integer, allocatable :: gw_grho_data_q0(:),  gw_geps_data_q0(:),gw_g_components_data_q0(:,:)
!!  integer, allocatable :: gw_nq_data_q0(:),gw_nmtx_max_data_q0(:),gw_fftgrid_data_q0(:),gw_qgrid_data_q0(:),gw_ng_data_q0(:)
!!!q0
!!
!!
!!!  integer(i8b), allocatable :: gw_nqi8(:)
!!
!!    real(DP),allocatable ::gw_qabs_q1(:)
!!    INTEGER :: gw_q_g_commonsubset_size_q1
!!    integer(DP),allocatable ::gw_q_g_commonsubset_indinrho_q1(:)
!!
!!!q0
!!    real(DP),allocatable ::gw_qabs_q0(:)
!!    INTEGER :: gw_q_g_commonsubset_size_q0
!!    integer(DP),allocatable ::gw_q_g_commonsubset_indinrho_q0(:)
!!!q0
!!
!!
!!!!!!!!!!!!!!!!!!!!
!!    integer(DP),allocatable ::gind_rho2psi_gw(:)
!!    real(DP) ::gvec_gw(3)
!!    integer(DP),allocatable ::gind_psi2rho_gw(:)
!!
!!    integer(DP),allocatable ::gind_rho2psi_gw_q0(:)
!!    real(DP) ::gvec_gw_q0(3)
!!    integer(DP),allocatable ::gind_psi2rho_gw_q0(:)
!!
!!    integer(DP),allocatable ::gind_rho2psi_gw_q1(:)
!!    real(DP) ::gvec_gw_q1(3)
!!    integer(DP),allocatable ::gind_psi2rho_gw_q1(:)
!!!!!!!!!!!!!!!!!!!!
!
!
!
!!
!!gw_ng_data 
!!gw_nmtx_max_data 
!!gw_nmtx_data 
!!gw_gind_eps2rho_data 
!!gw_gind_rho2eps_data 
!!gw_g_components_data 
!!gw_bvec_data 
!!gw_blat_data 
!!gw_qpts_data 
!!gw_nq_data 
!!gw_epsmat_diag_data 
!!gw_epsmat_full_data 
!!gw_q_g_commonsubset_indinrho 
!!gw_q_g_commonsubset_size
!
!
!
!  !real(dp) ,dimension(:,:), intent (inout) :: gw_vcoul_data,gw_qpts_data
!  !real(dp) ,allocatable, intent (inout) :: gw_vcoul_data(:,:),gw_qpts_data(:,:)
!  real(dp) ,allocatable, intent (inout) :: gw_qpts_data(:,:)
!!  real(dp) ,allocatable, intent (inout) :: gw_vcoul_data(:,:)
!  real(dp), allocatable,intent (inout)  :: gw_blat_data(:),gw_bvec_data(:,:)
!  integer, allocatable,intent (inout)  :: gw_gind_eps2rho_data(:,:), gw_gind_rho2eps_data(:,:),gw_nmtx_data(:)
!   integer, allocatable,intent (inout)  :: gw_g_components_data(:,:)
!   integer, allocatable  :: gw_grho_data(:),  gw_geps_data(:)
!   !integer, allocatable,intent (inout)  :: gw_grho_data(:),  gw_geps_data(:),gw_g_components_data(:,:)
!  integer, allocatable  :: gw_qgrid_data(:),gw_fftgrid_data(:)
!  integer, allocatable ,intent (inout) :: gw_nq_data(:),gw_nmtx_max_data(:),gw_ng_data(:)
!    
!    integer(DP),allocatable ::gw_q_g_commonsubset_indinrhotmp1(:)
!    real(DP),allocatable ,intent (inout) ::gw_qabs(:)
!    INTEGER ,intent (inout) :: gw_q_g_commonsubset_size
!    integer(DP),allocatable ,intent (inout) ::gw_q_g_commonsubset_indinrho(:)
!  real(dp), allocatable ,intent (inout) :: gw_epsmat_full_data(:,:,:,:,:,:)
!
!  real(dp), allocatable,intent (inout)  :: gw_epsmat_diag_data(:,:,:)
!
!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALL h5gn_members_f(file_id, "/mats", nmembers, error)
! write(*,*) "Number of root group member is " , nmembers
! do i = 0, nmembers - 1
!    CALL h5gget_obj_info_idx_f(file_id, "/mats", i, name_buffer, dtype, error)
! write(*,*) trim(name_buffer), dtype
! end do


h5filename=eps_filename_      ! Dataset name

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
 allocate(gw_%ng_data(h5dims1(1)))
 gw_%ng_data=reshape(h5dataset_data_integer,h5dims1)
 write(*,*)  'shape h5dataset',shape(h5dataset_data_integer)
 write(*,*)  'ng()',gw_%ng_data(:)
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
 allocate(gw_%nmtx_max_data(h5dims1(1)))
 gw_%nmtx_max_data=reshape(h5dataset_data_integer,h5dims1)
 write(*,*)  'shape h5dataset',shape(h5dataset_data_integer)
 write(*,*)  'nmtx_max()',gw_%nmtx_max_data(:)
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
 allocate(gw_%nmtx_data(h5dims1(1)))
 gw_%nmtx_data=reshape(h5dataset_data_integer,h5dims1)
 write(*,*)  'shape h5dataset',shape(h5dataset_data_integer)
 write(*,*)  'nmtx()',gw_%nmtx_data(:)
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
 allocate(gw_%gind_eps2rho_data(h5dims2(1),h5dims2(2)))
 gw_%gind_eps2rho_data=reshape(h5dataset_data_integer,h5dims2)
 write(*,*)  'shape h5dataset',shape(gw_%gind_eps2rho_data)
 write(*,*)  'gw_gind_eps2rho_data()',gw_%gind_eps2rho_data(1:100,1)
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
 allocate(gw_%gind_rho2eps_data(h5dims2(1),h5dims2(2)))
 gw_%gind_rho2eps_data=reshape(h5dataset_data_integer,h5dims2)
 write(*,*)  'shape h5dataset',shape(h5dataset_data_integer)
 write(*,*)  'gw_gind_rho2eps_data()',gw_%gind_rho2eps_data(1:100,1)
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
 allocate(gw_%g_components_data(h5dims2(1),h5dims2(2)))
 gw_%g_components_data=reshape(h5dataset_data_integer,h5dims2)
 write(*,*)  'shape h5dataset',shape(h5dataset_data_integer)
 write(*,*)  'gw_g_components_data()',gw_%g_components_data(:,1:7)
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
 allocate(gw_%bvec_data(h5dims2(1),h5dims2(2)))
 gw_%bvec_data=reshape(h5dataset_data_double,h5dims2)
 write(*,*)  'shape h5dataset',shape(h5dataset_data_double)
 write(*,*)  'gw_bvec_data()',gw_%bvec_data(:,:)
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
 allocate(gw_%blat_data(h5dims1(1)))
 gw_%blat_data=reshape(h5dataset_data_double,h5dims1)
 write(*,*)  'shape h5dataset',shape(h5dataset_data_double)
 write(*,*)  'gw_blat_data()',gw_%blat_data(:)
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
 allocate(gw_%qpts_data(h5dims2(1),h5dims2(2)))
 gw_%qpts_data=reshape(h5dataset_data_double,h5dims2)
 write(*,*)  'shape h5dataset',shape(h5dataset_data_double)
 write(*,*)  'gw_qpts_data()',gw_%qpts_data(:,:)
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
          !      write(*,*) 'sizeof(int(i4b)):',sizeof(gw_%nq)
!                write(*,*) 'sizeof(int(i8b)):',sizeof(gw_nqi8)
                write(*,*) 'sizeof(int):',sizeof(h5rank)
 h5dims1=h5dims
 allocate(gw_%nq_data(h5dims1(1)))
 gw_%nq_data=reshape(h5dataset_data_integer,h5dims1)
 write(*,*)  'shape h5dataset',shape(h5dataset_data_integer)
 write(*,*)  'gw_nq_data()',gw_%nq_data(:)
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
 allocate(gw_%epsmat_diag_data(h5dims3(1),h5dims3(2),h5dims3(3)))
 gw_%epsmat_diag_data=reshape(h5dataset_data_double,h5dims3)
 write(*,*)  'shape h5dataset',shape(h5dataset_data_double)
 write(*,*)  'gw_epsmat_diag_data(:,1,1)',gw_%epsmat_diag_data(:,1,:)
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
 allocate(gw_%epsmat_full_data(h5dims6(1),h5dims6(2),h5dims6(3),h5dims6(4),h5dims6(5),h5dims6(6)))
 gw_%epsmat_full_data=reshape(h5dataset_data_double,h5dims6)
 write(*,*)  'shape h5dataset',shape(h5dataset_data_double)
 write(*,*)  'gw_epsmat_full_data(:,1,1)diag',gw_%epsmat_full_data(:,1,1,1,1,1)
 write(*,*)  'gw_epsmat_full_data(:,1,1)diag',gw_%epsmat_full_data(:,1,1,1,1,2)
 write(*,*)  'gw_epsmat_full_data(:,1,1)diag',gw_%epsmat_full_data(:,1,1,1,1,3)
 write(*,*)  'gw_epsmat_full_data(:,1,1)wing',gw_%epsmat_full_data(:,1,1,1,1,1)
 write(*,*)  'gw_epsmat_full_data(:,1,1)wing',gw_%epsmat_full_data(:,1,2,1,1,1)
 write(*,*)  'gw_epsmat_full_data(:,1,1)wing',gw_%epsmat_full_data(:,1,3,1,1,1)
 write(*,*)  'gw_epsmat_full_data(:,1,1)wing',gw_%epsmat_full_data(:,1,4,1,1,1)
 write(*,*)  'gw_epsmat_full_data(:,1,1)wing',gw_%epsmat_full_data(:,1,5,1,1,1)
 write(*,*)  'gw_epsmat_full_data(:,1,1)wing',gw_%epsmat_full_data(:,1,6,1,1,1)
endif


!!!!!!!!!!!!!!!!!!!
! prep read gw h5 data

! qabs
write(*,*) gw_%bvec_data(:,1)
write(*,*) gw_%bvec_data(:,2)
write(*,*) gw_%bvec_data(:,3)
write(*,*) gw_%qpts_data(:,1)
write(*,*) gw_%qpts_data(:,2)
write(*,*) gw_%qpts_data(:,3)

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

!!!!!!!!!!!!!!
!  convert eps(q) g index to common gw-rho based g index
!     gw_%q_g_commonsubset_size
!    gw_%q_g_commonsubset2rho(:,:)
!    do ig = 1, gw_%ng_data(1)
!      do iq=1,gw_%nq_data(1)
!        gind_gw_%eps=gw_%gind_rho2eps_data(iq,ig)
!           if      (gindgw_%_eps<gw_%nmtx(iq))  then
!      enddo
!    enddo
!eps(gw_%gind_rho2eps_data(iq,1:gw_%nmtx_data(iq)))

allocate(gw_%q_g_commonsubset_indinrho(gw_%nmtx_max_data(1)))
gw_%q_g_commonsubset_indinrho(:)=0
gw_%q_g_commonsubset_indinrho(:)=gw_%gind_eps2rho_data(1:gw_%nmtx_data(1),1)

!write(*,*)  'gw_%q_g_commonsubset_indinrho',gw_%q_g_commonsubset_indinrho(1:10),shape(gw_%q_g_commonsubset_indinrho)

do iq=1,gw_%nq_data(1)
  do ig=1,gw_%nmtx_max_data(1)
    if(gw_%q_g_commonsubset_indinrho(ig)>0) then
      if (gw_%gind_rho2eps_data(gw_%q_g_commonsubset_indinrho(ig),iq)>gw_%nmtx_data(iq) ) then
         gw_%q_g_commonsubset_indinrho(ig)=0
       endif
    endif
  enddo
enddo
!write(*,*)  'gw_%q_g_commonsubset_indinrho',gw_%q_g_commonsubset_indinrho(:)
ig=0
  do ig1=1,gw_%nmtx_max_data(1)
    if(gw_%q_g_commonsubset_indinrho(ig1)>0) ig=ig+1
  enddo

write(*,*)  'gw_%q_g_commonsubset_indinrho',gw_%q_g_commonsubset_indinrho(:)
allocate(gw_%q_g_commonsubset_indinrhotmp1(ig))
ig1=1
do ig=1,gw_%nmtx_max_data(1)
  if(gw_%q_g_commonsubset_indinrho(ig)>0) then 
!     write(*,*) gw_q_g_commonsubset_indinrhotmp1(ig1),gw_%q_g_commonsubset_indinrho(ig)
     gw_%q_g_commonsubset_indinrhotmp1(ig1)=gw_%q_g_commonsubset_indinrho(ig) 
     ig1=ig1+1
  endif
enddo
deallocate(gw_%q_g_commonsubset_indinrho)
allocate(gw_%q_g_commonsubset_indinrho(size(gw_%q_g_commonsubset_indinrhotmp1)))
gw_%q_g_commonsubset_indinrho(:)=gw_%q_g_commonsubset_indinrhotmp1(:) 

write(*,*)  'gw_%q_g_commonsubset_indinrho',gw_%q_g_commonsubset_indinrho(:),shape(gw_%q_g_commonsubset_indinrho)
gw_%q_g_commonsubset_size=size(gw_%q_g_commonsubset_indinrho)
!  convert eps(q) g index to common gw-rho based g index
!!!!!!!!!!!!!


! prep read gw h5 data
!!!!!!!!!!!!!!!!!!!




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
    
    

contains
subroutine h5gw_read(h5filename,h5datasetname,h5dataset_data_double,h5dataset_Data_integer,h5dims,h5rank,h5error)
USE kinds, ONLY: DP,sgl
USE HDF5
  CHARACTER(LEN=256) :: h5groupname = "/mats"     ! Dataset name
  CHARACTER(LEN=256) :: h5name_buffer 
  INTEGER(HID_T) :: h5file_id       ! File identifier
  INTEGER(HID_T) :: h5group_id       ! Dataset identifier
  INTEGER(HID_T) :: h5dataset_id       ! Dataset identifier
  INTEGER(HID_T) :: h5datatype_id       ! Dataset identifier
  INTEGER(HID_T) :: h5dataspace_id

  INTEGER :: h5dataype       ! Dataset identifier
 
  CHARACTER(LEN=256), intent(in) :: h5filename      ! Dataset name
  CHARACTER(LEN=256) , intent(in) :: h5datasetname      ! Dataset name
  real(dp), allocatable , intent(inout) :: h5dataset_data_double(:)
  real(dp), allocatable :: data_out(:)
  integer, allocatable , intent(inout) :: h5dataset_data_integer(:)
  LOGICAL :: h5flag,h5flag_integer,h5flag_double           ! TRUE/FALSE flag to indicate 
  INTEGER     ::  h5nmembers,i,h5datasize
  INTEGER(HSIZE_T), allocatable :: h5maxdims(:)
  INTEGER(HSIZE_T), allocatable , intent(inout) :: h5dims(:)
  INTEGER  , intent(inout)    ::   h5rank
  INTEGER  , intent(inout)    ::   h5error ! Error flag
  INTEGER(HID_T) :: file_s1_t,h5_file_datatype 
  INTEGER(HID_T) :: mem_s1_t  ,h5_mem_datatype  
  INTEGER(HID_T) :: debugflag=00
  CALL h5open_f(h5error)
  if (h5error<debugflag) then
    write(*,*)  'h5error',h5error
  elseif (h5error<0) then 
    return(h5error)
  endif
  
    !h5 file
    CALL h5fopen_f (h5filename, H5F_ACC_RDWR_F, h5file_id, h5error)
    if (h5error<debugflag) then
      write(*,*)  'h5error',       h5error,trim(h5filename),'h5file_id', h5file_id
    elseif (h5error<0)  then
      return(h5error)
    endif
      !dataset
      CALL h5dopen_f(h5file_id,   trim(h5datasetname), h5dataset_id, h5error)
      if (h5error<debugflag) then
        write(*,*)  'h5error',       h5error, trim(h5datasetname),'h5dataset_id', h5dataset_id
      elseif (h5error<0)  then
        return(h5error)
      endif
        ! dataspace
        call h5dget_space_f(h5dataset_id, h5dataspace_id,  h5error) 
        if (h5error<debugflag) then
          write(*,*)  'h5error',       h5error,'h5dataspace_id',h5dataspace_id
        elseif (h5error<0)  then
          return(h5error)
        endif
          ! rank and shape
          call h5sget_simple_extent_ndims_f(h5dataspace_id, h5rank, h5error) 
          if (h5error<debugflag) then
            write(*,*)  'h5error',       h5error,'h5rank',h5rank
!            if (h5rank>0) write(*,*)   h5dims,h5maxdims
          elseif (h5error<0)  then
            return(h5error)
          endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111
! rank=0 scalar
          if(h5rank==0) then
            !write(*,*)  'h5error',       h5error,'h5rank',h5rank
            h5rank=1
!            write(*,*)  'h5error',       h5error,'h5rank',h5rank
            allocate(h5maxdims(h5rank))
!            write(*,*)  'h5error',       h5error,'h5dimssize',size(h5maxdims)
!            write(*,*)  'h5error',       h5error,'h5dimssize',size(h5dims)
            allocate(h5dims(h5rank))
            h5maxdims(1)=1
            h5dims(1)=1
            h5datasize=1
            do i =1,h5rank
              h5datasize=h5datasize*h5dims(i)
            enddo
            allocate(h5dataset_data_integer(1))
            allocate(h5dataset_data_double(1))
            !allocate(gw_nq(1))
            call H5Dget_type_f(h5dataset_id, h5_file_datatype, h5error);
            if (h5error<debugflag) then
              write(*,*)  'h5error',       h5error,'h5_file_datatype',h5_file_datatype
            elseif (h5error<0)  then
              return(h5error)
            endif
!!!!!!!!!!!!!!!
!debug comment out ok
            ! datatype of memory data, test datatype
            call H5Tget_native_type_f(h5_file_datatype,H5T_DIR_ASCEND_F, h5_mem_datatype,h5error)
              if (h5error<debugflag) then
                write(*,*)  'h5error',       h5error,'h5_mem_datatype',h5_mem_datatype
              elseif (h5error<0)  then
                return(h5error)
              endif
              call h5tequal_F(h5_mem_datatype,H5T_NATIVE_integer,h5flag,h5error)
              if (h5error<debugflag) then
                write(*,*)  'h5error',       h5error,'h5_mem_datatype',h5_mem_datatype,'H5T_NATIVE_integer'
              elseif (h5error<0)  then
                return(h5error)
              endif
              call h5tequal_F(h5_file_datatype,H5T_NATIVE_integer,h5flag,h5error)
              if (h5error<debugflag) then
                write(*,*)  'h5error',       h5error,'h5_file_datatype',h5_file_datatype,h5flag
              elseif (h5error<0)  then
                return(h5error)
              endif
!debug comment out ok
!!!!!!!!!!!!!!!

            call h5tequal_F(h5_file_datatype,H5T_NATIVE_integer,h5flag_integer,h5error)
            call h5tequal_F(h5_file_datatype,H5T_NATIVE_double,h5flag_double,h5error)
            if (h5flag_integer) then
              CALL h5dread_f(h5dataset_id,  h5_file_datatype, h5dataset_Data_integer(1), h5dims, h5error)
              if (h5error<debugflag) then
                write(*,*)  'h5data',h5error,       h5dataset_Data_integer
              elseif (h5error<0)  then
                return(h5error)
              endif
            elseif (h5flag_double) then
              CALL h5dread_f(h5dataset_id,  h5_file_datatype, h5dataset_Data_double(1), h5dims, h5error)
              if (h5error<debugflag) then
                write(*,*)  'h5data',h5error,       h5dataset_Data_double
              elseif (h5error<0)  then
                return(h5error)
              endif
            else
              write(*,*) 'h5 data type not supported'
            endif
            return(h5error)
          endif
! rank=0 scalar
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111

          allocate(h5dims(h5rank))
          allocate(h5maxdims(h5rank))
          if (h5error<debugflag) then
            write(*,*)  'h5error',       h5error,'h5rank'
          endif 
          call h5sget_simple_extent_dims_f(h5dataspace_id, h5dims, h5maxdims,h5error ) 
          if (h5error<debugflag) then
            write(*,*)  'h5error',       h5error, 'h5dims', h5dims,'h5maxdims',h5maxdims
          elseif (h5error<0)  then
            return(h5error)
          endif
          h5datasize=1
          do i =1,h5rank
            h5datasize=h5datasize*h5dims(i)
          enddo
          allocate(h5dataset_data_double(h5datasize))
          allocate(h5dataset_data_integer(h5datasize))
        ! datatype of dataset
        call H5Dget_type_f(h5dataset_id, h5_file_datatype, h5error);
          if (h5error<debugflag) then
            write(*,*)  'h5error',       h5error,'h5_file_datatype',h5_file_datatype
          elseif (h5error<0)  then
            return(h5error)
          endif
          ! datatype of memory data, test datatype
          call H5Tget_native_type_f(h5_file_datatype,H5T_DIR_ASCEND_F, h5_mem_datatype,h5error)
          if (h5error<debugflag) then
            write(*,*)  'h5error',       h5error,'h5_mem_datatype',h5_mem_datatype
          elseif (h5error<0)  then
            return(h5error)
          endif
          call h5tequal_F(h5_mem_datatype,H5T_NATIVE_DOUBLE,h5flag,h5error)
          if (h5error<debugflag) then
            write(*,*)  'h5error',       h5error,'h5_mem_datatype',h5_mem_datatype,h5flag
          elseif (h5error<0)  then
            return(h5error)
          endif
          call h5tequal_F(h5_file_datatype,H5T_NATIVE_DOUBLE,h5flag,h5error)
          if (h5error<debugflag) then
            write(*,*)  'h5error',       h5error,'h5_file_datatype',h5_file_datatype,h5flag
          elseif (h5error<0)  then
            return(h5error)
          endif
!!!!!!!!!!!!!!!!!!!!!!!!
!! read matrix old
!          call h5tequal_F(h5_file_datatype,H5T_NATIVE_integer,h5flag,h5error)
!          if (h5flag) then
!            CALL h5dread_f(h5dataset_id,  h5_file_datatype, h5dataset_Data_integer, h5dims, h5error)
!            if (h5error<debugflag) then
!              write(*,*)  'h5data',h5error,       h5dataset_Data_integer
!            elseif (h5error<0)  then
!              return(h5error)
!            endif
!          endif
!          call h5tequal_F(h5_file_datatype,H5T_NATIVE_double,h5flag,h5error)
!          if (h5flag) then
!            CALL h5dread_f(h5dataset_id,  h5_file_datatype, h5dataset_Data_double, h5dims, h5error)
!            if (h5error<debugflag) then
!              write(*,*)  'h5data',h5error,       h5dataset_Data_double
!            elseif (h5error<0)  then
!              return(h5error)
!            endif
!          endif
!! read matrix old
!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!
! read matrix 
            call h5tequal_F(h5_file_datatype,H5T_NATIVE_integer,h5flag_integer,h5error)
            call h5tequal_F(h5_file_datatype,H5T_NATIVE_double,h5flag_double,h5error)
            if (h5flag_integer) then
              CALL h5dread_f(h5dataset_id,  h5_file_datatype, h5dataset_Data_integer, h5dims, h5error)
              if (h5error<debugflag) then
                write(*,*)  'h5data',h5error,       h5dataset_Data_integer
              elseif (h5error<0)  then
                return(h5error)
              endif
            elseif (h5flag_double) then
              CALL h5dread_f(h5dataset_id,  h5_file_datatype, h5dataset_Data_double, h5dims, h5error)
              if (h5error<debugflag) then
                write(*,*)  'h5data',h5error,       h5dataset_Data_double
              elseif (h5error<0)  then
                return(h5error)
              endif
            else
              write(*,*) 'h5 data type not supported'
            endif
 
! read matrix 
!!!!!!!!!!!!!!!!!!!!!!!!


!        CALL h5dread_f(h5dataset_id,  h5_file_datatype, h5dataset_data, h5dims, h5error)
!        if (h5error<debugflag) then
!          write(*,*)  'h5error',       h5error
!        elseif (h5error<0)  then
!          return(h5error)
!        endif
      CALL h5dclose_f(h5dataset_id, h5error)
    CALL h5fclose_f(h5file_id, h5error)
  CALL h5close_f(h5error)
end subroutine h5gw_read



end subroutine gw_eps_init
    
