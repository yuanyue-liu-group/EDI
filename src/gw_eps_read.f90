subroutine gw_eps_read(eps_filename_,gw_)
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
 
  h5filename=trim(eps_filename_)      ! Dataset name

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
  ! write(*,*) 'rank,h5dims',p_rank,h5dims(:), allocated(h5dims)
  call h5gw_read(h5filename,h5datasetname,h5dataset_data_double,h5dataset_Data_integer,h5dims,h5rank,h5error)
  if (h5error<0)  write(*,*)  'h5error',h5error
  if (h5rank/=1) then
      write(*,*)  'h5rank error(should be 1)',h5rank 
  else
      h5dims1=h5dims
  
      if ( .not. allocated(gw_%ng_data)) then
          allocate(gw_%ng_data(h5dims1(1)))
      else
          deallocate(gw_%ng_data)
      allocate(gw_%ng_data(h5dims1(1)))
  endif
  
     
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
      
      if (  allocated(gw_%nmtx_max_data)) then
          deallocate(gw_%nmtx_max_data)
      endif
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
      if (  allocated(gw_%nmtx_data)) then
          deallocate(gw_%nmtx_data)
      endif
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
      if (  allocated(gw_%gind_eps2rho_data)) then
          deallocate(gw_%gind_eps2rho_data)
      endif
      allocate(gw_%gind_eps2rho_data(h5dims2(1),h5dims2(2)))
      gw_%gind_eps2rho_data=reshape(h5dataset_data_integer,h5dims2)
      write(*,*)  'shape h5dataset',shape(gw_%gind_eps2rho_data)
      write(*,*)  'gw_gind_eps2rho_data()',gw_%gind_eps2rho_data(1:40,1)
      write(*,*)  'gw_gind_eps2rho_data()',gw_%gind_eps2rho_data(1:40,2)
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
      if (  allocated(gw_%gind_rho2eps_data)) then
          deallocate(gw_%gind_rho2eps_data)
      endif
      allocate(gw_%gind_rho2eps_data(h5dims2(1),h5dims2(2)))
      gw_%gind_rho2eps_data=reshape(h5dataset_data_integer,h5dims2)
      write(*,*)  'shape h5dataset',shape(h5dataset_data_integer)
      write(*,*)  'gw_gind_rho2eps_data()',gw_%gind_rho2eps_data(1:40,1)
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
      if (  allocated(gw_%g_components_data)) then
          deallocate(gw_%g_components_data)
      endif
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
      if (  allocated(gw_%bvec_data)) then
          deallocate(gw_%bvec_data)
      endif
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
      if (  allocated(gw_%blat_data)) then
          deallocate(gw_%blat_data)
      endif
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
      if (  allocated(gw_%qpts_data)) then
          deallocate(gw_%qpts_data)
      endif
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
      !write(*,*) 'sizeof(int(i4b)):',sizeof(gw_%nq)
      !write(*,*) 'sizeof(int(i8b)):',sizeof(gw_nqi8)
      write(*,*) 'sizeof(int):',sizeof(h5rank)
      h5dims1=h5dims
      if (  allocated(gw_%nq_data)) then
          deallocate(gw_%nq_data)
      endif
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
      if (  allocated(gw_%epsmat_diag_data)) then
         deallocate(gw_%epsmat_diag_data)
      endif
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
      if (  allocated(gw_%epsmat_full_data)) then
          deallocate(gw_%epsmat_full_data)
      endif
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
      deallocate(h5dims)
      deallocate(h5dataset_Data_integer)
      deallocate(h5dataset_Data_double)
  endif
 
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
    INTEGER(HID_T) :: debugflag=01
    ! if debugflag<=10, not print epsilon data, else, print
    INTEGER(HID_T)                               :: loc_id, attr_id, data_type, mem_type
    integer:: p_rank,p_size,ik
    !  CALL h5open_f(h5error)
    
    
    !call  mpi_comm_rank(mpi_comm_world,p_rank,ik)
    !call  mpi_comm_size(mpi_comm_world,p_size,ik)
    ! write(*,*) 'rank,h5dims',p_rank,h5dims(:), allocated(h5dims)


    if (h5error<debugflag) then
        write(*,*)  'h5error',h5error
    elseif (h5error<0) then 
        return(h5error)
    endif
  
    !h5 file
    CALL h5fopen_f (h5filename, H5F_ACC_RDWR_F, h5file_id, h5error)
    if (h5error<debugflag) then
        write(*,*)  'h5error',       h5error,'h5filename',trim(h5filename),'h5file_id', h5file_id
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
            !if (h5rank>0) write(*,*)   h5dims,h5maxdims
          elseif (h5error<0)  then
            return(h5error)
          endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1111
! rank=0 scalar
          if(h5rank==0) then
            !write(*,*)  'h5error',       h5error,'h5rank',h5rank
            h5rank=1
            !write(*,*)  'h5error',       h5error,'h5rank',h5rank
            allocate(h5maxdims(h5rank))
            !write(*,*)  'h5error',       h5error,'h5dimssize',size(h5maxdims)
            !write(*,*)  'h5error',       h5error,'h5dimssize',size(h5dims)
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
!!!!!!!!!!!!!!!!
!!debug comment out ok
!            ! datatype of memory data, test datatype
!            call H5Tget_native_type_f(h5_file_datatype,H5T_DIR_ASCEND_F, h5_mem_datatype,h5error)
!              if (h5error<debugflag) then
!                write(*,*)  'h5error',       h5error,'h5_mem_datatype',h5_mem_datatype
!              elseif (h5error<0)  then
!                return(h5error)
!              endif
!!              call h5tequal_F(h5_mem_datatype,H5T_NATIVE_integer,h5flag,h5error)
!              if (h5error<debugflag) then
!                write(*,*)  'h5error',       h5error,'h5_mem_datatype',h5_mem_datatype,'H5T_NATIVE_integer'
!              elseif (h5error<0)  then
!                return(h5error)
!              endif
!!              call h5tequal_F(h5_file_datatype,H5T_NATIVE_integer,h5flag,h5error)
!              if (h5error<debugflag) then
!                write(*,*)  'h5error',       h5error,'h5_file_datatype',h5_file_datatype,h5flag
!              elseif (h5error<0)  then
!                return(h5error)
!              endif
!
!! qeh5_module bug
!        CALL H5Tcopy_f( H5T_NATIVE_INTEGER, mem_type, ierr )
!write(*,*) 'ierr        H5T_NATIVE_INTEGER',ierr,H5T_NATIVE_INTEGER, mem_type,sizeof(H5T_NATIVE_INTEGER), sizeof(mem_type)
!! qeh5_module bug
!
!!debug comment out ok
!!!!!!!!!!!!!!!!

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

! qeh5_module bug
!        CALL H5Tcopy_f( H5T_NATIVE_INTEGER, mem_type, ierr )
!write(*,*) 'ierr        H5T_NATIVE_INTEGER',ierr,H5T_NATIVE_INTEGER, mem_type,sizeof(H5T_NATIVE_INTEGER), sizeof(mem_type)
! qeh5_module bug


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
            write(*,*)  'h5error',       h5error,'h5_mem_datatype',h5_mem_datatype,h5flag
          if (h5error<debugflag) then
            write(*,*)  'h5error',       h5error,'h5_mem_datatype',h5_mem_datatype,h5flag
          elseif (h5error<0)  then
            return(h5error)
          endif
          call h5tequal_F(h5_file_datatype,H5T_NATIVE_DOUBLE,h5flag,h5error)
            write(*,*)  'h5error',       h5error,'h5_file_datatype',h5_file_datatype,h5flag
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
              if (h5error<debugflag-10) then
                write(*,*)  'h5data',h5error,       h5dataset_Data_integer
              elseif (h5error<0)  then
                return(h5error)
              endif
            elseif (h5flag_double) then
              CALL h5dread_f(h5dataset_id,  h5_file_datatype, h5dataset_Data_double, h5dims, h5error)
              if (h5error<debugflag-10) then
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



  ! leads to bug in later read_wfc in qe
  !  CALL h5close_f(h5error)
  ! if uncomment leads to bug in later read_wfc in qe
  end subroutine h5gw_read
end subroutine gw_eps_read


