

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
  INTEGER(HSIZE_T), allocatable :: h5maxdims(:)
  INTEGER     ::  h5nmembers,i,h5datasize
  INTEGER(HSIZE_T), allocatable , intent(inout) :: h5dims(:)
  INTEGER  , intent(inout)    ::   h5rank
  INTEGER  , intent(inout)    ::   h5error ! Error flag
  INTEGER(HID_T) :: file_s1_t,h5_file_datatype 
  INTEGER(HID_T) :: mem_s1_t  ,h5_mem_datatype  
  INTEGER(HID_T) :: debugflag=02
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
            h5rank=1
            allocate(h5dims(h5rank))
            allocate(h5maxdims(h5rank))
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



