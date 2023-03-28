
CHARACTER(LEN=256) :: eps_filename_
type(gw_eps_data),intent (inout) ,target:: gw_
  CHARACTER(LEN=256) :: h5filename      ! Dataset name
  CHARACTER(LEN=256) :: h5datasetname = "matrix-diagonal"     ! Dataset name
 INTEGER     ::   h5rank,h5error ! Error flag
  real(dp), allocatable :: h5dataset_data_double(:), data_out(:)
  integer, allocatable :: h5dataset_data_integer(:)
  INTEGER(HSIZE_T), allocatable :: h5dims(:),h5maxdims(:)

  integer :: h5dims1(1),h5dims2(2),h5dims3(3),h5dims4(4),h5dims5(5),h5dims6(6)
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


 


  INTEGER(HID_T) :: h5file_id       ! File identifier
  INTEGER(HID_T) :: h5group_id       ! Dataset identifier
  INTEGER(HID_T) :: h5dataset_id       ! Dataset identifier
  INTEGER(HID_T) :: h5datatype_id       ! Dataset identifier
  INTEGER(HID_T) :: h5dataspace_id


 
