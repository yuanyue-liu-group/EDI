    SUBROUTINE readevc( iuni, filename, root_in_group, intra_group_comm,  &
         ik, xk, ispin, npol, wfc, ngw, gamma_only, nbnd, igl, ngwl, &
         b1, b2, b3, mill_k, scalef, ierr )
      !
      !! Processor "root_in_group" reads wfc and related information from file 
      !! "filename.*" (* = dat if fortran binary, * = hdf5 if HDF5),
      !! distributes wfc on "intra_group_comm"
      !! if ierr is present, return 0 if everything is ok, /= 0 if not
      !------------------------------------------------------------------------
      !
      USE mp_wave,     ONLY : splitwf, splitkg
      USE mp,          ONLY : mp_bcast, mp_size, mp_rank, mp_max
      !
      USE  qeh5_base_module

      IMPLICIT NONE
      !
      INTEGER,            INTENT(IN)    :: iuni
      CHARACTER(LEN=*),   INTENT(IN)    :: filename
      INTEGER,            INTENT(IN)    :: root_in_group, intra_group_comm
      INTEGER,            INTENT(IN)    :: ik
      INTEGER,            INTENT(IN)    :: ngwl
      INTEGER,            INTENT(INOUT) :: ngw, nbnd, ispin, npol
      COMPLEX(DP),        INTENT(OUT)   :: wfc(:,:)
      INTEGER,            INTENT(IN)    :: igl(:)
      REAL(DP),           INTENT(OUT)   :: scalef
      REAL(DP),           INTENT(OUT)   :: xk(3)
      REAL(DP),           INTENT(OUT)   :: b1(3), b2(3), b3(3)
      INTEGER,            INTENT(OUT)   :: mill_k(:,:)
      LOGICAL,            INTENT(OUT)   :: gamma_only
      INTEGER, OPTIONAL,  INTENT(OUT)   :: ierr
      !
      INTEGER                           :: j
      INTEGER, ALLOCATABLE              :: itmp(:,:)
      COMPLEX(DP), ALLOCATABLE, TARGET  :: wtmp(:)
      COMPLEX(DP), POINTER              :: wtmp2(:)
      INTEGER                           :: ierr_
      INTEGER                           :: igwx, igwx_, npwx, ik_, nbnd_
      INTEGER                           :: me_in_group, nproc_in_group
      LOGICAL                           :: ionode_in_group
      TYPE (qeh5_file)    ::   h5file
      TYPE (qeh5_dataset) ::   h5dset_wfc, h5dset_mill
      CHARACTER(LEN=8)    ::   char_buf 
      !
      !
      igwx = MAXVAL( igl(1:ngwl) )
      CALL mp_max( igwx, intra_group_comm )
      !
         CALL qeh5_openfile( h5file, TRIM(filename)//'.hdf5', ACTION = 'read', ERROR = ierr_)
write(*,*) '1.1 ierr_ openfile        ', ierr_
      END IF
      CALL mp_bcast( ierr_, root_in_group, intra_group_comm )
      IF ( PRESENT(ierr) ) THEN
         ierr = ierr_
         IF ( ierr /= 0 ) RETURN
      ELSE
         CALL errore( 'read_wfc ', &
              'cannot open restart file ' // TRIM(filename) //' for reading', ierr_ )
      END IF
      !
      IF ( ionode_in_group ) THEN
write(*,*)  ' 1.1 h5, ik_',  ik_
         CALL qeh5_read_attribute (h5file%id, "ik", ik_)
write(*,*)  ' 1.1 h5, ik_',  ik_
         CALL qeh5_read_attribute (h5file%id, "xk",xk, RANK =1, DIMENSIONS = [3])
         CALL qeh5_read_attribute (h5file%id, "ispin", ispin)
         CALL qeh5_read_attribute (h5file%id, "gamma_only", char_buf, MAXLEN = len(char_buf) )
         IF (TRIM(char_buf) =='.TRUE.' .OR. TRIM(char_buf)=='.true.') THEN 
            gamma_only = .TRUE. 
         ELSE 
            gamma_only = .FALSE.
         END IF
         CALL qeh5_read_attribute (h5file%id, "scale_factor",scalef)
write(*,*)  ' 1.1 h5, ngw',    ngw
         CALL qeh5_read_attribute (h5file%id, "ngw", ngw)
write(*,*)  ' 1.1 h5, ngw',    ngw
         CALL qeh5_read_attribute (h5file%id, "nbnd", nbnd_)
write(*,*)  ' 1.1 h5, nbnd_',    nbnd_
         CALL qeh5_read_attribute (h5file%id, "npol",npol)
         CALL qeh5_read_attribute (h5file%id, "igwx",igwx_)
write(*,*)  ' 1.1 h5, npol',    npol
      END IF
      !
      CALL mp_bcast( ik_,    root_in_group, intra_group_comm )
      CALL mp_bcast( xk,     root_in_group, intra_group_comm )
      CALL mp_bcast( ispin,  root_in_group, intra_group_comm )
      CALL mp_bcast( gamma_only, root_in_group, intra_group_comm )
      CALL mp_bcast( scalef, root_in_group, intra_group_comm )
      CALL mp_bcast( ngw,    root_in_group, intra_group_comm )
      CALL mp_bcast( igwx_,  root_in_group, intra_group_comm )
      CALL mp_bcast( npol,   root_in_group, intra_group_comm )
      CALL mp_bcast( nbnd_,   root_in_group, intra_group_comm )
      !
      npwx = SIZE( wfc, 1 ) / npol
write(*,*)  '    1.2  npwx , SIZE( wfc, 1 ) , npol',      npwx , SIZE( wfc, 1 ) , npol
      !
      IF ( ionode_in_group ) THEN 
         ALLOCATE( itmp( 3,MAX( igwx_, igwx ) ) )
       CALL qeh5_open_dataset(h5file, h5dset_mill, ACTION = 'read', NAME = 'MillerIndices')
       IF ( h5dset_mill%filespace%dims(2) .GT. MAX(igwx_, igwx)  ) &
          CALL errore ( 'read_wfc', 'real dimensions of Miller Indices dataset do not  match with igwx attribute', 8) 
       ! no reading of b1, b2, and b3 from file. They should be already set. 
       CALL qeh5_read_dataset ( itmp(:,1), h5dset_mill) 
       CALL qeh5_close ( h5dset_mill) 
         IF ( igwx > igwx_ ) itmp(1:3,igwx_+1:igwx) = 0
      ELSE
         ALLOCATE( itmp( 3, 1 ) )
      END IF
      CALL splitkg( mill_k(:,:), itmp, ngwl, igl, me_in_group, &
           nproc_in_group, root_in_group, intra_group_comm )
      DEALLOCATE (itmp)
      !
      IF ( ionode_in_group ) THEN 
         ALLOCATE( wtmp( npol*MAX( igwx_, igwx ) ) )
         IF ( npol == 2 ) wtmp2 => wtmp(igwx_+1:2*igwx_)
         CALL qeh5_open_dataset( h5file, h5dset_wfc, ACTION = 'read', NAME = 'evc')
         CALL qeh5_set_space ( h5dset_wfc, wtmp(1), RANK = 1, DIMENSIONS = [npol*igwx_], MODE = 'm') 
      ELSE
         ALLOCATE( wtmp(1) )
         IF ( npol == 2 ) wtmp2 => wtmp( 1:1 )
      ENDIF
      nbnd = nbnd_ 
      DO j = 1, nbnd_ 
         !
         IF ( j <= SIZE( wfc, 2 ) ) THEN
            !
            IF ( ionode_in_group ) THEN 
               
               CALL qeh5_set_file_hyperslab (h5dset_wfc, OFFSET = [0,j-1], COUNT = [2*npol*igwx_,1] )
               CALL qeh5_read_dataset (wtmp, h5dset_wfc )
               IF ( igwx > igwx_ ) wtmp((npol*igwx_+1):npol*igwx) = 0.0_DP
               !
            END IF
            !
            IF ( npol == 2 ) THEN
               !
               ! Quick-and-dirty noncolinear case - mergewf should be modified
               ! Collect into wtmp(1:igwx_) first set of plane wave components
               !
               CALL splitwf( wfc(1:npwx,       j), wtmp ,   &
                    ngwl, igl, me_in_group, nproc_in_group, root_in_group, &
                    intra_group_comm )
               !
               ! Collect into wtmp(igwx_+1:2*igwx_) the second set of plane wave
               ! components - instead of wtmp(igwx_+1:2*igwx_), pointer wtmp2
               ! is used, in order to prevent a bogus out-of-bound error
               !
               CALL splitwf( wfc(npwx+1:2*npwx,j), wtmp2,  &
                    ngwl, igl, me_in_group, nproc_in_group, root_in_group, &
                    intra_group_comm )
            ELSE
               CALL splitwf( wfc(:,j), wtmp, ngwl, igl, me_in_group, &
                    nproc_in_group, root_in_group, intra_group_comm )
            END IF
            !
         END IF
         !
      END DO
      !
      IF ( ionode_in_group ) THEN
         CALL qeh5_close(h5dset_wfc) 
         CALL qeh5_close(h5file)
      END IF
      !
      IF ( npol == 2 ) NULLIFY ( wtmp2 )
      DEALLOCATE( wtmp )
      !
      RETURN
      !
    END SUBROUTINE readevc
    !

