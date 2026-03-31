MODULE edi_read_hr
  USE kinds, ONLY : dp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_hr_file

CONTAINS

  SUBROUTINE read_hr_file(seedname, nwann, nrpts, ndegen, irvec, hr)
    USE io_global, ONLY : ionode, ionode_id, stdout
    USE mp, ONLY : mp_bcast
    USE mp_world, ONLY : world_comm
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: seedname
    INTEGER, INTENT(OUT) :: nwann, nrpts
    INTEGER, ALLOCATABLE, INTENT(OUT) :: ndegen(:), irvec(:,:)
    COMPLEX(dp), ALLOCATABLE, INTENT(OUT) :: hr(:,:,:)

    CHARACTER(LEN=256) :: fname
    INTEGER :: iunit, ios, ir, ii, jj, r1, r2, r3, idx_i, idx_j
    REAL(dp) :: re_part, im_part
    CHARACTER(LEN=256) :: line

    iunit = 88
    fname = TRIM(seedname) // '_hr.dat'

    nwann = 0
    nrpts = 0

    IF (ionode) THEN
       OPEN(iunit, FILE=TRIM(fname), STATUS='old', FORM='formatted', IOSTAT=ios)
       IF (ios /= 0) CALL errore('read_hr_file', 'Cannot open ' // TRIM(fname), 1)

       READ(iunit, '(A)') line
       READ(iunit, *) nwann
       READ(iunit, *) nrpts

       ALLOCATE(ndegen(nrpts))
       ir = 0
       DO WHILE (ir < nrpts)
          READ(iunit, *, IOSTAT=ios) ndegen(ir+1:MIN(ir+15, nrpts))
          IF (ios /= 0) EXIT
          ir = MIN(ir + 15, nrpts)
       ENDDO

       ALLOCATE(irvec(3, nrpts))
       ALLOCATE(hr(nwann, nwann, nrpts))
       hr = (0.0_dp, 0.0_dp)

       DO ir = 1, nrpts
          DO jj = 1, nwann
             DO ii = 1, nwann
                READ(iunit, *) r1, r2, r3, idx_i, idx_j, re_part, im_part
                IF (ii == 1 .AND. jj == 1) irvec(:, ir) = (/r1, r2, r3/)
                hr(idx_i, idx_j, ir) = CMPLX(re_part, im_part, dp)
             ENDDO
          ENDDO
       ENDDO

       CLOSE(iunit)
       WRITE(stdout, '(5X,A,A)') 'Read H(R) from ', TRIM(fname)
       WRITE(stdout, '(5X,A,I4,A,I6)') 'nwann = ', nwann, ', nrpts = ', nrpts
    ENDIF

    CALL mp_bcast(nwann, ionode_id, world_comm)
    CALL mp_bcast(nrpts, ionode_id, world_comm)

    IF (.NOT. ionode) THEN
       ALLOCATE(ndegen(nrpts))
       ALLOCATE(irvec(3, nrpts))
       ALLOCATE(hr(nwann, nwann, nrpts))
    ENDIF

    CALL mp_bcast(ndegen, ionode_id, world_comm)
    CALL mp_bcast(irvec, ionode_id, world_comm)
    CALL mp_bcast(hr, ionode_id, world_comm)
  END SUBROUTINE read_hr_file

END MODULE edi_read_hr
