MODULE io
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: readwfc

CONTAINS

  SUBROUTINE readwfc(ipool, recn, evc0)
    USE kinds, ONLY : DP
    USE io_files, ONLY : prefix, tmp_dir
    USE wvfct, ONLY : npwx
    USE wvfct, ONLY : nbnd
    USE noncollin_module, ONLY : npol
    USE mp_global, ONLY : nproc_pool, me_pool, npool
    USE low_lvl, ONLY : set_ndnmbr
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: recn, ipool
    COMPLEX(KIND=DP), INTENT(OUT) :: evc0(npwx * npol, nbnd)

    CHARACTER(LEN=256) :: tempfile
    CHARACTER(LEN=4) :: nd_nmbr0
    INTEGER :: unf_recl, ios, lrwfc
    REAL(KIND=DP) :: dummy

    lrwfc = 2 * nbnd * npwx * npol

    CALL set_ndnmbr(ipool, me_pool, nproc_pool, npool, nd_nmbr0)
#if defined(__MPI)
    tempfile = TRIM(tmp_dir) // TRIM(prefix) // '.wfc' // nd_nmbr0
#else
    tempfile = TRIM(tmp_dir) // TRIM(prefix) // '.wfc'
#endif
    INQUIRE(IOLENGTH=unf_recl) dummy
    unf_recl = unf_recl * lrwfc

    OPEN(20, FILE=tempfile, FORM='unformatted', ACCESS='direct', &
         IOSTAT=ios, RECL=unf_recl)
    IF (ios /= 0) CALL errore('readwfc', 'error opening wfc file', 20)
    READ(20, REC=recn) evc0
    CLOSE(20, STATUS='keep')
  END SUBROUTINE readwfc

END MODULE io
