MODULE wigner_seitz_edi
  USE kinds, ONLY : dp
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC :: wigner_seitz_k

  INTEGER, PARAMETER :: NSEARCH = 3
  REAL(dp), PARAMETER :: WS_DIST_TOL = 1.0d-6

CONTAINS

  SUBROUTINE wigner_seitz_k(nk1, nk2, nk3, at, irvec, ndegen, wslen, nrr)
    USE constants, ONLY : eps8
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nk1, nk2, nk3
    REAL(dp), INTENT(IN) :: at(3, 3)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: irvec(:, :)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: ndegen(:)
    REAL(dp), ALLOCATABLE, INTENT(OUT) :: wslen(:)
    INTEGER, INTENT(OUT) :: nrr

    INTEGER :: n1, n2, n3, i1, i2, i3, ndiff(3)
    INTEGER :: nrr_max, ir, nfound
    REAL(dp) :: dist(((2*NSEARCH+1)**3)), dmin, rvec(3), tot(3)
    INTEGER :: irvec_tmp(3, 20*nk1*nk2*nk3)
    INTEGER :: ndegen_tmp(20*nk1*nk2*nk3)
    REAL(dp) :: wslen_tmp(20*nk1*nk2*nk3)
    INTEGER :: icnt, i

    nrr_max = 20 * nk1 * nk2 * nk3
    nrr = 0

    DO n1 = -nk1, nk1
      DO n2 = -nk2, nk2
        DO n3 = -nk3, nk3
          icnt = 0
          DO i1 = -NSEARCH, NSEARCH
            DO i2 = -NSEARCH, NSEARCH
              DO i3 = -NSEARCH, NSEARCH
                icnt = icnt + 1
                ndiff(1) = n1 - i1 * nk1
                ndiff(2) = n2 - i2 * nk2
                ndiff(3) = n3 - i3 * nk3
                tot(:) = DBLE(ndiff(1)) * at(:, 1) + &
                         DBLE(ndiff(2)) * at(:, 2) + &
                         DBLE(ndiff(3)) * at(:, 3)
                dist(icnt) = SQRT(DOT_PRODUCT(tot, tot))
              ENDDO
            ENDDO
          ENDDO

          dmin = MINVAL(dist)
          rvec(:) = DBLE(n1) * at(:, 1) + DBLE(n2) * at(:, 2) + DBLE(n3) * at(:, 3)

          nfound = 0
          DO i = 1, icnt
            IF (ABS(dist(i) - dmin) < WS_DIST_TOL) nfound = nfound + 1
          ENDDO

          IF (nfound > 0) THEN
            nrr = nrr + 1
            IF (nrr > nrr_max) CALL errore('wigner_seitz_k', 'nrr exceeds nrr_max', 1)
            irvec_tmp(1, nrr) = n1
            irvec_tmp(2, nrr) = n2
            irvec_tmp(3, nrr) = n3
            ndegen_tmp(nrr) = nfound
            wslen_tmp(nrr) = SQRT(DOT_PRODUCT(rvec, rvec))
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    ALLOCATE(irvec(3, nrr))
    ALLOCATE(ndegen(nrr))
    ALLOCATE(wslen(nrr))
    irvec(:, 1:nrr) = irvec_tmp(:, 1:nrr)
    ndegen(1:nrr) = ndegen_tmp(1:nrr)
    wslen(1:nrr) = wslen_tmp(1:nrr)
  END SUBROUTINE wigner_seitz_k

END MODULE wigner_seitz_edi
