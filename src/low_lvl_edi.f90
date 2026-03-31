MODULE low_lvl
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: set_ndnmbr, s, px, py, p_z, dz2, dxz, dyz, dx2my2, dxy
  PUBLIC :: fz3, fxz2, fyz2, fzx2my2, fxyz, fxx2m3y2, fy3x2my2
  PUBLIC :: utility_zdotu
CONTAINS
    PURE FUNCTION s()
    !-----------------------------------------------------------------------
    !!
    !! s-orbital
    !!
    USE kinds,         ONLY : DP
    USE ep_constants, ONLY : fpi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP) :: s

    s = 1.d0 / DSQRT(fpi)
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION s
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    PURE FUNCTION px(cost, phi)
    !-----------------------------------------------------------------------
    !!
    !! p-orbital
    !!
    USE kinds,         ONLY : DP
    USE ep_constants, ONLY : fpi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: cost
    !! cos(t)
    REAL(KIND = DP), INTENT(in) :: phi
    !! Phi
    REAL(KIND = DP) :: px
    !! Output
    !
    ! Local variable
    REAL(KIND = DP) :: sint
    !! sin(t)
    !
    sint = DSQRT(ABS(1.d0 - cost * cost))
    px =  DSQRT(3.d0 / fpi) * sint * cos(phi)
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION px
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    PURE FUNCTION py(cost, phi)
    !-----------------------------------------------------------------------
    !!
    !! p-orbital
    !!
    USE kinds,         ONLY : DP
    USE ep_constants, ONLY : fpi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: cost
    !! Cos
    REAL(KIND = DP), INTENT(in) :: phi
    !! Phi
    REAL(KIND = DP) :: py
    !! Output
    !
    ! Local variables
    REAL(KIND = DP) :: sint
    !! Sin(t)
    !
    sint = DSQRT(ABS(1.d0 - cost * cost))
    py =  DSQRT(3.d0 / fpi) * sint * sin(phi)
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION py
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    PURE FUNCTION p_z(cost)
    !-----------------------------------------------------------------------
    !!
    !! p-orbital
    !!
    USE kinds,         ONLY : DP
    USE ep_constants, ONLY : fpi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: cost
    !! cos(t)
    REAL(KIND = DP) :: p_z
    !! Output
    !
    p_z =  DSQRT(3.d0 / fpi) * cost
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION p_z
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    PURE FUNCTION dz2(cost)
    !-----------------------------------------------------------------------
    !!
    !! d-orbital
    !!
    USE kinds, ONLY : DP
    USE ep_constants, ONLY : fpi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: cost
    !! cos(t)
    REAL(KIND = DP) :: dz2
    !! Output
    !
    dz2 =  DSQRT(1.25d0 / fpi) * (3.d0 * cost * cost - 1.d0)
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION dz2
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    PURE FUNCTION dxz(cost, phi)
    !-----------------------------------------------------------------------
    !!
    !! d-orbital
    !!
    USE kinds,         ONLY : DP
    USE ep_constants, ONLY : fpi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: cost
    !! cos(t)
    REAL(KIND = DP), INTENT(in) :: phi
    !! phi
    REAL(KIND = DP) :: dxz
    !! Output
    REAL(KIND = DP) :: sint
    !! sin(t)
    !
    sint = DSQRT(ABS(1.d0 - cost * cost))
    dxz =  DSQRT(15.d0 / fpi) * sint * cost * COS(phi)
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION dxz
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    PURE FUNCTION dyz(cost, phi)
    !-----------------------------------------------------------------------
    !!
    !! d-orbital
    !!
    USE kinds,         ONLY : DP
    USE ep_constants, ONLY : fpi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: cost
    !! cos(t)
    REAL(KIND = DP), INTENT(in) :: phi
    !! phi
    REAL(KIND = DP) :: dyz
    !! output
    !
    ! Local variables
    REAL(KIND = DP) :: sint
    !! sin(t)
    !
    sint = DSQRT(ABS(1.d0 - cost * cost))
    dyz =  DSQRT(15.d0 / fpi) * sint * cost * SIN(phi)
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION dyz
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    PURE FUNCTION dx2my2(cost, phi)
    !-----------------------------------------------------------------------
    !!
    !! d-orbital
    !!
    USE kinds,         ONLY : DP
    USE ep_constants, ONLY : fpi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: cost
    !! cos(t)
    REAL(KIND = DP), INTENT(in) :: phi
    !! phi
    REAL(KIND = DP) :: dx2my2
    !! Output
    !
    ! Local variables
    REAL(KIND = DP) :: sint
    !! sin(t)
    !
    sint = DSQRT(ABS(1.d0 - cost * cost))
    dx2my2 =  DSQRT(3.75d0 / fpi) * sint * sint * COS(2.d0 * phi)
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION dx2my2
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    PURE FUNCTION dxy(cost, phi)
    !-----------------------------------------------------------------------
    !!
    !! d-orbital
    !!
    USE kinds,         ONLY : DP
    USE ep_constants, ONLY : fpi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: cost
    !! cos(t)
    REAL(KIND = DP), INTENT(in) :: phi
    !! phi
    REAL(KIND = DP) :: dxy
    !! Output
    !
    ! Local variables
    REAL(KIND = DP) :: sint
    !! sin(t)
    !
    sint = DSQRT(ABS(1.d0 - cost * cost))
    dxy =  DSQRT(3.75d0 / fpi) * sint * sint * SIN(2.d0 * phi)
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION dxy
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    PURE FUNCTION fz3(cost)
    !-----------------------------------------------------------------------
    !!
    !! f-orbital
    !!
    USE kinds,         ONLY : DP
    USE ep_constants, ONLY : pi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: cost
    !! cos(t)
    REAL(KIND = DP) :: fz3
    !! Output
    !
    fz3 = 0.25d0 * DSQRT(7.d0 / pi) * (5.d0 * cost * cost - 3.d0) * cost
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION fz3
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    PURE FUNCTION fxz2(cost, phi)
    !-----------------------------------------------------------------------
    !!
    !! f-orbital
    !!
    USE kinds,         ONLY : DP
    USE ep_constants, ONLY : pi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: cost
    !! cos(t)
    REAL(KIND = DP), INTENT(in) :: phi
    !! phi
    REAL(KIND = DP) :: fxz2
    !! Output
    !
    ! Local variables
    REAL(KIND = DP) :: sint
    !! sin(t)
    !
    sint = DSQRT(ABS(1.d0 - cost * cost))
    fxz2 = 0.25d0 * DSQRT(10.5d0 / pi) * (5.d0 * cost * cost - 1.d0) * sint * COS(phi)
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION fxz2
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    PURE FUNCTION fyz2(cost, phi)
    !-----------------------------------------------------------------------
    !!
    !! f-orbital
    !!
    USE kinds,         ONLY : DP
    USE ep_constants, ONLY : pi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: cost
    !! cos(t)
    REAL(KIND = DP), INTENT(in) :: phi
    !! phi
    REAL(KIND = DP) :: fyz2
    !! Output
    !
    ! Local variables
    REAL(KIND = DP) :: sint
    !! sin(t)
    !
    sint = DSQRT(ABS(1.d0 - cost * cost))
    fyz2 = 0.25d0 * DSQRT(10.5d0 / pi) * (5.d0 * cost * cost - 1.d0) * sint * SIN(phi)
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION fyz2
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    PURE FUNCTION fzx2my2(cost, phi)
    !-----------------------------------------------------------------------
    !!
    !! f-orbital
    !!
    USE kinds,         ONLY : DP
    USE ep_constants, ONLY : pi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: cost
    !! cos(t)
    REAL(KIND = DP), INTENT(in) :: phi
    !! phi
    REAL(KIND = DP) :: fzx2my2
    !! Output
    !
    ! Local variables
    REAL(KIND = DP) :: sint
    !! sin(t)
    !
    sint = DSQRT(ABS(1.d0 - cost * cost))
    fzx2my2 = 0.25d0 * DSQRT(105d0/pi) * sint * sint * cost * COS(2.d0 * phi)
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION fzx2my2
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    PURE FUNCTION fxyz(cost, phi)
    !-----------------------------------------------------------------------
    !!
    !! f-orbital
    !!
    USE kinds,         ONLY : DP
    USE ep_constants, ONLY : pi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: cost
    !! cos(t)
    REAL(KIND = DP), INTENT(in) :: phi
    !! phi
    REAL(KIND = DP) :: fxyz
    !! Output
    !
    ! Local variables
    REAL(KIND = DP) :: sint
    !! sin(t)
    !
    sint = DSQRT(ABS(1.d0 - cost * cost))
    fxyz = 0.25d0 * DSQRT(105d0 / pi) * sint * sint * cost * SIN(2.d0 * phi)
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION fxyz
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    PURE FUNCTION fxx2m3y2(cost, phi)
    !-----------------------------------------------------------------------
    !!
    !! f-orbital
    !!
    USE kinds,         ONLY : DP
    USE ep_constants, ONLY : pi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: cost
    !! cos(t)
    REAL(KIND = DP), INTENT(in) :: phi
    !! phi
    REAL(KIND = DP) :: fxx2m3y2
    !! Output
    !
    ! Local variables
    REAL(KIND = DP) :: sint
    !! sin(t)
    !
    sint = DSQRT(ABS(1.d0 - cost * cost))
    fxx2m3y2 = 0.25d0 * DSQRT(17.5d0 / pi) * sint * sint * sint * COS(3.d0 * phi)
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION fxx2m3y2
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    FUNCTION fy3x2my2(cost, phi)
    !-----------------------------------------------------------------------
    !!
    !! f-orbital
    !!
    USE kinds,         ONLY : DP
    USE ep_constants, ONLY : pi
    !
    IMPLICIT NONE
    !
    REAL(KIND = DP), INTENT(in) :: cost
    !! cos(t)
    REAL(KIND = DP), INTENT(in) :: phi
    !! phi
    REAL(KIND = DP) :: fy3x2my2
    !! Output
    !
    ! Local variables
    REAL(KIND = DP) :: sint
    !! sin(t)
    !
    sint = DSQRT(ABS(1.d0 - cost * cost))
    fy3x2my2 = 0.25d0 * DSQRT(17.5d0 / pi) * sint * sint * sint * SIN(3.d0 * phi)
    !
    RETURN
    !-----------------------------------------------------------------------
    END FUNCTION fy3x2my2
    SUBROUTINE set_ndnmbr(pool, proc, procp, npool, ndlab)
    !----------------------------------------------------------------
    !!
    !!  create ndlab label from pool and proc numbers
    !!
    !!  The rule for deciding the node number is based on
    !!  the restriction set in startup.f90 that every pool
    !!  has the same number of procs.
    !!
    IMPLICIT NONE
    !
    CHARACTER(LEN = 4), INTENT(out) :: ndlab
    !! Label
    INTEGER, INTENT(in) :: pool
    !! Pool = 1,..., npool
    INTEGER, INTENT(in) :: proc
    !! Processor = 0,..., nproc_pool-1
    INTEGER, INTENT(in) :: procp
    !!
    INTEGER, INTENT(in) :: npool
    !!
    ! Local variables
    INTEGER :: node
    !! Number of nodes
    INTEGER :: nprocs
    !!
    !
    nprocs = npool * procp
    !
    node = (pool - 1) * procp + proc + 1
    !
    ndlab = '    '
    IF (nprocs < 10) THEN
      WRITE(ndlab(1:1), '(i1)') node
    ELSEIF (nprocs < 100 ) then
      IF (node < 10) THEN
        WRITE(ndlab(1:1), '(i1)') node
      ELSE
        WRITE(ndlab(1:2), '(i2)') node
      ENDIF
    ELSEIF (nprocs < 100) THEN
      IF (node < 10) THEN
        WRITE(ndlab(1:1), '(i1)') node
      ELSEIF (node < 100) THEN
        WRITE(ndlab(1:2), '(i2)') node
      ELSE
        WRITE(ndlab(1:3), '(i3)') node
      ENDIF
    ELSE
      IF (node < 10) THEN
        WRITE(ndlab(1:1), '(i1)') node
      ELSEIF (node < 100) THEN
        WRITE(ndlab(1:2), '(i2)') node
      ELSEIF (node < 1000) THEN
        WRITE(ndlab(1:3), '(i3)') node
      ELSE
        WRITE(ndlab(1:4), '(i4)') node
      ENDIF
    ENDIF
    !
    !----------------------------------------------------------------
    END SUBROUTINE set_ndnmbr

    FUNCTION utility_zdotu(n, zx, incx, zy, incy)
      USE kinds, ONLY : DP
      IMPLICIT NONE
      INTEGER, INTENT(in) :: n, incx, incy
      COMPLEX(KIND=DP), INTENT(in) :: zx(*), zy(*)
      COMPLEX(KIND=DP) :: utility_zdotu
      COMPLEX(KIND=DP), EXTERNAL :: ZDOTU
      utility_zdotu = ZDOTU(n, zx, incx, zy, incy)
    END FUNCTION utility_zdotu

END MODULE low_lvl
