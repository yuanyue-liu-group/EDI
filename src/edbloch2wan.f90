MODULE edbloch2wan
  USE kinds, ONLY : dp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: edbloch2wane, edbloch2wanr

CONTAINS

  SUBROUTINE edbloch2wane(nbnd, nbndsub, nks, nkstot, xk, cu, cuq, &
                           edmatkq, nrr, irvec, wslen, edmatw)
    USE cell_base, ONLY : at, bg, alat
    USE ep_constants, ONLY : twopi, ci, czero, cone, bohr2ang, zero
    USE io_global, ONLY : ionode_id, ionode, stdout
    USE mp_global, ONLY : inter_pool_comm
    USE mp, ONLY : mp_sum
    USE mp_world, ONLY : mpime
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nbnd, nbndsub, nks, nkstot, nrr
    REAL(dp), INTENT(IN) :: xk(3, nks)
    COMPLEX(dp), INTENT(IN) :: cu(nbnd, nbndsub, nks)
    COMPLEX(dp), INTENT(IN) :: cuq(nbnd, nbndsub, nks)
    COMPLEX(dp), INTENT(IN) :: edmatkq(nbnd, nbnd, nks)
    INTEGER, INTENT(IN) :: irvec(3, nrr)
    REAL(dp), INTENT(IN) :: wslen(nrr)
    COMPLEX(dp), INTENT(OUT) :: edmatw(nbndsub, nbndsub, nrr)

    INTEGER :: ik, ir, ibnd, jbnd, mbnd
    REAL(dp) :: rdotk, tmp
    COMPLEX(dp) :: cfac
    COMPLEX(dp), ALLOCATABLE :: edms(:,:,:), eptmp(:,:)

    ALLOCATE(edms(nbndsub, nbndsub, nks))
    ALLOCATE(eptmp(nbnd, nbndsub))

    edms = czero
    DO ik = 1, nks
       eptmp = czero
       CALL ZGEMM('N', 'N', nbnd, nbndsub, nbnd, cone, &
                   edmatkq(:,:,ik), nbnd, cuq(:,:,ik), nbnd, czero, eptmp, nbnd)
       CALL ZGEMM('C', 'N', nbndsub, nbndsub, nbnd, cone, &
                   cu(:,:,ik), nbnd, eptmp, nbnd, czero, edms(:,:,ik), nbndsub)
    ENDDO

    CALL cryst_to_cart(nks, xk, at, -1)

    edmatw = czero
    DO ir = 1, nrr
       DO ik = 1, nks
          rdotk = twopi * DOT_PRODUCT(xk(:, ik), DBLE(irvec(:, ir)))
          cfac = EXP(-ci * rdotk) / DBLE(nkstot)
          edmatw(:, :, ir) = edmatw(:, :, ir) + cfac * edms(:, :, ik)
       ENDDO
    ENDDO
    CALL mp_sum(edmatw, inter_pool_comm)

    CALL cryst_to_cart(nks, xk, bg, 1)

    IF (mpime == ionode_id) THEN
       OPEN(UNIT=99, FILE='decay.M')
       WRITE(99, '(a)') '# Spatial decay of e-d matrix element in Wannier basis'
       DO ir = 1, nrr
          tmp = zero
          DO jbnd = 1, nbndsub
             DO ibnd = 1, nbndsub
                tmp = MAX(tmp, ABS(edmatw(ibnd, jbnd, ir)))
             ENDDO
          ENDDO
          WRITE(99, '(5x, f15.8, 2x, E22.14)') wslen(ir) * alat * bohr2ang, tmp
       ENDDO
       CLOSE(99)
       WRITE(stdout, '(5X,A)') 'Wrote decay.M (spatial decay of M in Wannier basis)'
    ENDIF

    DEALLOCATE(edms, eptmp)
  END SUBROUTINE edbloch2wane

  SUBROUTINE edbloch2wanr(nbndsub, nks, nkstot, xk, edmatkq_wan, &
                            nrr, irvec, edmatw)
    USE cell_base, ONLY : at, bg
    USE ep_constants, ONLY : twopi, ci, czero
    USE mp_global, ONLY : inter_pool_comm
    USE mp, ONLY : mp_sum
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nbndsub, nks, nkstot, nrr
    REAL(dp), INTENT(IN) :: xk(3, nks)
    COMPLEX(dp), INTENT(IN) :: edmatkq_wan(nbndsub, nbndsub, nks)
    INTEGER, INTENT(IN) :: irvec(3, nrr)
    COMPLEX(dp), INTENT(OUT) :: edmatw(nbndsub, nbndsub, nrr)

    INTEGER :: ik, ir
    REAL(dp) :: rdotk
    COMPLEX(dp) :: cfac

    CALL cryst_to_cart(nks, xk, at, -1)

    edmatw = czero
    DO ir = 1, nrr
       DO ik = 1, nks
          rdotk = twopi * DOT_PRODUCT(xk(:, ik), DBLE(irvec(:, ir)))
          cfac = EXP(-ci * rdotk) / DBLE(nkstot)
          edmatw(:, :, ir) = edmatw(:, :, ir) + cfac * edmatkq_wan(:, :, ik)
       ENDDO
    ENDDO
    CALL mp_sum(edmatw, inter_pool_comm)

    CALL cryst_to_cart(nks, xk, bg, 1)
  END SUBROUTINE edbloch2wanr

END MODULE edbloch2wan
