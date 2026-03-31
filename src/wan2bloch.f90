MODULE wan2bloch_edi
  USE kinds, ONLY : dp
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC :: hambloch2wan_edi
  PUBLIC :: hamwan2bloch_edi
  PUBLIC :: hamwan2bloch_with_evec
  PUBLIC :: edmatwan2bloch
  PUBLIC :: edmatwan2bloch_2d
  PUBLIC :: get_cfac

  REAL(dp), PARAMETER :: twopi = 6.283185307179586_dp
  COMPLEX(dp), PARAMETER :: ci = (0.0_dp, 1.0_dp)
  COMPLEX(dp), PARAMETER :: czero = (0.0_dp, 0.0_dp)
  COMPLEX(dp), PARAMETER :: cone = (1.0_dp, 0.0_dp)

CONTAINS

  SUBROUTINE hambloch2wan_edi(nbndsub, nks, nkstot, et_sub, xk_cryst, cu, &
                               nrr, irvec, chw)
    USE mp_global, ONLY : inter_pool_comm
    USE mp, ONLY : mp_sum
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nbndsub, nks, nkstot, nrr
    REAL(dp), INTENT(IN) :: et_sub(nbndsub, nks)
    REAL(dp), INTENT(IN) :: xk_cryst(3, nks)
    COMPLEX(dp), INTENT(IN) :: cu(nbndsub, nbndsub, nks)
    INTEGER, INTENT(IN) :: irvec(3, nrr)
    COMPLEX(dp), INTENT(OUT) :: chw(nbndsub, nbndsub, nrr)

    INTEGER :: ik, ir, ibnd, jbnd, mbnd
    REAL(dp) :: rdotk
    COMPLEX(dp) :: cfac, ctmp
    COMPLEX(dp), ALLOCATABLE :: chs(:, :, :)

    ALLOCATE(chs(nbndsub, nbndsub, nks))

    chs = czero
    DO ik = 1, nks
       DO jbnd = 1, nbndsub
          DO ibnd = 1, jbnd
             ctmp = czero
             DO mbnd = 1, nbndsub
                ctmp = ctmp + CONJG(cu(mbnd, ibnd, ik)) * et_sub(mbnd, ik) * cu(mbnd, jbnd, ik)
             ENDDO
             chs(ibnd, jbnd, ik) = ctmp
             chs(jbnd, ibnd, ik) = CONJG(ctmp)
          ENDDO
       ENDDO
    ENDDO

    chw = czero
    DO ir = 1, nrr
       DO ik = 1, nks
          rdotk = twopi * DOT_PRODUCT(xk_cryst(:, ik), DBLE(irvec(:, ir)))
          cfac = EXP(-ci * rdotk) / DBLE(nkstot)
          chw(:, :, ir) = chw(:, :, ir) + cfac * chs(:, :, ik)
       ENDDO
    ENDDO
    CALL mp_sum(chw, inter_pool_comm)

    DEALLOCATE(chs)
  END SUBROUTINE hambloch2wan_edi

  SUBROUTINE hamwan2bloch_edi(nbnd, nrr, ndegen, eig, chw, cfac)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nbnd, nrr
    INTEGER, INTENT(IN) :: ndegen(nrr)
    REAL(dp), INTENT(OUT) :: eig(nbnd)
    COMPLEX(dp), INTENT(IN) :: cfac(nrr)
    COMPLEX(dp), INTENT(IN) :: chw(nbnd, nbnd, nrr)

    INTEGER :: ibnd, jbnd, ir, info, lwork
    COMPLEX(dp) :: chf(nbnd, nbnd)
    COMPLEX(dp) :: cfac_degen(nrr)
    REAL(dp) :: w(nbnd)
    REAL(dp), ALLOCATABLE :: rwork(:)
    COMPLEX(dp), ALLOCATABLE :: work(:)

    DO ir = 1, nrr
       IF (ndegen(ir) > 0) THEN
          cfac_degen(ir) = cfac(ir) / DBLE(ndegen(ir))
       ELSE
          cfac_degen(ir) = czero
       ENDIF
    ENDDO

    chf = czero
    CALL ZGEMV('n', nbnd*nbnd, nrr, cone, chw, nbnd*nbnd, cfac_degen, 1, cone, chf, 1)

    DO jbnd = 1, nbnd
       DO ibnd = 1, jbnd
          chf(ibnd, jbnd) = (chf(ibnd, jbnd) + CONJG(chf(jbnd, ibnd))) * 0.5_dp
          chf(jbnd, ibnd) = CONJG(chf(ibnd, jbnd))
       ENDDO
    ENDDO

    lwork = MAX(1, 2*nbnd - 1)
    ALLOCATE(work(lwork))
    ALLOCATE(rwork(MAX(1, 3*nbnd - 2)))

    CALL ZHEEV('V', 'U', nbnd, chf, nbnd, w, work, lwork, rwork, info)
    IF (info /= 0) CALL errore('hamwan2bloch_edi', 'ZHEEV diagonalization failed', info)

    eig = w

    DEALLOCATE(work, rwork)
  END SUBROUTINE hamwan2bloch_edi

  SUBROUTINE get_cfac(nrr, irvec, xk_cryst, cfac)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nrr
    INTEGER, INTENT(IN) :: irvec(3, nrr)
    REAL(dp), INTENT(IN) :: xk_cryst(3)
    COMPLEX(dp), INTENT(OUT) :: cfac(nrr)

    INTEGER :: ir
    REAL(dp) :: rdotk

    DO ir = 1, nrr
       rdotk = twopi * DOT_PRODUCT(xk_cryst, DBLE(irvec(:, ir)))
       cfac(ir) = EXP(ci * rdotk)
    ENDDO
  END SUBROUTINE get_cfac

  SUBROUTINE hamwan2bloch_with_evec(nbnd, nrr, ndegen, eig, evec, chw, cfac)
    !-----------------------------------------------------------------------
    ! Same as hamwan2bloch_edi but ALSO returns eigenvectors U(k).
    ! eig(nbnd)           : eigenvalues
    ! evec(nbnd, nbnd)    : eigenvectors (columns = Bloch states in Wannier basis)
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nbnd, nrr
    INTEGER, INTENT(IN) :: ndegen(nrr)
    REAL(dp), INTENT(OUT) :: eig(nbnd)
    COMPLEX(dp), INTENT(OUT) :: evec(nbnd, nbnd)
    COMPLEX(dp), INTENT(IN) :: cfac(nrr)
    COMPLEX(dp), INTENT(IN) :: chw(nbnd, nbnd, nrr)

    INTEGER :: ibnd, jbnd, ir, info, lwork
    COMPLEX(dp) :: chf(nbnd, nbnd)
    COMPLEX(dp) :: cfac_degen(nrr)
    REAL(dp) :: w(nbnd)
    REAL(dp), ALLOCATABLE :: rwork(:)
    COMPLEX(dp), ALLOCATABLE :: work(:)

    DO ir = 1, nrr
       IF (ndegen(ir) > 0) THEN
          cfac_degen(ir) = cfac(ir) / DBLE(ndegen(ir))
       ELSE
          cfac_degen(ir) = czero
       ENDIF
    ENDDO

    chf = czero
    CALL ZGEMV('n', nbnd*nbnd, nrr, cone, chw, nbnd*nbnd, cfac_degen, 1, cone, chf, 1)

    ! Hermitianize
    DO jbnd = 1, nbnd
       DO ibnd = 1, jbnd
          chf(ibnd, jbnd) = (chf(ibnd, jbnd) + CONJG(chf(jbnd, ibnd))) * 0.5_dp
          chf(jbnd, ibnd) = CONJG(chf(ibnd, jbnd))
       ENDDO
    ENDDO

    lwork = MAX(1, 2*nbnd - 1)
    ALLOCATE(work(lwork))
    ALLOCATE(rwork(MAX(1, 3*nbnd - 2)))

    CALL ZHEEV('V', 'U', nbnd, chf, nbnd, w, work, lwork, rwork, info)
    IF (info /= 0) CALL errore('hamwan2bloch_with_evec', 'ZHEEV failed', info)

    eig = w
    evec = chf   ! columns of chf are the eigenvectors after ZHEEV

    DEALLOCATE(work, rwork)
  END SUBROUTINE hamwan2bloch_with_evec

  SUBROUTINE edmatwan2bloch(nbndsub, nrr, ndegen, edmatw, cfac, edmatf)
    !-----------------------------------------------------------------------
    ! Wannier → Bloch interpolation for electron-defect matrix elements.
    !
    ! Input:  edmatw(nbndsub, nbndsub, nrr) — M(R) in Wannier basis
    !         cfac(nrr) — phase factors exp(ik·R) for the target k-point
    !         ndegen(nrr) — Wigner-Seitz degeneracy
    !
    ! Output: edmatf(nbndsub, nbndsub) — M_W(k) in Wannier basis at k
    !
    ! To get M in the Bloch basis, the caller must rotate:
    !   M_B = U† · edmatf · U
    ! where U comes from diagonalizing H(k) via hamwan2bloch_with_evec.
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nbndsub, nrr
    INTEGER, INTENT(IN) :: ndegen(nrr)
    COMPLEX(dp), INTENT(IN) :: edmatw(nbndsub, nbndsub, nrr)
    COMPLEX(dp), INTENT(IN) :: cfac(nrr)
    COMPLEX(dp), INTENT(OUT) :: edmatf(nbndsub, nbndsub)

    INTEGER :: ir
    COMPLEX(dp) :: cfac_degen(nrr)

    DO ir = 1, nrr
       IF (ndegen(ir) > 0) THEN
          cfac_degen(ir) = cfac(ir) / DBLE(ndegen(ir))
       ELSE
          cfac_degen(ir) = czero
       ENDIF
    ENDDO

    edmatf = czero
    ! M_W(k) = Σ_R exp(ik·R)/ndegen(R) · M(R)
    CALL ZGEMV('n', nbndsub*nbndsub, nrr, cone, edmatw, nbndsub*nbndsub, &
                cfac_degen, 1, cone, edmatf, 1)

  END SUBROUTINE edmatwan2bloch

  SUBROUTINE edmatwan2bloch_2d(nbndsub, nrr, ndegen, edmatw_2d, cfac_ki, cfac_kf, edmatf)
    !-----------------------------------------------------------------------
    ! Full double-FT Wannier → Bloch interpolation for M(R, R').
    !
    ! M_W(k_i, k_f) = Σ_{R, R'} exp(ik_i·R) exp(-ik_f·R')
    !                   × M(R, R') / [ndegen(R_e) · ndegen(R_p)]
    !
    ! cfac_ki(ir) = exp(+i k_i · R_ir)  — phase for initial k
    ! cfac_kf(ir) = exp(+i k_f · R_ir)  — phase for final k
    ! Paper Eq.6 inverse: exp(-ik_i·R) × exp(+ik_f·R')
    !   → uses CONJG(cfac_ki) for bra, cfac_kf for ket
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nbndsub, nrr
    INTEGER, INTENT(IN) :: ndegen(nrr)
    COMPLEX(dp), INTENT(IN) :: edmatw_2d(nbndsub, nbndsub, nrr, nrr)
    COMPLEX(dp), INTENT(IN) :: cfac_ki(nrr), cfac_kf(nrr)
    COMPLEX(dp), INTENT(OUT) :: edmatf(nbndsub, nbndsub)

    INTEGER :: ir_e, ir_p
    COMPLEX(dp) :: phase_e, phase_p, weight
    REAL(dp) :: dege, degp

    edmatf = czero
    DO ir_p = 1, nrr
       degp = DBLE(ndegen(ir_p))
       IF (degp < 0.5_dp) CYCLE
       phase_p = cfac_kf(ir_p) / degp   ! exp(+ik_f·R') / ndeg_p  [paper Eq.6: ket = +]
       DO ir_e = 1, nrr
          dege = DBLE(ndegen(ir_e))
          IF (dege < 0.5_dp) CYCLE
          weight = CONJG(cfac_ki(ir_e)) / dege * phase_p  ! exp(-ik_i·R)/ndeg_e × exp(+ik_f·R')/ndeg_p  [paper Eq.6: bra = -]
          edmatf(:,:) = edmatf(:,:) + weight * edmatw_2d(:,:,ir_e,ir_p)
       ENDDO
    ENDDO

  END SUBROUTINE edmatwan2bloch_2d

END MODULE wan2bloch_edi
