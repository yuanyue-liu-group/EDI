SUBROUTINE get_betavkb( npw_, igk_, q_, vkb_, &
                         nat_perturb, ityp_perturb, tau_perturb, nkb_perturb )
    USE kinds,        ONLY : DP
    USE ions_base,    ONLY : nat, ntyp => nsp, ityp, tau
    USE cell_base,    ONLY : tpiba, bg, omega
    USE constants,    ONLY : tpi
    USE gvect,        ONLY : eigts1, eigts2, eigts3, mill, g, ngm
    Use edic_mod,     Only : eigts1_perturb, eigts2_perturb, eigts3_perturb
    USE wvfct,        ONLY : npwx
    USE vlocal,       ONLY : strf
    USE fft_base,     ONLY : dfftp
    USE beta_mod,     ONLY : interp_beta
    USE uspp,         ONLY : nkb, nhtol, nhtolm, indv
    USE uspp_param,   ONLY : upf, lmaxkb, nhm, nh, nbetam
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: npw_
    INTEGER, INTENT(IN) :: igk_(npw_)
    INTEGER, INTENT(IN) :: nat_perturb, ityp_perturb(nat_perturb), nkb_perturb
    REAL(DP), INTENT(IN) :: q_(3)
    REAL(DP), INTENT(IN) :: tau_perturb(3,nat_perturb)
    COMPLEX(DP), INTENT(OUT) :: vkb_(npwx,nkb_perturb)
    !
    INTEGER :: ig, lm, na, nt, nb, ih, jkb
    REAL(DP) :: arg
    REAL(DP), ALLOCATABLE :: gk(:,:), qg(:), ylm(:,:), vq(:,:), vkb1(:,:)
    COMPLEX(DP) :: phase, pref
    COMPLEX(DP), ALLOCATABLE :: sk(:)
    !
    IF (lmaxkb < 0) RETURN
    !
    CALL start_clock( 'get_betavkb' )
    !
    ALLOCATE( eigts1_perturb(-dfftp%nr1:dfftp%nr1, nat_perturb) )
    ALLOCATE( eigts2_perturb(-dfftp%nr2:dfftp%nr2, nat_perturb) )
    ALLOCATE( eigts3_perturb(-dfftp%nr3:dfftp%nr3, nat_perturb) )
    !
    CALL struc_fact( nat_perturb, tau_perturb, ntyp, ityp_perturb, ngm, g, bg, &
                     dfftp%nr1, dfftp%nr2, dfftp%nr3, &
                     strf, eigts1_perturb, eigts2_perturb, eigts3_perturb )
    !
    ALLOCATE( gk(3, npw_) )
    ALLOCATE( qg(npw_) )
    ALLOCATE( ylm(npw_, (lmaxkb+1)**2) )
    ALLOCATE( vq(npw_, nbetam) )
    ALLOCATE( vkb1(npw_, nhm) )
    ALLOCATE( sk(npw_) )
    !
    vkb_(:,:) = (0.0_DP, 0.0_DP)
    !
    DO ig = 1, npw_
       gk(1,ig) = q_(1) + g(1, igk_(ig))
       gk(2,ig) = q_(2) + g(2, igk_(ig))
       gk(3,ig) = q_(3) + g(3, igk_(ig))
       qg(ig) = gk(1,ig)**2 + gk(2,ig)**2 + gk(3,ig)**2
    ENDDO
    !
    CALL ylmr2( (lmaxkb+1)**2, npw_, gk, qg, ylm )
    !
    DO ig = 1, npw_
       qg(ig) = SQRT(qg(ig)) * tpiba
    ENDDO
    !
    jkb = 0
    DO nt = 1, ntyp
       !
       CALL interp_beta( nt, npw_, qg, vq )
       !
       DO ih = 1, nh(nt)
          nb = indv(ih, nt)
          lm = nhtolm(ih, nt)
          DO ig = 1, npw_
             vkb1(ig, ih) = ylm(ig, lm) * vq(ig, nb)
          ENDDO
       ENDDO
       !
       DO na = 1, nat_perturb
          IF (ityp_perturb(na) == nt) THEN
             !
             arg = ( q_(1) * tau_perturb(1,na) + &
                     q_(2) * tau_perturb(2,na) + &
                     q_(3) * tau_perturb(3,na) ) * tpi
             phase = CMPLX( COS(arg), -SIN(arg), KIND=DP )
             !
             DO ig = 1, npw_
                sk(ig) = eigts1_perturb(mill(1, igk_(ig)), na) * &
                         eigts2_perturb(mill(2, igk_(ig)), na) * &
                         eigts3_perturb(mill(3, igk_(ig)), na)
             ENDDO
             !
             DO ih = 1, nh(nt)
                jkb = jkb + 1
                pref = (0.d0, -1.d0)**nhtol(ih, nt) * phase
                DO ig = 1, npw_
                   vkb_(ig, jkb) = vkb1(ig, ih) * sk(ig) * pref
                ENDDO
                DO ig = npw_+1, npwx
                   vkb_(ig, jkb) = (0.0_DP, 0.0_DP)
                ENDDO
             ENDDO
             !
          ENDIF
       ENDDO
    ENDDO
    !
    DEALLOCATE( gk, ylm, vq, qg, sk, vkb1 )
    DEALLOCATE( eigts1_perturb, eigts2_perturb, eigts3_perturb )
    !
    CALL stop_clock( 'get_betavkb' )
    !
    RETURN
    !
END SUBROUTINE get_betavkb
