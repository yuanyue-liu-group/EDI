SUBROUTINE init_us_2_perturb( npw_,  igk_, q_, vkb_ ,nat_perturb,ityp_perturb,tau_perturb,nkb_perturb)
    !----------------------------------------------------------------------
    !! Calculates beta functions (Kleinman-Bylander projectors), with
    !! structure factor, for all atoms, in reciprocal space.
    !
    USE kinds,        ONLY : DP
    USE ions_base,    ONLY : nat, ntyp => nsp, ityp, tau
    USE cell_base,    ONLY : tpiba, bg
    USE constants,    ONLY : tpi
    USE gvect,        ONLY : eigts1, eigts2, eigts3, mill, g, ngm
    Use gvect_perturb, Only : eigts1_perturb, eigts2_perturb, eigts3_perturb
    USE wvfct,        ONLY : npwx
    USE vlocal,       ONLY : strf
    USE fft_base,         ONLY : dfftp
    USE uspp_data,           ONLY : nqx, dq, tab, tab_d2y, spline_ps
    USE m_gth,        ONLY : mk_ffnl_gth
    USE splinelib
    USE uspp,         ONLY : nkb, nhtol, nhtolm, indv
    USE uspp_param,   ONLY : upf, lmaxkb, nhm, nh
  
  
  
    USE kinds,        ONLY : DP
    USE ions_base,    ONLY : nat, ityp, tau
    USE cell_base,    ONLY : tpiba, omega
    USE gvect,        ONLY : eigts1, eigts2, eigts3, mill, g
    USE wvfct,        ONLY : npwx
    USE uspp,         ONLY : nkb
    USE fft_base ,    ONLY : dfftp
    
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: npw_
    !! number of PWs 
    INTEGER, INTENT(IN) :: igk_(npw_)
    INTEGER, INTENT(IN) :: nat_perturb,ityp_perturb(nat_perturb),nkb_perturb
    !! indices of G in the list of q+G vectors
    REAL(DP), INTENT(IN) :: q_(3)
    REAL(DP), INTENT(IN) :: tau_perturb(3,nat_perturb)
    !! q vector (2pi/a units)
    COMPLEX(DP), INTENT(OUT) :: vkb_(npwx,nkb_perturb)
    !! beta functions (npw_ <= npwx)
    !
    ! ... Local variables
    !
    INTEGER :: i0, i1, i2, i3, ig, ig_orig, lm, na, nt, nb, ih, jkb
    REAL(DP) :: px, ux, vx, wx, arg
    REAL(DP), ALLOCATABLE :: gk(:,:), qg(:), vq(:), ylm(:,:), vkb1(:,:)
    COMPLEX(DP) :: phase, pref
    COMPLEX(DP), ALLOCATABLE :: sk(:)
    REAL(DP), ALLOCATABLE :: xdata(:)
    INTEGER :: iq
    ! cache blocking parameters
    INTEGER, PARAMETER :: blocksize = 256
    INTEGER :: iblock, numblock, realblocksize
    !
    !
    IF (lmaxkb < 0) RETURN
    !
    CALL start_clock( 'init_us_2_perturb' )
    !
    ! write(*,'(3i4,i5,3f10.5)') size(tab,1), size(tab,2), size(tab,3), size(vq), q_
    !
    ! setting cache blocking size
    numblock = (npw_+blocksize-1)/blocksize
    ALLOCATE( eigts1_perturb(-dfftp%nr1:dfftp%nr1,nat_perturb) )
    ALLOCATE( eigts2_perturb(-dfftp%nr2:dfftp%nr2,nat_perturb) )
    ALLOCATE( eigts3_perturb(-dfftp%nr3:dfftp%nr3,nat_perturb) )
  
  !write(*,*) '1',eigts1_perturb
  !write(*,*) 'shape(eigts1_perturb)',shape(eigts1_perturb)
  !write(*,*) 'shape(eigts1)',shape(eigts1)
    CALL struc_fact( nat_perturb, tau_perturb, ntyp, ityp_perturb, ngm, g, bg, &
                     dfftp%nr1, dfftp%nr2, dfftp%nr3, &
                     strf, eigts1_perturb, eigts2_perturb, eigts3_perturb )
  !write(*,*) 'shape(eigts1_perturb)',shape(eigts1_perturb)
    !
    IF (spline_ps) THEN
      ALLOCATE( xdata(nqx) )
      DO iq = 1, nqx
        xdata(iq) = (iq - 1) * dq
      ENDDO
    ENDIF
    !
  !$omp parallel private(vkb1, sk, qg, vq, ylm, gk, ig_orig, &
  !$omp                  realblocksize, jkb, px, ux, vx, wx, &
  !$omp                  i0, i1, i2, i3, lm, arg, phase, pref)
    !
    ALLOCATE( vkb1(blocksize,nhm) )
    ALLOCATE( sk(blocksize) )
    ALLOCATE( qg(blocksize) )
    ALLOCATE( vq(blocksize) )
    ALLOCATE( ylm(blocksize,(lmaxkb+1)**2) )
    ALLOCATE( gk(3,blocksize) )
    !
  !$omp do
    DO iblock = 1, numblock
       !
       realblocksize = MIN(npw_-(iblock-1)*blocksize,blocksize)
       !
       DO ig = 1, realblocksize
          ig_orig = (iblock-1)*blocksize+ig
          gk(1,ig) = q_(1) + g(1,igk_(ig_orig) )
          gk(2,ig) = q_(2) + g(2,igk_(ig_orig) )
          gk(3,ig) = q_(3) + g(3,igk_(ig_orig) )
          qg(ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
       ENDDO
       !
       CALL ylmr2( (lmaxkb+1)**2, realblocksize, gk, qg, ylm(1:realblocksize,:) )
       !
       ! set now qg=|q+G| in atomic units
       !
       DO ig = 1, realblocksize
          qg(ig) = SQRT(qg(ig))*tpiba
       ENDDO
       !
       ! |beta_lm(q)> = (4pi/omega).Y_lm(q).f_l(q).(i^l).S(q)
       jkb = 0
       DO nt = 1, ntyp
          ! ... calculate beta in G-space using an interpolation table:
          !     f_l(q)=\int _0 ^\infty dr r^2 f_l(r) j_l(q.r)
          DO nb = 1, upf(nt)%nbeta
             !
             IF ( upf(nt)%is_gth ) THEN
                CALL mk_ffnl_gth( nt, nb, realblocksize, omega, qg, vq )
             ELSE
                DO ig = 1, realblocksize
                   IF (spline_ps) THEN
                      vq(ig) = splint(xdata, tab(:,nb,nt), tab_d2y(:,nb,nt), qg(ig))
                   ELSE
                      px = qg(ig) / dq - INT( qg(ig)/dq )
                      ux = 1.d0 - px
                      vx = 2.d0 - px
                      wx = 3.d0 - px
                      i0 = INT( qg(ig)/dq ) + 1
                      i1 = i0 + 1
                      i2 = i0 + 2
                      i3 = i0 + 3
                      vq(ig) = tab(i0,nb,nt) * ux * vx * wx / 6.d0 + &
                               tab(i1,nb,nt) * px * vx * wx / 2.d0 - &
                               tab(i2,nb,nt) * px * ux * wx / 2.d0 + &
                               tab(i3,nb,nt) * px * ux * vx / 6.d0
                   ENDIF
                ENDDO
             ENDIF
             ! add spherical harmonic part  (Y_lm(q)*f_l(q)) 
             DO ih = 1, nh(nt)
                IF (nb == indv(ih,nt) ) THEN
                   !l = nhtol(ih, nt)
                   lm = nhtolm(ih,nt)
                   DO ig = 1, realblocksize
                      vkb1(ig,ih) = ylm(ig,lm) * vq(ig)
                   ENDDO
                ENDIF
             ENDDO
             !
          ENDDO
          !
          ! vkb1 contains all betas including angular part for type nt
          ! now add the structure factor and factor (-i)^l
          !
          DO na = 1, nat_perturb
             ! ordering: first all betas for atoms of type 1
             !           then  all betas for atoms of type 2  and so on
             IF (ityp_perturb(na) == nt) THEN
                !
                arg = ( q_(1) * tau_perturb(1,na) + &
                        q_(2) * tau_perturb(2,na) + &
                        q_(3) * tau_perturb(3,na) ) * tpi
                phase = CMPLX( COS(arg), -SIN(arg) ,KIND=DP)
                !
                DO ig = 1, realblocksize
                   ig_orig = (iblock-1)*blocksize+ig
                   sk(ig) = eigts1_perturb(mill(1,igk_(ig_orig)), na) * &
                            eigts2_perturb(mill(2,igk_(ig_orig)), na) * &
                            eigts3_perturb(mill(3,igk_(ig_orig)), na)
                ENDDO
                !            write(*,*),'eigts3',eigts3_perturb(mill(3,igk_(ig_orig)), na)
                !
                DO ih = 1, nh(nt)
                   jkb = jkb + 1
                   pref = (0.d0, -1.d0)**nhtol(ih, nt) * phase
                   DO ig = 1, realblocksize
                      vkb_((iblock-1)*blocksize+ig, jkb) = vkb1(ig,ih) * sk(ig) * pref
                   ENDDO
                   ! clean up garbage in the last block
                   IF (iblock == numblock) THEN
                      DO ig = npw_+1, npwx
                         vkb_(ig, jkb) = (0.0_DP, 0.0_DP)
                      ENDDO
                   ENDIF
                ENDDO
                !
             ENDIF
             !
          ENDDO
          !
       ENDDO
       !
    ENDDO
  !$omp end do nowait
    DEALLOCATE( gk )
    DEALLOCATE( ylm )
    DEALLOCATE( vq )
    DEALLOCATE( qg )
    DEALLOCATE( sk )
    DEALLOCATE( vkb1 )
    deALLOCATE( eigts1_perturb)
    deALLOCATE( eigts2_perturb)
    deALLOCATE( eigts3_perturb)
  !$omp end parallel
    !
    IF (spline_ps) DEALLOCATE( xdata )
    !
    CALL stop_clock( 'init_us_2_perturb' )
    !
    !
    RETURN
    !
  END SUBROUTINE init_us_2_perturb