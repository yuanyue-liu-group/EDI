subroutine calcmdefect_ml_rs_noncolin(ibnd0,ibnd,ik0,ik)
    Use kinds,          Only : dp
    USE cell_base,       ONLY : alat, ibrav, omega, at, bg, celldm, wmass
    USE constants, ONLY: tpi, e2, eps6,pi
    Use fft_base,       Only : dfftp, dffts
    USE fft_interfaces, ONLY : fwfft, invfft
    USE klist , ONLY: nks, nelec, xk, wk, degauss, ngauss, igk_k, ngk
    Use wavefunctions_calcmdefect, Only : evc1,evc2,evc3,evc4,&
                               psic1, psic2, psic3, psic4
    USE wvfct, ONLY: npwx, nbnd, wg, et, g2kin
    USE v_type, Only :  V_0, V_loc
    Use input_parameters_calcmdefect, Only: m_loc
    Implicit None 

    complex(dp) :: ml_up, ml_down
    real(dp) :: d1, d2, d3

    INTEGER :: ibnd, ik, ik0,ibnd0, ig

    REAL(dp)::arg,argt,argt2
    COMPLEX(DP)::phase
    INTEGER :: npw,  ispin, nbndup, nbnddown, &
                nk , ikk,ikk0,  inr, ig1, ig2
    INTEGER:: irx,iry,irz
    INTEGER:: irx2,iry2,irz2
    INTEGER:: irx1,iry1,irz1
    
    INTEGER :: ix0,ix1,ix2
    INTEGER :: iy0,iy1,iy2
    INTEGER :: iz0,iz1,iz2, ikpsi0, ikpsi1, ikpsi2
    integer :: ir1mod,ir2mod,ir3mod,irnmod

    !write(*,*) V_loc(20:40,1)
    !write(*,*) V_loc(20:40,2)
    !write(*,*) V_0%plot(20:40)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!! vl in real super2prim, module
    !auxr(:) =  vrs(:,1)
    !psiprod(:)=0.00
    !vgk_perturb(:)=0.00
    ml_up=0
    ml_down=0
    d1=((1.0/dffts%nr1*at(1,1))*(xk(1,ik)-xk(1,ik0)) +&
        (1.0/dffts%nr1*at(2,1))*(xk(2,ik)-xk(2,ik0)) +&
        (1.0/dffts%nr1*at(3,1))*(xk(3,ik)-xk(3,ik0)) )*tpi 
    d2=((1.0/dffts%nr2*at(1,2))*(xk(1,ik)-xk(1,ik0)) +&
        (1.0/dffts%nr2*at(2,2))*(xk(2,ik)-xk(2,ik0)) +&
        (1.0/dffts%nr2*at(3,2))*(xk(3,ik)-xk(3,ik0)) )*tpi 
    d3=((1.0/dffts%nr3*at(1,3))*(xk(1,ik)-xk(1,ik0)) +&
        (1.0/dffts%nr3*at(2,3))*(xk(2,ik)-xk(2,ik0)) +&
        (1.0/dffts%nr3*at(3,3))*(xk(3,ik)-xk(3,ik0)) )*tpi 
    
    psic2(1:dffts%nnr) = (0.d0,0.d0)
    psic1(1:dffts%nnr) = (0.d0,0.d0)
    psic4(1:dffts%nnr) = (0.d0,0.d0)
    psic3(1:dffts%nnr) = (0.d0,0.d0)

    DO ig = 1, ngk(ik)
       psic2 (dffts%nl (igk_k(ig,ik) ) ) = evc2 (ig, ibnd)
    ENDDO
    DO ig = 1, ngk(ik0)
       psic1 (dffts%nl (igk_k(ig,ik0) ) ) = evc1 (ig, ibnd0)
    ENDDO

    DO ig = 1, ngk(ik)
       psic4 (dffts%nl (igk_k(ig,ik) ) ) = evc2 (ig+npwx, ibnd)
    ENDDO
    DO ig = 1, ngk(ik0)
       psic3 (dffts%nl (igk_k(ig,ik0) ) ) = evc1 (ig+npwx, ibnd0)
    ENDDO

    

    CALL invfft ('Wave', psic2, dffts)
    CALL invfft ('Wave', psic1, dffts)
    CALL invfft ('Wave', psic4, dffts)
    CALL invfft ('Wave', psic3, dffts)
    
    !write(*,*) psic1(1:10)
    !write(*,*) psic2(1:10)
    !write(*,*) psic3(1:10)
    !write(*,*) psic4(1:10)
  
    arg=0
    inr=0
    do irz =0, V_0%nr3-1
    ir3mod=irz-(irz/(dffts%nr3))*dffts%nr3
    do iry =0, V_0%nr2-1
    ir2mod=iry-(iry/(dffts%nr2))*dffts%nr2
    do irx =0, V_0%nr1-1
    ir1mod=irx-(irx/(dffts%nr1))*dffts%nr1
    
    arg=irz*d3+iry*d2+irx*d1
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !move vloc center 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    arg=tpi*(real(irx)/dffts%nr1*at(1,1)+real(iry)/dffts%nr2*at(1,2)+real(irz)/dffts%nr3*at(1,3))*(xk(1,ik)-xk(1,ik0)) +&
        tpi*(real(irx)/dffts%nr1*at(2,1)+real(iry)/dffts%nr2*at(2,2)+real(irz)/dffts%nr3*at(2,3))*(xk(2,ik)-xk(2,ik0)) +&
        tpi*(real(irx)/dffts%nr1*at(3,1)+real(iry)/dffts%nr2*at(3,2)+real(irz)/dffts%nr3*at(3,3))*(xk(3,ik)-xk(3,ik0))   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! shift arg center
    arg=irz*d3+(iry-iry/(V_0%nr2/2+1)*V_0%nr2)*d2+(irx-irx/(V_0%nr1/2+1)*V_0%nr1)*d1
    !arg=irz*d3+(iry-iry/(dffts%nr2/2+1)*dffts%nr1)*d2+(irx-irx/(dffts%nr1/2+1)*dffts%nr1)*d1
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !arg=irz*d3+iry*d2+irx*d1
    
    phase=CMPLX(COS(arg),SIN(arg),kind=dp)
    inr=inr+1
    irnmod=(ir3mod)*dffts%nr1*dffts%nr2+(ir2mod)*dffts%nr1+ir1mod+1
    ml_up=ml_up+CONJG(psic1(irnmod))*psic2(irnmod)*phase*V_loc(inr, 1)
    ml_down=ml_down+CONJG(psic3(irnmod))*psic4(irnmod)*phase*V_loc(inr, 2)

    if ( irnmod<0 .or. irnmod>dffts%nnr ) then
       write (*,*) 'grid mismatch', irnmod, dffts%nnr 
    endif
    
    
    if (irz==dffts%nr3/2) then
            argt= atan2(real(CONJG(psic1(irnmod))*psic2(irnmod)*phase),aimag(CONJG(psic1(irnmod))*psic2(irnmod)*phase))
            argt2= atan2(real(CONJG(psic1(irnmod))*psic2(irnmod)),aimag(CONJG(psic1(irnmod))*psic2(irnmod)))
            if (argt<0) argt=argt+tpi
    endif

       
    enddo
    enddo
    enddo
    ml_up=ml_up/dffts%nnr
    ml_down=ml_down/dffts%nnr
    m_loc = ml_up+ml_down
    write (*,1001) 'Ml_up ki->kf ',ik0,ik, ml_up, abs(ml_up)
    write (*,1001) 'Ml_down ki->kf ',ik0,ik,  ml_down, abs(ml_down)
    write (*,1002) 'Ml ki->kf ',ik0,ik,  ml_up+ml_down, abs(ml_up+ml_down)
1001 format(A16,I9,I9," ( ",e17.9," , ",e17.9," ) ",e17.9)
1002 format(A16,I9,I9," ( ",e17.9," , ",e17.9," ) ",e17.9/)
   
    
    
 end subroutine calcmdefect_ml_rs_noncolin

 SUBROUTINE calcmdefect_mnl_ks_soc(ibnd0,ibnd,ik0,ik)

    Use kinds,    Only : dp
    USE klist , ONLY: nks, nelec, xk, wk, degauss, ngauss, igk_k, ngk
    Use becmod,          Only : calbec, allocate_bec_type, deallocate_bec_type
    Use becomd_perturb, Only : becp1_perturb,becp2_perturb
    USE v_type, Only :  V_0, V_p
    USE uspp_param, ONLY: nh
    USE wvfct, ONLY: npwx, nbnd, wg, et, g2kin
    USE uspp, ONLY: nkb, vkb, dvan, dvan_so
    Use wavefunctions_calcmdefect, Only : evc1,evc2,evc3,evc4
    Use input_parameters_calcmdefect, Only: m_nloc

    Implicit None

    INTEGER :: ibnd, ik, ik0,ibnd0
    integer :: nr1x_perturb, nr2x_perturb, nr3x_perturb, nr1_perturb, nr2_perturb, nr3_perturb, &
                nat_perturb, ntyp_perturb, ibrav_perturb, plot_num_perturb,  i_perturb,nkb_perturb
    integer :: iunplot_perturb, ios_perturb, ipol_perturb, na_perturb, nt_perturb, &
                ir_perturb, ndum_perturb
    integer :: ijkb0, ih, jh, ikb, jkb
    Complex(dp) :: mnl, mnl_d, mnl_p, mnltot
    Complex(dp), allocatable :: vkb_perturb(:,:)
 

    nkb_perturb=0

    DO nt_perturb = 1, V_0%ntyp
        DO na_perturb = 1, V_0%nat
           IF(V_0%ityp (na_perturb) == nt_perturb)THEN
               nkb_perturb = nkb_perturb + nh (nt_perturb)
           ENDIF
        ENDDO
     ENDDO

    CALL allocate_bec_type ( nkb_perturb, nbnd, becp1_perturb )
    CALL allocate_bec_type ( nkb_perturb, nbnd, becp2_perturb )
    ALLOCATE(vkb_perturb(npwx,nkb_perturb))

    Call init_us_2_perturb (ngk(ik0), igk_k(1,ik0), xk (1, ik0), vkb_perturb,V_0%nat,V_0%ityp,V_0%tau,nkb_perturb)
    CALL calbec ( ngk(ik0), vkb_perturb, evc1, becp1_perturb )
         
    CALL init_us_2_perturb (ngk(ik), igk_k(1,ik), xk (1, ik), vkb_perturb,V_0%nat,V_0%ityp,V_0%tau,nkb_perturb)
    CALL calbec ( ngk(ik), vkb_perturb, evc2, becp2_perturb )
         
         

    ijkb0 = 0
    mnl_d=0
    mnltot=0

    DO nt_perturb = 1, V_0%ntyp
      DO na_perturb = 1, V_0%nat
         IF(V_0%ityp (na_perturb) == nt_perturb)THEN
            DO ih = 1, nh (nt_perturb)
               ikb = ijkb0 + ih
               
               mnl_d=mnl_d+&
               conjg(becp1_perturb%nc(ikb,1,ibnd0))*becp2_perturb%nc(ikb,1,ibnd) * dvan_so(ih,ih,1,nt_perturb)&
               + conjg(becp1_perturb%nc(ikb,1,ibnd0))*becp2_perturb%nc(ikb,2,ibnd) * dvan_so(ih,ih,2,nt_perturb)&
               + conjg(becp1_perturb%nc(ikb,2,ibnd0))*becp2_perturb%nc(ikb,1,ibnd) * dvan_so(ih,ih,3,nt_perturb)&
               + conjg(becp1_perturb%nc(ikb,2,ibnd0))*becp2_perturb%nc(ikb,2,ibnd) * dvan_so(ih,ih,4,nt_perturb)
               DO jh = ( ih + 1 ), nh(nt_perturb)
                  jkb = ijkb0 + jh

                  mnl_d=mnl_d + &
                     conjg(becp1_perturb%nc(ikb, 1, ibnd0))*becp2_perturb%nc(jkb, 1, ibnd) &
                     * dvan_so(ih,jh,1,nt_perturb) &
                     + &
                     conjg(becp1_perturb%nc(jkb, 1, ibnd0))*becp2_perturb%nc(ikb, 1, ibnd) &
                     * dvan_so(jh,ih,1,nt_perturb) &

                     + &

                     conjg(becp1_perturb%nc(ikb, 1, ibnd0))*becp2_perturb%nc(jkb, 2, ibnd) &
                     * dvan_so(ih,jh,2,nt_perturb) &
                     + &
                     conjg(becp1_perturb%nc(jkb, 1, ibnd0))*becp2_perturb%nc(ikb, 2, ibnd) &
                     * dvan_so(jh,ih,2,nt_perturb) &

                     + &

                     conjg(becp1_perturb%nc(ikb, 2, ibnd0))*becp2_perturb%nc(jkb, 1, ibnd) &
                     * dvan_so(ih,jh,3,nt_perturb) &
                     + &
                     conjg(becp1_perturb%nc(jkb, 2, ibnd0))*becp2_perturb%nc(ikb, 1, ibnd) &
                     * dvan_so(jh,ih,3,nt_perturb) &

                     + &

                     conjg(becp1_perturb%nc(ikb, 2, ibnd0))*becp2_perturb%nc(jkb, 2, ibnd) &
                     * dvan_so(ih,jh,4,nt_perturb) &
                     + &
                     conjg(becp1_perturb%nc(jkb, 2, ibnd0))*becp2_perturb%nc(ikb, 2, ibnd) &
                     * dvan_so(jh,ih,4,nt_perturb) 
               End Do
            End Do
            ijkb0 = ijkb0 + nh (nt_perturb)
         End If
      End Do
   End Do
   mnltot=mnltot+mnl_d*wg(ibnd,ik)!
   CALL deallocate_bec_type (  becp1_perturb )
   CALL deallocate_bec_type (  becp2_perturb )
   DEALLOCATE(vkb_perturb)

   nkb_perturb=0
  
  
      DO nt_perturb = 1, V_p%ntyp
         DO na_perturb = 1, V_p%nat
            IF(V_p%ityp (na_perturb) == nt_perturb)THEN
                nkb_perturb = nkb_perturb + nh (nt_perturb)
            ENDIF
         ENDDO
      ENDDO
      
      
      CALL allocate_bec_type ( nkb_perturb, nbnd, becp1_perturb )
      CALL allocate_bec_type ( nkb_perturb, nbnd, becp2_perturb )
      ALLOCATE(vkb_perturb(npwx,nkb_perturb))
      
      CALL init_us_2_perturb (ngk(ik0), igk_k(1,ik0), xk (1, ik0), vkb_perturb,V_p%nat,V_p%ityp,V_p%tau,nkb_perturb)
      CALL calbec ( ngk(ik0), vkb_perturb, evc1, becp1_perturb )
      
      CALL init_us_2_perturb (ngk(ik), igk_k(1,ik), xk (1, ik), vkb_perturb,V_p%nat,V_p%ityp,V_p%tau,nkb_perturb)
      CALL calbec ( ngk(ik), vkb_perturb, evc2, becp2_perturb )
      
      
      ijkb0 = 0
      mnl_p=0
      mnltot=0
      DO nt_perturb = 1, V_p%ntyp
         DO na_perturb = 1, V_p%nat
            IF(V_p%ityp (na_perturb) == nt_perturb)THEN
               DO ih = 1, nh (nt_perturb)
                  ikb = ijkb0 + ih
                  
                     mnl_p=mnl_p&
                     + conjg(becp1_perturb%nc(ikb,1,ibnd0))*becp2_perturb%nc(ikb,1,ibnd) * dvan_so(ih,ih,1,nt_perturb) &
                     + conjg(becp1_perturb%nc(ikb,1,ibnd0))*becp2_perturb%nc(ikb,2,ibnd) * dvan_so(ih,ih,2,nt_perturb) &
                     + conjg(becp1_perturb%nc(ikb,2,ibnd0))*becp2_perturb%nc(ikb,1,ibnd) * dvan_so(ih,ih,3,nt_perturb) &
                     + conjg(becp1_perturb%nc(ikb,2,ibnd0))*becp2_perturb%nc(ikb,2,ibnd) * dvan_so(ih,ih,4,nt_perturb)
      
                  DO jh = ( ih + 1 ), nh(nt_perturb)
                     jkb = ijkb0 + jh
               
                        mnl_p=mnl_p + &
                           conjg(becp1_perturb%nc(ikb, 1, ibnd0))*becp2_perturb%nc(jkb, 1, ibnd) &
                           * dvan_so(ih,jh,1,nt_perturb) &
                           + &
                           conjg(becp1_perturb%nc(jkb, 1, ibnd0))*becp2_perturb%nc(ikb, 1, ibnd) &
                           * dvan_so(jh,ih,1,nt_perturb) &

                           + &

                           conjg(becp1_perturb%nc(ikb, 1, ibnd0))*becp2_perturb%nc(jkb, 2, ibnd) &
                           * dvan_so(ih,jh,2,nt_perturb) &
                           + &
                           conjg(becp1_perturb%nc(jkb, 1, ibnd0))*becp2_perturb%nc(ikb, 2, ibnd) &
                           * dvan_so(jh,ih,2,nt_perturb) &

                           + &

                           conjg(becp1_perturb%nc(ikb, 2, ibnd0))*becp2_perturb%nc(jkb, 1, ibnd) &
                           * dvan_so(ih,jh,3,nt_perturb) &
                           + &
                           conjg(becp1_perturb%nc(jkb, 2, ibnd0))*becp2_perturb%nc(ikb, 1, ibnd) &
                           * dvan_so(jh,ih,3,nt_perturb) &

                           + &

                           conjg(becp1_perturb%nc(ikb, 2, ibnd0))*becp2_perturb%nc(jkb, 2, ibnd) &
                           * dvan_so(ih,jh,4,nt_perturb) &
                           + &
                           conjg(becp1_perturb%nc(jkb, 2, ibnd0))*becp2_perturb%nc(ikb, 2, ibnd) &
                           * dvan_so(jh,ih,4,nt_perturb) 

      
                  ENDDO
      
               ENDDO
               ijkb0 = ijkb0 + nh (nt_perturb)
            ENDIF
         ENDDO
      ENDDO
      mnltot=mnltot+mnl_p*wg(ibnd,ik)!
      m_nloc=mnl_d-mnl_p
   
      CALL deallocate_bec_type (  becp1_perturb )
      CALL deallocate_bec_type (  becp2_perturb )
      DEALLOCATE(vkb_perturb)
  1001 format(A16,I9,I9, " ( ",e17.9," , ",e17.9," ) ",e17.9)
  1002 format(A16,I9,I9," ( ",e17.9," , ",e17.9," ) ",e17.9/)
      write (*,1001) 'Mnl_d ki->kf ', ik0,ik, mnl_d, abs(mnl_d)
      write (*,1001) 'Mnl_p ki->kf ', ik0,ik, mnl_p, abs(mnl_p)
      write (*,1002) 'Mnl ki->kf ', ik0,ik, mnl_d-mnl_p, abs(mnl_d-mnl_p)
 End Subroutine