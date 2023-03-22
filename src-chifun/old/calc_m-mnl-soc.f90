 SUBROUTINE calcmdefect_mnl_ks_soc(ibnd0,ibnd,ik0,ik)
    Use kinds,    Only : dp
    USE klist , ONLY: nks, nelec, xk, wk, degauss, ngauss, igk_k, ngk
    Use becmod,          Only : calbec, allocate_bec_type, deallocate_bec_type
    Use edic_mod, Only : becp1_perturb,becp2_perturb
    USE edic_mod, Only :  V_0, V_p
    USE uspp_param, ONLY: nh
    USE wvfct, ONLY: npwx, nbnd, wg, et, g2kin
    USE uspp, ONLY: nkb, vkb, dvan, dvan_so
    Use edic_mod,only : evc1,evc2,evc3,evc4
    Use edic_mod, Only: m_nloc

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

    Call init_us_2_sc (ngk(ik0), igk_k(1,ik0), xk (1, ik0), vkb_perturb,V_0%nat,V_0%ityp,V_0%tau,nkb_perturb)
    CALL calbec ( ngk(ik0), vkb_perturb, evc1, becp1_perturb )
         
    CALL init_us_2_sc (ngk(ik), igk_k(1,ik), xk (1, ik), vkb_perturb,V_0%nat,V_0%ityp,V_0%tau,nkb_perturb)
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
      
      CALL init_us_2_sc (ngk(ik0), igk_k(1,ik0), xk (1, ik0), vkb_perturb,V_p%nat,V_p%ityp,V_p%tau,nkb_perturb)
      CALL calbec ( ngk(ik0), vkb_perturb, evc1, becp1_perturb )
      
      CALL init_us_2_sc (ngk(ik), igk_k(1,ik), xk (1, ik), vkb_perturb,V_p%nat,V_p%ityp,V_p%tau,nkb_perturb)
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
