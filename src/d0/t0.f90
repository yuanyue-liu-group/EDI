
      nkb_perturb=0
  
    DO nt_perturb = 1, V_p%ntyp
        DO na_perturb = 1, V_p%nat
           IF(V_p%ityp (na_perturb) == nt_perturb)THEN
               nkb_perturb = nkb_perturb + nh (nt_perturb)
           ENDIF
        ENDDO
     ENDDO
      
    CALL allocate_bec_type ( nkb, nbnd, becp )
    CALL allocate_bec_type ( nkb, nbnd, becp1 )
    CALL allocate_bec_type ( nkb, nbnd, becp2 )
    CALL allocate_bec_type ( nkb_perturb, nbnd, becp_perturb )
    CALL allocate_bec_type ( nkb_perturb, nbnd, becp1_perturb )
    CALL allocate_bec_type ( nkb_perturb, nbnd, becp2_perturb )
    ALLOCATE(vkb_perturb(npwx,nkb_perturb))
    
    !CALL init_us_2 (ngk(ik), igk_k(1,ik), xk (1, ik), vkb)
    !CALL calbec ( ngk(ik), vkb, evc, becp )
    
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
!                ENDIF
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
       
      CALL deallocate_bec_type (  becp )
      CALL deallocate_bec_type (  becp1 )
      CALL deallocate_bec_type (  becp2 )
      CALL deallocate_bec_type (  becp_perturb )
      CALL deallocate_bec_type (  becp1_perturb )
      CALL deallocate_bec_type (  becp2_perturb )
      DEALLOCATE(vkb_perturb)
  
