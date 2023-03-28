      SUBROUTINE calcmdefect_mnl_ks_soc(ibnd0,ibnd,ik0,ik, V_d, V_p)
    
         USE becmod, ONLY: becp, calbec, allocate_bec_type, deallocate_bec_type
         USE becmod, ONLY: becp1,becp2,becp_perturb,becp1_perturb,becp2_perturb 
          
         INTEGER :: ibnd, ik, ik0,ibnd0
         type(V_file) :: V_d, V_p
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!! initialization
         nkb_perturb=0
     
     
         DO nt_perturb = 1, V_d%ntyp
            DO na_perturb = 1, V_d%nat
               IF(V_d%ityp (na_perturb) == nt_perturb)THEN
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
         
         !!!!!! initialization
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         CALL init_us_2 (ngk(ik), igk_k(1,ik), xk (1, ik), vkb)
         CALL calbec ( ngk(ik), vkb, evc, becp )
         
         CALL init_us_2_perturb (ngk(ik0), igk_k(1,ik0), xk (1, ik0), vkb_perturb,V_d%nat,V_d%ityp,V_d%tau,nkb_perturb)
         CALL calbec ( ngk(ik0), vkb_perturb, evc1, becp1_perturb )
         
         CALL init_us_2_perturb (ngk(ik), igk_k(1,ik), xk (1, ik), vkb_perturb,V_d%nat,V_d%ityp,V_d%tau,nkb_perturb)
         CALL calbec ( ngk(ik), vkb_perturb, evc2, becp2_perturb )
         
         ijkb0 = 0
         mnl_d=0
         mnltot=0
   
         DO nt_perturb = 1, V_d%ntyp
            DO na_perturb = 1, V_d%nat
               IF(V_d%ityp (na_perturb) == nt_perturb)THEN
                  DO ih = 1, nh (nt_perturb)
                     ikb = ijkb0 + ih
                     IF(gamma_only)THEN
                        mnl_d=mnl_d+becp1%r(ikb,ibnd0)*becp2%r(ikb,ibnd) &
                           * dvan_so(ih,ih,1,nt_perturb)
                     ELSE
                        !write(*,*) dvan_so(ih,ih,1:4,nt_perturb)
                        mnl_d=mnl_d+&
                           conjg(becp1_perturb%nc(ikb,1,ibnd0))*becp2_perturb%nc(ikb,1,ibnd) * dvan_so(ih,ih,1,nt_perturb)&
                           + conjg(becp1_perturb%nc(ikb,1,ibnd0))*becp2_perturb%nc(ikb,2,ibnd) * dvan_so(ih,ih,2,nt_perturb)&
                           + conjg(becp1_perturb%nc(ikb,2,ibnd0))*becp2_perturb%nc(ikb,1,ibnd) * dvan_so(ih,ih,3,nt_perturb)&
                           + conjg(becp1_perturb%nc(ikb,2,ibnd0))*becp2_perturb%nc(ikb,2,ibnd) * dvan_so(ih,ih,4,nt_perturb)
                     ENDIF
                     DO jh = ( ih + 1 ), nh(nt_perturb)
                        jkb = ijkb0 + jh
                        IF(gamma_only)THEN
                           mnl_d=mnl_d + &
                              (becp1%r(ikb,ibnd0)*becp2%r(jkb,ibnd)+&
                                 becp1%r(jkb,ibnd0)*becp2%r(ikb,ibnd))&
                               * dvan_so(ih,jh,1,nt_perturb)
                        ELSE
                           !write(*,*) 'test'
                           !write(*,*) dvan_so(ih,jh,1:4,nt_perturb)
                           !write(*,*) dvan_so(jh,ih,1:4,nt_perturb)
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
                              
                        ENDIF
         
                     ENDDO
                     
                  enddo
                  ijkb0 = ijkb0 + nh(nt_perturb)
               endif
            enddo
         enddo
         mnltot=mnltot+mnl_d*wg(ibnd,ik)


      CALL deallocate_bec_type (  becp )
      CALL deallocate_bec_type (  becp1 )
      CALL deallocate_bec_type (  becp2 )
      CALL deallocate_bec_type (  becp_perturb )
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
      
      
      CALL allocate_bec_type ( nkb, nbnd, becp )
      CALL allocate_bec_type ( nkb, nbnd, becp1 )
      CALL allocate_bec_type ( nkb, nbnd, becp2 )
      CALL allocate_bec_type ( nkb_perturb, nbnd, becp_perturb )
      CALL allocate_bec_type ( nkb_perturb, nbnd, becp1_perturb )
      CALL allocate_bec_type ( nkb_perturb, nbnd, becp2_perturb )
      ALLOCATE(vkb_perturb(npwx,nkb_perturb))
      
      !!!!!! initialization
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL init_us_2 (ngk(ik), igk_k(1,ik), xk (1, ik), vkb)
      CALL calbec ( ngk(ik), vkb, evc, becp )
      
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
                  IF(gamma_only)THEN
                     mnl_p=mnl_p+becp1%r(ikb,ibnd0)*becp2%r(ikb,ibnd) &
                        * dvan_so(ih,ih,1,nt_perturb)
                  ELSE
                     mnl_p=mnl_p&
                     + conjg(becp1_perturb%nc(ikb,1,ibnd0))*becp2_perturb%nc(ikb,1,ibnd) * dvan_so(ih,ih,1,nt_perturb) &
                     + conjg(becp1_perturb%nc(ikb,1,ibnd0))*becp2_perturb%nc(ikb,2,ibnd) * dvan_so(ih,ih,2,nt_perturb) &
                     + conjg(becp1_perturb%nc(ikb,2,ibnd0))*becp2_perturb%nc(ikb,1,ibnd) * dvan_so(ih,ih,3,nt_perturb) &
                     + conjg(becp1_perturb%nc(ikb,2,ibnd0))*becp2_perturb%nc(ikb,2,ibnd) * dvan_so(ih,ih,4,nt_perturb)
                  ENDIF
                  DO jh = ( ih + 1 ), nh(nt_perturb)
                     jkb = ijkb0 + jh
                     IF(gamma_only)THEN
                        mnl_p=mnl_p + &
                           (becp1%r(ikb,ibnd0)*becp2%r(jkb,ibnd)+&
                              becp1%r(jkb,ibnd0)*becp2%r(ikb,ibnd))&
                            * dvan_so(ih,jh,1,nt_perturb)
                     ELSE
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

                     ENDIF
      
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
  1001 format(A16,I9,I9, " ( ",e17.9," , ",e17.9," ) ",e17.9)
  1002 format(A16,I9,I9," ( ",e17.9," , ",e17.9," ) ",e17.9/)
      write (stdout,1001) 'Mnl_d ki->kf ', ik0,ik, mnl_d, abs(mnl_d)
      write (stdout,1001) 'Mnl_p ki->kf ', ik0,ik, mnl_p, abs(mnl_p)
      write (stdout,1002) 'Mnl ki->kf ', ik0,ik, mnl_d-mnl_p, abs(mnl_d-mnl_p)
         
         END SUBROUTINE calcmdefect_mnl_ks_soc
    
