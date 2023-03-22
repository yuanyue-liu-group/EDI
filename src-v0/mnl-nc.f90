SUBROUTINE calcmdefect_mnl_ks_noncolin(ibnd,ibnd0,ik,ik0,v_mnl,mnonlocal)
  USE kinds, ONLY: DP
  USE io_global, ONLY: stdout

  USE wvfct, ONLY: npwx, nbnd
  USE klist , ONLY: nks, nelec, xk, wk, degauss, ngauss, igk_k, ngk

  USE uspp_param, ONLY: nh
  USE uspp, ONLY: nkb, vkb, dvan

  Use edic_mod, Only : evc1,evc2,V_file
  Use edic_mod, Only : becp1_perturb,becp2_perturb

  USE becmod, ONLY: becp, calbec, allocate_bec_type, deallocate_bec_type
     
  COMPLEX(DP), ALLOCATABLE ::vkb_perturb(:,:)
  COMPLEX(DP) :: mnl
  COMPLEX(DP) ,intent(inout):: mnonlocal
  type(V_file),intent(in):: v_mnl
  INTEGER :: ibnd, ik, ik0,ibnd0, na_perturb, nt_perturb
  integer ::nkb_perturb
  integer :: ijkb0, ih, jh, ikb, jkb
  
  !INTEGER,intent(in) :: bnd_idx_i, kp_idx_i, kp_idx_f,bnd_idx_f

  nkb_perturb=0
  DO nt_perturb = 1, v_mnl%ntyp
     DO na_perturb = 1, v_mnl%nat
        IF(v_mnl%ityp (na_perturb) == nt_perturb)THEN
            nkb_perturb = nkb_perturb + nh (nt_perturb)
        ENDIF
     ENDDO
  ENDDO
  
  
  CALL allocate_bec_type ( nkb_perturb, nbnd, becp1_perturb )
  CALL allocate_bec_type ( nkb_perturb, nbnd, becp2_perturb )
  ALLOCATE(vkb_perturb(npwx,nkb_perturb))
  CALL init_us_2_sc (ngk(ik0), igk_k(1,ik0), xk (1, ik0), vkb_perturb,v_mnl%nat,v_mnl%ityp,v_mnl%tau,nkb_perturb)
  CALL calbec ( ngk(ik0), vkb_perturb, evc1, becp1_perturb )
  
  CALL init_us_2_sc (ngk(ik), igk_k(1,ik), xk (1, ik), vkb_perturb,v_mnl%nat,v_mnl%ityp,v_mnl%tau,nkb_perturb)
  CALL calbec ( ngk(ik), vkb_perturb, evc2, becp2_perturb )
  
  
  ijkb0 = 0
  mnl=0
  DO nt_perturb = 1, v_mnl%ntyp
    DO na_perturb = 1, v_mnl%nat
       IF(v_mnl%ityp (na_perturb) == nt_perturb)THEN
           DO ih = 1, nh (nt_perturb)
              ikb = ijkb0 + ih
              mnl=mnl+conjg(becp1_perturb%nc(ikb,1,ibnd0))*becp2_perturb%nc(ikb,1,ibnd) * dvan(ih,ih,nt_perturb) &
                     +conjg(becp1_perturb%nc(ikb,2,ibnd0))*becp2_perturb%nc(ikb,2,ibnd) * dvan(ih,ih,nt_perturb)
              DO jh = ( ih + 1 ), nh(nt_perturb)
                 jkb = ijkb0 + jh
                 mnl=mnl + &
                     (conjg(becp1_perturb%nc(ikb, 1, ibnd0))*becp2_perturb%nc(jkb, 1, ibnd)+&
                      conjg(becp1_perturb%nc(jkb, 1, ibnd0))*becp2_perturb%nc(ikb, 1, ibnd))&
                        * dvan(ih,jh,nt_perturb) &
                     + &
                     (conjg(becp1_perturb%nc(ikb, 2, ibnd0))*becp2_perturb%nc(jkb, 2, ibnd)+&
                      conjg(becp1_perturb%nc(jkb, 2, ibnd0))*becp2_perturb%nc(ikb, 2, ibnd))&
                        * dvan(ih,jh,nt_perturb)
              ENDDO
  
           ENDDO
           ijkb0 = ijkb0 + nh (nt_perturb)
        ENDIF
     ENDDO
  ENDDO
   
  CALL deallocate_bec_type (  becp1_perturb )
  CALL deallocate_bec_type (  becp2_perturb )
  DEALLOCATE(vkb_perturb)

  write (stdout,1002) 'Mnl ki->kf ', ik0,ik, mnl, abs(mnl)
1001 format(A16,I9,I9, " ( ",e17.9," , ",e17.9," ) ",e17.9)
1002 format(A16,I9,I9," ( ",e17.9," , ",e17.9," ) ",e17.9/)
  mnonlocal=mnl
END SUBROUTINE calcmdefect_mnl_ks_noncolin


