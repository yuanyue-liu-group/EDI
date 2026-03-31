Program test_sum_rule
  !-----------------------------------------------------------------------
  ! Verify the energy sum rule:
  !   T(n,k) + V_loc(n,k) + V_NL(n,k) = ε(n,k)
  !
  ! Uses the SAME code path as EDI matrix element computation:
  !   - invfft for ψ → u(r) on primitive cell grid
  !   - (1/N_r) Σ_r  u*(r) V_KS(r) u(r)  for V_loc  [same as M_local]
  !   - calbec + becp* D becp                for V_NL  [same as M_nonlocal]
  !
  ! If the sum rule holds, the local and nonlocal matrix element code
  ! is unitarily correct and the matrix elements are in Ry.
  !
  ! Input: &test_nml  prefix_in, outdir_in /
  !   (same as QE NSCF output directory)
  !
  ! Usage: srun -n N edi_test_sr.x -nk N -i test.in
  !-----------------------------------------------------------------------
  USE kinds, ONLY : dp
  USE io_files, ONLY : prefix, tmp_dir, restart_dir
  USE wvfct, ONLY : npwx, nbnd, et
  USE klist, ONLY : nkstot, nks, xk, igk_k, ngk
  USE noncollin_module, ONLY : npol, noncolin, lspinorb
  USE fft_base, ONLY : dffts
  USE fft_interfaces, ONLY : invfft
  USE pw_restart_new, ONLY : read_collected_wfc
  USE scf, ONLY : v, vltot
  USE uspp_param, ONLY : nh
  USE uspp, ONLY : nkb, vkb, dvan, dvan_so
  USE becmod, ONLY : bec_type, calbec, allocate_bec_type, deallocate_bec_type
  USE ions_base, ONLY : nat, ntyp => nsp, ityp
  USE gvect, ONLY : g, ngm, igtongl
  USE gvect, ONLY : eigts1, eigts2, eigts3, mill
  USE cell_base, ONLY : tpiba2, omega
  USE environment, ONLY : environment_start, environment_end
  USE mp_global, ONLY : mp_startup, mp_global_end
  USE mp_bands, ONLY : nbgrp
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp, ONLY : mp_bcast
  USE mp_world, ONLY : world_comm
  USE constants, ONLY : rytoev, tpi
  USE lsda_mod, ONLY : nspin
  USE uspp_init, ONLY : init_us_2

  IMPLICIT NONE

  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  INTEGER :: ik, ibnd, ig, ir, ipol, ios
  INTEGER :: na, nt, ih, jh, ikb, jkb, ijkb0
  REAL(dp) :: T_nk, V_loc_nk, V_NL_nk, eps_nk, err, sum_nk
  REAL(dp) :: max_err, avg_err, gk2
  INTEGER :: ntot
  COMPLEX(dp), ALLOCATABLE :: evc(:,:), psic_tmp(:)
  REAL(dp), ALLOCATABLE :: V_KS(:)
  TYPE(bec_type) :: becp_work
  LOGICAL :: need_wf
  CHARACTER(LEN=256) :: prefix_in, outdir_in

  NAMELIST /test_nml/ prefix_in, outdir_in

  CALL mp_startup(start_images=.TRUE.)
  CALL environment_start('TEST_SUM_RULE')

  IF (nbgrp > 1) CALL errore('test_sum_rule', 'band groups not supported', nbgrp)

  prefix_in = 'pwscf'
  outdir_in = './'

  IF (ionode) THEN
     CALL input_from_file()
     READ(5, NML=test_nml, IOSTAT=ios)
     IF (ios < 0) ios = 0
  ENDIF
  CALL mp_bcast(ios, ionode_id, world_comm)
  IF (ios > 0) CALL errore('test_sum_rule', 'error reading test_nml', ios)
  CALL mp_bcast(prefix_in, ionode_id, world_comm)
  CALL mp_bcast(outdir_in, ionode_id, world_comm)

  prefix = TRIM(prefix_in)
  tmp_dir = trimcheck(TRIM(outdir_in))

  ! Read primitive cell data (wfc + charge → potential reconstruction)
  need_wf = .TRUE.
  CALL read_file_new(need_wf)

  IF (ionode) THEN
     WRITE(stdout, '(/,5X,A)') REPEAT('=', 60)
     WRITE(stdout, '(5X,A)')   'Energy sum rule test'
     WRITE(stdout, '(5X,A)')   'T(n,k) + V_loc(n,k) + V_NL(n,k) = eps(n,k)'
     WRITE(stdout, '(5X,A)')   REPEAT('=', 60)
     WRITE(stdout, '(5X,A,I6)') 'nkstot = ', nkstot
     WRITE(stdout, '(5X,A,I6)') 'nbnd   = ', nbnd
     WRITE(stdout, '(5X,A,I6)') 'npol   = ', npol
     WRITE(stdout, '(5X,A,I6)') 'nat    = ', nat
     WRITE(stdout, '(5X,A,I6)') 'nkb    = ', nkb
     FLUSH(stdout)
  ENDIF

  ! Set up local pseudopotential on real-space grid: vltot = Vloc_pseudo + Ewald
  CALL setlocal()

  ! Construct total KS local potential: V_KS(r) = V_H(r) + V_xc(r) + V_loc_pseudo(r)
  ! v%of_r(:,1) = Hartree + XC (from v_of_rho, called by read_file_new)
  ! vltot(:) = local pseudopotential + Ewald (from setlocal)
  ALLOCATE(V_KS(dffts%nnr))
  V_KS(:) = v%of_r(1:dffts%nnr, 1) + vltot(1:dffts%nnr)
  IF (ionode) THEN
     WRITE(stdout, '(5X,A,2ES14.6)') 'V_KS range: ', MINVAL(V_KS), MAXVAL(V_KS)
     FLUSH(stdout)
  ENDIF

  ! Allocate workspace
  ALLOCATE(evc(npwx * npol, nbnd))
  ALLOCATE(psic_tmp(dffts%nnr))

  max_err = 0.0_dp
  avg_err = 0.0_dp
  ntot = 0

  IF (ionode) THEN
     WRITE(stdout, '(/,5X,A)') &
          '  ik  ib     T(eV)      V_loc(eV)    V_NL(eV)     Sum(eV)      eps(eV)      err(meV)'
     WRITE(stdout, '(5X,A)') REPEAT('-', 90)
  ENDIF

  DO ik = 1, nks
     ! Read wavefunctions
     CALL read_collected_wfc(restart_dir(), ik, evc)

     ! Set up beta projectors for this k-point (QE standard routine)
     IF (nkb > 0) CALL init_us_2(ngk(ik), igk_k(1,ik), xk(1,ik), vkb)
     CALL allocate_bec_type(nkb, nbnd, becp_work)
     IF (nkb > 0) CALL calbec(ngk(ik), vkb, evc, becp_work)

     DO ibnd = 1, nbnd
        ! ============================================================
        ! 1. KINETIC ENERGY: T = Σ_{σ,G} |c_{σ}(G)|² × |k+G|²
        ! ============================================================
        T_nk = 0.0_dp
        DO ipol = 1, npol
           DO ig = 1, ngk(ik)
              gk2 = ( (xk(1,ik) + g(1, igk_k(ig,ik)))**2 + &
                      (xk(2,ik) + g(2, igk_k(ig,ik)))**2 + &
                      (xk(3,ik) + g(3, igk_k(ig,ik)))**2 ) * tpiba2
              T_nk = T_nk + gk2 * ABS(evc(ig + (ipol-1)*npwx, ibnd))**2
           ENDDO
        ENDDO

        ! ============================================================
        ! 2. LOCAL POTENTIAL: (1/N_r) Σ_{σ,r} |u_σ(r)|² × V_KS(r)
        !    [Same code path as EDI M_local with V_folded = V_KS]
        ! ============================================================
        V_loc_nk = 0.0_dp
        DO ipol = 1, npol
           psic_tmp(:) = (0.0_dp, 0.0_dp)
           DO ig = 1, ngk(ik)
              psic_tmp(dffts%nl(igk_k(ig, ik))) = evc(ig + (ipol-1)*npwx, ibnd)
           ENDDO
           CALL invfft('Wave', psic_tmp, dffts)
           DO ir = 1, dffts%nnr
              V_loc_nk = V_loc_nk + DBLE(CONJG(psic_tmp(ir)) * psic_tmp(ir)) * V_KS(ir)
           ENDDO
        ENDDO
        V_loc_nk = V_loc_nk / DBLE(dffts%nnr)

        ! ============================================================
        ! 3. NONLOCAL: Σ_{a,i,j,σ,σ'} becp*(i,σ,n) D(i,j,σσ') becp(j,σ',n)
        !    [Same code path as EDI M_nonlocal]
        ! ============================================================
        V_NL_nk = 0.0_dp
        ijkb0 = 0
        DO nt = 1, ntyp
           DO na = 1, nat
              IF (ityp(na) == nt) THEN
                 IF (lspinorb) THEN
                    ! SOC: 4 spin channels
                    DO ih = 1, nh(nt)
                       ikb = ijkb0 + ih
                       DO jh = 1, nh(nt)
                          jkb = ijkb0 + jh
                          V_NL_nk = V_NL_nk + DBLE( &
                               CONJG(becp_work%nc(ikb,1,ibnd)) * becp_work%nc(jkb,1,ibnd) * dvan_so(ih,jh,1,nt) + &
                               CONJG(becp_work%nc(ikb,1,ibnd)) * becp_work%nc(jkb,2,ibnd) * dvan_so(ih,jh,2,nt) + &
                               CONJG(becp_work%nc(ikb,2,ibnd)) * becp_work%nc(jkb,1,ibnd) * dvan_so(ih,jh,3,nt) + &
                               CONJG(becp_work%nc(ikb,2,ibnd)) * becp_work%nc(jkb,2,ibnd) * dvan_so(ih,jh,4,nt) )
                       ENDDO
                    ENDDO
                 ELSE
                    ! Scalar-relativistic
                    DO ih = 1, nh(nt)
                       ikb = ijkb0 + ih
                       DO jh = 1, nh(nt)
                          jkb = ijkb0 + jh
                          V_NL_nk = V_NL_nk + DBLE( &
                               CONJG(becp_work%k(ikb, ibnd)) * becp_work%k(jkb, ibnd) ) * dvan(ih, jh, nt)
                       ENDDO
                    ENDDO
                 ENDIF
                 ijkb0 = ijkb0 + nh(nt)
              ENDIF
           ENDDO
        ENDDO

        ! ============================================================
        ! 4. Compare with eigenvalue
        ! ============================================================
        eps_nk = et(ibnd, ik)
        sum_nk = T_nk + V_loc_nk + V_NL_nk
        err = ABS(sum_nk - eps_nk)
        max_err = MAX(max_err, err)
        avg_err = avg_err + err
        ntot = ntot + 1

        IF (ionode) THEN
           WRITE(stdout, '(5X,I4,I4,6F13.6)') &
                ik, ibnd, T_nk*rytoev, V_loc_nk*rytoev, V_NL_nk*rytoev, &
                sum_nk*rytoev, eps_nk*rytoev, err*rytoev*1000.0_dp
        ENDIF
     ENDDO

     CALL deallocate_bec_type(becp_work)
  ENDDO

  IF (ionode) THEN
     WRITE(stdout, '(/,5X,A)') REPEAT('=', 60)
     WRITE(stdout, '(5X,A,ES12.4,A)') 'Max error:     ', max_err * rytoev * 1000.0_dp, ' meV'
     WRITE(stdout, '(5X,A,ES12.4,A)') 'Average error: ', (avg_err/ntot) * rytoev * 1000.0_dp, ' meV'
     IF (max_err * rytoev * 1000.0_dp < 1.0_dp) THEN
        WRITE(stdout, '(5X,A)') 'PASSED: Sum rule verified (< 1 meV)'
     ELSE
        WRITE(stdout, '(5X,A)') 'WARNING: Sum rule error > 1 meV'
     ENDIF
     WRITE(stdout, '(5X,A)') REPEAT('=', 60)
  ENDIF

  DEALLOCATE(evc, psic_tmp, V_KS)

  CALL environment_end('TEST_SUM_RULE')
  CALL mp_global_end()

End Program test_sum_rule
