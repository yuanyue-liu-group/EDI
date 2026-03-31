Program test_potential
  !-----------------------------------------------------------------------
  ! Diagnostic: read a QE save directory and print the KS potential
  ! components separately: V_H+V_xc (v%of_r), V_loc_ionic (vltot),
  ! and their sum. Also writes a cube file for comparison with pp.x.
  !
  ! This runs the SAME code path as get_vloc_onthefly:
  !   read_file_new → post_xml_init → setlocal + v_of_rho
  !
  ! Input: &test_pot_nml  prefix_in, outdir_in /
  ! Usage: srun -n 1 edi_test_pot.x -i test_pot.in
  !-----------------------------------------------------------------------
  USE kinds, ONLY : dp
  USE io_files, ONLY : prefix, tmp_dir
  USE scf, ONLY : v, vltot
  USE fft_base, ONLY : dfftp, dffts
  USE ions_base, ONLY : nat, ntyp => nsp, ityp, tau, atm
  USE cell_base, ONLY : at, alat, omega
  USE environment, ONLY : environment_start, environment_end
  USE mp_global, ONLY : mp_startup, mp_global_end
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp, ONLY : mp_bcast
  USE mp_world, ONLY : world_comm

  IMPLICIT NONE

  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  INTEGER :: ios, ir, irx, iry, irz, inr, nval, ia, iunit
  REAL(dp) :: v_hxc_min, v_hxc_max, v_hxc_avg
  REAL(dp) :: vlt_min, vlt_max, vlt_avg
  REAL(dp) :: vtot_min, vtot_max, vtot_avg
  REAL(dp) :: voxel(3)
  REAL(dp), ALLOCATABLE :: vtot(:), v_hxc(:), v_ionic(:)
  REAL(dp), ALLOCATABLE :: z_avg_hxc(:), z_avg_ionic(:), z_avg_tot(:)
  CHARACTER(LEN=256) :: prefix_in, outdir_in, fname
  LOGICAL :: wfc_is_collected

  NAMELIST /test_pot_nml/ prefix_in, outdir_in

  CALL mp_startup(start_images=.TRUE.)
  CALL environment_start('TEST_POT')

  prefix_in = 'pwscf'
  outdir_in = './'

  IF (ionode) THEN
     CALL input_from_file()
     READ(5, NML=test_pot_nml, IOSTAT=ios)
     IF (ios < 0) ios = 0
  ENDIF
  CALL mp_bcast(ios, ionode_id, world_comm)
  IF (ios > 0) CALL errore('test_potential', 'error reading namelist', ios)
  CALL mp_bcast(prefix_in, ionode_id, world_comm)
  CALL mp_bcast(outdir_in, ionode_id, world_comm)

  prefix = TRIM(prefix_in)
  tmp_dir = trimcheck(TRIM(outdir_in))

  ! Read QE data — same as get_vloc_onthefly
  wfc_is_collected = .FALSE.
  CALL read_file_new(wfc_is_collected)

  IF (ionode) THEN
     WRITE(stdout, '(/,5X,A)') REPEAT('=', 60)
     WRITE(stdout, '(5X,A)')   'Potential diagnostic'
     WRITE(stdout, '(5X,A)')   REPEAT('=', 60)
     WRITE(stdout, '(5X,A,3I6)') 'dfftp grid: ', dfftp%nr1, dfftp%nr2, dfftp%nr3
     WRITE(stdout, '(5X,A,3I6)') 'dfftp padded: ', dfftp%nr1x, dfftp%nr2x, dfftp%nr3x
     WRITE(stdout, '(5X,A,I10)') 'dfftp%nnr = ', dfftp%nnr
     WRITE(stdout, '(5X,A,3I6)') 'dffts grid: ', dffts%nr1, dffts%nr2, dffts%nr3
     WRITE(stdout, '(5X,A,3I6)') 'dffts padded: ', dffts%nr1x, dffts%nr2x, dffts%nr3x
     WRITE(stdout, '(5X,A,I10)') 'dffts%nnr = ', dffts%nnr
     WRITE(stdout, '(5X,A,I6)')  'nat = ', nat
     WRITE(stdout, '(5X,A,F12.6)') 'alat = ', alat
     WRITE(stdout, '(5X,A,F12.6)') 'omega = ', omega
     FLUSH(stdout)

     ! Extract components
     ALLOCATE(v_hxc(dfftp%nnr))
     ALLOCATE(v_ionic(dfftp%nnr))
     ALLOCATE(vtot(dfftp%nnr))

     v_hxc(:) = v%of_r(:, 1)        ! V_Hartree + V_xc
     v_ionic(:) = vltot(:)           ! V_loc_ionic (local pseudopotential)
     vtot(:) = v_hxc(:) + v_ionic(:) ! total = what EDI uses

     ! Statistics (Ry)
     WRITE(stdout, '(/,5X,A)') 'Potential statistics (Ry):'
     WRITE(stdout, '(5X,A,3ES14.6)') '  V_Hxc   min/max/avg: ', &
          MINVAL(v_hxc(1:dfftp%nnr)), MAXVAL(v_hxc(1:dfftp%nnr)), &
          SUM(v_hxc(1:dfftp%nnr))/dfftp%nnr
     WRITE(stdout, '(5X,A,3ES14.6)') '  V_ionic min/max/avg: ', &
          MINVAL(v_ionic(1:dfftp%nnr)), MAXVAL(v_ionic(1:dfftp%nnr)), &
          SUM(v_ionic(1:dfftp%nnr))/dfftp%nnr
     WRITE(stdout, '(5X,A,3ES14.6)') '  V_total min/max/avg: ', &
          MINVAL(vtot(1:dfftp%nnr)), MAXVAL(vtot(1:dfftp%nnr)), &
          SUM(vtot(1:dfftp%nnr))/dfftp%nnr

     ! Z-averaged profiles
     ALLOCATE(z_avg_hxc(dfftp%nr3))
     ALLOCATE(z_avg_ionic(dfftp%nr3))
     ALLOCATE(z_avg_tot(dfftp%nr3))
     z_avg_hxc = 0.0_dp
     z_avg_ionic = 0.0_dp
     z_avg_tot = 0.0_dp

     DO irz = 0, dfftp%nr3 - 1
        DO iry = 0, dfftp%nr2 - 1
           DO irx = 0, dfftp%nr1 - 1
              ! QE linear index: x-fastest
              inr = irx + iry * dfftp%nr1x + irz * dfftp%nr1x * dfftp%nr2x + 1
              z_avg_hxc(irz+1) = z_avg_hxc(irz+1) + v_hxc(inr)
              z_avg_ionic(irz+1) = z_avg_ionic(irz+1) + v_ionic(inr)
              z_avg_tot(irz+1) = z_avg_tot(irz+1) + vtot(inr)
           ENDDO
        ENDDO
     ENDDO
     z_avg_hxc = z_avg_hxc / DBLE(dfftp%nr1 * dfftp%nr2)
     z_avg_ionic = z_avg_ionic / DBLE(dfftp%nr1 * dfftp%nr2)
     z_avg_tot = z_avg_tot / DBLE(dfftp%nr1 * dfftp%nr2)

     ! Write z-average data
     fname = TRIM(prefix_in) // '_pot_zavg.dat'
     iunit = 99
     OPEN(iunit, FILE=TRIM(fname), STATUS='REPLACE')
     WRITE(iunit, '(A)') '# z(bohr)  V_Hxc(Ry)  V_ionic(Ry)  V_total(Ry)'
     DO irz = 1, dfftp%nr3
        WRITE(iunit, '(4ES16.8)') &
             DBLE(irz-1) / DBLE(dfftp%nr3) * at(3,3) * alat, &
             z_avg_hxc(irz), z_avg_ionic(irz), z_avg_tot(irz)
     ENDDO
     CLOSE(iunit)
     WRITE(stdout, '(5X,A,A)') 'Z-average written: ', TRIM(fname)

     ! Write total potential cube file (same method as pp.x plot_num=1)
     fname = TRIM(prefix_in) // '_vtot_test.cube'
     iunit = 98
     OPEN(iunit, FILE=TRIM(fname), STATUS='REPLACE')
     WRITE(iunit, '(A)') 'V_total = V_Hxc + V_ionic [Ry] (test_potential)'
     WRITE(iunit, '(A)') 'Same as pp.x plot_num=1'
     WRITE(iunit, '(I5, 3F12.6)') nat, 0.0_dp, 0.0_dp, 0.0_dp
     voxel(:) = at(:, 1) * alat / DBLE(dfftp%nr1)
     WRITE(iunit, '(I5, 3F12.6)') dfftp%nr1, voxel
     voxel(:) = at(:, 2) * alat / DBLE(dfftp%nr2)
     WRITE(iunit, '(I5, 3F12.6)') dfftp%nr2, voxel
     voxel(:) = at(:, 3) * alat / DBLE(dfftp%nr3)
     WRITE(iunit, '(I5, 3F12.6)') dfftp%nr3, voxel
     DO ia = 1, nat
        WRITE(iunit, '(I5, 4F12.6)') ityp(ia), 0.0_dp, &
             tau(1, ia) * alat, tau(2, ia) * alat, tau(3, ia) * alat
     ENDDO
     ! Cube format: x outermost, z innermost
     DO irx = 0, dfftp%nr1 - 1
        DO iry = 0, dfftp%nr2 - 1
           nval = 0
           DO irz = 0, dfftp%nr3 - 1
              inr = irx + iry * dfftp%nr1x + irz * dfftp%nr1x * dfftp%nr2x + 1
              WRITE(iunit, '(ES13.5)', ADVANCE='NO') vtot(inr)
              nval = nval + 1
              IF (MOD(nval, 6) == 0) WRITE(iunit, *)
           ENDDO
           IF (MOD(nval, 6) /= 0) WRITE(iunit, *)
        ENDDO
     ENDDO
     CLOSE(iunit)
     WRITE(stdout, '(5X,A,A)') 'Cube file written: ', TRIM(fname)

     ! Also write V_ionic only cube for comparison
     fname = TRIM(prefix_in) // '_vionic_test.cube'
     OPEN(iunit, FILE=TRIM(fname), STATUS='REPLACE')
     WRITE(iunit, '(A)') 'V_ionic = vltot [Ry] (local pseudopotential only)'
     WRITE(iunit, '(A)') 'Generated by test_potential'
     WRITE(iunit, '(I5, 3F12.6)') nat, 0.0_dp, 0.0_dp, 0.0_dp
     voxel(:) = at(:, 1) * alat / DBLE(dfftp%nr1)
     WRITE(iunit, '(I5, 3F12.6)') dfftp%nr1, voxel
     voxel(:) = at(:, 2) * alat / DBLE(dfftp%nr2)
     WRITE(iunit, '(I5, 3F12.6)') dfftp%nr2, voxel
     voxel(:) = at(:, 3) * alat / DBLE(dfftp%nr3)
     WRITE(iunit, '(I5, 3F12.6)') dfftp%nr3, voxel
     DO ia = 1, nat
        WRITE(iunit, '(I5, 4F12.6)') ityp(ia), 0.0_dp, &
             tau(1, ia) * alat, tau(2, ia) * alat, tau(3, ia) * alat
     ENDDO
     DO irx = 0, dfftp%nr1 - 1
        DO iry = 0, dfftp%nr2 - 1
           nval = 0
           DO irz = 0, dfftp%nr3 - 1
              inr = irx + iry * dfftp%nr1x + irz * dfftp%nr1x * dfftp%nr2x + 1
              WRITE(iunit, '(ES13.5)', ADVANCE='NO') v_ionic(inr)
              nval = nval + 1
              IF (MOD(nval, 6) == 0) WRITE(iunit, *)
           ENDDO
           IF (MOD(nval, 6) /= 0) WRITE(iunit, *)
        ENDDO
     ENDDO
     CLOSE(iunit)
     WRITE(stdout, '(5X,A,A)') 'Cube file written: ', TRIM(fname)

     WRITE(stdout, '(5X,A)') REPEAT('=', 60)

     DEALLOCATE(v_hxc, v_ionic, vtot)
     DEALLOCATE(z_avg_hxc, z_avg_ionic, z_avg_tot)
  ENDIF

  CALL environment_end('TEST_POT')
  CALL mp_global_end()

End Program test_potential
