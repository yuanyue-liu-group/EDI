MODULE bz_symmetry
  USE kinds, ONLY : dp
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC :: read_symmetry_info
  PUBLIC :: gen_qirr_2d, gen_qirr_3d
  PUBLIC :: get_kindex_round
  PUBLIC :: get_kfrac_2d, get_kfrac_3d

CONTAINS

  SUBROUTINE read_symmetry_info(filename, sym_mat, nsym, ndim)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    REAL(dp), ALLOCATABLE, INTENT(OUT) :: sym_mat(:,:,:)
    INTEGER, INTENT(OUT) :: nsym
    INTEGER, INTENT(IN) :: ndim

    INTEGER :: iunit, ios, i, j, tmp_nsym, ngroups, remainder
    CHARACTER(LEN=256) :: line_buf
    CHARACTER(LEN=256), ALLOCATABLE :: words(:)
    INTEGER, ALLOCATABLE :: sym_int(:,:,:)
    REAL(dp), ALLOCATABLE :: sym_full(:,:,:)
    INTEGER :: row_vals(18)

    iunit = 99
    OPEN(UNIT=iunit, FILE=filename, STATUS='old', IOSTAT=ios)
    IF (ios /= 0) THEN
       WRITE(*,*) 'Error opening symmetry file: ', TRIM(filename)
       nsym = 0
       RETURN
    ENDIF

    nsym = 0
    DO
       READ(iunit, '(A)', IOSTAT=ios) line_buf
       IF (ios /= 0) EXIT
       IF (INDEX(line_buf, 'symmetry') > 0) THEN
          READ(line_buf, *, IOSTAT=ios) nsym
          IF (ios /= 0) CYCLE
          IF (nsym > 0) EXIT
       ENDIF
    ENDDO

    IF (nsym == 0) THEN
       CLOSE(iunit)
       RETURN
    ENDIF

    ALLOCATE(sym_full(nsym, 3, 3))
    sym_full = 0.0_dp

    ngroups = nsym / 6
    remainder = MOD(nsym, 6)

    DO i = 1, ngroups
       READ(iunit, '(A)', IOSTAT=ios) line_buf
       READ(iunit, *) row_vals(1:18)
       DO j = 0, 5
          sym_full((i-1)*6+j+1, 1, 1:3) = DBLE(row_vals(j*3+1:j*3+3))
       ENDDO
       READ(iunit, *) row_vals(1:18)
       DO j = 0, 5
          sym_full((i-1)*6+j+1, 2, 1:3) = DBLE(row_vals(j*3+1:j*3+3))
       ENDDO
       READ(iunit, *) row_vals(1:18)
       DO j = 0, 5
          sym_full((i-1)*6+j+1, 3, 1:3) = DBLE(row_vals(j*3+1:j*3+3))
       ENDDO
    ENDDO

    IF (remainder > 0) THEN
       READ(iunit, '(A)', IOSTAT=ios) line_buf
       READ(iunit, *) row_vals(1:remainder*3)
       DO j = 0, remainder-1
          sym_full(ngroups*6+j+1, 1, 1:3) = DBLE(row_vals(j*3+1:j*3+3))
       ENDDO
       READ(iunit, *) row_vals(1:remainder*3)
       DO j = 0, remainder-1
          sym_full(ngroups*6+j+1, 2, 1:3) = DBLE(row_vals(j*3+1:j*3+3))
       ENDDO
       READ(iunit, *) row_vals(1:remainder*3)
       DO j = 0, remainder-1
          sym_full(ngroups*6+j+1, 3, 1:3) = DBLE(row_vals(j*3+1:j*3+3))
       ENDDO
    ENDIF

    CLOSE(iunit)

    IF (ndim == 2) THEN
       ALLOCATE(sym_mat(nsym, 2, 2))
       sym_mat(:, 1:2, 1:2) = sym_full(:, 1:2, 1:2)
    ELSE
       ALLOCATE(sym_mat(nsym, 3, 3))
       sym_mat = sym_full
    ENDIF

    DEALLOCATE(sym_full)
  END SUBROUTINE read_symmetry_info

  INTEGER FUNCTION get_kindex_round(kfx, kfy, nqf1, nqf2)
    REAL(dp), INTENT(IN) :: kfx, kfy
    INTEGER, INTENT(IN) :: nqf1, nqf2
    INTEGER :: ikx, iky
    ikx = MOD(NINT(kfx * nqf1), nqf1)
    IF (ikx < 0) ikx = ikx + nqf1
    iky = MOD(NINT(kfy * nqf2), nqf2)
    IF (iky < 0) iky = iky + nqf2
    get_kindex_round = iky + ikx * nqf2
  END FUNCTION get_kindex_round

  SUBROUTINE get_kfrac_2d(ik, nqf1, nqf2, kfx, kfy)
    INTEGER, INTENT(IN) :: ik, nqf1, nqf2
    REAL(dp), INTENT(OUT) :: kfx, kfy
    kfx = DBLE(ik / nqf2) / DBLE(nqf1)
    kfy = DBLE(MOD(ik, nqf2)) / DBLE(nqf2)
  END SUBROUTINE get_kfrac_2d

  SUBROUTINE get_kfrac_3d(ik, nk1, nk2, nk3, kfx, kfy, kfz)
    INTEGER, INTENT(IN) :: ik, nk1, nk2, nk3
    REAL(dp), INTENT(OUT) :: kfx, kfy, kfz
    INTEGER :: nk23
    nk23 = nk2 * nk3
    kfx = DBLE(ik / nk23) / DBLE(nk1)
    kfy = DBLE(MOD(ik, nk23) / nk3) / DBLE(nk2)
    kfz = DBLE(MOD(ik, nk3)) / DBLE(nk3)
  END SUBROUTINE get_kfrac_3d

  SUBROUTINE gen_qirr_2d(nqf1, nqf2, sym_mat, nsym, ibz2bz, nibz, bz2ibz)
    INTEGER, INTENT(IN) :: nqf1, nqf2, nsym
    REAL(dp), INTENT(IN) :: sym_mat(nsym, 2, 2)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: ibz2bz(:)
    INTEGER, INTENT(OUT) :: nibz
    INTEGER, INTENT(OUT) :: bz2ibz(0:nqf1*nqf2-1)

    INTEGER :: iq, isym, nqtot, ikrot
    REAL(dp) :: qx_bz(2), qrot(2)
    INTEGER :: bz_done(0:nqf1*nqf2-1)
    INTEGER :: ibz_tmp(nqf1*nqf2)

    nqtot = nqf1 * nqf2
    bz_done = -1
    bz2ibz = -1
    nibz = 0

    DO iq = 0, nqtot - 1
       IF (bz_done(iq) >= 0) CYCLE

       nibz = nibz + 1
       ibz_tmp(nibz) = iq

       CALL get_kfrac_2d(iq, nqf1, nqf2, qx_bz(1), qx_bz(2))

       DO isym = 1, nsym
          qrot(1) = sym_mat(isym, 1, 1) * qx_bz(1) + sym_mat(isym, 1, 2) * qx_bz(2)
          qrot(2) = sym_mat(isym, 2, 1) * qx_bz(1) + sym_mat(isym, 2, 2) * qx_bz(2)
          ikrot = get_kindex_round(qrot(1), qrot(2), nqf1, nqf2)
          IF (bz_done(ikrot) < 0) THEN
             bz_done(ikrot) = isym
             bz2ibz(ikrot) = nibz - 1
          ENDIF
       ENDDO
    ENDDO

    ALLOCATE(ibz2bz(nibz))
    ibz2bz(1:nibz) = ibz_tmp(1:nibz)
  END SUBROUTINE gen_qirr_2d

  SUBROUTINE gen_qirr_3d(nk1, nk2, nk3, sym_mat, nsym, ibz2bz, nibz, bz2ibz)
    INTEGER, INTENT(IN) :: nk1, nk2, nk3, nsym
    REAL(dp), INTENT(IN) :: sym_mat(nsym, 3, 3)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: ibz2bz(:)
    INTEGER, INTENT(OUT) :: nibz
    INTEGER, INTENT(OUT) :: bz2ibz(0:nk1*nk2*nk3-1)

    INTEGER :: iq, isym, nqtot, ikx, iky, ikz, ikrot
    REAL(dp) :: qx_bz(3), qrot(3)
    INTEGER :: bz_done(0:nk1*nk2*nk3-1)
    INTEGER :: ibz_tmp(nk1*nk2*nk3)

    nqtot = nk1 * nk2 * nk3
    bz_done = -1
    bz2ibz = -1
    nibz = 0

    DO iq = 0, nqtot - 1
       IF (bz_done(iq) >= 0) CYCLE

       nibz = nibz + 1
       ibz_tmp(nibz) = iq

       CALL get_kfrac_3d(iq, nk1, nk2, nk3, qx_bz(1), qx_bz(2), qx_bz(3))

       DO isym = 1, nsym
          qrot(1) = sym_mat(isym, 1, 1)*qx_bz(1) + sym_mat(isym, 1, 2)*qx_bz(2) + sym_mat(isym, 1, 3)*qx_bz(3)
          qrot(2) = sym_mat(isym, 2, 1)*qx_bz(1) + sym_mat(isym, 2, 2)*qx_bz(2) + sym_mat(isym, 2, 3)*qx_bz(3)
          qrot(3) = sym_mat(isym, 3, 1)*qx_bz(1) + sym_mat(isym, 3, 2)*qx_bz(2) + sym_mat(isym, 3, 3)*qx_bz(3)

          ikx = MOD(NINT(qrot(1)*nk1), nk1); IF (ikx < 0) ikx = ikx + nk1
          iky = MOD(NINT(qrot(2)*nk2), nk2); IF (iky < 0) iky = iky + nk2
          ikz = MOD(NINT(qrot(3)*nk3), nk3); IF (ikz < 0) ikz = ikz + nk3
          ikrot = ikx * nk2 * nk3 + iky * nk3 + ikz

          IF (bz_done(ikrot) < 0) THEN
             bz_done(ikrot) = isym
             bz2ibz(ikrot) = nibz - 1
          ENDIF
       ENDDO
    ENDDO

    ALLOCATE(ibz2bz(nibz))
    ibz2bz(1:nibz) = ibz_tmp(1:nibz)
  END SUBROUTINE gen_qirr_3d

END MODULE bz_symmetry
