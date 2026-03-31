MODULE wannierize_edi
  USE kinds, ONLY : dp
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC :: load_umat_ukk
  PUBLIC :: load_umat_w90

CONTAINS

  SUBROUTINE load_umat_ukk(filukk, nbndep, nbndsub, nkstot, cu, lwin)
    USE io_global, ONLY : ionode, ionode_id
    USE mp, ONLY : mp_bcast
    USE mp_images, ONLY : intra_image_comm
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: filukk
    INTEGER, INTENT(IN) :: nbndep, nbndsub, nkstot
    COMPLEX(dp), INTENT(OUT) :: cu(nbndep, nbndsub, nkstot)
    LOGICAL, INTENT(OUT) :: lwin(nbndep, nkstot)

    INTEGER :: ik, ibnd, jbnd, ios, iunit, dum1, dum2
    iunit = 98

    cu = (0.0_dp, 0.0_dp)
    lwin = .TRUE.

    IF (ionode) THEN
       OPEN(iunit, FILE=TRIM(filukk), STATUS='old', FORM='formatted', IOSTAT=ios)
       IF (ios /= 0) CALL errore('load_umat_ukk', 'error opening ukk file: '//TRIM(filukk), 1)

       READ(iunit, *) dum1, dum2

       DO ibnd = 1, nbndep
          READ(iunit, *) dum1
       ENDDO

       DO ik = 1, nkstot
          DO ibnd = 1, nbndep
             DO jbnd = 1, nbndsub
                READ(iunit, *) cu(ibnd, jbnd, ik)
             ENDDO
          ENDDO
       ENDDO

       DO ik = 1, nkstot
          DO ibnd = 1, nbndep
             READ(iunit, *) lwin(ibnd, ik)
          ENDDO
       ENDDO

       CLOSE(iunit)
    ENDIF

    CALL mp_bcast(cu, ionode_id, intra_image_comm)
    CALL mp_bcast(lwin, ionode_id, intra_image_comm)
  END SUBROUTINE load_umat_ukk

  SUBROUTINE load_umat_w90(seedname, nbndsub, nkstot, cu)
    USE io_global, ONLY : ionode, ionode_id
    USE mp, ONLY : mp_bcast
    USE mp_images, ONLY : intra_image_comm
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: seedname
    INTEGER, INTENT(IN) :: nbndsub, nkstot
    COMPLEX(dp), INTENT(OUT) :: cu(nbndsub, nbndsub, nkstot)

    CHARACTER(LEN=512) :: fname
    INTEGER :: ik, ibnd, jbnd, ios, iunit, nk_file, nb_file
    INTEGER :: dum_ik
    REAL(dp) :: re_part, im_part
    REAL(dp) :: xk_dum(3)

    iunit = 98
    cu = (0.0_dp, 0.0_dp)

    IF (ionode) THEN
       fname = TRIM(seedname) // '_u.mat'
       OPEN(iunit, FILE=TRIM(fname), STATUS='old', FORM='formatted', IOSTAT=ios)
       IF (ios /= 0) CALL errore('load_umat_w90', 'error opening ' // TRIM(fname), 1)

       READ(iunit, '(A)') fname
       READ(iunit, *) nk_file, nb_file, dum_ik
       IF (nk_file /= nkstot) &
          CALL errore('load_umat_w90', 'nk mismatch in _u.mat', nk_file)
       IF (nb_file /= nbndsub) &
          CALL errore('load_umat_w90', 'nbnd mismatch in _u.mat', nb_file)

       DO ik = 1, nkstot
          READ(iunit, *)
          READ(iunit, *) xk_dum
          DO jbnd = 1, nbndsub
             DO ibnd = 1, nbndsub
                READ(iunit, *) re_part, im_part
                cu(ibnd, jbnd, ik) = CMPLX(re_part, im_part, dp)
             ENDDO
          ENDDO
       ENDDO

       CLOSE(iunit)
    ENDIF

    CALL mp_bcast(cu, ionode_id, intra_image_comm)
  END SUBROUTINE load_umat_w90

END MODULE wannierize_edi
