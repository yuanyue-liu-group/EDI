MODULE delta_weights
  USE kinds, ONLY : dp
  IMPLICIT NONE
  PRIVATE
  SAVE

  PUBLIC :: get_triangular_weight
  PUBLIC :: get_qweight_2d
  PUBLIC :: get_qweight_wrap_2d
  PUBLIC :: kindex_add_2d
  PUBLIC :: kindex_sub_2d
  PUBLIC :: kindex_add_3d
  PUBLIC :: get_tetrahedral_weight
  PUBLIC :: get_qweight_3d
  PUBLIC :: get_delta_weight

  REAL(dp), PARAMETER :: pi_val = 3.14159265358979323846_dp
  REAL(dp), PARAMETER :: kB_eV = 8.621738d-5

CONTAINS

  INTEGER FUNCTION kindex_add_2d(ik, iq, nqf1, nqf2)
    INTEGER, INTENT(IN) :: ik, iq, nqf1, nqf2
    INTEGER :: ikqx, ikqy
    ikqx = MOD(ik/nqf2 + iq/nqf2, nqf1)
    ikqy = MOD(MOD(ik, nqf2) + MOD(iq, nqf2), nqf2)
    kindex_add_2d = ikqx * nqf2 + ikqy
  END FUNCTION kindex_add_2d

  INTEGER FUNCTION kindex_sub_2d(ikf, iki, nqf1, nqf2)
    ! Compute q = kf - ki on the 2D grid (0-indexed), with periodic wrapping
    INTEGER, INTENT(IN) :: ikf, iki, nqf1, nqf2
    INTEGER :: iqx, iqy
    iqx = MOD(ikf/nqf2 - iki/nqf2 + nqf1, nqf1)
    iqy = MOD(MOD(ikf, nqf2) - MOD(iki, nqf2) + nqf2, nqf2)
    kindex_sub_2d = iqx * nqf2 + iqy
  END FUNCTION kindex_sub_2d

  INTEGER FUNCTION kindex_add_3d(ik, iq, nk1, nk2, nk3)
    INTEGER, INTENT(IN) :: ik, iq, nk1, nk2, nk3
    INTEGER :: ix_k, iy_k, iz_k, ix_q, iy_q, iz_q
    INTEGER :: nk23
    nk23 = nk2 * nk3
    ix_k = ik / nk23
    iy_k = MOD(ik, nk23) / nk3
    iz_k = MOD(ik, nk3)
    ix_q = iq / nk23
    iy_q = MOD(iq, nk23) / nk3
    iz_q = MOD(iq, nk3)
    kindex_add_3d = MOD(ix_k + ix_q, nk1) * nk23 + &
                    MOD(iy_k + iy_q, nk2) * nk3 + &
                    MOD(iz_k + iz_q, nk3)
  END FUNCTION kindex_add_3d

  REAL(dp) FUNCTION get_triangular_weight(ea, eb, ec)
    REAL(dp), INTENT(IN) :: ea, eb, ec
    REAL(dp) :: E0, E1, E2, w_add, sorted(3)
    INTEGER :: idx(3)

    IF (ea >= 0.0_dp .AND. eb >= 0.0_dp .AND. ec >= 0.0_dp) THEN
       get_triangular_weight = 0.0_dp
       RETURN
    ENDIF
    IF (ea <= 0.0_dp .AND. eb <= 0.0_dp .AND. ec <= 0.0_dp) THEN
       get_triangular_weight = 0.0_dp
       RETURN
    ENDIF

    sorted(1) = ea; sorted(2) = eb; sorted(3) = ec
    CALL sort3(sorted)
    E0 = sorted(1); E1 = sorted(2); E2 = sorted(3)

    IF (ea < eb .AND. eb < ec) THEN
       IF (eb > 0.0_dp) THEN
          w_add = -1.0_dp * E0 / (E1 - E0) / (E2 - E0)
          w_add = w_add * (2.0_dp + E0/(E1-E0) + E0/(E2-E0))
       ELSE
          w_add = E2*E2 / (E2-E1) / (E2-E0) / (E2-E0)
       ENDIF
    ELSE IF (ea < ec .AND. ec < eb) THEN
       IF (ec > 0.0_dp) THEN
          w_add = -1.0_dp * E0 / (E1 - E0) / (E2 - E0)
          w_add = w_add * (2.0_dp + E0/(E1-E0) + E0/(E2-E0))
       ELSE
          w_add = E2*E2 / (E2-E1) / (E2-E0) / (E2-E0)
       ENDIF
    ELSE IF (eb < ea .AND. ea < ec) THEN
       IF (ea > 0.0_dp) THEN
          w_add = E0*E0 / (E1-E0) / (E2-E0) / (E1-E0)
       ELSE
          w_add = E2*E2 / (E2-E1) / (E2-E0) / (E2-E1)
       ENDIF
    ELSE IF (ec < ea .AND. ea < eb) THEN
       IF (ea > 0.0_dp) THEN
          w_add = E0*E0 / (E1-E0) / (E2-E0) / (E1-E0)
       ELSE
          w_add = E2*E2 / (E2-E1) / (E2-E0) / (E2-E1)
       ENDIF
    ELSE IF (eb < ec .AND. ec < ea) THEN
       IF (ec > 0.0_dp) THEN
          w_add = E0*E0 / (E1-E0) / (E2-E0) / (E2-E0)
       ELSE
          w_add = E2 / (E2-E1) / (E2-E0)
          w_add = w_add * (2.0_dp - E2/(E2-E1) - E2/(E2-E0))
       ENDIF
    ELSE IF (ec < eb .AND. eb < ea) THEN
       IF (eb > 0.0_dp) THEN
          w_add = E0*E0 / (E1-E0) / (E2-E0) / (E2-E0)
       ELSE
          w_add = E2 / (E2-E1) / (E2-E0)
          w_add = w_add * (2.0_dp - E2/(E2-E1) - E2/(E2-E0))
       ENDIF
    ELSE
       get_triangular_weight = 0.0_dp
       RETURN
    ENDIF

    get_triangular_weight = w_add / 2.0_dp
  END FUNCTION get_triangular_weight

  REAL(dp) FUNCTION get_qweight_2d(iq, ik, ibnd, jbnd, imode, &
                                    freq, bande, nqf, F1, F2, nqtot, nbnd, nmode)
    INTEGER, INTENT(IN) :: iq, ik, ibnd, jbnd, imode, nqf, nqtot, nbnd, nmode
    REAL(dp), INTENT(IN) :: freq(0:nqtot-1, nmode), bande(0:nqtot-1, nbnd)
    REAL(dp), INTENT(IN) :: F1, F2

    INTEGER, PARAMETER :: num_nest = 6
    INTEGER :: nest_ivec(num_nest), iq_nest(num_nest), ikq_nest(num_nest)
    INTEGER :: ikq, i, iv1, iv2
    REAL(dp) :: E0_f1, E0_f2, E_nest_f1(num_nest), E_nest_f2(num_nest)
    REAL(dp) :: weight_f1, weight_f2, emin_f1, emax_f1, emin_f2, emax_f2

    nest_ivec(1) = (nqf-1)*nqf
    nest_ivec(2) = (nqf-1)*nqf + 1
    nest_ivec(3) = 1
    nest_ivec(4) = nqf
    nest_ivec(5) = nqf*2 - 1
    nest_ivec(6) = nqf - 1

    ikq = kindex_add_2d(ik, iq, nqf, nqf)
    DO i = 1, num_nest
       iq_nest(i) = kindex_add_2d(iq, nest_ivec(i), nqf, nqf)
       ikq_nest(i) = kindex_add_2d(ik, iq_nest(i), nqf, nqf)
    ENDDO

    E0_f1 = bande(ik, ibnd) - bande(ikq, jbnd) - freq(iq, imode)
    E0_f2 = bande(ik, ibnd) - bande(ikq, jbnd) + freq(iq, imode)
    DO i = 1, num_nest
       E_nest_f1(i) = bande(ik, ibnd) - bande(ikq_nest(i), jbnd) - freq(iq_nest(i), imode)
       E_nest_f2(i) = bande(ik, ibnd) - bande(ikq_nest(i), jbnd) + freq(iq_nest(i), imode)
    ENDDO

    weight_f1 = 0.0_dp
    emin_f1 = MIN(E0_f1, MINVAL(E_nest_f1))
    emax_f1 = MAX(E0_f1, MAXVAL(E_nest_f1))
    IF (emin_f1 < 0.0_dp .AND. emax_f1 > 0.0_dp) THEN
       DO i = 1, num_nest
          iv1 = i
          iv2 = MOD(i, num_nest) + 1
          weight_f1 = weight_f1 + get_triangular_weight(E0_f1, E_nest_f1(iv1), E_nest_f1(iv2))
       ENDDO
       weight_f1 = weight_f1 * F1
    ENDIF

    weight_f2 = 0.0_dp
    emin_f2 = MIN(E0_f2, MINVAL(E_nest_f2))
    emax_f2 = MAX(E0_f2, MAXVAL(E_nest_f2))
    IF (emin_f2 < 0.0_dp .AND. emax_f2 > 0.0_dp) THEN
       DO i = 1, num_nest
          iv1 = i
          iv2 = MOD(i, num_nest) + 1
          weight_f2 = weight_f2 + get_triangular_weight(E0_f2, E_nest_f2(iv1), E_nest_f2(iv2))
       ENDDO
       weight_f2 = weight_f2 * F2
    ENDIF

    get_qweight_2d = pi_val * (1.0_dp / DBLE(nqf**2)) * (weight_f1 + weight_f2)
  END FUNCTION get_qweight_2d

  SUBROUTINE get_qweight_wrap_2d(ik, ibnd, freq, bande, nqf, efermi, &
                                  nqtot, nbnd, nmode, weight_iq)
    INTEGER, INTENT(IN) :: ik, ibnd, nqf, nqtot, nbnd, nmode
    REAL(dp), INTENT(IN) :: freq(0:nqtot-1, nmode), bande(0:nqtot-1, nbnd)
    REAL(dp), INTENT(IN) :: efermi
    REAL(dp), INTENT(OUT) :: weight_iq(0:nqtot-1, nmode, nbnd)

    REAL(dp) :: bande_shifted(0:nqtot-1, nbnd)
    REAL(dp) :: kT, ikibe, freq_cut, tmpfreq, wgq, wgkq, F1, F2
    INTEGER :: iq, im, jbnd, ikq

    kT = 300.0_dp * kB_eV
    bande_shifted = bande - efermi
    ikibe = bande_shifted(ik, ibnd)
    freq_cut = MAXVAL(freq) + 100.0_dp / DBLE(nqf)
    weight_iq = 0.0_dp

    DO iq = 0, nqtot - 1
       DO im = 1, nmode
          ikq = kindex_add_2d(ik, iq, nqf, nqf)
          tmpfreq = freq(iq, im)
          IF (tmpfreq > 1.0d-50) THEN
             wgq = 1.0_dp / (EXP(tmpfreq / kT) - 1.0_dp)
          ELSE
             wgq = 0.0_dp
          ENDIF
          DO jbnd = 1, nbnd
             IF (ABS(bande_shifted(ikq, jbnd) - ikibe) > freq_cut) CYCLE
             wgkq = 1.0_dp / (EXP(bande_shifted(ikq, jbnd) / kT) + 1.0_dp)
             F1 = 1.0_dp - wgkq + wgq
             F2 = wgkq + wgq
             weight_iq(iq, im, jbnd) = get_qweight_2d(iq, ik, ibnd, jbnd, im, &
                                        freq, bande_shifted, nqf, F1, F2, nqtot, nbnd, nmode)
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE get_qweight_wrap_2d

  REAL(dp) FUNCTION get_tetrahedral_weight(ea, eb, ec, ed)
    REAL(dp), INTENT(IN) :: ea, eb, ec, ed
    REAL(dp) :: e(4), w_add
    INTEGER :: npos, nneg

    e(1) = ea; e(2) = eb; e(3) = ec; e(4) = ed
    CALL sort4(e)

    npos = 0; nneg = 0
    IF (e(1) >= 0.0_dp) npos = npos + 1
    IF (e(2) >= 0.0_dp) npos = npos + 1
    IF (e(3) >= 0.0_dp) npos = npos + 1
    IF (e(4) >= 0.0_dp) npos = npos + 1
    IF (e(1) <= 0.0_dp) nneg = nneg + 1
    IF (e(2) <= 0.0_dp) nneg = nneg + 1
    IF (e(3) <= 0.0_dp) nneg = nneg + 1
    IF (e(4) <= 0.0_dp) nneg = nneg + 1

    IF (npos == 4 .OR. nneg == 4) THEN
       get_tetrahedral_weight = 0.0_dp
       RETURN
    ENDIF

    IF (e(1) < 0.0_dp .AND. e(2) < 0.0_dp .AND. e(3) < 0.0_dp .AND. e(4) > 0.0_dp) THEN
       w_add = e(4)**2 / (e(4)-e(1)) / (e(4)-e(2)) / (e(4)-e(3))
    ELSE IF (e(1) < 0.0_dp .AND. e(2) < 0.0_dp .AND. e(3) > 0.0_dp .AND. e(4) > 0.0_dp) THEN
       w_add = 1.0_dp / (e(4)-e(1)) / (e(3)-e(1)) * &
               (e(3)*e(4) - e(1)*(e(3)+e(4)-e(1)) + e(2)*(e(3)+e(4)-2.0_dp*e(1))) / &
               (e(4)-e(2)) * (-e(1)) / (e(3)-e(2))
       w_add = w_add + (-e(1))**2 * e(4) / (e(4)-e(1))**2 / (e(3)-e(1)) / (e(4)-e(2))
       w_add = w_add + (-e(2))**2 * e(3) / (e(3)-e(2))**2 / (e(4)-e(2)) / (e(3)-e(1))
    ELSE IF (e(1) < 0.0_dp .AND. e(2) > 0.0_dp .AND. e(3) > 0.0_dp .AND. e(4) > 0.0_dp) THEN
       w_add = (-e(1))**2 / (e(2)-e(1)) / (e(3)-e(1)) / (e(4)-e(1))
    ELSE
       get_tetrahedral_weight = 0.0_dp
       RETURN
    ENDIF

    get_tetrahedral_weight = w_add / 6.0_dp
  END FUNCTION get_tetrahedral_weight

  REAL(dp) FUNCTION get_qweight_3d(iq, ik, ibnd, jbnd, imode, &
                                    freq, bande, nk1, nk2, nk3, &
                                    F1, F2, nqtot, nbnd, nmode)
    INTEGER, INTENT(IN) :: iq, ik, ibnd, jbnd, imode
    INTEGER, INTENT(IN) :: nk1, nk2, nk3, nqtot, nbnd, nmode
    REAL(dp), INTENT(IN) :: freq(0:nqtot-1, nmode), bande(0:nqtot-1, nbnd)
    REAL(dp), INTENT(IN) :: F1, F2

    INTEGER, PARAMETER :: num_nest = 6
    INTEGER :: nest_ivec(num_nest), iq_nest(num_nest), ikq_nest(num_nest)
    INTEGER :: ikq, i, nk23
    REAL(dp) :: E0_f1, E0_f2, E_nest_f1(num_nest), E_nest_f2(num_nest)
    REAL(dp) :: weight_f1, weight_f2, emin_f1, emax_f1, emin_f2, emax_f2

    nk23 = nk2 * nk3
    nest_ivec(1) = nk23
    nest_ivec(2) = nk3
    nest_ivec(3) = 1
    nest_ivec(4) = nk23 + nk3
    nest_ivec(5) = nk23 + 1
    nest_ivec(6) = nk3 + 1

    ikq = kindex_add_3d(ik, iq, nk1, nk2, nk3)
    DO i = 1, num_nest
       iq_nest(i) = kindex_add_3d(iq, nest_ivec(i), nk1, nk2, nk3)
       ikq_nest(i) = kindex_add_3d(ik, iq_nest(i), nk1, nk2, nk3)
    ENDDO

    E0_f1 = bande(ik, ibnd) - bande(ikq, jbnd) - freq(iq, imode)
    E0_f2 = bande(ik, ibnd) - bande(ikq, jbnd) + freq(iq, imode)
    DO i = 1, num_nest
       E_nest_f1(i) = bande(ik, ibnd) - bande(ikq_nest(i), jbnd) - freq(iq_nest(i), imode)
       E_nest_f2(i) = bande(ik, ibnd) - bande(ikq_nest(i), jbnd) + freq(iq_nest(i), imode)
    ENDDO

    weight_f1 = 0.0_dp
    emin_f1 = MIN(E0_f1, MINVAL(E_nest_f1))
    emax_f1 = MAX(E0_f1, MAXVAL(E_nest_f1))
    IF (emin_f1 < 0.0_dp .AND. emax_f1 > 0.0_dp) THEN
       DO i = 1, num_nest - 1
          weight_f1 = weight_f1 + get_tetrahedral_weight(E0_f1, E_nest_f1(i), &
                                   E_nest_f1(i+1), E_nest_f1(MOD(i+1, num_nest)+1))
       ENDDO
       weight_f1 = weight_f1 * F1
    ENDIF

    weight_f2 = 0.0_dp
    emin_f2 = MIN(E0_f2, MINVAL(E_nest_f2))
    emax_f2 = MAX(E0_f2, MAXVAL(E_nest_f2))
    IF (emin_f2 < 0.0_dp .AND. emax_f2 > 0.0_dp) THEN
       DO i = 1, num_nest - 1
          weight_f2 = weight_f2 + get_tetrahedral_weight(E0_f2, E_nest_f2(i), &
                                   E_nest_f2(i+1), E_nest_f2(MOD(i+1, num_nest)+1))
       ENDDO
       weight_f2 = weight_f2 * F2
    ENDIF

    get_qweight_3d = pi_val * (1.0_dp / DBLE(nk1*nk2*nk3)) * (weight_f1 + weight_f2)
  END FUNCTION get_qweight_3d

  REAL(dp) FUNCTION get_delta_weight(iq, ik, ibnd, jbnd, imode, &
                                      freq, bande, nk1, nk2, nk3, &
                                      F1, F2, nqtot, nbnd, nmode, ndim)
    INTEGER, INTENT(IN) :: iq, ik, ibnd, jbnd, imode
    INTEGER, INTENT(IN) :: nk1, nk2, nk3, nqtot, nbnd, nmode, ndim
    REAL(dp), INTENT(IN) :: freq(0:nqtot-1, nmode), bande(0:nqtot-1, nbnd)
    REAL(dp), INTENT(IN) :: F1, F2

    IF (ndim == 2) THEN
       get_delta_weight = get_qweight_2d(iq, ik, ibnd, jbnd, imode, &
                           freq, bande, nk1, F1, F2, nqtot, nbnd, nmode)
    ELSE
       get_delta_weight = get_qweight_3d(iq, ik, ibnd, jbnd, imode, &
                           freq, bande, nk1, nk2, nk3, F1, F2, nqtot, nbnd, nmode)
    ENDIF
  END FUNCTION get_delta_weight

  SUBROUTINE sort3(a)
    REAL(dp), INTENT(INOUT) :: a(3)
    REAL(dp) :: tmp
    IF (a(1) > a(2)) THEN; tmp = a(1); a(1) = a(2); a(2) = tmp; ENDIF
    IF (a(2) > a(3)) THEN; tmp = a(2); a(2) = a(3); a(3) = tmp; ENDIF
    IF (a(1) > a(2)) THEN; tmp = a(1); a(1) = a(2); a(2) = tmp; ENDIF
  END SUBROUTINE sort3

  SUBROUTINE sort4(a)
    REAL(dp), INTENT(INOUT) :: a(4)
    REAL(dp) :: tmp
    INTEGER :: i, j
    DO i = 1, 3
       DO j = 1, 4-i
          IF (a(j) > a(j+1)) THEN
             tmp = a(j); a(j) = a(j+1); a(j+1) = tmp
          ENDIF
       ENDDO
    ENDDO
  END SUBROUTINE sort4

END MODULE delta_weights
