MODULE parallelism
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: fkbounds, poolgather

CONTAINS

  SUBROUTINE fkbounds(nktot, lower_bnd, upper_bnd)
    USE mp_global, ONLY : my_pool_id, npool
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nktot
    INTEGER, INTENT(OUT) :: lower_bnd, upper_bnd
    INTEGER :: nkl, nkr
#if defined(__MPI)
    nkl = nktot / npool
    nkr = nktot - nkl * npool
    IF (my_pool_id < nkr) nkl = nkl + 1
    lower_bnd = my_pool_id * nkl + 1
    IF (my_pool_id >= nkr) lower_bnd = my_pool_id * nkl + 1 + nkr
    upper_bnd = lower_bnd + nkl - 1
#else
    lower_bnd = 1
    upper_bnd = nktot
#endif
  END SUBROUTINE fkbounds

  SUBROUTINE poolgather(ndim, nkstot, nks, arr_loc, arr_all)
    USE kinds, ONLY : DP
    USE mp_global, ONLY : inter_pool_comm
    USE mp, ONLY : mp_sum
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ndim, nkstot, nks
    REAL(KIND=DP), INTENT(IN) :: arr_loc(ndim, nks)
    REAL(KIND=DP), INTENT(OUT) :: arr_all(ndim, nkstot)
    arr_all = 0.0_DP
    arr_all(:, 1:nks) = arr_loc(:, 1:nks)
    CALL mp_sum(arr_all, inter_pool_comm)
  END SUBROUTINE poolgather

END MODULE parallelism
