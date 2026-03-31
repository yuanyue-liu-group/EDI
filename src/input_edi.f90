MODULE input
  USE kinds, ONLY : DP
  IMPLICIT NONE
  SAVE

  REAL(KIND=DP), ALLOCATABLE :: xk_loc(:,:)
  REAL(KIND=DP), ALLOCATABLE :: xk_all(:,:)
  REAL(KIND=DP), ALLOCATABLE :: xk_cryst(:,:)

  INTEGER :: nkc1 = 0, nkc2 = 0, nkc3 = 0

  CHARACTER(LEN=256) :: filukk = 'filukk'
  CHARACTER(LEN=256) :: kmap = ' '

  LOGICAL :: scdm_proj = .FALSE.
  CHARACTER(LEN=15) :: scdm_entanglement = 'isolated'
  REAL(KIND=DP) :: scdm_mu = 0.0_DP
  REAL(KIND=DP) :: scdm_sigma = 1.0_DP

  LOGICAL :: eig_read = .FALSE.
  LOGICAL :: wannier_plot = .FALSE.
  REAL(KIND=DP) :: wannier_plot_radius = 3.5_DP
  REAL(KIND=DP) :: wannier_plot_scale = 1.0_DP
  INTEGER :: wannier_plot_supercell(3) = (/3, 3, 3/)
  LOGICAL :: reduce_unk = .FALSE.
  LOGICAL :: exciton = .FALSE.

  INTEGER, PARAMETER :: nwanxx = 200
  CHARACTER(LEN=256) :: proj(200) = ' '
  CHARACTER(LEN=256) :: wdata(200) = ' '
  CHARACTER(LEN=256) :: bands_skipped = ' '
  INTEGER :: nbndsub = 0
  INTEGER :: iprint = 2
  REAL(KIND=DP) :: dis_win_min = -9999.0_DP
  REAL(KIND=DP) :: dis_win_max =  9999.0_DP
  REAL(KIND=DP) :: dis_froz_min = -9999.0_DP
  REAL(KIND=DP) :: dis_froz_max =  9999.0_DP
  INTEGER :: num_iter = 200
  LOGICAL :: auto_projections = .TRUE.

END MODULE input
