MODULE global_var
  USE kinds, ONLY : DP
  IMPLICIT NONE
  SAVE

  INTEGER :: nbndep = 0
  INTEGER :: nbndskip = 0
  INTEGER, ALLOCATABLE :: ibndkept(:)
  COMPLEX(KIND=DP), ALLOCATABLE :: umat(:,:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: umat_all(:,:,:)
  REAL(KIND=DP), ALLOCATABLE :: xkq(:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: dmec(:,:,:,:)
  COMPLEX(KIND=DP), ALLOCATABLE :: pmec(:,:,:,:)
  INTEGER, ALLOCATABLE :: wanplotlist(:)
  INTEGER :: num_wannier_plot = 0

END MODULE global_var
