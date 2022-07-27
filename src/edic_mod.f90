Module edic_mod
    USE kinds,            ONLY : DP
    USE control_flags,    ONLY : gamma_only, smallmem
    USE gvect,            ONLY : gstart
    USE noncollin_module, ONLY : noncolin, npol
    Use becmod,           Only : bec_type
  
    Implicit None
    Save

! <beta|psi> with proper supercell handling 
    TYPE (bec_type) :: becp1_perturb  ! <beta|psi>
    TYPE (bec_type) :: becp2_perturb  ! <beta|psi>

! factor with proper supercell handling 
    Complex(dp), allocatable :: eigts1_perturb(:,:), eigts2_perturb(:,:), eigts3_perturb(:,:)


! primitive wavefunctions 
    complex(dp), allocatable, target :: &
        evc1(: , :), &
        evc2(: , :), &
        evc3(: , :), &
        evc4(: , :) 
    complex(dp), allocatable, target :: &
        psic1(:), &
        psic2(:), &
        psic3(:), &
        psic4(:) 

! input parameter
   CHARACTER(LEN=256) :: vperturb_filename_p='vloc-p.dat'
   CHARACTER(LEN=256) :: vperturb_filename_d='vloc-d.dat'
  CHARACTER(LEN=256) :: eps_filename='eps.dat'
          character(len=256) :: V_0_filename = 'none', Bxc_1_filename='none', Bxc_2_filename='none', Bxc_3_filename='none'
          character(len=256) :: V_p_filename='none'
          character(len=256) :: V_up_filename='none', V_down_filename='none'
          INTEGER :: kpoint_initial
          INTEGER :: kpoint_final
          INTEGER :: bnd_initial
          INTEGER :: bnd_final
          LOGICAL :: calcmlocal =.false.
          LOGICAL :: calcmnonlocal =.false.
          LOGICAL :: calcmcharge =.false.
          LOGICAL :: mcharge_dolfa =.false.
          REAL :: k0screen_read=0.0

  NAMELIST / calcmcontrol / vperturb_filename_p,vperturb_filename_d,eps_filename, kpoint_initial, kpoint_final, &
            bnd_initial, bnd_final, calcmlocal,calcmnonlocal,calcmcharge, mcharge_dolfa,k0screen_read,&
            V_0_filename, Bxc_1_filename, Bxc_2_filename, Bxc_3_filename, V_p_filename,&
            V_up_filename, V_down_filename

            Complex(dp) :: m_loc, m_nloc

! potential type
        type :: V_file
                integer :: unit
                character (len=75) :: filename
                character (len=75) :: title
                character (len=3) ,allocatable :: atm(:)
                integer :: nr1x, nr2x, nr3x, nr1, nr2, nr3, &
                        nat, ntyp, ibrav, plot_num,  i,nkb
                integer :: iunplot, ios, ipol, na, nt, &
                        ir, ndum
                real(dp) :: celldm(6), gcutm, dual, ecut,  at(3,3), omega, alat
                integer, allocatable:: ityp(:), intyp(:)
                real(dp),allocatable:: zv(:), tau(:, :)  , plot(:)
        end type V_file

        type(V_file), target :: V_0, Bxc_1, Bxc_2, Bxc_3, V_p

        real(dp), allocatable, target :: &
                V_loc(:, :)

End Module edic_mod

