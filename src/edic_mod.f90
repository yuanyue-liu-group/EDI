Module edic_mod
    USE kinds,            ONLY : DP
    USE control_flags,    ONLY : gamma_only, smallmem
    USE gvect,            ONLY : gstart
    !USE noncollin_module, ONLY : noncolin, npol
    Use becmod,           Only : bec_type
      USE io_files,  ONLY : prefix, tmp_dir, nwordwfc, iunwfc, restart_dir
  
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
          character(len=256) :: V_d_filename = 'none', Bxc_1_filename='none', Bxc_2_filename='none', Bxc_3_filename='none'
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
  CHARACTER (len=256) :: filband, filp, outdir
  LOGICAL :: lsigma(4), lsym, lp, no_overlap, plot_2d, wfc_is_collected, exst
  LOGICAL :: lspinorb, noncolin 
  LOGICAL :: lvacalign=.true.
  LOGICAL :: lcorealign=.false.
  INTEGER :: spin_component, firstk, lastk
  INTEGER :: nspin
                integer :: vac_idx=0
          REAL :: core_v_d=0.0
          REAL :: core_v_p=0.0

  NAMELIST / calcmcontrol / vperturb_filename_p,vperturb_filename_d,eps_filename, kpoint_initial, kpoint_final, &
            bnd_initial, bnd_final, calcmlocal,calcmnonlocal,calcmcharge, mcharge_dolfa,k0screen_read,&
            V_d_filename, Bxc_1_filename, Bxc_2_filename, Bxc_3_filename, V_p_filename,&
            V_up_filename, V_down_filename,&
   outdir, prefix, filband, filp, spin_component, lsigma,&
                       lsym, lp, filp, firstk, lastk, no_overlap, plot_2d,&
noncolin , lspinorb  ,nspin,lvacalign,lcorealign,vac_idx,core_v_d,core_v_p

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
!                integer :: vac_direction
                real(dp) :: celldm(6), gcutm, dual, ecut,  at(3,3), omega, alat
                integer, allocatable:: ityp(:), intyp(:)
                real(dp),allocatable:: zv(:), tau(:, :)  , plot(:)
        end type V_file

        type(V_file), target :: V_d, Bxc_1, Bxc_2, Bxc_3, V_p

                real(dp) :: v_p_shift,v_d_shift
        real(dp), allocatable, target :: &
                V_nc(:, :),V_colin(:)

End Module edic_mod

