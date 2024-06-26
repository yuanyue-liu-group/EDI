Module edic_mod
  USE kinds,            ONLY : DP
  USE control_flags,    ONLY : gamma_only, smallmem
  USE gvect,            ONLY : gstart
  Use becmod,           Only : bec_type
  USE io_files,  ONLY : prefix, tmp_dir, nwordwfc, iunwfc, restart_dir
 
  Implicit None
  Save

  !scattering weight sum delta(Ei-Ef)
  type ::s_wt
    integer,allocatable::kp_idx(:,:)!i_pair, ki and kf
    integer,allocatable::bnd_idx(:,:)!i_pair, ibnd and fbnd
    real,allocatable::k_coord(:,:,:)! i_pair, kxyz,kif
    real,allocatable::e_pair(:,:)! i_pair, ei and ef
    real,allocatable::v_pair(:,:,:)! i_pair, vxyz,kif
    real,allocatable::wt(:)!i_pair
    complex,allocatable::m(:)!i_pair
    complex,allocatable::mc(:)!i_pair

    complex,allocatable::ml(:)!i_pair
    complex,allocatable::mnl_p(:)!i_pair
    complex,allocatable::mnl_d(:)!i_pair
    complex,allocatable::mc_lfa(:)!i_pair
    complex,allocatable::mc_nolfa(:)!i_pair

    complex,allocatable::mc_lfa_2dns(:)!i_pair
    complex,allocatable::mc_lfa_2dtf(:)!i_pair
    complex,allocatable::mc_lfa_2dqeh(:)!i_pair
    complex,allocatable::mc_lfa_2dgw(:)!i_pair

    complex,allocatable::mc_lfa_3dns(:)!i_pair
    complex,allocatable::mc_lfa_3dtf(:)!i_pair
    complex,allocatable::mc_lfa_3dqeh(:)!i_pair
    complex,allocatable::mc_lfa_3dgw(:)!i_pair

    complex,allocatable::mc_lfa_3dcutns(:)!i_pair
    complex,allocatable::mc_lfa_3dcuttf(:)!i_pair
    complex,allocatable::mc_lfa_3dcutqeh(:)!i_pair
    complex,allocatable::mc_lfa_3dcutgw(:)!i_pair

    complex,allocatable::mc_nolfa_2dns(:)!i_pair
    complex,allocatable::mc_nolfa_2dtf(:)!i_pair
    complex,allocatable::mc_nolfa_2dqeh(:)!i_pair
    complex,allocatable::mc_nolfa_2dgw(:)!i_pair

    complex,allocatable::mc_nolfa_3dns(:)!i_pair
    complex,allocatable::mc_nolfa_3dtf(:)!i_pair
    complex,allocatable::mc_nolfa_3dqeh(:)!i_pair
    complex,allocatable::mc_nolfa_3dgw(:)!i_pair

    complex,allocatable::mc_nolfa_3dcutns(:)!i_pair
    complex,allocatable::mc_nolfa_3dcuttf(:)!i_pair
    complex,allocatable::mc_nolfa_3dcutqeh(:)!i_pair
    complex,allocatable::mc_nolfa_3dcutgw(:)!i_pair
    integer:: npairs
  end type 
  
  type(S_wt)::bndkp_pair
  
  !scattering weight sum delta(Ei-Ef)
  type ::s_wt_rj
  integer::k_initial!i_pair, kif
  integer::k_final! i_pair, kxyz,kif
  integer::bnd_initial!i_pair, kif
  integer::bnd_final! i_pair, kxyz,kif
  real::wt!i_pair
  end type 
  
  type(S_wt_rj),allocatable::wfc_pair(:)
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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! control file variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! epsilon defect potential files
  real(DP),allocatable:: qeh_eps_data (:,:)
  real(DP),allocatable:: chi_data (:,:)
  CHARACTER(LEN=256) :: qeh_eps_filename='eps.dat'
  CHARACTER(LEN=256) :: chi_filename='chi.dat'
  integer :: nchilines, nqxofchi,nqyofchi,nqzofchi 
  CHARACTER(LEN=256) :: gw_epsmat_filename='epsmat.h5'
  CHARACTER(LEN=256) :: gw_eps0mat_filename='eps0mat.h5'
  character(len=256) :: V_d_filename = 'none', Bxc_1_d_filename='none', Bxc_2_d_filename='none', Bxc_3_d_filename='none'
  character(len=256) :: V_p_filename = 'none', Bxc_1_p_filename='none', Bxc_2_p_filename='none', Bxc_3_p_filename='none'
  character(len=256) :: V_up_filename='none', V_down_filename='none'

  ! TF info
  REAL :: k0screen_read=-10.0
  
  ! control: dummy keys 
  INTEGER :: kpoint_initial
  INTEGER :: kpoint_final
  INTEGER :: bnd_initial
  INTEGER :: bnd_final

  ! keys
  LOGICAL :: calcmlocal =.false.
  LOGICAL :: calcmnonlocal =.false.
  LOGICAL :: calcmcharge =.false.
  LOGICAL :: mcharge_dolfa =.false.

  CHARACTER (len=256) ::  outdir
  CHARACTER (len=256) ::  eps_type
  CHARACTER (len=256) ::  wt_filename,klist_filename,ev_filename
  LOGICAL ::  wfc_is_collected, exst
  LOGICAL :: lspinorb=.false., noncolin =.false.
  LOGICAL :: lvacalign=.true.
  LOGICAL :: lcorealign=.false.
  LOGICAL :: doqeh=.false.
  LOGICAL :: do3d=.false.
  LOGICAL :: do2d=.false.
  LOGICAL :: dogwfull=.false.
  LOGICAL :: dogwdiag=.false.
  INTEGER :: spin_component, firstk, lastk
  INTEGER :: nspin
  integer :: vac_idx=0
  REAL :: core_v_d=0.0
  REAL :: core_v_p=0.0
  REAL(dp), PARAMETER :: machine_eps = 1.d-4


  NAMELIST / calcmcontrol / &
            wt_filename,klist_filename,ev_filename,   outdir, prefix, &
            calcmlocal,calcmnonlocal,calcmcharge, mcharge_dolfa,k0screen_read,&
            qeh_eps_filename,gw_epsmat_filename,gw_eps0mat_filename, &
            eps_type,doqeh,dogwfull,dogwdiag,do2d,do3d,&
            V_d_filename, Bxc_1_d_filename, Bxc_2_d_filename, Bxc_3_d_filename,&
            V_p_filename, Bxc_1_p_filename, Bxc_2_p_filename, Bxc_3_p_filename,&
            noncolin , lspinorb  ,nspin,lvacalign,lcorealign,vac_idx,core_v_d,core_v_p,lconstalign,Eshift_d,Eshift_p

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
          real(dp),allocatable:: zv(:), tau(:, :)  , pot(:)
  end type V_file

  type(V_file), target :: V_d, Bxc_1_d, Bxc_2_d, Bxc_3_d, V_p, Bxc_1_p, Bxc_2_p, Bxc_3_p


  real(dp) :: v_p_shift,v_d_shift
  real(dp), allocatable, target :: &
                V_nc(:, :),V_nc1(:, :),V_nc2(:, :),V_nc3(:, :),V_nc4(:, :),V_colin(:)



  
  type :: gw_eps_data
    real(dp), allocatable :: epsmat_diag_data(:,:,:)
    real(dp), allocatable :: epsmat_full_data(:,:,:,:,:,:)
    real(dp), allocatable :: vcoul_data(:,:),qpts_data(:,:)
    real(dp), allocatable :: blat_data(:),bvec_data(:,:)
  
    integer, allocatable :: nq_data(:),nmtx_max_data(:),fftgrid_data(:),qgrid_data(:),ng_data(:)
    integer, allocatable :: gind_rho2eps_data(:,:)
    integer, allocatable :: gind_eps2rho_data(:,:)
    integer, allocatable :: nmtx_data(:)
    integer, allocatable :: grho_data(:),  geps_data(:),g_components_data(:,:)
  
    real(DP),allocatable ::qabs(:)
  
    integer(DP),allocatable ::q_g_commonsubset_indinrhotmp1(:)
    INTEGER :: q_g_commonsubset_size
    integer(DP),allocatable ::q_g_commonsubset_indinrho(:)
  
    integer(DP),allocatable ::gind_rho2psi(:)
    real(DP) ::gvec(3)
    integer(DP),allocatable ::gind_psi2rho(:)
  end type  gw_eps_data
  type(gw_eps_data),target:: gw_epsq1_data
  type(gw_eps_data),target:: gw_epsq0_data

End Module edic_mod

