Module edic_mod
    USE kinds,            ONLY : DP
    USE control_flags,    ONLY : gamma_only, smallmem
    USE gvect,            ONLY : gstart
    !USE noncollin_module, ONLY : noncolin, npol
    Use becmod,           Only : bec_type
      USE io_files,  ONLY : prefix, tmp_dir, nwordwfc, iunwfc, restart_dir
  
    Implicit None
    Save
!scattering weight sum delta(Ei-Ef)
type ::s_wt

integer,allocatable::kp_idx(:,:)!i_pair, kif
integer,allocatable::bnd_idx(:,:)!i_pair, kif
real,allocatable::k_coord(:,:,:)! i_pair, kxyz,kif
real,allocatable::wt(:)!i_pair
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

real(DP),allocatable:: eps_data (:,:)
! input parameter
!   CHARACTER(LEN=256) :: vperturb_filename_p='vloc-p.dat'
!   CHARACTER(LEN=256) :: vperturb_filename_d='vloc-d.dat'
  CHARACTER(LEN=256) :: eps_filename='eps.dat'
  CHARACTER(LEN=256) :: gw_epsmat_filename='epsmat.h5'
  CHARACTER(LEN=256) :: gw_eps0mat_filename='epsmat.h5'
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
  !CHARACTER (len=256) :: filband, filp, outdir
  CHARACTER (len=256) ::  outdir
  CHARACTER (len=256) ::  wt_filename,klist_filename,ev_filename
  !LOGICAL :: lsigma(4), lsym, lp, no_overlap, plot_2d, wfc_is_collected, exst
  LOGICAL ::  wfc_is_collected, exst
  LOGICAL :: lspinorb=.false., noncolin =.false.
  LOGICAL :: lvacalign=.true.
  LOGICAL :: lcorealign=.false.
  INTEGER :: spin_component, firstk, lastk
  INTEGER :: nspin
                integer :: vac_idx=0
          REAL :: core_v_d=0.0
          REAL :: core_v_p=0.0
REAL(dp), PARAMETER :: machine_eps = 1.d-4


!  NAMELIST / calcmcontrol / vperturb_filename_p,vperturb_filename_d,eps_filename, kpoint_initial, kpoint_final, &
  NAMELIST / calcmcontrol / eps_filename,gw_epsmat_filename,gw_eps0mat_filename, kpoint_initial, kpoint_final, &
            bnd_initial, bnd_final, calcmlocal,calcmnonlocal,calcmcharge, mcharge_dolfa,k0screen_read,&
            V_d_filename, Bxc_1_d_filename, Bxc_2_d_filename, Bxc_3_d_filename, Bxc_1_p_filename, Bxc_2_p_filename, Bxc_3_p_filename, V_p_filename,&
            V_up_filename, V_down_filename,&
wt_filename,klist_filename,ev_filename,&
   outdir, prefix, &!filband, filp, spin_component, lsigma,&
                    !   lsym, lp, filp, firstk, lastk, no_overlap, plot_2d,&
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
                real(dp),allocatable:: zv(:), tau(:, :)  , pot(:)
        end type V_file

        type(V_file), target :: V_d, Bxc_1, Bxc_2, Bxc_3, V_p

                real(dp) :: v_p_shift,v_d_shift
        real(dp), allocatable, target :: &
                V_nc(:, :),V_colin(:)



!!!!!!!!!!!!!!!!!
!type :: gw_eps_data
!
!!  real(dp), allocatable :: gw_epsmat_diag_data(:,:,2),  gw_eps0mat_diag_data(:,:,2)
!  real(dp), allocatable :: gw_epsmat_diag_data_q1(:,:,:),  gw_epsmat_diag_data_q0(:,:,:)
!  !complex(dp), allocatable :: gw_epsmat_diag_data(:,:,:),  gw_eps0mat_diag_data(:,:,:)
!!  real(dp), allocatable :: gw_epsmat_full_data(:,1,1,:,:,2),  gw_eps0mat_full_data(:,1,1,:,:,2)
!  real(dp), allocatable :: gw_epsmat_full_data_q1(:,:,:,:,:,:),  gw_epsmat_full_data_q0(:,:,:,:,:,:)
!!  real(dp), allocatable :: gw_epsallmat_full_data(:,1,1,:,:,2)
!  real(dp), allocatable :: gw_epsmat_full_data_qall(:,:,:,:,:,:)
!
!  real(dp), allocatable :: gw_vcoul_data_q1(:,:),gw_qpts_Data_q1(:,:)
!  real(dp), allocatable :: gw_blat_data_q1(:),gw_bvec_Data_q1(:,:)
!  integer, allocatable :: gw_gind_eps2rho_data_q1(:,:), gw_gind_rho2eps_data_q1(:,:),gw_nmtx_data_q1(:)
!
!!q0
!  real(dp), allocatable :: gw_vcoul_data_q0(:,:),gw_qpts_Data_q0(:,:)
!  real(dp), allocatable :: gw_blat_data_q0(:),gw_bvec_Data_q0(:,:)
!  integer, allocatable :: gw_gind_eps2rho_data_q0(:,:), gw_gind_rho2eps_data_q0(:,:),gw_nmtx_data_q0(:)
!!q0
!
!
!   integer, allocatable :: gw_grho_data_q1(:),  gw_geps_data_q1(:),gw_g_components_data_q1(:,:)
!  integer, allocatable :: gw_nq_data_q1(:),gw_nmtx_max_data_q1(:),gw_fftgrid_data_q1(:),gw_qgrid_data_q1(:),gw_ng_data_q1(:)
!
!!q0
!   integer, allocatable :: gw_grho_data_q0(:),  gw_geps_data_q0(:),gw_g_components_data_q0(:,:)
!  integer, allocatable :: gw_nq_data_q0(:),gw_nmtx_max_data_q0(:),gw_fftgrid_data_q0(:),gw_qgrid_data_q0(:),gw_ng_data_q0(:)
!!q0
!
!
!!  integer(i8b), allocatable :: gw_nqi8(:)
!
!    real(DP),allocatable ::gw_qabs_q1(:)
!    INTEGER :: gw_q_g_commonsubset_size_q1
!    integer(DP),allocatable ::gw_q_g_commonsubset_indinrho_q1(:)
!
!!q0
!    real(DP),allocatable ::gw_qabs_q0(:)
!    INTEGER :: gw_q_g_commonsubset_size_q0
!    integer(DP),allocatable ::gw_q_g_commonsubset_indinrho_q0(:)
!!q0
!!!!!!!!!!!!!!!!!!!
!    integer(DP),allocatable ::gind_rho2psi_gw(:)
!    real(DP) ::gvec_gw(3)
!    integer(DP),allocatable ::gind_psi2rho_gw(:)
!
!    integer(DP),allocatable ::gind_rho2psi_gw_q0(:)
!    real(DP) ::gvec_gw_q0(3)
!    integer(DP),allocatable ::gind_psi2rho_gw_q0(:)
!
!    integer(DP),allocatable ::gind_rho2psi_gw_q1(:)
!    real(DP) ::gvec_gw_q1(3)
!    integer(DP),allocatable ::gind_psi2rho_gw_q1(:)
!!!!!!!!!!!!!!!!!!!


!end type  gw_eps_data
!!!!!!!!!!!!!!!!


type :: gw_eps_data
  real(dp), allocatable :: epsmat_diag_data(:,:,:)
  real(dp), allocatable :: epsmat_full_data(:,:,:,:,:,:)
  real(dp), allocatable :: vcoul_data(:,:),qpts_data(:,:)
  real(dp), allocatable :: blat_data(:),bvec_data(:,:)
  integer, allocatable :: gind_eps2rho_data(:,:), gind_rho2eps_data(:,:),nmtx_data(:)
   integer, allocatable :: grho_data(:),  geps_data(:),g_components_data(:,:)
  integer, allocatable :: nq_data(:),nmtx_max_data(:),fftgrid_data(:),qgrid_data(:),ng_data(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  integer(i8b), allocatable :: gw_nqi8(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer(DP),allocatable ::q_g_commonsubset_indinrhotmp1(:)
    real(DP),allocatable ::qabs(:)
    INTEGER :: q_g_commonsubset_size
    integer(DP),allocatable ::q_g_commonsubset_indinrho(:)

    integer(DP),allocatable ::gind_rho2psi(:)
    real(DP) ::gvec(3)
    integer(DP),allocatable ::gind_psi2rho(:)
end type  gw_eps_data
type(gw_eps_data),target:: gw_epsq1_data
type(gw_eps_data),target:: gw_epsq0_data

End Module edic_mod

