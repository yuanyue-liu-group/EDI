Module becomd_perturb
    USE kinds,            ONLY : DP
    USE control_flags,    ONLY : gamma_only, smallmem
    USE gvect,            ONLY : gstart
    USE noncollin_module, ONLY : noncolin, npol
    Use becmod,           Only : bec_type
  
    Save
    TYPE (bec_type) :: becp1_perturb  ! <beta|psi>
    TYPE (bec_type) :: becp2_perturb  ! <beta|psi>
End Module becomd_perturb  