Module input_parameters_calcmdefect
      Use kinds,    Only : dp

      Implicit None

      Save

   CHARACTER(LEN=256) :: vperturb_filename='vloc.dat'
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

  NAMELIST / calcmcontrol / vperturb_filename,eps_filename, kpoint_initial, kpoint_final, &
            bnd_initial, bnd_final, calcmlocal,calcmnonlocal,calcmcharge, mcharge_dolfa,k0screen_read,&
            V_0_filename, Bxc_1_filename, Bxc_2_filename, Bxc_3_filename, V_p_filename,&
            V_up_filename, V_down_filename

            Complex(dp) :: m_loc, m_nloc

End Module input_parameters_calcmdefect
