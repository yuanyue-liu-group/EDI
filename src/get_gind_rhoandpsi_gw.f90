!subroutine get_gind_rhoandpsi_gw(gw_,ik0,ik)
subroutine get_gind_rhoandpsi_gw(gw_)
  USE kinds, ONLY: DP,sgl
  USE gvect, ONLY: ngm, gstart, g, gg, gcutm, igtongl
  USE klist , ONLY: nks, nelec, xk, wk, degauss, ngauss, igk_k, ngk
  use edic_mod, only: gw_eps_data
  use edic_mod, only: machine_eps
  type(gw_eps_data),intent (inout) ,target:: gw_
  real(dp)::dgtmp

  allocate(gw_%gind_rho2psi(gw_%q_g_commonsubset_size))
  allocate(gw_%gind_psi2rho(ngm))
  gw_%gind_rho2psi(:)=0
  gw_%gind_psi2rho(:)=0
  do ig1 = 1, gw_%q_g_commonsubset_size
    do ig2=1,ngm

      gw_%gvec= gw_%g_components_data(1,gw_%q_g_commonsubset_indinrho(ig1))*gw_%bvec_data(:,1)+ &
                gw_%g_components_data(2,gw_%q_g_commonsubset_indinrho(ig1))*gw_%bvec_data(:,2)+ &
                gw_%g_components_data(3,gw_%q_g_commonsubset_indinrho(ig1))*gw_%bvec_data(:,3)

      dgtmp=norm2(g(1:3,ig2)-gw_%gvec)
      if (dgtmp<machine_eps)then
        gw_%gind_rho2psi(ig1)=ig2
        gw_%gind_psi2rho(ig2)=ig1
      endif
    enddo
  enddo

  write(*,*)  'Routine get_gind_rho2eps done',shape(gw_%gind_rho2psi),shape(gw_%gind_psi2rho)

end subroutine get_gind_rhoandpsi_gw




