
subroutine get_gind_rhoandpsi_gw(gind_rho2psi_gw,gind_psi2rho_gw,gw_ng_data,&
gw_q_g_commonsubset_size,gvec_gw,gw_bvec_data,gw_q_g_commonsubset_indinrho)

    integer(DP),allocatable ,intent(inout)::gind_rho2psi_gw(:)
    integer(DP),allocatable ,intent(inout)::gind_psi2rho_gw(:)
    real(DP) ,intent(inout)::gvec_gw(3)
    real(DP) ,intent(inout)::gw_bvec_data(3,3)
  integer, allocatable ,intent (inout) :: gw_ng_data(:)
  integer, allocatable  :: gw_qgrid_data(:),gw_nq_data(:),gw_fftgrid_data(:),gw_nmtx_max_data(:)
    INTEGER ,intent(in) :: gw_q_g_commonsubset_size
    integer(DP),allocatable ,intent(in) ::gw_q_g_commonsubset_indinrho(:)
real(dp)::dgtmp
   integer, allocatable  :: gw_g_components_data(:,:)
   !integer, allocatable,intent (inout)  :: gw_g_components_data(:,:)

      allocate(gind_rho2psi_gw(gw_ng_data(1)))
      allocate(gind_psi2rho_gw(ngm))
!write(*,*) 'gw1'
      gind_rho2psi_gw(:)=0
      gind_psi2rho_gw(:)=0
         npw = ngk(ik)
      write(*,*)  gw_ng_data(1),ngm
      write(*,*)  gw_q_g_commonsubset_size,npwx,npw
      do ig1 = 1, gw_q_g_commonsubset_size
        do ig2=1,npw

          gvec_gw= gw_g_components_data(1,gw_q_g_commonsubset_indinrho(ig1))*gw_bvec_data(:,1)+ &
                   gw_g_components_data(2,gw_q_g_commonsubset_indinrho(ig1))*gw_bvec_data(:,2)+ &
                   gw_g_components_data(3,gw_q_g_commonsubset_indinrho(ig1))*gw_bvec_data(:,3)

          dgtmp=norm2(g(1:3,igk_k(ig2,ik0))-gvec_gw)
!write(*,*) 'gw2'
          if (dgtmp<machine_eps)then
      !      write(*,*) 'gw rho gind to psi gind: ig1,ig2', ig1, ig2,dgtmp,gvec_gw
      !      write(*,*) gw_g_components_data(:,ig1),g(1:3,igk_k(ig2,ik0)),gvec_gw
            write(*,*) 'gw rho gind to psi gind: ig1,ig2', ig1, ig2,dgtmp,gvec_gw
            write(*,*) 'ig2',ig2,ik0,igk_k(ig2,ik0),g(1:3,igk_k(ig2,ik0))
            gind_rho2psi_gw(ig1)=ig2
            gind_psi2rho_gw(ig2)=ig1
          endif
        enddo
      enddo

      write(*,*)  gw_ng_data(1),ngm

end subroutine get_gind_rhoandpsi_gw




