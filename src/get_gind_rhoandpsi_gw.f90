
!subroutine get_gind_rhoandpsi_gw(gind_rho2psi_gw,gind_psi2rho_gw,gw_ng_data,&
!gw_q_g_commonsubset_size,gvec_gw,gw_bvec_data,gw_q_g_commonsubset_indinrho)
subroutine get_gind_rhoandpsi_gw(gw_)

USE kinds, ONLY: DP,sgl
USE gvect, ONLY: ngm, gstart, g, gg, gcutm, igtongl
USE klist , ONLY: nks, nelec, xk, wk, degauss, ngauss, igk_k, ngk
use edic_mod, only: gw_eps_data
use edic_mod, only: machine_eps
!type(gw_eps_data),target:: gw_
type(gw_eps_data),intent (inout) ,target:: gw_


!    integer(DP),allocatable ,intent(inout)::gind_rho2psi_gw(:)
!    integer(DP),allocatable ,intent(inout)::gind_psi2rho_gw(:)
!    real(DP) ,intent(inout)::gvec_gw(3)
!    real(DP) ,intent(inout)::gw_bvec_data(3,3)
!  integer, allocatable ,intent (inout) :: gw_ng_data(:)
!  integer, allocatable  :: gw_qgrid_data(:),gw_nq_data(:),gw_fftgrid_data(:),gw_nmtx_max_data(:)
!    INTEGER ,intent(in) :: gw_q_g_commonsubset_size
!    integer(DP),allocatable ,intent(in) ::gw_q_g_commonsubset_indinrho(:)
real(dp)::dgtmp
!   integer, allocatable  :: gw_g_components_data(:,:)
!   !integer, allocatable,intent (inout)  :: gw_g_components_data(:,:)

      allocate(gw_%gind_rho2psi(gw_%ng_data(1)))
      allocate(gw_%gind_psi2rho(ngm))
!write(*,*) 'gw1'
      gw_%gind_rho2psi(:)=0
      gw_%gind_psi2rho(:)=0
         npw = ngk(ik)
      write(*,*)  gw_%ng_data(1),ngm
      write(*,*)  gw_%q_g_commonsubset_size,npwx,npw
      do ig1 = 1, gw_%q_g_commonsubset_size
        do ig2=1,npw

          gw_%gvec= gw_%g_components_data(1,gw_%q_g_commonsubset_indinrho(ig1))*gw_%bvec_data(:,1)+ &
                   gw_%g_components_data(2,gw_%q_g_commonsubset_indinrho(ig1))*gw_%bvec_data(:,2)+ &
                   gw_%g_components_data(3,gw_%q_g_commonsubset_indinrho(ig1))*gw_%bvec_data(:,3)

          dgtmp=norm2(g(1:3,igk_k(ig2,ik0))-gw_%gvec)
!write(*,*) 'gw2'
          if (dgtmp<machine_eps)then
      !      write(*,*) 'gw rho gind to psi gind: ig1,ig2', ig1, ig2,dgtmp,gw_%gvec
      !      write(*,*) gw_%g_components_data(:,ig1),g(1:3,igk_k(ig2,ik0)),gw_%gvec
            write(*,*) 'gw rho gind to psi gind: ig1,ig2', ig1, ig2,dgtmp,gw_%gvec
            write(*,*) 'ig2',ig2,ik0,igk_k(ig2,ik0),g(1:3,igk_k(ig2,ik0))
            gw_%gind_rho2psi(ig1)=ig2
            gw_%gind_psi2rho(ig2)=ig1
          endif
        enddo
      enddo

      write(*,*)  gw_%ng_data(1),ngm

end subroutine get_gind_rhoandpsi_gw




