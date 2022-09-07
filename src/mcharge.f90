
    SUBROUTINE calcmdefect_charge_nolfa(ibnd,ibnd0,ik,ik0)
USE kinds, ONLY: DP,sgl
USE fft_base,  ONLY: dfftp, dffts
USE gvect, ONLY: ngm, gstart, g, gg, gcutm, igtongl
USE klist , ONLY: nks, nelec, xk, wk, degauss, ngauss, igk_k, ngk
!      Use edic_mod,   only: V_file, V_loc, V_0, Bxc_1, Bxc_2, Bxc_3, V_p
      Use edic_mod, Only : evc1,evc2,evc3,evc4,&
                               psic1, psic2, psic3, psic4
    use splinelib, only: dosplineint,spline,splint

      Use edic_mod,   only: gw_epsq1_data,gw_epsq0_data
USE HDF5
use edic_mod, only: machine_eps
    USE cell_base, ONLY: omega, alat, tpiba2, at, bg, tpiba
 USE constants, ONLY: tpi, e2, eps6,pi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! charge
COMPLEX(DP), ALLOCATABLE ::  mlat1(:),mlat2(:)
INTEGER :: iscx, iscy,nscx,nscy
REAL(dp)::k0screen, kbT,deltak,deltakG0,deltakG, qxy,qz,lzcutoff
INTEGER:: icount,jcount,kcount
real(DP):: mscreen,mcharge, rmod
INTEGER:: Nlzcutoff,iNlzcutoff,flag1,flag2, nNlzcutoff,Ngzcutoff
!!!!! eps data file 
integer :: nepslines
real(DP),allocatable:: eps_data (:,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!
!!
!!!!!!!!!!!!hdf5
!!  CHARACTER(LEN=256) :: h5filename      ! Dataset name
!!  CHARACTER(LEN=256) :: h5datasetname = "matrix-diagonal"     ! Dataset name
!! INTEGER     ::   h5rank,h5error ! Error flag
!!  !INTEGER     ::  i, j
!!!  real(dp), DIMENSION(3,1465,2) :: h5dataset_data, data_out ! Data buffers
!!!  real(dp), DIMENSION(3,1465,2) :: h5dataset_data, data_out ! Data buffers
!!  real(dp), allocatable :: h5dataset_data_double(:), data_out(:)
!!  integer, allocatable :: h5dataset_data_integer(:)
!!  INTEGER(HSIZE_T), allocatable :: h5dims(:),h5maxdims(:)
!! integer :: h5dims1(1),h5dims2(2),h5dims3(3),h5dims4(4),h5dims5(5),h5dims6(6)


!
!!  real(dp), allocatable :: gw_epsmat_diag_data(:,:,2),  gw_eps0mat_diag_data(:,:,2)
!  real(dp), allocatable :: gw_epsmat_diag_data_q1(:,:,:),  gw_epsmat_diag_data_q0(:,:,:)
!  !complex(dp), allocatable :: gw_epsmat_diag_data(:,:,:),  gw_eps0mat_diag_data(:,:,:)
!!  real(dp), allocatable :: gw_epsmat_full_data(:,1,1,:,:,2),  gw_eps0mat_full_data(:,1,1,:,:,2)
!  real(dp), allocatable :: gw_epsq1_data%epsmat_full_data(:,:,:,:,:,:),  gw_epsq0_data%epsmat_full_data(:,:,:,:,:,:)
!!  real(dp), allocatable :: gw_epsallmat_full_data(:,1,1,:,:,2)
!  real(dp), allocatable :: gw_epsmat_full_data_qall(:,:,:,:,:,:)
!
!  real(dp), allocatable :: gw_vcoul_data_q1(:,:),gw_epsq0_data%qpts_data(:,:)
!  real(dp), allocatable :: gw_blat_data_q1(:),gw_epsq1_data%bvec_data(:,:)
!  integer, allocatable :: gw_gind_eps2rho_data_q1(:,:), gw_epsq1_data%gind_rho2eps_data(:,:),gw_epsq1_data%nmtx_data(:)
!
!!q0
!  real(dp), allocatable :: gw_vcoul_data_q0(:,:),gw_qpts_Data_q0(:,:)
!  real(dp), allocatable :: gw_blat_data_q0(:),gw_bvec_Data_q0(:,:)
!  integer, allocatable :: gw_gind_eps2rho_data_q0(:,:), gw_epsq0_data%gind_rho2eps_data(:,:),gw_epsq0_data%nmtx_data(:)
!!q0
!
!
!   integer, allocatable :: gw_grho_data_q1(:),  gw_geps_data_q1(:),gw_g_components_data_q1(:,:)
!  integer, allocatable :: gw_epsq1_data%nq_data(:),gw_nmtx_max_data_q1(:),gw_fftgrid_data_q1(:),gw_qgrid_data_q1(:),gw_ng_data_q1(:)
!
!!q0
!   integer, allocatable :: gw_grho_data_q0(:),  gw_geps_data_q0(:),gw_g_components_data_q0(:,:)
!  integer, allocatable :: gw_epsq0_data%nq_data(:),gw_nmtx_max_data_q0(:),gw_fftgrid_data_q0(:),gw_qgrid_data_q0(:),gw_ng_data_q0(:)
!!q0
!
!
!!  integer(i8b), allocatable :: gw_nqi8(:)
!
!    real(DP),allocatable ::gw_epsq1_data%qabs(:)
!    INTEGER :: gw_q_g_commonsubset_size_q1
!    integer(DP),allocatable ::gw_epsq1_data%q_g_commonsubset_indinrho(:)
!
!!q0
!    real(DP),allocatable ::gw_epsq0_data%qabs(:)
!    INTEGER :: gw_q_g_commonsubset_size_q0
!    integer(DP),allocatable ::gw_epsq0_data%q_g_commonsubset_indinrho(:)
!!q0
!!!!!!!!!!!!!!!!!!!
!    integer(DP),allocatable ::gind_rho2psi_gw(:)
!    real(DP) ::gvec_gw(3)
    integer(DP),allocatable ::gind_psi2rho_gw(:)
!
!    integer(DP),allocatable ::gind_rho2psi_gw_q0(:)
!    real(DP) ::gvec_gw_q0(3)
!    integer(DP),allocatable ::gind_psi2rho_gw_q0(:)
!
!    integer(DP),allocatable ::gind_rho2psi_gw_q1(:)
!    real(DP) ::gvec_gw_q1(3)
!    integer(DP),allocatable ::gind_psi2rho_gw_q1(:)
!!!!!!!!!!!!!!!!!!!







    COMPLEX(DP) ::  mcharge0,mcharge1,mcharge2,mcharge3,mcharge4,mcharge5,mcharge6
    INTEGER :: ibnd, ik, ik0,ibnd0
    real(DP) , allocatable::  eps_data_dy(:)
    real(DP) :: epsk, deltakG_para,q2d_coeff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!GW
!|M|=int  u1(r) u2(r) [int eps(r,r')V(r')]
!   =int  u1(G')u2(G'+G) V(q+G)
!                        V(q+G)=int eps(q,G,G')V(q+G')
!GW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    COMPLEX(DP),allocatable ::w_gw(:)
    real(DP),allocatable ::epsint_q1_tmp1(:)
    real(DP),allocatable ::epsint_q1_tmp2(:)
    real(DP),allocatable ::epsint_q1_tmp3(:)
    real(DP),allocatable ::epsint_q1_tmp4(:)
    real(DP),allocatable ::epsint_q0_tmp1(:)
    real(DP),allocatable ::epsint_q0_tmp2(:)
    real(DP),allocatable ::epsint_q0_tmp3(:)
    real(DP),allocatable ::epsint_q0_tmp4(:)
    real(DP),allocatable ::w1(:)
    real(DP) ::epsinttmp1s
    real(DP) ::epsinttmp2s
    real(DP) ::epsinttmp3s
    real(DP) ::epsinttmp4s
    complex(DP),allocatable ::epsmat_inted(:,:) ! interpolated eps matrix
    COMPLEX(DP) :: epstmp1,epstmp2
    INTEGER :: iq1,iq2,nqgrid_gw=48!fixme
!    INTEGER :: gw_ng
    INTEGER :: gind_gw_eps1,gind_gw_eps2
    INTEGER :: gind_psi1,gind_psi2
    INTEGER :: gind_gwrho1,gind_gwrho2
real(dp)::q1(3)
logical:: interpolate_2d,interpolate_smallq1d=.false.

    INTEGER :: gw_q_g_commonsubset_size
    COMPLEX(DP) ::  mcharge0gw,mcharge1gw,mcharge2gw,mcharge3gw,mcharge4gw,mcharge5gw,mcharge6gw


 Nlzcutoff=dffts%nr3/2
    lzcutoff=Nlzcutoff*alat/dffts%nr1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get interpolated eps matrix
! eps(0,0) fine after check
    icount=0
    allocate(epsint_q0_tmp1(gw_epsq0_data%nq_data(1)))
    allocate(epsint_q0_tmp2(gw_epsq0_data%nq_data(1)))
    allocate(epsint_q0_tmp3(gw_epsq0_data%nq_data(1)))
    allocate(epsint_q0_tmp4(gw_epsq0_data%nq_data(1)))


    allocate(epsmat_inted(gw_q_g_commonsubset_size,gw_q_g_commonsubset_size))
    epsmat_inted(:,:)=0.0
    

interpolate_2d=.false.
interpolate_smallq1d=.false.
       if(abs(norm2((xk(1:3,ik0)-xk(1:3,ik))*tpiba))<tpiba*(2*3**.5/3.0)*2.0/nqgrid_gw) then
             interpolate_smallq1d=.true.
       else
             interpolate_2d=.true.
       endif

if(interpolate_2d) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2d simple interpolate prepare fixme
    deltakG=norm2(xk(1:3,ik0)-xk(1:3,ik))*tpiba
             deltakG_para=deltakG
    allocate(w1(gw_epsq1_data%nq_data(1)))
    w1(:)=0.0
    do iq1 = 1, gw_epsq1_data%nq_data(1)
        q1(:)=   gw_epsq0_data%qpts_data(1,iq1)*gw_epsq1_data%bvec_data(:,1)+ &
              gw_epsq0_data%qpts_data(2,iq1)*gw_epsq1_data%bvec_data(:,2)+ &
              gw_epsq0_data%qpts_data(3,iq1)*gw_epsq1_data%bvec_data(:,3)
       if(abs(norm2((xk(1:3,ik0)-xk(1:3,ik))*tpiba)-norm2(q1(:)))<tpiba*(2*3**.5/3.0)*1.0/nqgrid_gw) then
         w1(iq1)=1/abs(norm2((xk(1:3,ik0)-xk(1:3,ik))*tpiba)-norm2(q1(:)))
       endif
    enddo
    write(*,*) 'gw_debug w1',w1(:)
    !do iq1 = 1, gw_nq_data(1)
    !   if(abs(norm2((xk(1:3,ik0)-xk(1:3,ik))*tpiba)-norm2(q1))<machine_eps*1e-6) 
    !enddo
! 2d simple interpolate prepare fixme
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
endif



    do ig1 = 1, gw_q_g_commonsubset_size
      do ig2 = 1, gw_q_g_commonsubset_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! simple interpolate
        !epstmp1=(gw_epsmat_full_data(1,ig1,ig2,1,1,iq1),&
        !         gw_epsmat_full_data(2,ig1,ig2,1,1,iq1))
        !epstmp2=(gw_epsmat_full_data(1,ig1,ig2,1,1,iq2),&
        !         gw_epsmat_full_data(2,ig1,ig2,1,1,iq2))

        !      eps_gw=gw_epsmat_full_data(:,ig1,ig2,1,1,iq1)*wq1+&
        !       gw_epsmat_full_data(:,ig1,ig2,1,1,iq2)*wq2
! simple interpolate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


             if( interpolate_smallq1d) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!spline interpolate fixme 
!! change to 2d space matrix interpolate
!! change boundary condition

        !*gw_blat_data(1)
        !write(*,*) 'gw3.1 gind',ig1,ig2,gw_gind_rho2eps_data(1,1:5)
        !write(*,*) 'gw3.1 gind',ig1,ig2,gw_gind_rho2eps_data(2,1:5)
        !write(*,*) 'gw3.1 gind',ig1,ig2,gw_gind_rho2eps_data(3,1:5)
        !write(*,*) 'gw3.1',gw_q_g_commonsubset_indinrho(1:5)
        !write(*,*) 'gw3.1',gw_gind_rho2eps_data(1,1:5)
        !write(*,*) 'gw3.1',ig1,ig2
                do iq1=1,gw_epsq0_data%nq_data(1)
        !write(*,*) 'gw3.1.1',iq1, gind_gw_eps2,gw_q_g_commonsubset_indinrho(ig2),gw_gind_rho2eps_data(gw_q_g_commonsubset_indinrho(ig2),iq1)
                   gind_gw_eps1=gw_epsq0_data%gind_rho2eps_data(gw_epsq0_data%q_g_commonsubset_indinrho(ig1),iq1)
                   gind_gw_eps2=gw_epsq0_data%gind_rho2eps_data(gw_epsq0_data%q_g_commonsubset_indinrho(ig2),iq1)
        !write(*,*) 'gw3.1.2',iq1
                   if  (gind_gw_eps2>gw_epsq0_data%nmtx_data(iq1).or. gind_gw_eps2>gw_epsq0_data%nmtx_data(iq1)  )  &
write(*,*) 'gindex of eps qpts messedup'
        !write(*,*) 'gw3.1.3',1,gind_gw_eps1,gind_gw_eps2,1,1,iq1
        !write(*,*) 'gw3.1.3',gw_epsmat_full_data(1,gind_gw_eps1,gind_gw_eps2,1,1,iq1)
                   epsint_q0_tmp1(iq1)=gw_epsq0_data%epsmat_full_data(1,gind_gw_eps1,gind_gw_eps2,1,1,iq1)
        !write(*,*) 'gw3.1.3',iq1
                   epsint_q0_tmp2(iq1)=gw_epsq0_data%epsmat_full_data(2,gind_gw_eps1,gind_gw_eps2,1,1,iq1)
        !write(*,*) 'gw3.1.4',iq1
                enddo
        !write(*,*) 'gw3.2',ig1,ig2
                call  spline(gw_epsq0_data%qabs(:),epsint_q0_tmp1(:),0.0_DP,0.0_DP,epsint_q0_tmp3(:))
                epsinttmp1s= splint(gw_epsq0_data%qabs(:),epsint_q0_tmp1(:),epsint_q0_tmp3(:),deltakG_para)
                if (deltakG_para>maxval(gw_epsq0_data%qabs(:)))  epsinttmp1s=minval(epsint_q0_tmp1(:))
        
                call  spline(gw_epsq0_data%qabs(:),epsint_q0_tmp2(:),0.0_DP,0.0_DP,epsint_q0_tmp4(:))
                epsinttmp2s= splint(gw_epsq0_data%qabs(:),epsint_q0_tmp2(:),epsint_q0_tmp4(:),deltakG_para)
                if (deltakG_para>maxval(gw_epsq0_data%qabs(:)))  epsinttmp2s=minval(epsint_q0_tmp2(:))
        !        epsmat_inted(gw_gind_rho2eps_data(ig1,iq1),gw_gind_rho2eps_data(ig2,iq1))=complex(epsinttmp1s,epsinttmp2s)
                epsmat_inted(ig1,ig2)=complex(epsinttmp1s,epsinttmp2s)
        
!!spline interpolate fixme 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
endif

  if (interpolate_2d) then
if(sum(w1(:))<machine_eps) then
write(*,*) 'eps 2d interpolation error'
stop -1
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2d simple interpolate  fixme
                 do iq1 = 1, gw_epsq1_data%nq_data(1)
                   gind_gw_eps1=gw_epsq1_data%gind_rho2eps_data(gw_epsq1_data%q_g_commonsubset_indinrho(ig1),iq1)
                   gind_gw_eps2=gw_epsq1_data%gind_rho2eps_data(gw_epsq1_data%q_g_commonsubset_indinrho(ig2),iq1)
                   epsinttmp1s=gw_epsq1_data%epsmat_full_data(1,gind_gw_eps1,gind_gw_eps2,1,1,iq1)
                   epsinttmp2s=gw_epsq1_data%epsmat_full_data(2,gind_gw_eps1,gind_gw_eps2,1,1,iq1)
                   epsmat_inted(ig1,ig2)=epsmat_inted(ig1,ig2)+complex(epsinttmp1s,epsinttmp2s)*w1(iq1)
                 enddo
                 epsmat_inted(ig1,ig2)=epsmat_inted(ig1,ig2)/sum(w1(:))
! 2d simple interpolate  fixme
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 endif
 
        write(*,*) 'gw_debug epsmat_inted ig1,ig2,q',epsmat_inted(ig1,ig2),'ig1',ig1,'ig2',ig2,deltakG_para
               
              enddo
            enddo

! get interpolated eps matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get w(g)
!write(*,*) 'gw4'
    allocate(w_gw(ngk(ik0)))
    w_gw(:)=0.0

    qxy=norm2(xk(1:2,ik0)-xk(1:2,ik))*tpiba
    qz= (( xk(3,ik0)-xk(3,ik))**2)**0.5*tpiba
q2d_coeff=(1-(cos(qz*lzcutoff)-sin(qz*lzcutoff)*qz/qxy)*exp(-(qxy*lzcutoff)))
    DO ig1 = 1, ngk(ik0)
      Do ig2=1, ngk(ik)
           icount=icount+1
           deltakG=norm2(g(1:3,igk_k(ig1,ik0))&
                      -g(1:3,igk_k(ig2,ik))&
                      +xk(1:3,ik0)-xk(1:3,ik))*tpiba
           !w_gw(ig1)=w_gw(ig1)+epsmat_inted(gind_psi2rho_gw(ig1),gind_psi2rho_gw(ig2))*(tpi/(deltakG))
           w_gw(ig1)=w_gw(ig1)+epsmat_inted(gind_psi2rho_gw(ig1),gind_psi2rho_gw(ig2))*4*pi/(deltakG**2)*q2d_coeff
      Enddo

write(*,*) 'gw_debug W_gw vs q, ig1, g, w',ig1,norm2(g(1:3,igk_k(ig1,ik0))),w_gw(ig1) ,abs(w_gw(ig1) )
    Enddo
 
! get w(g)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!GW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 
    allocate(eps_data_dy(size(eps_data(1,:))))
    call spline(eps_data(1,:),eps_data(2,:),0.0_DP,0.0_DP,eps_data_dy(:))
    mcharge1=0
    mcharge2=0.00
    mcharge3=0.00
    icount=0
    DO ig1 = 1, ngk(ik0)
      Do ig2=1, ngk(ik)
           mcharge0=conjg(evc1(ig1,ibnd0))*evc2(ig2,ibnd)
           icount=icount+1
           deltakG=norm2(g(1:2,igk_k(ig1,ik0))&
                      -g(1:2,igk_k(ig2,ik))&
                      +xk(1:2,ik0)-xk(1:2,ik))*tpiba

             deltakG_para=deltakG
             epsk= splint(eps_data(1,:),eps_data(2,:),eps_data_dy(:),deltakG_para)
             if (deltak>maxval(eps_data(1,:)))      epsk=minval(eps_data(2,:))
           mcharge1=mcharge1+mcharge0*tpi/deltakG
           mcharge2=mcharge2+mcharge0*tpi/(deltakG**2+k0screen**2)**0.5
           mcharge3=mcharge3+mcharge0*tpi/(deltakG)*epsk
      Enddo
    Enddo
    mcharge1=mcharge1/dffts%nnr
    mcharge2=mcharge2/dffts%nnr
    mcharge3=mcharge3/dffts%nnr
    write(*,*)  'Mcharge2DnoLFAns noki->kf ',ik0,ik, mcharge1, abs(mcharge1),icount
    write(*,*)  'Mcharge2DnoLFAs  noki->kf ',ik0,ik, mcharge2, abs(mcharge2),icount , 'k0screen', k0screen
    write(*,*)  'Mcharge2DnoLFAes noki->kf ',ik0,ik, mcharge3, abs(mcharge3),icount , 'epsk', epsk
    
    
    
    mcharge1=0
    mcharge2=0
    mcharge3=0
    mcharge4=0
    mcharge5=0
    mcharge6=0
    mcharge1gw=0
    mcharge2gw=0
    mcharge3gw=0
    mcharge4gw=0
    mcharge5gw=0
    mcharge6gw=0
    DO ig1 = 1, ngk(ik0)
      Do ig2=1, ngk(ik)
    
         mcharge0=conjg(evc1(ig1,ibnd0))*evc2(ig2,ibnd)
         deltakG=norm2(g(:,igk_k(ig1,ik0))&
                    -g(:,igk_k(ig2,ik))&
                    +xk(:,ik0)-xk(:,ik))*tpiba
    

         qxy=norm2(g(1:2,igk_k(ig1,ik0))&
                    -g(1:2,igk_k(ig2,ik))&
                    +xk(1:2,ik0)-xk(1:2,ik))*tpiba
    
         qz= ((g(3,igk_k(ig1,ik0))-g(3,igk_k(ig2,ik))+ &
              xk(3,ik0)-xk(3,ik))**2)**0.5*tpiba

             epsk= splint(eps_data(1,:),eps_data(2,:),eps_data_dy(:),qxy)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! bug, should be epsk(q)=epsk(q_//)
             !epsk= splint(eps_data(1,:),eps_data(2,:),eps_data_dy(:),deltakG)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             if (deltak>maxval(eps_data(1,:)))      epsk=minval(eps_data(2,:))
q2d_coeff=(1-(cos(qz*lzcutoff)-sin(qz*lzcutoff)*qz/qxy)*exp(-(qxy*lzcutoff)))

             mcharge1=mcharge1+mcharge0*4*pi/(deltakG**2)
             mcharge2=mcharge2+mcharge0*4*pi/(deltakG**2+k0screen**2)
             mcharge3=mcharge3+mcharge0*4*pi/(deltakG**2)*epsk
! write(*,*) 'mcharge3',mcharge3,mcharge0,4*pi,(deltakG**2),epsk
             mcharge4=mcharge4+mcharge0*4*pi/(deltakG**2)            *q2d_coeff
             mcharge5=mcharge5+mcharge0*4*pi/(deltakG**2+k0screen**2)*q2d_coeff
             mcharge6=mcharge6+mcharge0*4*pi/(deltakG**2)*epsk       *q2d_coeff



             do iq = 1,ngk(ik0) 
               if (norm2(g(1:3,igk_k(iq,ik0))-(g(:,igk_k(ig1,ik0))-g(:,igk_k(ig2,ik))))<machine_eps) then
                 mcharge1gw=mcharge1gw+mcharge0*w_gw(iq)
                 mcharge2gw=mcharge2gw+mcharge0*w_gw(iq)            *q2d_coeff
                 write(*,*) 'gw_debug W in M, ig1,ig2,iq,g1,g2,q,w_gw(iq)',&
           ig1,ig2,iq,g(:,igk_k(ig1,ik0)),g(:,igk_k(ig2,ik)) ,g(1:3,igk_k(iq,ik0)),w_gw(iq) 
               endif
            
             Enddo
      Enddo
    Enddo
    mcharge1gw=mcharge1/dffts%nnr
    mcharge2gw=mcharge2/dffts%nnr
    mcharge1=mcharge1/dffts%nnr
    mcharge2=mcharge2/dffts%nnr
    mcharge3=mcharge3/dffts%nnr
    mcharge4=mcharge4/dffts%nnr
    mcharge5=mcharge5/dffts%nnr
    mcharge6=mcharge6/dffts%nnr
    write(*,*)  'Mcharge3DnoLFAgw    0ki->kf ',ik0,ik,    mcharge1gw, abs(mcharge1gw)
    write(*,*)  'Mcharge3DnoLFAns    0ki->kf ',ik0,ik,    mcharge1, abs(mcharge1)
    write(*,*)  'Mcharge3DnoLFAs     0ki->kf ',ik0,ik,    mcharge2, abs(mcharge2) , 'k0screen', k0screen
    write(*,*)  'Mcharge3DnoLFAes    0ki->kf ',ik0,ik,    mcharge3, abs(mcharge3) , 'epsk', epsk
    write(*,*)  'Mcharge3DcutnoLFAgw 0ki->kf ',ik0,ik,    mcharge2gw, abs(mcharge2gw)
    write(*,*)  'Mcharge3DcutnoLFAns 0ki->kf ',ik0,ik,    mcharge4, abs(mcharge4)
    write(*,*)  'Mcharge3DcutnoLFAs  0ki->kf ',ik0,ik,    mcharge5, abs(mcharge5) , 'k0screen', k0screen
    write(*,*)  'Mcharge3DcutnoLFAes 0ki->kf ',ik0,ik,    mcharge6, abs(mcharge6) , 'epsk', epsk
    
    END SUBROUTINE calcmdefect_charge_nolfa
    


