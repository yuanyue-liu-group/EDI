SUBROUTINE calcmdefect_charge_nolfa(ibnd,ibnd0,ik,ik0,noncolin,mcharge)
  USE kinds, ONLY: DP
  Use edic_mod, Only : qeh_eps_data,eps_type
  Use edic_mod, Only : dogwfull,dogwdiag,doqeh,do2d,do3d
  Use edic_mod, Only : evc1,evc2
  USE fft_base,  ONLY: dfftp, dffts
  USE gvect, ONLY: g
  USE gvect, ONLY: ngm
  USE klist , ONLY:  xk, igk_k, ngk
  use splinelib, only: spline,splint
  USE cell_base, ONLY:  alat, tpiba
  USE constants, ONLY: tpi, pi
  use edic_mod, only: machine_eps,k0screen_read

  Use edic_mod,   only: gw_epsq1_data,gw_epsq0_data
  USE HDF5
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! charge
  REAL(dp) ::k0screen
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! k0screen not initilized correctly
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !REAL(dp) ,intent(in)::k0screen
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(dp):: kbT,deltak,deltakG0,deltakG, qxy,qz,lzcutoff
  INTEGER:: icount,jcount,kcount
  real(DP):: mscreen, rmod
  complex(DP),intent(inout):: mcharge
  INTEGER:: Nlzcutoff,iNlzcutoff,flag1,flag2, nNlzcutoff,Ngzcutoff
  !!!!! eps data file 
  integer :: nepslines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(DP),allocatable ::gind_psi2rho_gw(:)
  INTEGER :: gw_q_g_commonsubset_size
  COMPLEX(DP) ::  mcharge00,mcharge0,mcharge1,mcharge2,mcharge3,mcharge4,mcharge5,mcharge6
  INTEGER :: ibnd, ik, ik0,ibnd0
  real(DP) , allocatable::  eps_data_dy(:)
  real(DP) :: epsk, deltakG_para,q2d_coeff
  logical :: noncolin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!GW
!|M|=int  u1(r) u2(r) [int eps(r,r')V(r')]
!   =int  u1(G')u2(G'+G) V(q+G)
!                        V(q+G)=int eps(q,G,G')V(q+G')
!GW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  COMPLEX(DP),allocatable ::w_gw(:)
  COMPLEX(DP),allocatable ::w_gw_non0(:)
  real(dp),allocatable ::g_of_w_gw_non0(:,:)
  COMPLEX(DP),allocatable ::w_gw_tmp1(:)
  COMPLEX(DP),allocatable ::w_gw_tmp2(:)
  complex(DP),allocatable ::epsmat_inted(:,:) ! interpolated eps matrix(epsilon^-1)
  complex(DP),allocatable ::epsmat_inv(:,:) ! inverted interpolated eps matrix(epsilon^+1)
  complex(DP),allocatable ::epsmat_lindhard(:,:) ! interpolated eps matrix(epsilon^-1)
  INTEGER :: iq1,iq2,nqgrid_gw=48!fixme
  INTEGER :: dg(3)
  real(dp)::q1(3)
  logical:: interpolate_2d,interpolate_smallq1d=.false.

  COMPLEX(DP) ::  mcharge0gw,mcharge1gw,mcharge2gw,mcharge3gw,mcharge4gw,mcharge5gw,mcharge6gw


  write(*,*) 'Start Mcharge Calculation'
  !if(eps_type=='gw')then
  if(dogwfull .or. dogwdiag) then
    Nlzcutoff=dffts%nr3/2
    lzcutoff=Nlzcutoff*alat/dffts%nr1
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! get interpolated eps matrix
    ! eps(0,0) fine after check
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    !! no need to use this judge, since q0 q1 have different common g subset 
    ! if (abs(norm2(real(gw_epsq1_data%q_g_commonsubset_indinrho-gw_epsq0_data%q_g_commonsubset_indinrho)))>machine_eps) then
    !      write(*,*) gw_epsq1_data%q_g_commonsubset_indinrho
    !      write(*,*) gw_epsq0_data%q_g_commonsubset_indinrho
    !      write(*,*) gw_epsq1_data%q_g_commonsubset_indinrho-gw_epsq0_data%q_g_commonsubset_indinrho
    !      stop ('epsq0 and epsq1 file q_g_commonsubset_indinrho not matching')
    ! endif
    !
    ! if (abs(norm2(real(gw_epsq1_data%gind_psi2rho-gw_epsq0_data%gind_psi2rho)))>machine_eps) then
    !      write(*,*) gw_epsq1_data%gind_psi2rho
    !      write(*,*) gw_epsq0_data%gind_psi2rho
    !      write(*,*) gw_epsq1_data%gind_psi2rho-gw_epsq0_data%gind_psi2rho
    !      stop ('epsq0 and epsq1 file gind_psi2rho not matching')
    ! endif
    k0screen=k0screen_read
    write(*,*) 'k0sc',k0screen

    interpolate_2d=.false.
    interpolate_smallq1d=.false.
    if(abs(norm2((xk(1:3,ik0)-xk(1:3,ik))*tpiba))<tpiba*(2*3**.5/3.0)*2.0/nqgrid_gw*0.5) then
          interpolate_smallq1d=.true.

          if (allocated(gind_psi2rho_gw))   deallocate(gind_psi2rho_gw)
          allocate(gind_psi2rho_gw(size(gw_epsq0_data%gind_psi2rho)))
          gind_psi2rho_gw(:)=gw_epsq0_data%gind_psi2rho(:)
          gw_q_g_commonsubset_size=gw_epsq0_data%q_g_commonsubset_size

          if (allocated(epsmat_inted)) deallocate(epsmat_inted)
          allocate(epsmat_inted(gw_q_g_commonsubset_size,gw_q_g_commonsubset_size))
          epsmat_inted(:,:)=(0.0,0.0)
          write(*,*) allocated(epsmat_inted)

          write(*,*) 'interp 1d, common g subset size',gw_q_g_commonsubset_size
          call interp_eps_1d(epsmat_inted,gw_q_g_commonsubset_size,ik0,ik)
    else
          interpolate_2d=.true.

          if (allocated(gind_psi2rho_gw))   deallocate(gind_psi2rho_gw)
          allocate(gind_psi2rho_gw(size(gw_epsq1_data%gind_psi2rho)))
          gind_psi2rho_gw(:)=gw_epsq1_data%gind_psi2rho(:)
          gw_q_g_commonsubset_size=gw_epsq1_data%q_g_commonsubset_size

          if (allocated(epsmat_inted)) deallocate(epsmat_inted)
          allocate(epsmat_inted(gw_q_g_commonsubset_size,gw_q_g_commonsubset_size))
          epsmat_inted(:,:)=(0.0,0.0)
          write(*,*) allocated(epsmat_inted)


          write(*,*) 'interp 2d, common g subset size',gw_q_g_commonsubset_size
          call interp_eps_2d(epsmat_inted,gw_q_g_commonsubset_size,gind_psi2rho_gw,ik0,ik)
    endif

    ! get interpolated eps matrix
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (allocated(epsmat_inv)) deallocate(epsmat_inv)
    allocate(epsmat_inv(gw_q_g_commonsubset_size,gw_q_g_commonsubset_size))
    if (allocated(epsmat_lindhard)) deallocate(epsmat_lindhard)
    allocate(epsmat_lindhard(gw_q_g_commonsubset_size,gw_q_g_commonsubset_size))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! lindhard model eps
    write(*,*) 'gw-lin1 inted'
    write(*,*) 'gw-lin1',shape(epsmat_inted),epsmat_inted(1,1)
    write(*,*) 'gw-lin1',shape(epsmat_inted),epsmat_inted(1,2)
    write(*,*) 'gw-lin1',shape(epsmat_inted),epsmat_inted(1,3)
    write(*,*) 'gw-lin1',shape(epsmat_inted),epsmat_inted(1,4)
    write(*,*) 'gw-lin1',shape(epsmat_inted),epsmat_inted(1,1)
    write(*,*) 'gw-lin1',shape(epsmat_inted),epsmat_inted(2,2)
    write(*,*) 'gw-lin1',shape(epsmat_inted),epsmat_inted(3,3)
    write(*,*) 'gw-lin1',shape(epsmat_inted),epsmat_inted(4,4)
    !call  mpi_barrier(gid)
    !call flush(6)
    !epsmat_inv(:,:)=epsmat_inted(:,:)
    call mat_inv(epsmat_inted,epsmat_inv)
    write(*,*) 'gw-lin2 inv'
    write(*,*) 'gw-lin2',shape(epsmat_inv),epsmat_inv(1,1)
    write(*,*) 'gw-lin2',shape(epsmat_inv),epsmat_inv(1,2)
    write(*,*) 'gw-lin2',shape(epsmat_inv),epsmat_inv(1,3)
    write(*,*) 'gw-lin2',shape(epsmat_inv),epsmat_inv(1,4)
    write(*,*) 'gw-lin2',shape(epsmat_inv),epsmat_inv(1,1)
    write(*,*) 'gw-lin2',shape(epsmat_inv),epsmat_inv(2,2)
    write(*,*) 'gw-lin2',shape(epsmat_inv),epsmat_inv(3,3)
    write(*,*) 'gw-lin2',shape(epsmat_inv),epsmat_inv(4,4)
    !call  mpi_barrier(gid)
    !call flush(6)
    !epsmat_lindhard(:,:)=epsmat_inv(:,:)

    !deltakG=norm2(g(:,igk_k(ig1,ik0))&
    !             -g(:,igk_k(ig2,ik))&
    !           +xk(:,ik0)-xk(:,ik))*tpiba
    
    epsmat_lindhard(:,:)=(0.0,0.0)

    qxy=norm2(xk(1:2,ik0)-xk(1:2,ik))*tpiba
    qz= (( xk(3,ik0)-xk(3,ik))**2)**0.5*tpiba
    q2d_coeff=(1-(cos(qz*lzcutoff)-sin(qz*lzcutoff)*qz/qxy)*exp(-(qxy*lzcutoff)))

    do ig1 = 1, ngm
      do ig2 = 1, ngm
    !DO ig1 = 1, ngk(ik0)
    !  Do ig2=1, ngk(ik)
    !       icount=icount+1
         if(gind_psi2rho_gw(ig1)>0 .and. gind_psi2rho_gw(ig2)>0)then
           if(norm2(g(1:2,ig1)-g(1:2,ig2))<machine_eps) then
               deltakG=norm2(g(1:3,ig1)&
                            -g(1:3,ig2)&
                          +xk(1:3,ik0)-xk(1:3,ik))*tpiba
               write(*,*) deltakG,4*pi/(deltakG**2)*q2d_coeff*k0screen/(lzcutoff*2)
               write(*,*) 4*pi,(deltakG**2),q2d_coeff,k0screen,(lzcutoff*2)
               epsmat_inv(gind_psi2rho_gw(ig1),gind_psi2rho_gw(ig2))=&
               epsmat_inv(gind_psi2rho_gw(ig1),gind_psi2rho_gw(ig2))+&
                  4*pi/(deltakG**2)*q2d_coeff*k0screen/(lzcutoff*2)
           endif
         endif
      enddo
    enddo
    write(*,*) 'gw-lin3 inv lin'
    write(*,*) 'gw-lin3',shape(epsmat_inv),epsmat_inv(1,1)
    write(*,*) 'gw-lin3',shape(epsmat_inv),epsmat_inv(1,2)
    write(*,*) 'gw-lin3',shape(epsmat_inv),epsmat_inv(1,3)
    write(*,*) 'gw-lin3',shape(epsmat_inv),epsmat_inv(1,4)
    write(*,*) 'gw-lin3',shape(epsmat_inv),epsmat_inv(1,1)
    write(*,*) 'gw-lin3',shape(epsmat_inv),epsmat_inv(2,2)
    write(*,*) 'gw-lin3',shape(epsmat_inv),epsmat_inv(3,3)
    write(*,*) 'gw-lin3',shape(epsmat_inv),epsmat_inv(4,4)
    call mat_inv(epsmat_inv,epsmat_lindhard)
    write(*,*) 'gw-lin4 lin'
    write(*,*) 'gw-lin4',shape(epsmat_lindhard),epsmat_lindhard(1,1)
    write(*,*) 'gw-lin4',shape(epsmat_lindhard),epsmat_lindhard(1,2)
    write(*,*) 'gw-lin4',shape(epsmat_lindhard),epsmat_lindhard(1,3)
    write(*,*) 'gw-lin4',shape(epsmat_lindhard),epsmat_lindhard(1,4)
    write(*,*) 'gw-lin4',shape(epsmat_lindhard),epsmat_lindhard(1,1)
    write(*,*) 'gw-lin4',shape(epsmat_lindhard),epsmat_lindhard(2,2)
    write(*,*) 'gw-lin4',shape(epsmat_lindhard),epsmat_lindhard(3,3)
    write(*,*) 'gw-lin4',shape(epsmat_lindhard),epsmat_lindhard(4,4)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get w(g)
    !write(*,*) 'gw4'
    write(*,*) 'gw-lin5 w'

!    allocate(w_gw(ngk(ik0)))
    allocate(w_gw(ngm))
    write(*,*)shape(w_gw)
    w_gw(:)=0.0

    qxy=norm2(xk(1:2,ik0)-xk(1:2,ik))*tpiba
    qz= (( xk(3,ik0)-xk(3,ik))**2)**0.5*tpiba
    q2d_coeff=(1-(cos(qz*lzcutoff)-sin(qz*lzcutoff)*qz/qxy)*exp(-(qxy*lzcutoff)))
    write(*,*) 'epsmat_lindhard(gind_psi2rho_gw(ig1),gind_psi2rho_gw(ig2))',maxval(gind_psi2rho_gw(:))
!    DO ig1 = 1, ngk(ik0)
!      Do ig2=1, ngk(ik)
    icount=0
    DO ig1 = 1, ngm
      Do ig2=1, ngm
        if(gind_psi2rho_gw(ig1)>0 .and. gind_psi2rho_gw(ig2)>0)then
           
           !write(*,*) 'deltakG',g(1:3,igk_k(ig1,ik0)),g(1:3,igk_k(ig2,ik)),xk(1:3,ik0)-xk(1:3,ik)
           deltakG=norm2(g(1:3,ig1)-g(1:3,ig2) +xk(1:3,ik0)-xk(1:3,ik))*tpiba
           deltakG=norm2(g(1:3,ig1)*0-g(1:3,ig2) +xk(1:3,ik0)-xk(1:3,ik))*tpiba
           !w_gw(ig1)=w_gw(ig1)+epsmat_inted(gind_psi2rho_gw(ig1),gind_psi2rho_gw(ig2))*(tpi/(deltakG))
!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!
! start here assign gind_psi2rho
!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!

           !if(deltakG>machine_eps)then
           !write(*,*) 'epsmat_lindhard(gind_psi2rho_gw(ig1),gind_psi2rho_gw(ig2))',  ig1,ig2
           !write(*,*) '2',w_gw(ig1),gind_psi2rho_gw(ig1),gind_psi2rho_gw(ig2)
           !write(*,*) '2',epsmat_lindhard(gind_psi2rho_gw(ig1),gind_psi2rho_gw(ig2))
           !write(*,*) '2',4*pi/(deltakG**2)*q2d_coeff
           !    w_gw(ig1)=w_gw(ig1)+epsmat_lindhard(gind_psi2rho_gw(ig1),gind_psi2rho_gw(ig2))*4*pi/(deltakG**2)*q2d_coeff
           !else
 
           !    w_gw(ig1)=w_gw(ig1)+1.0/machine_eps
           w_gw(ig1)=w_gw(ig1)+epsmat_lindhard(gind_psi2rho_gw(ig1),gind_psi2rho_gw(ig2))*4*pi/(deltakG**2)*q2d_coeff/(lzcutoff*2)
           !endif
         write(*,*) 'gw_debug W_gw vs q:dk, ig1, g, w',xk(:,ik0)-xk(:,ik),ig1,g(:,ig1),deltakG,w_gw(ig1) ,abs(w_gw(ig1) )
        endif
      Enddo
      if(    abs(w_gw(ig1))>machine_eps) then
         icount=icount+1
         !write(*,*) 'gw_debug W_gw vs q:dk, ig1, g, w',xk(:,ik0)-xk(:,ik),ig1,g(:,ig1),w_gw(ig1) ,abs(w_gw(ig1) )
      endif
    Enddo
    write(*,*) 'w nonzero part number', icount

    if (allocated(w_gw_non0)) deallocate(w_gw_non0)
    allocate(w_gw_non0(icount))
    w_gw_non0(:)=(0.0,0.0)

    if (allocated(g_of_w_gw_non0)) deallocate(g_of_w_gw_non0)
    allocate(g_of_w_gw_non0(3,icount))
    g_of_w_gw_non0(:,:)=0.0

    icount=0
    DO ig1 = 1, ngm
      if(    abs(w_gw(ig1))>machine_eps) then
         icount=icount+1
         w_gw_non0(icount)=w_gw(ig1)
         g_of_w_gw_non0(:,icount)=g(:,ig1)
      endif
    Enddo
    write(*,*)'nonzero w',w_gw_non0
    write(*,*)'g_of_w_gw_non0',g_of_w_gw_non0
 
 
! get w(g)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!GW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  endif







  if(dogwdiag ) then
  endif



  !elseif(eps_type=='qeh')then
  if(doqeh ) then
 
    allocate(eps_data_dy(size(qeh_eps_data(1,:))))
    call spline(qeh_eps_data(1,:),qeh_eps_data(2,:),0.0_DP,0.0_DP,eps_data_dy(:))

  endif

  !else
  !  stop ('eps_type incorrect')
  !endif

  if(do2d ) then
    mcharge1=0
    mcharge2=0.00
    mcharge3=0.00
    icount=0
    mcharge0=0
    DO ig1 = 1, ngk(ik0)
      Do ig2=1, ngk(ik)
           !mcharge0=conjg(evc1(ig1,ibnd0))*evc2(ig2,ibnd)
           if (.not. noncolin )then
              mcharge0=conjg(evc1(ig1,ibnd0))*evc2(ig2,ibnd)
           else
              mcharge0=conjg(evc1(ig1,ibnd0))*evc2(ig2,ibnd) &
                               +conjg(evc1(ig1+npwx,ibnd0))*evc2(ig2+npwx,ibnd)
           endif
           icount=icount+1
           deltakG=norm2(g(1:2,igk_k(ig1,ik0))&
                      -g(1:2,igk_k(ig2,ik))&
                      +xk(1:2,ik0)-xk(1:2,ik))*tpiba

           deltakG_para=deltakG
           !if(eps_type=='qeh')then
           !  epsk= splint(qeh_eps_data(1,:),qeh_eps_data(2,:),eps_data_dy(:),deltakG_para)
           !  if (deltak>maxval(qeh_eps_data(1,:)))      epsk=minval(qeh_eps_data(2,:))
           !  mcharge1=mcharge1+mcharge0*tpi/(deltakG)*epsk
           !elseif(eps_type=='tf')then
           !  mcharge1=mcharge1+mcharge0*tpi/(deltakG**2+k0screen**2)**0.5
           !endif
        
           mcharge1=mcharge1+mcharge0*tpi/deltakG
           mcharge2=mcharge2+mcharge0*tpi/(deltakG**2+k0screen**2)**0.5
           if (doqeh)  then
             epsk= splint(qeh_eps_data(1,:),qeh_eps_data(2,:),eps_data_dy(:),deltakG_para)
             if (deltak>maxval(qeh_eps_data(1,:)))      epsk=minval(qeh_eps_data(2,:))
             mcharge3=mcharge3+mcharge0*tpi/(deltakG)*epsk
           endif
      Enddo
      !write(*,*)  'mcharge ig1',ig1
    Enddo
    !mcharge1=mcharge1/dffts%nnr
    !mcharge2=mcharge2/dffts%nnr
    !mcharge3=mcharge3/dffts%nnr
    write(*,*)  'mcharge start ',ik0,ik, mcharge0, abs(mcharge0),icount
    write(*,*)  'Mcharge2DnoLFAns noki->kf ',ik0,ik, mcharge1, abs(mcharge1),icount
    write(*,*)  'Mcharge2DnoLFAs  noki->kf ',ik0,ik, mcharge2, abs(mcharge2),icount , 'k0screen', k0screen
    write(*,*)  'Mcharge2DnoLFAes noki->kf ',ik0,ik, mcharge3, abs(mcharge3),icount , 'epsk', epsk
    !write(*,*)  '1Mcharge2DnoLFAes ki->kf ',ik0,ik, mcharge1, abs(mcharge1),icount , 'epsk', epsk
  
  endif

  if(do3d ) then
    
    
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

    mcharge0=0
    mcharge00=0
    DO ig1 = 1, ngk(ik0)
      Do ig2=1, ngk(ik)
    
         !mcharge0=conjg(evc1(ig1,ibnd0))*evc2(ig2,ibnd)
         if (.not. noncolin )then
            mcharge0=conjg(evc1(ig1,ibnd0))*evc2(ig2,ibnd)
            if (norm2(g(:,igk_k(ig1,ik0))-g(:,igk_k(ig2,ik)))<machine_eps) then
                mcharge00=mcharge00+conjg(evc1(ig1,ibnd0))*evc2(ig2,ibnd)
            endif
         else
            mcharge0=conjg(evc1(ig1,ibnd0))*evc2(ig2,ibnd) &
                             +conjg(evc1(ig1+npwx,ibnd0))*evc2(ig2+npwx,ibnd)
            if (norm2(g(:,igk_k(ig1,ik0))-g(:,igk_k(ig2,ik)))<machine_eps) then
                 mcharge00=mcharge00+conjg(evc1(ig1,ibnd0))*evc2(ig2,ibnd) &
                             +conjg(evc1(ig1+npwx,ibnd0))*evc2(ig2+npwx,ibnd)
            endif
         endif

         deltakG=norm2(g(:,igk_k(ig1,ik0))&
                      -g(:,igk_k(ig2,ik))&
                  +xk(:,ik0)-xk(:,ik))*tpiba
    

         qxy=norm2(g(1:2,igk_k(ig1,ik0))&
                    -g(1:2,igk_k(ig2,ik))&
                    +xk(1:2,ik0)-xk(1:2,ik))*tpiba
    
         qz= ((g(3,igk_k(ig1,ik0))-g(3,igk_k(ig2,ik))+ &
              xk(3,ik0)-xk(3,ik))**2)**0.5*tpiba


         !deltakG=norm2(&
         !    (g(1,igk_k(ig1,ik0))-g(1,igk_k(ig2,ik))+xk(1,ik0)-xk(1,ik))*gw_epsq0_data%bvec_data(:,1)+&
         !    (g(2,igk_k(ig1,ik0))-g(2,igk_k(ig2,ik))+xk(2,ik0)-xk(2,ik))*gw_epsq0_data%bvec_data(:,2)+&
         !    (g(3,igk_k(ig1,ik0))-g(3,igk_k(ig3,ik))+xk(3,ik0)-xk(3,ik))*gw_epsq0_data%bvec_data(:,3)&
         !           )*tpiba
    

         !qxy=norm2(&
         !    (g(1,igk_k(ig1,ik0))-g(1,igk_k(ig2,ik))+xk(1,ik0)-xk(1,ik))*gw_epsq0_data%bvec_data(:,1)+&
         !    (g(2,igk_k(ig1,ik0))-g(2,igk_k(ig2,ik))+xk(2,ik0)-xk(2,ik))*gw_epsq0_data%bvec_data(:,2)&
         !           )*tpiba
         !qz=norm2(&
         !    (g(3,igk_k(ig1,ik0))-g(3,igk_k(ig3,ik))+xk(3,ik0)-xk(3,ik))*gw_epsq0_data%bvec_data(:,3)&
         !           )*tpiba



         !q2d_coeff=(1-(cos(qz*lzcutoff)-sin(qz*lzcutoff)*qz/qxy)*exp(-(qxy*lzcutoff)))
         q2d_coeff=(1-(cos(qz*lzcutoff))*exp(-(qxy*lzcutoff)))

         !if(eps_type=='qeh')then
         if(doqeh)then
             epsk= splint(qeh_eps_data(1,:),qeh_eps_data(2,:),eps_data_dy(:),qxy)
 !  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !   bug, should be epsk(q)=epsk(q_//)
               !epsk= splint(qeh_eps_data(1,:),qeh_eps_data(2,:),eps_data_dy(:),deltakG)
 !  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             if (deltak>maxval(qeh_eps_data(1,:)))      epsk=minval(qeh_eps_data(2,:))

             mcharge3=mcharge3+mcharge0*4*pi/(deltakG**2)*epsk
             !write(*,*) 'mcharge3',mcharge3,mcharge0,4*pi,(deltakG**2),epsk
             mcharge1=mcharge1+mcharge0*4*pi/(deltakG**2)
             mcharge2=mcharge2+mcharge0*4*pi/(deltakG**2+k0screen**2)
             mcharge4=mcharge4+mcharge0*4*pi/(deltakG**2)            *q2d_coeff
             mcharge5=mcharge5+mcharge0*4*pi/(deltakG**2+k0screen**2)*q2d_coeff
             mcharge6=mcharge6+mcharge0*4*pi/(deltakG**2)*epsk       *q2d_coeff
         endif



!         if(dogwfull)then
!         !elseif(eps_type=='gw')then
!             do iq = 1,ngk(ik0) 
!               if (norm2(g(1:3,igk_k(iq,ik0))-(g(:,igk_k(ig1,ik0))-g(:,igk_k(ig2,ik))))<machine_eps) then
!                 mcharge1gw=mcharge1gw+mcharge0*w_gw(iq)
!                 mcharge2gw=mcharge2gw+mcharge0*w_gw(iq)            *q2d_coeff
!                 !write(*,*) 'gw_debug W in M, ig1,ig2,iq,g1,g2,q,w_gw(iq)',&
!                 !          ig1,ig2,iq,g(:,igk_k(ig1,ik0)),g(:,igk_k(ig2,ik)) ,g(1:3,igk_k(iq,ik0)),w_gw(iq) 
!               endif
!            
!             Enddo
!         endif
         mcharge1=mcharge1+mcharge0*4*pi/(deltakG**2)
         mcharge2=mcharge2+mcharge0*4*pi/(deltakG**2+k0screen**2)
         mcharge4=mcharge4+mcharge0*4*pi/(deltakG**2)            *q2d_coeff
         mcharge5=mcharge5+mcharge0*4*pi/(deltakG**2+k0screen**2)*q2d_coeff

         dg=(g(:,igk_k(ig1,ik0))-g(:,igk_k(ig2,ik)))
         if(dogwfull .or. dogwdiag)then
             !do iq = 1,ngk(ik0) 
             do iq = 1,size(w_gw_non0)
               if (norm2(g_of_w_gw_non0(:,iq)-dg)<machine_eps) then
                 mcharge1gw=mcharge1gw+mcharge0*w_gw_non0(iq)
                 mcharge2gw=mcharge2gw+mcharge0*w_gw_non0(iq)            *q2d_coeff
!                 write(*,*) 'gw_debug W in M, ig1,ig2,iq,g1,g2,q,w_gw(iq)',&
!                           ig1,ig2,iq,g(:,igk_k(ig1,ik0)),g(:,igk_k(ig2,ik)) ,g(1:3,igk_k(iq,ik0)),w_gw_non0(iq) ,mcharge0,mcharge1gw,mcharge2gw
               endif
             Enddo
         endif
 
      Enddo
      write(*,*) 'mcharge gw_debug W in M, ig1,ig2,iq,g1,g2,q,w_gw(iq)',ig1,mcharge1gw,mcharge2gw
    Enddo
    !mcharge1gw=mcharge1gw/dffts%nnr
    !mcharge2gw=mcharge2gw/dffts%nnr
    !mcharge1=mcharge1/dffts%nnr
    !mcharge2=mcharge2/dffts%nnr
    !mcharge3=mcharge3/dffts%nnr
    !mcharge4=mcharge4/dffts%nnr
    !mcharge5=mcharge5/dffts%nnr
    !mcharge6=mcharge6/dffts%nnr
    write(*,*)  'Mcharge0 ',ik0,ik,    mcharge00, abs(mcharge00)
    write(*,*)  'Mcharge3DnoLFAgw    0ki->kf ',ik0,ik,    mcharge1gw, abs(mcharge1gw)
    write(*,*)  'Mcharge3DnoLFAns    0ki->kf ',ik0,ik,    mcharge1, abs(mcharge1)
    write(*,*)  'Mcharge3DnoLFAs     0ki->kf ',ik0,ik,    mcharge2, abs(mcharge2) , 'k0screen', k0screen
    write(*,*)  'Mcharge3DnoLFAes    0ki->kf ',ik0,ik,    mcharge3, abs(mcharge3) , 'epsk', epsk
    write(*,*)  'Mcharge3DcutnoLFAgw 0ki->kf ',ik0,ik,    mcharge2gw, abs(mcharge2gw)
    write(*,*)  'Mcharge3DcutnoLFAns 0ki->kf ',ik0,ik,    mcharge4, abs(mcharge4)
    write(*,*)  'Mcharge3DcutnoLFAs  0ki->kf ',ik0,ik,    mcharge5, abs(mcharge5) , 'k0screen', k0screen
    write(*,*)  'Mcharge3DcutnoLFAes 0ki->kf ',ik0,ik,    mcharge6, abs(mcharge6) , 'epsk', epsk
    mcharge=mcharge2gw

  endif
  
contains
  subroutine mat_inv(A,Ainv)
    complex(dp), dimension(:,:), intent(in)    :: A
    complex(dp), dimension(:,:), intent(inout) :: Ainv
  
    complex(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info
  
    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI
  
    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    !write(*,*) 'inv 1',info
    Ainv = A
    n = size(A,1)
  
    !DGETRF computes an LU factorization of a general M-by-N matrix A
    !using partial pivoting with row interchanges.
    !write(*,*) 'inv 1',info
    !call  mpi_barrier(gid)
    call flush(6)
    call DGETRF(n, n, Ainv, n, ipiv, info)
    !write(*,*) 'inv 1',info
  
    !call  mpi_barrier(gid)
    call flush(6)
    if (info /= 0) then
       stop 'Matrix is numerically singular!'
    end if
  
    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call DGETRI(n, Ainv, n, ipiv, work, n, info)
  
    if (info /= 0) then
       stop 'Matrix inversion failed!'
    end if
  end subroutine mat_inv

subroutine interp_eps_1d(epsmat_inted,gw_q_g_commonsubset_size,ik0,ik)
  USE kinds, ONLY: DP
  Use edic_mod,   only: gw_epsq1_data,gw_epsq0_data
  USE klist , ONLY:  xk
  use splinelib, only: spline,splint
  real(DP),allocatable ::epsint_q0_tmp1(:)
  real(DP),allocatable ::epsint_q0_tmp2(:)
  real(DP),allocatable ::epsint_q0_tmp3(:)
  real(DP),allocatable ::epsint_q0_tmp4(:)
  real(DP) ::epsinttmp1s
  real(DP) ::epsinttmp2s
  real(DP) ::epsinttmp3s
  real(DP) ::epsinttmp4s
  COMPLEX(DP) :: epstmp1,epstmp2
  INTEGER :: gind_gw_eps1,gind_gw_eps2
  real(DP) ::  deltak_para

  complex(DP),intent(inout),allocatable ::epsmat_inted(:,:) ! interpolated eps matrix(epsilon^-1)
  INTEGER ,intent(in):: gw_q_g_commonsubset_size
  !integer(DP),intent(inout),allocatable ::gind_psi2rho_gw(:)
  integer ,intent(in):: ik0,ik


    write(*,*) 'enter interp1d'
    !write(*,*) allocated(gind_psi2rho_gw)
    write(*,*) allocated(epsmat_inted)
!    if (allocated(gind_psi2rho_gw)) then
!        write(*,*) 'enter interp1d'
!    else
!        write(*,*) 'enter interp1d'
!    endif
!
!    deallocate(gind_psi2rho_gw)
    !if (allocated(gind_psi2rho_gw))   deallocate(gind_psi2rho_gw)
    !allocate(gind_psi2rho_gw(size(gw_epsq0_data%gind_psi2rho)))
    !gind_psi2rho_gw(:)=gw_epsq0_data%gind_psi2rho(:)
    !write(*,*) 'gind_psi2rho_gw',gind_psi2rho_gw(1:10)


    if (allocated(epsmat_inted)) deallocate(epsmat_inted)
    allocate(epsmat_inted(gw_q_g_commonsubset_size,gw_q_g_commonsubset_size))
    epsmat_inted(:,:)=(0.0,0.0)

    deltak_para=norm2(xk(1:3,ik0)-xk(1:3,ik))*tpiba ! fixme



    if (allocated(epsint_q0_tmp1)) deallocate(epsint_q0_tmp1)
    if (allocated(epsint_q0_tmp2)) deallocate(epsint_q0_tmp2)
    if (allocated(epsint_q0_tmp3)) deallocate(epsint_q0_tmp3)
    if (allocated(epsint_q0_tmp4)) deallocate(epsint_q0_tmp4)
    allocate(epsint_q0_tmp1(gw_epsq0_data%nq_data(1)))
    allocate(epsint_q0_tmp2(gw_epsq0_data%nq_data(1)))
    allocate(epsint_q0_tmp3(gw_epsq0_data%nq_data(1)))
    allocate(epsint_q0_tmp4(gw_epsq0_data%nq_data(1)))

    write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(1,1)
    write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(1,2)
    write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(1,3)
    write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(1,4)
    write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(1,1)
    write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(2,2)
    write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(3,3)
    write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(4,4)
 
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


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!spline interpolate fixme 
        !! change to 2d space matrix interpolate
        !! change boundary condition
        !if( interpolate_smallq1d) then

            !icount=0


            !write(*,*) 'interp 1d'
        
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
               if  (gind_gw_eps1>gw_epsq0_data%nmtx_data(iq1).or. gind_gw_eps2>gw_epsq0_data%nmtx_data(iq1)  )  &
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
            epsinttmp1s= splint(gw_epsq0_data%qabs(:),epsint_q0_tmp1(:),epsint_q0_tmp3(:),deltak_para)
            !write(*,*) 'gw3.2',epsinttmp1s
            if (deltak_para>maxval(gw_epsq0_data%qabs(:)))  epsinttmp1s=minval(epsint_q0_tmp1(:))
            if (deltak_para<minval(gw_epsq0_data%qabs(:)))  epsinttmp1s=epsint_q0_tmp2(size(epsint_q0_tmp2))
                
            call  spline(gw_epsq0_data%qabs(:),epsint_q0_tmp2(:),0.0_DP,0.0_DP,epsint_q0_tmp4(:))
            epsinttmp2s= splint(gw_epsq0_data%qabs(:),epsint_q0_tmp2(:),epsint_q0_tmp4(:),deltak_para)
            !write(*,*) 'gw3.2',epsinttmp2s
            if (deltak_para>maxval(gw_epsq0_data%qabs(:)))  epsinttmp2s=minval(epsint_q0_tmp2(:))
            if (deltak_para<minval(gw_epsq0_data%qabs(:)))  epsinttmp2s=epsint_q0_tmp2(size(epsint_q0_tmp2))
            !epsmat_inted(gw_gind_rho2eps_data(ig1,iq1),gw_gind_rho2eps_data(ig2,iq1))=complex(epsinttmp1s,epsinttmp2s)
            epsmat_inted(ig1,ig2)=complex(epsinttmp1s,epsinttmp2s)
            !write(*,*) 'gw3.2',ig1,ig2,epsinttmp1s,epsinttmp2s,epsmat_inted(ig1,ig2)
                
        !endif
        !!spline interpolate fixme 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      enddo
    enddo

    !write(*,*) 'gw-lin1 inted'
    !write(*,*) 'gw-lin1',shape(epsmat_inted),epsmat_inted(1,1)
    !write(*,*) 'gw-lin1',shape(epsmat_inted),epsmat_inted(1,2)
    !write(*,*) 'gw-lin1',shape(epsmat_inted),epsmat_inted(1,3)
    !write(*,*) 'gw-lin1',shape(epsmat_inted),epsmat_inted(1,4)
    !write(*,*) 'gw-lin1',shape(epsmat_inted),epsmat_inted(1,1)
    !write(*,*) 'gw-lin1',shape(epsmat_inted),epsmat_inted(2,2)
    !write(*,*) 'gw-lin1',shape(epsmat_inted),epsmat_inted(3,3)
    !write(*,*) 'gw-lin1',shape(epsmat_inted),epsmat_inted(4,4)

end subroutine  interp_eps_1d


subroutine interp_eps_2d(epsmat_inted,gw_q_g_commonsubset_size,gind_psi2rho_gw,ik0,ik)
  USE kinds, ONLY: DP
  Use edic_mod,   only: gw_epsq1_data,gw_epsq0_data
  USE klist , ONLY:  xk
  real(DP) ::epsinttmp1s
  real(DP) ::epsinttmp2s
  real(DP) ::epsinttmp3s
  real(DP) ::epsinttmp4s
  real(DP),allocatable ::w1(:)
  INTEGER :: gind_gw_eps1,gind_gw_eps2
  !real(DP) ::  deltak_para

  complex(DP),intent(inout),allocatable ::epsmat_inted(:,:) ! interpolated eps matrix(epsilon^-1)
  INTEGER ,intent(inout):: gw_q_g_commonsubset_size
  integer(DP),intent(inout),allocatable ::gind_psi2rho_gw(:)
  integer ,intent(in):: ik0,ik
  real(dp)::q1(3)


    if (allocated(gind_psi2rho_gw)) deallocate(gind_psi2rho_gw)
    allocate(gind_psi2rho_gw(size(gw_epsq1_data%gind_psi2rho)))
    gind_psi2rho_gw(:)=gw_epsq1_data%gind_psi2rho(:)
    gw_q_g_commonsubset_size=gw_epsq1_data%q_g_commonsubset_size


    if (allocated(epsmat_inted)) deallocate(epsmat_inted)
    allocate(epsmat_inted(gw_q_g_commonsubset_size,gw_q_g_commonsubset_size))
    epsmat_inted(:,:)=(0.0,0.0)

    !deltak_para=norm2(xk(1:3,ik0)-xk(1:3,ik))*tpiba ! fixme

    write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(1,1)
    write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(1,2)
    write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(1,3)
    write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(1,4)
    write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(1,1)
    write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(2,2)
    write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(3,3)
    write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(4,4)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 2d simple interpolate prepare fixme
    !if(interpolate_2d) then
        !deltakG=norm2(xk(1:3,ik0)-xk(1:3,ik))*tpiba
        !deltakG_para=deltakG
        allocate(w1(gw_epsq1_data%nq_data(1)))
        w1(:)=0.0
        do iq1 = 1, gw_epsq1_data%nq_data(1)
           q1(:)= gw_epsq1_data%qpts_data(1,iq1)*gw_epsq1_data%bvec_data(:,1)+ &
                  gw_epsq1_data%qpts_data(2,iq1)*gw_epsq1_data%bvec_data(:,2)+ &
                  gw_epsq1_data%qpts_data(3,iq1)*gw_epsq1_data%bvec_data(:,3)
           write(*,*)'q1 dk',(xk(1:3,ik0)-xk(1:3,ik)),q1
           if(abs(norm2((xk(1:3,ik0)-xk(1:3,ik))*tpiba)-norm2(q1(:)))<tpiba*(2*3**.5/3.0)*8.0/nqgrid_gw) then
             w1(iq1)=1/abs(norm2((xk(1:3,ik0)-xk(1:3,ik))*tpiba)-norm2(q1(:)))
           endif
        enddo

        write(*,*) 'gw_debug w1',w1(:)

        if(sum(w1(:))<machine_eps) then
            write(*,*) 'eps 2d interpolation error'
            stop -1
        endif
        !do iq1 = 1, gw_nq_data(1)
        !   if(abs(norm2((xk(1:3,ik0)-xk(1:3,ik))*tpiba)-norm2(q1))<machine_eps*1e-6) 
        !enddo
    !endif
    ! 2d simple interpolate prepare fixme
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do ig1 = 1, gw_q_g_commonsubset_size
      do ig2 = 1, gw_q_g_commonsubset_size
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 2d simple interpolate  fixme
        !if (interpolate_2d) then
            !write(*,*) 'interp 2d'
            do iq1 = 1, gw_epsq1_data%nq_data(1)
                gind_gw_eps1=gw_epsq1_data%gind_rho2eps_data(gw_epsq1_data%q_g_commonsubset_indinrho(ig1),iq1)
                gind_gw_eps2=gw_epsq1_data%gind_rho2eps_data(gw_epsq1_data%q_g_commonsubset_indinrho(ig2),iq1)
                epsinttmp1s=gw_epsq1_data%epsmat_full_data(1,gind_gw_eps1,gind_gw_eps2,1,1,iq1)
                epsinttmp2s=gw_epsq1_data%epsmat_full_data(2,gind_gw_eps1,gind_gw_eps2,1,1,iq1)
                epsmat_inted(ig1,ig2)=epsmat_inted(ig1,ig2)+complex(epsinttmp1s,epsinttmp2s)*w1(iq1)
            enddo
            epsmat_inted(ig1,ig2)=epsmat_inted(ig1,ig2)/sum(w1(:))
        !endif
        ! 2d simple interpolate  fixme
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
        !write(*,*) 'gw_debug epsmat_inted ig1,ig2,q',epsmat_inted(ig1,ig2),'ig1',ig1,'ig2',ig2,deltakG_para
       
      enddo
    enddo
! get interpolated eps matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end subroutine  interp_eps_2d

END SUBROUTINE calcmdefect_charge_nolfa
