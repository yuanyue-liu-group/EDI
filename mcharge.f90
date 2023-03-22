SUBROUTINE calcmdefect_charge_nolfa(ibnd,ibnd0,ik,ik0,noncolin,mcharge)
  USE kinds, ONLY: DP
  Use edic_mod, Only : qeh_eps_data,eps_type
  Use edic_mod, Only : dogwfull,dogwdiag,doqeh,do2d,do3d
  Use edic_mod, Only : evc1,evc2
  Use edic_mod, Only : psic1, psic2
  USE fft_base,  ONLY: dfftp, dffts
  USE wvfct, ONLY: npwx, nbnd, wg, et, g2kin
  USE gvect, ONLY: g
  USE gvect, ONLY: ngm
  USE klist , ONLY:  xk, igk_k, ngk
  use splinelib, only: spline,splint
  USE cell_base, ONLY:  alat, tpiba,omega
  USE constants, ONLY: tpi, pi
  use edic_mod, only: machine_eps,k0screen_read

  USE fft_interfaces, ONLY : fwfft, invfft
  USE ions_base,        ONLY : nat, ntyp => nsp, ityp, tau, zv, atm
  USE cell_base,        ONLY : at, bg, omega, alat, celldm, ibrav
  USE gvecw,            ONLY : ecutwfc
  USE gvect,            ONLY : gcutm
  USE gvecs,            ONLY : dual


  Use edic_mod,   only: gw_epsq1_data,gw_epsq0_data
  USE HDF5
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! charge
  REAL(dp) ::k0screen

  REAL(dp):: kbT,deltak,deltakG0,deltakG, qxy,qz,lzcutoff
  INTEGER:: icount,jcount,kcount
  real(DP):: mscreen, rmod
  complex(DP),intent(inout):: mcharge
  INTEGER:: Nlzcutoff,iNlzcutoff,flag1,flag2, nNlzcutoff,Ngzcutoff
  integer :: nepslines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(DP),allocatable ::gind_psi2rho_gw(:)
  INTEGER :: gw_q_g_commonsubset_size
  COMPLEX(DP) ::  mcharge00,mcharge0,mcharge1,mcharge2,mcharge3,mcharge4,mcharge5,mcharge6
  COMPLEX(DP) ::  mcharge01,mcharge02,mcharge03,mcharge04,mcharge05,mcharge06,mcharge07,mcharge08
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

  REAL(dp)::arg,argt,argt2,rs(3)
  COMPLEX(DP)::phase

  write(*,*) 'Start Mcharge Calculation'
  if(dogwfull .or. dogwdiag) then
    Nlzcutoff=dffts%nr3/2
    lzcutoff=Nlzcutoff*alat/dffts%nr1
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! get interpolated eps matrix
    ! eps(0,0) fine after check
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

    call mat_inv(epsmat_inted,epsmat_inv)
   
    epsmat_lindhard(:,:)=(0.0,0.0)

    qxy=norm2(xk(1:2,ik0)-xk(1:2,ik))*tpiba
    qz= (( xk(3,ik0)-xk(3,ik))**2)**0.5*tpiba
    q2d_coeff=(1-(cos(qz*lzcutoff)-sin(qz*lzcutoff)*qz/qxy)*exp(-(qxy*lzcutoff)))

    do ig1 = 1, ngm
      do ig2 = 1, ngm
         if(gind_psi2rho_gw(ig1)>0 .and. gind_psi2rho_gw(ig2)>0)then
           if(norm2(g(1:2,ig1)-g(1:2,ig2))<machine_eps) then
               deltakG=norm2(g(1:3,ig1)&
                            -g(1:3,ig2)&
                          +xk(1:3,ik0)-xk(1:3,ik))*tpiba
               epsmat_inv(gind_psi2rho_gw(ig1),gind_psi2rho_gw(ig2))=&
               epsmat_inv(gind_psi2rho_gw(ig1),gind_psi2rho_gw(ig2))+&
                  4*pi/(deltakG**2)*q2d_coeff*k0screen/(lzcutoff*2)
           endif
         endif
      enddo
    enddo
    call mat_inv(epsmat_inv,epsmat_lindhard)
    allocate(w_gw(ngm))
    w_gw(:)=0.0

    qxy=norm2(xk(1:2,ik0)-xk(1:2,ik))*tpiba
    qz= (( xk(3,ik0)-xk(3,ik))**2)**0.5*tpiba
    q2d_coeff=(1-(cos(qz*lzcutoff)-sin(qz*lzcutoff)*qz/qxy)*exp(-(qxy*lzcutoff)))
    icount=0
    DO ig1 = 1, ngm
      Do ig2=1, ngm
        if(gind_psi2rho_gw(ig1)>0 .and. gind_psi2rho_gw(ig2)>0)then
           
           deltakG=norm2(g(1:3,ig1)-g(1:3,ig2) +xk(1:3,ik0)-xk(1:3,ik))*tpiba
           deltakG=norm2(g(1:3,ig1)*0-g(1:3,ig2) +xk(1:3,ik0)-xk(1:3,ik))*tpiba
           rs(1)=1.0/0.529* 0.5*3.168 
           rs(2)=1.0/0.529* 0.28867513459*3.168 
           rs(3)=1.0/0.529* 25*0.06265727744
           arg=(-g(1,ig2) +xk(1,ik0)-xk(1,ik))*tpiba*rs(1)&
              +(-g(2,ig2) +xk(2,ik0)-xk(2,ik))*tpiba*rs(2)&
              +(-g(3,ig2) +xk(3,ik0)-xk(3,ik))*tpiba*rs(3)
           phase=CMPLX(COS(-arg),SIN(-arg),kind=dp)
           w_gw(ig1)=w_gw(ig1)+epsmat_lindhard(gind_psi2rho_gw(ig1),gind_psi2rho_gw(ig2))*4*pi/(deltakG**2)*q2d_coeff*phase
        endif
      Enddo
      if(    abs(w_gw(ig1))>machine_eps) then
         icount=icount+1
      endif
    Enddo

    psic1(1:dffts%nnr) = (0.d0,0.d0)
    psic1 (dffts%nl (1:ngm ) ) = w_gw (1:ngm )
    CALL plot_io ('Vck.dat', 'V from w_gw invfft', dfftp%nr1, dfftp%nr2, dfftp%nr3,&
            dfftp%nr1, dfftp%nr2, dfftp%nr3, nat,    ntyp, ibrav, celldm, at,&
            gcutm, dual, ecutwfc, 0, atm, ityp, zv, tau, psic1, + 1) 
    CALL invfft ('Rho', psic1, dffts)
    CALL plot_io ('Vcr.dat', 'V from w_gw invfft', dfftp%nr1, dfftp%nr2, dfftp%nr3,&
            dfftp%nr1, dfftp%nr2, dfftp%nr3, nat,    ntyp, ibrav, celldm, at,&
            gcutm, dual, ecutwfc, 0, atm, ityp, zv, tau, psic1, + 1) 


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
          deltakG=norm2(g(1:2,igk_k(ig1,ik0))&
                      -g(1:2,igk_k(ig2,ik))&
                      +xk(1:2,ik0)-xk(1:2,ik))*tpiba
 
  endif

  if(doqeh ) then
 
    allocate(eps_data_dy(size(qeh_eps_data(1,:))))
    call spline(qeh_eps_data(1,:),qeh_eps_data(2,:),0.0_DP,0.0_DP,eps_data_dy(:))

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
    mcharge01=0
    mcharge02=0
    DO ig1 = 1, ngk(ik0)
      Do ig2=1, ngk(ik)
    
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
                 mcharge01=mcharge01+conjg(evc1(ig1,ibnd0))*evc2(ig2,ibnd) 
                 mcharge02=mcharge02+conjg(evc1(ig1+npwx,ibnd0))*evc2(ig2+npwx,ibnd) 
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

         q2d_coeff=(1-(cos(qz*lzcutoff))*exp(-(qxy*lzcutoff)))
         dg=(g(:,igk_k(ig1,ik0))-g(:,igk_k(ig2,ik)))
         if(dogwfull .or. dogwdiag)then
             !do iq = 1,ngk(ik0) 
             do iq = 1,size(w_gw_non0)
               if (norm2(g_of_w_gw_non0(:,iq)-dg)<machine_eps) then
           rs(1)=1.0/0.529* 0.5*3.168 
           rs(2)=1.0/0.529* 0.28867513459*3.168 
           rs(3)=1.0/0.529* 25*0.06265727744
           arg=(dg(1) +xk(1,ik0)-xk(1,ik))*tpiba*rs(1)&
              +(dg(2) +xk(2,ik0)-xk(2,ik))*tpiba*rs(2)&
              +(dg(3) +xk(3,ik0)-xk(3,ik))*tpiba*rs(3)
           phase=CMPLX(COS(arg),SIN(arg),kind=dp)
                 mcharge1gw=mcharge1gw+mcharge0*w_gw_non0(iq)  *phase
                 mcharge2gw=mcharge2gw+mcharge0*w_gw_non0(iq)  *phase          *q2d_coeff
               endif
             Enddo
         endif
 
      Enddo
    Enddo
    mcharge1gw=mcharge1gw/omega
    mcharge2gw=mcharge2gw/omega
    write(*,*)  'Mcharge3DnoLFAgw    0ki->kf ',ik0,ik,    mcharge1gw, abs(mcharge1gw)
    write(*,*)  'Mcharge3DcutnoLFAgw 0ki->kf ',ik0,ik,    mcharge2gw, abs(mcharge2gw)
    mcharge=mcharge2gw

  endif
  
contains
  subroutine mat_inv(A,Ainv)
    complex(dp), dimension(:,:), intent(in)    :: A
    complex(dp), dimension(:,:), intent(inout) :: Ainv
  
    complex(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info
  
    external DGETRF
    external DGETRI
  
    Ainv = A
    n = size(A,1)
  
    call DGETRF(n, n, Ainv, n, ipiv, info)
    if (info /= 0) then
       write(*,*)'info',info,Ainv(info,info),Ainv(1,1),Ainv(2,2),Ainv(3,3),Ainv(4,4),Ainv(5,5)
       write(*,*)'info',info,Ainv(info,info),Ainv(1,1),Ainv(2,1),Ainv(3,1),Ainv(4,1),Ainv(5,1)
       stop 'Matrix is numerically singular!'
    end if
  
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
  real(DP),allocatable ::qabs_tmp(:)
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

    if (allocated(epsmat_inted)) deallocate(epsmat_inted)
    allocate(epsmat_inted(gw_q_g_commonsubset_size,gw_q_g_commonsubset_size))
    epsmat_inted(:,:)=(0.0,0.0)

    deltak_para=norm2(xk(1:3,ik0)-xk(1:3,ik))*tpiba ! fixme



    if (allocated(qabs_tmp))       deallocate(qabs_tmp)
    if (allocated(epsint_q0_tmp1)) deallocate(epsint_q0_tmp1)
    if (allocated(epsint_q0_tmp2)) deallocate(epsint_q0_tmp2)
    if (allocated(epsint_q0_tmp3)) deallocate(epsint_q0_tmp3)
    if (allocated(epsint_q0_tmp4)) deallocate(epsint_q0_tmp4)
    allocate(      qabs_tmp(0:gw_epsq0_data%nq_data(1)))
    allocate(epsint_q0_tmp1(0:gw_epsq0_data%nq_data(1)))
    allocate(epsint_q0_tmp2(0:gw_epsq0_data%nq_data(1)))
    allocate(epsint_q0_tmp3(0:gw_epsq0_data%nq_data(1)))
    allocate(epsint_q0_tmp4(0:gw_epsq0_data%nq_data(1)))
    do ig1 = 1, gw_q_g_commonsubset_size
      do ig2 = 1, gw_q_g_commonsubset_size
            epsint_q0_tmp1(0)=0.0
            epsint_q0_tmp2(0)=0.0
            qabs_tmp(0)=0.0
            do iq1=1,gw_epsq0_data%nq_data(1)
               gind_gw_eps1=gw_epsq0_data%gind_rho2eps_data(gw_epsq0_data%q_g_commonsubset_indinrho(ig1),iq1)
               gind_gw_eps2=gw_epsq0_data%gind_rho2eps_data(gw_epsq0_data%q_g_commonsubset_indinrho(ig2),iq1)
               if  (gind_gw_eps1>gw_epsq0_data%nmtx_data(iq1).or. gind_gw_eps2>gw_epsq0_data%nmtx_data(iq1)  )  &
                    write(*,*) 'gindex of eps qpts messedup'
               epsint_q0_tmp1(iq1)=gw_epsq0_data%epsmat_full_data(1,gind_gw_eps1,gind_gw_eps2,1,1,iq1)
               epsint_q0_tmp2(iq1)=gw_epsq0_data%epsmat_full_data(2,gind_gw_eps1,gind_gw_eps2,1,1,iq1)
               qabs_tmp(iq1)=gw_epsq0_data%qabs(iq1)
            enddo
            if (ig1==ig2) epsint_q0_tmp1(0)=0.9
            epsint_q0_tmp1(0)=gw_epsq0_data%epsmat_full_data(1,gind_gw_eps1,gind_gw_eps2,1,1,1)
            epsint_q0_tmp2(0)=gw_epsq0_data%epsmat_full_data(2,gind_gw_eps1,gind_gw_eps2,1,1,1)


            call  spline(qabs_tmp(:),epsint_q0_tmp1(:),0.0_DP,0.0_DP,epsint_q0_tmp3(:))
            epsinttmp1s= splint(qabs_tmp(:),epsint_q0_tmp1(:),epsint_q0_tmp3(:),deltak_para)
                
            call  spline(qabs_tmp(:),epsint_q0_tmp2(:),0.0_DP,0.0_DP,epsint_q0_tmp4(:))
            epsinttmp2s= splint(qabs_tmp(:),epsint_q0_tmp2(:),epsint_q0_tmp4(:),deltak_para)
            
            epsmat_inted(ig1,ig2)=complex(epsinttmp1s,epsinttmp2s)
      enddo
    enddo
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

  complex(DP),intent(inout),allocatable ::epsmat_inted(:,:) ! interpolated eps matrix(epsilon^-1)
  INTEGER ,intent(inout):: gw_q_g_commonsubset_size
  integer(DP),intent(inout),allocatable ::gind_psi2rho_gw(:)
  integer ,intent(in):: ik0,ik
  real(dp)::q1(3)
  integer::symop(24,3,3),nsym
  nsym=24
symop(1,1,:)=(/  1 , 0 , 0 /)
symop(1,2,:)=(/  0 , 1 , 0 /)
symop(1,3,:)=(/  0 , 0 , 1 /)
       
symop(2,1,:)=(/  0 , 1 , 0 /)
symop(2,2,:)=(/ -1 ,-1 , 0 /)
symop(2,3,:)=(/  0 , 0 , 1 /)
       
symop(3,1,:)=(/ -1 , 0 , 0 /)
symop(3,2,:)=(/  0 ,-1 , 0 /)
symop(3,3,:)=(/  0 , 0 ,-1 /)
       
symop(4,1,:)=(/  0 ,-1 , 0 /)
symop(4,2,:)=(/  1 , 1 , 0 /)
symop(4,3,:)=(/  0 , 0 ,-1 /)
       
symop(5,1,:)=(/  1 , 0 , 0 /)
symop(5,2,:)=(/ -1 ,-1 , 0 /)
symop(5,3,:)=(/  0 , 0 ,-1 /)
       
symop(6,1,:)=(/  1 , 1 , 0 /)
symop(6,2,:)=(/  0 ,-1 , 0 /)
symop(6,3,:)=(/  0 , 0 ,-1 /)
       
symop(7,1,:)=(/ -1 , 0 , 0 /)
symop(7,2,:)=(/  1 , 1 , 0 /)
symop(7,3,:)=(/  0 , 0 , 1 /)
       
symop(8,1,:)=(/ -1 ,-1 , 0 /)
symop(8,2,:)=(/  0 , 1 , 0 /)
symop(8,3,:)=(/  0 , 0 , 1 /)
       
symop(9,1,:)=(/  -1,  0,  0/)  
symop(9,2,:)=(/   0, -1,  0/)  
symop(9,3,:)=(/   0,  0,  1/)  

symop(10,1,:)=(/  -1, -1,  0/)  
symop(10,2,:)=(/   1,  0,  0/)  
symop(10,3,:)=(/   0,  0,  1/)  

symop(11,1,:)=(/   1,  0,  0/)  
symop(11,2,:)=(/   0,  1,  0/)  
symop(11,3,:)=(/   0,  0, -1/)  

symop(12,1,:)=(/   1,  1,  0/)  
symop(12,2,:)=(/  -1,  0,  0/)  
symop(12,3,:)=(/   0,  0, -1/)  

symop(13,1,:)=(/   1,  1,  0/)  
symop(13,2,:)=(/  -1,  0,  0/)  
symop(13,3,:)=(/   0,  0,  1/)  

symop(14,1,:)=(/  -1, -1,  0/)  
symop(14,2,:)=(/   0,  1,  0/)  
symop(14,3,:)=(/   0,  0, -1/)  

symop(15,1,:)=(/  -1, -1,  0/)  
symop(15,2,:)=(/   1,  0,  0/)  
symop(15,3,:)=(/   0,  0, -1/)  

symop(16,1,:)=(/   1,  1,  0/)  
symop(16,2,:)=(/   0, -1,  0/)  
symop(16,3,:)=(/   0,  0,  1/)  

symop(17,1,:)=(/ -1 , 0 , 0 /) 
symop(17,2,:)=(/  1 , 1 , 0 /) 
symop(17,3,:)=(/  0 , 0 ,-1 /) 

symop(18,1,:)=(/  0 ,-1 , 0 /) 
symop(18,2,:)=(/ -1 , 0 , 0 /) 
symop(18,3,:)=(/  0 , 0 ,-1 /) 

symop(19,1,:)=(/  1 , 0 , 0 /) 
symop(19,2,:)=(/ -1 ,-1 , 0 /) 
symop(19,3,:)=(/  0 , 0 , 1 /) 

symop(20,1,:)=(/  0 , 1 , 0 /) 
symop(20,2,:)=(/  1 , 0 , 0 /) 
symop(20,3,:)=(/  0 , 0 , 1 /) 

symop(21,1,:)=(/  0 ,-1 , 0 /)             
symop(21,2,:)=(/  1 , 1 , 0 /)             
symop(21,3,:)=(/  0 , 0 , 1 /)             

symop(22,1,:)=(/  0 , 1 , 0 /)             
symop(22,2,:)=(/  1 , 0 , 0 /)             
symop(22,3,:)=(/  0 , 0 ,-1 /)             

symop(23,1,:)=(/  0 , 1 , 0 /)             
symop(23,2,:)=(/ -1 ,-1 , 0 /)             
symop(23,3,:)=(/  0 , 0 ,-1 /)             

symop(24,1,:)=(/  0 ,-1 , 0 /)             
symop(24,2,:)=(/ -1 , 0 , 0 /)             
symop(24,3,:)=(/  0 , 0 , 1 /)             

    if (allocated(gind_psi2rho_gw)) deallocate(gind_psi2rho_gw)
    allocate(gind_psi2rho_gw(size(gw_epsq1_data%gind_psi2rho)))
    gind_psi2rho_gw(:)=gw_epsq1_data%gind_psi2rho(:)
    gw_q_g_commonsubset_size=gw_epsq1_data%q_g_commonsubset_size


    if (allocated(epsmat_inted)) deallocate(epsmat_inted)
    allocate(epsmat_inted(gw_q_g_commonsubset_size,gw_q_g_commonsubset_size))
    epsmat_inted(:,:)=(0.0,0.0)
   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 2d simple interpolate prepare fixme
        allocate(w1(gw_epsq1_data%nq_data(1)))
        w1(:)=0.0
        do iq1 = 1, gw_epsq1_data%nq_data(1)
           q1(:)= gw_epsq1_data%qpts_data(1,iq1)*gw_epsq1_data%bvec_data(:,1)+ &
                  gw_epsq1_data%qpts_data(2,iq1)*gw_epsq1_data%bvec_data(:,2)+ &
                  gw_epsq1_data%qpts_data(3,iq1)*gw_epsq1_data%bvec_data(:,3)

           do ig1=1,nsym
               q1(:)= (gw_epsq1_data%qpts_data(:,iq1)*symop(ig1,1,:))*gw_epsq1_data%bvec_data(:,1)+ &
                      (gw_epsq1_data%qpts_data(:,iq1)*symop(ig1,2,:))*gw_epsq1_data%bvec_data(:,2)+ &
                      (gw_epsq1_data%qpts_data(:,iq1)*symop(ig1,3,:))*gw_epsq1_data%bvec_data(:,3)

               if(abs(norm2((xk(1:3,ik0)-xk(1:3,ik))*tpiba)-norm2(q1(:)+gw_epsq1_data%bvec_data(:,1)))<&  
                       tpiba*(2*3**.5/3.0)*8.0/nqgrid_gw .or. &
                  abs(norm2((xk(1:3,ik0)-xk(1:3,ik))*tpiba)-norm2(q1(:)+gw_epsq1_data%bvec_data(:,2)))<&  
                       tpiba*(2*3**.5/3.0)*8.0/nqgrid_gw .or. &
                  abs(norm2((xk(1:3,ik0)-xk(1:3,ik))*tpiba)-norm2(q1(:)-gw_epsq1_data%bvec_data(:,1)))<&  
                       tpiba*(2*3**.5/3.0)*8.0/nqgrid_gw .or. &
                  abs(norm2((xk(1:3,ik0)-xk(1:3,ik))*tpiba)-norm2(q1(:)-gw_epsq1_data%bvec_data(:,2)))<&  
                       tpiba*(2*3**.5/3.0)*8.0/nqgrid_gw) then
                 w1(iq1)=w1(iq1)+1/abs(norm2((xk(1:3,ik0)-xk(1:3,ik))*tpiba)-norm2(q1(:)))
               endif
           enddo
        enddo

        if(sum(w1(:))<machine_eps) then
            write(*,*) 'eps 2d interpolation error'
            stop 'eps 2d interpolation error'
        endif
    ! 2d simple interpolate prepare fixme
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    do ig1 = 1, gw_q_g_commonsubset_size
      do ig2 = 1, gw_q_g_commonsubset_size
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 2d simple interpolate  fixme
            do iq1 = 1, gw_epsq1_data%nq_data(1)

                gind_gw_eps1=gw_epsq1_data%gind_rho2eps_data(gw_epsq1_data%q_g_commonsubset_indinrho(ig1),iq1)
                gind_gw_eps2=gw_epsq1_data%gind_rho2eps_data(gw_epsq1_data%q_g_commonsubset_indinrho(ig2),iq1)

                if  (gind_gw_eps1<gw_epsq1_data%nmtx_data(iq1).and. gind_gw_eps2<gw_epsq1_data%nmtx_data(iq1)  )  then
                    epsinttmp1s=gw_epsq1_data%epsmat_full_data(1,gind_gw_eps1,gind_gw_eps2,1,1,iq1)
                    epsinttmp2s=gw_epsq1_data%epsmat_full_data(2,gind_gw_eps1,gind_gw_eps2,1,1,iq1)
                    epsmat_inted(ig1,ig2)=epsmat_inted(ig1,ig2)+complex(epsinttmp1s,epsinttmp2s)*w1(iq1)
                 endif
            enddo
            epsmat_inted(ig1,ig2)=epsmat_inted(ig1,ig2)/sum(w1(:))
        ! 2d simple interpolate  fixme
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
       
      enddo
    enddo
! get interpolated eps matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end subroutine  interp_eps_2d

END SUBROUTINE calcmdefect_charge_nolfa
