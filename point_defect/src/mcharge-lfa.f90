SUBROUTINE calcmdefect_charge_lfa(ibnd,ibnd0,ik,ik0,noncolin,mcharge)
  USE kinds, ONLY: DP
  Use edic_mod, Only : evc1,evc2,qeh_eps_data
  USE fft_base,  ONLY: dfftp, dffts
  USE gvect, ONLY: g
  USE klist , ONLY:  xk, igk_k, ngk
  use splinelib, only: spline,splint
  USE cell_base, ONLY:  alat, tpiba
  USE constants, ONLY: tpi, pi
  use edic_mod, only: machine_eps
  USE cell_base, ONLY:  alat, tpiba,omega
  COMPLEX(DP) ::  mcharge0,mcharge1,mcharge2,mcharge3,mcharge4,mcharge5,mcharge6
  COMPLEX(DP) ,intent(inout)::  mcharge
  INTEGER :: ibnd, ik, ik0,ibnd0

  REAL(dp) ::k0screen
  REAL(dp):: kbT,deltak,deltakG0,deltakG, qxy,qz,lzcutoff
  INTEGER:: icount,jcount,kcount
  real(DP):: mscreen, rmod
  INTEGER:: Nlzcutoff,iNlzcutoff,flag1,flag2, nNlzcutoff,Ngzcutoff
  integer :: nepslines
  real(DP) , allocatable::  eps_data_dy(:)
  real(DP) :: epsk, deltak_para,q2d_coeff
  logical :: noncolin

  Nlzcutoff=dffts%nr3/2
  lzcutoff=Nlzcutoff*alat/dffts%nr1

  k0screen=k0screen_read
  allocate(eps_data_dy(size(qeh_eps_data(1,:))))
  call spline(qeh_eps_data(1,:),qeh_eps_data(2,:),0.0_DP,0.0_DP,eps_data_dy(:))

  mcharge0=0
  DO ig1 = 1, ngk(ik0)
    Do ig2=1, ngk(ik)
      if (sum(abs(g(:,igk_k(ig1,ik0))-g(:,igk_k(ig2,ik))))<machine_eps) then
        if (.not. noncolin )then
           mcharge0=mcharge0+conjg(evc1(ig1,ibnd0))*evc2(ig2,ibnd)
        else
           mcharge0=mcharge0+conjg(evc1(ig1,ibnd0))*evc2(ig2,ibnd) &
                            +conjg(evc1(ig1+npwx,ibnd0))*evc2(ig2+npwx,ibnd)
        endif
      endif
    Enddo
  Enddo
  deltak=norm2(xk(:,ik0)-xk(:,ik))*tpiba
  deltak_para=norm2(xk(1:2,ik0)-xk(1:2,ik))*tpiba
  epsk=1.0

  qxy=norm2(xk(1:2,ik0)-xk(1:2,ik))*tpiba
  qz= (( xk(3,ik0)-xk(3,ik))**2)**0.5*tpiba
  q2d_coeff= (1-(cos(qz*lzcutoff)-sin(qz*lzcutoff)*qz/qxy)*exp(-(qxy*lzcutoff)))
  mcharge1=mcharge0*4*pi/(deltak**2)
  mcharge2=mcharge0*4*pi/(deltak**2+k0screen**2)
  mcharge3=mcharge0*4*pi/(deltak**2)*epsk
  mcharge4=mcharge0*4*pi/(deltak**2)*            q2d_coeff
  mcharge5=mcharge0*4*pi/(deltak**2+k0screen**2)*q2d_coeff
  mcharge6=mcharge0*4*pi/(deltak**2)*epsk*       q2d_coeff
   
  mcharge1=mcharge1/omega
  mcharge2=mcharge2/omega
  mcharge3=mcharge3/omega
  mcharge4=mcharge4/omega
  mcharge5=mcharge5/omega
  mcharge6=mcharge6/omega
  write(*,*)  'Mcharge3DLFAns     ki->kf ',ik0,ik,    mcharge1, abs(mcharge1)
  write(*,*)  'Mcharge3DLFAs      ki->kf ',ik0,ik,    mcharge2, abs(mcharge2) , 'k0screen', k0screen
  write(*,*)  'Mcharge3DLFAes     ki->kf ',ik0,ik,    mcharge3, abs(mcharge3) , 'epsk',epsk 
  write(*,*)  'Mcharge3DcutLFAns  ki->kf ',ik0,ik,    mcharge4, abs(mcharge4)
  write(*,*)  'Mcharge3DcutLFAs   ki->kf ',ik0,ik,    mcharge5, abs(mcharge5) , 'k0screen', k0screen
  write(*,*)  'Mcharge3DcutLFAes  ki->kf ',ik0,ik,    mcharge6, abs(mcharge6) , 'epsk',epsk
   
  mcharge0=0
  icount=0
  DO ig1 = 1, ngk(ik0)
    Do ig2=1, ngk(ik)
      if (sum(abs(g(1:2,igk_k(ig1,ik0))-g(1:2,igk_k(ig2,ik))))<machine_eps) then
         icount=icount+1
        if (.not. noncolin )then
           mcharge0=mcharge0+conjg(evc1(ig1,ibnd0))*evc2(ig2,ibnd)
        else
           mcharge0=mcharge0+conjg(evc1(ig1,ibnd0))*evc2(ig2,ibnd) &
                            +conjg(evc1(ig1+npwx,ibnd0))*evc2(ig2+npwx,ibnd)
        endif
      endif
    Enddo
  Enddo
  deltak=((xk(1,ik)-xk(1,ik0))**2&
         +(xk(2,ik)-xk(2,ik0))**2)**0.5*tpiba
  
  mcharge1=mcharge0*tpi/deltak
  mcharge2=mcharge0*tpi/(deltak**2+k0screen**2)**0.5
  mcharge3=mcharge0*tpi/(deltak**2)**0.5*epsk
  mcharge1=mcharge1/omega
  mcharge2=mcharge2/omega
  mcharge3=mcharge3/omega
  write(*,*)   'mcharge0       0ki->kf ',ik0,ik, mcharge0, abs(mcharge0)
  write(*,*)   'Mcharge2DLFAns 0ki->kf ',ik0,ik, mcharge1, abs(mcharge1)
  write(*,*)   'Mcharge2DLFAs  0ki->kf ',ik0,ik, mcharge2, abs(mcharge2) , 'k0screen', k0screen
  write(*,*)   'Mcharge2DLFAes 0ki->kf ',ik0,ik, mcharge3, abs(mcharge3) , 'epsk', epsk
  mcharge=mcharge3
END SUBROUTINE calcmdefect_charge_lfa
