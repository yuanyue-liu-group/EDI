  
    SUBROUTINE calcmdefect_charge_lfa(ibnd0,ibnd,ik0,ik)
USE kinds, ONLY: DP,sgl
USE fft_base,  ONLY: dfftp, dffts
USE gvect, ONLY: ngm, gstart, g, gg, gcutm, igtongl
USE klist , ONLY: nks, nelec, xk, wk, degauss, ngauss, igk_k, ngk
 !     Use edic_mod,   only: V_file, V_loc, V_0, Bxc_1, Bxc_2, Bxc_3, V_p
      Use edic_mod, Only : evc1,evc2,evc3,evc4,&
                               psic1, psic2, psic3, psic4
    use splinelib, only: dosplineint,spline,splint
    COMPLEX(DP) ::  mcharge0,mcharge1,mcharge2,mcharge3,mcharge4,mcharge5,mcharge6
    INTEGER :: ibnd, ik, ik0,ibnd0

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



    real(DP) , allocatable::  eps_data_dy(:)
    real(DP) :: epsk, deltak_para,q2d_coeff
    !k0screen=tpiba*0.01
    allocate(eps_data_dy(size(eps_data(1,:))))
    call spline(eps_data(1,:),eps_data(2,:),0.0_DP,0.0_DP,eps_data_dy(:))

    !k0screen=tpiba*0.01
    mcharge0=0
    DO ig1 = 1, ngk(ik0)
      Do ig2=1, ngk(ik)
        if (sum(abs(g(:,igk_k(ig1,ik0))-g(:,igk_k(ig2,ik))))<machine_eps) then
             mcharge0=mcharge0+conjg(evc1(ig1,ibnd0))*evc2(ig2,ibnd)
        endif
      Enddo
    Enddo
    deltak=norm2(xk(:,ik0)-xk(:,ik))*tpiba
    deltak_para=norm2(xk(1:2,ik0)-xk(1:2,ik))*tpiba
    epsk= splint(eps_data(1,:),eps_data(2,:),eps_data_dy(:),deltak_para)
    if (deltak>maxval(eps_data(1,:)))      epsk=minval(eps_data(2,:))

    qxy=norm2(xk(1:2,ik0)-xk(1:2,ik))*tpiba
    qz= (( xk(3,ik0)-xk(3,ik))**2)**0.5*tpiba
    q2d_coeff= (1-(cos(qz*lzcutoff)-sin(qz*lzcutoff)*qz/qxy)*exp(-(qxy*lzcutoff)))
    mcharge1=mcharge0*4*pi/(deltak**2)
    mcharge2=mcharge0*4*pi/(deltak**2+k0screen**2)
    mcharge3=mcharge0*4*pi/(deltak**2)*epsk
    mcharge4=mcharge0*4*pi/(deltak**2)*            q2d_coeff!(1-(cos(qz*lzcutoff)-sin(qz*lzcutoff)*qz/qxy)*exp(-(qxy*lzcutoff)))
    mcharge5=mcharge0*4*pi/(deltak**2+k0screen**2)*q2d_coeff!(1-(cos(qz*lzcutoff)-sin(qz*lzcutoff)*qz/qxy)*exp(-(qxy*lzcutoff)))
    mcharge6=mcharge0*4*pi/(deltak**2)*epsk*       q2d_coeff!(1-(cos(qz*lzcutoff)-sin(qz*lzcutoff)*qz/qxy)*exp(-(qxy*lzcutoff)))
     
    mcharge1=mcharge1/dffts%nnr
    mcharge2=mcharge2/dffts%nnr
    mcharge3=mcharge3/dffts%nnr
    mcharge4=mcharge4/dffts%nnr
    mcharge5=mcharge5/dffts%nnr
    mcharge6=mcharge6/dffts%nnr
    write(*,*)  'mcharge0           ki->kf ',ik0,ik,    mcharge0, abs(mcharge0)
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
           mcharge0=mcharge0+conjg(evc1(ig1,ibnd0))*evc2(ig2,ibnd)
        endif
      Enddo
    Enddo
    deltak=((xk(1,ik)-xk(1,ik0))**2&
           +(xk(2,ik)-xk(2,ik0))**2)**0.5*tpiba
    
    mcharge1=mcharge0*tpi/deltak
    mcharge2=mcharge0*tpi/(deltak**2+k0screen**2)**0.5
    mcharge3=mcharge0*tpi/(deltak**2)**0.5*epsk
    mcharge1=mcharge1/dffts%nnr
    mcharge2=mcharge2/dffts%nnr
    mcharge3=mcharge3/dffts%nnr
    write(*,*)   'mcharge0       0ki->kf ',ik0,ik, mcharge0, abs(mcharge0)
    write(*,*)   'Mcharge2DLFAns 0ki->kf ',ik0,ik, mcharge1, abs(mcharge1)
    write(*,*)   'Mcharge2DLFAs  0ki->kf ',ik0,ik, mcharge2, abs(mcharge2) , 'k0screen', k0screen
    write(*,*)   'Mcharge2DLFAes 0ki->kf ',ik0,ik, mcharge3, abs(mcharge3) , 'epsk', epsk
 
     
 
    


    END SUBROUTINE calcmdefect_charge_lfa
 
