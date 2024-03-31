SUBROUTINE calcmdefect_ml_rs(bnd_idx_f,bnd_idx_i,kp_idx_f,kp_idx_i,mlocal)
  USE kinds, ONLY: DP,sgl
  USE cell_base,       ONLY : alat, ibrav, omega, at, bg, celldm, wmass
  USE scf, ONLY: rho, rho_core, rhog_core, v, vltot, vrs
  Use edic_mod,   only: V_colin,v_d
  USE fft_base,  ONLY: dfftp, dffts
  USE fft_interfaces, ONLY : fwfft, invfft
  Use edic_mod, Only : evc1,evc2, psic1, psic2
  USE klist , ONLY: nks, nelec, xk, wk, degauss, ngauss, igk_k, ngk
  use edic_mod, only: bndkp_pair
  INTEGER,intent(in) :: bnd_idx_i, kp_idx_i, kp_idx_f,bnd_idx_f
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! intermediate data
  COMPLEX(DP), ALLOCATABLE ::   psiprod(:),vgk(:),vgk_perturb(:),vkb_perturb(:,:)
  COMPLEX(DP) :: mnl, ml,mltot,mltot1,mltot2,mnltot,psicnorm,psicprod,enl1,phaseft,psicprod1
  LOGICAL :: offrange
  REAL(dp)::arg,argt,argt2
  COMPLEX(DP)::phase
  INTEGER:: irx,iry,irz
  
  COMPLEX(DP) ::  ml0,ml1,ml2, ml3,ml4,ml5,ml6,ml7
  
  COMPLEX(DP) ,intent(inout)::  mlocal

  ALLOCATE(psiprod(dfftp%nnr))
  ALLOCATE(vgk(dfftp%nnr))
  ALLOCATE(vgk_perturb(dfftp%nnr))
  psiprod(:)=0.00
  vgk_perturb(:)=0.00
  ml=0
  psicprod=0
  psicprod1=0
  d1=((1.0/dffts%nr1*at(1,1))*(xk(1,kp_idx_f)-xk(1,kp_idx_i)) +&
      (1.0/dffts%nr1*at(2,1))*(xk(2,kp_idx_f)-xk(2,kp_idx_i)) +&
      (1.0/dffts%nr1*at(3,1))*(xk(3,kp_idx_f)-xk(3,kp_idx_i)) )*tpi 
  d2=((1.0/dffts%nr2*at(1,2))*(xk(1,kp_idx_f)-xk(1,kp_idx_i)) +&
      (1.0/dffts%nr2*at(2,2))*(xk(2,kp_idx_f)-xk(2,kp_idx_i)) +&
      (1.0/dffts%nr2*at(3,2))*(xk(3,kp_idx_f)-xk(3,kp_idx_i)) )*tpi 
  d3=((1.0/dffts%nr3*at(1,3))*(xk(1,kp_idx_f)-xk(1,kp_idx_i)) +&
      (1.0/dffts%nr3*at(2,3))*(xk(2,kp_idx_f)-xk(2,kp_idx_i)) +&
      (1.0/dffts%nr3*at(3,3))*(xk(3,kp_idx_f)-xk(3,kp_idx_i)) )*tpi 
  
  psic2(1:dffts%nnr) = (0.d0,0.d0)
  psic1(1:dffts%nnr) = (0.d0,0.d0)
  ig=1
  DO ig = 1, ngk(kp_idx_i)
     psic1 (dffts%nl (igk_k(ig,kp_idx_i) ) ) = evc1 (ig, bnd_idx_i)
  ENDDO
  DO ig = 1, ngk(kp_idx_f)
     psic2 (dffts%nl (igk_k(ig,kp_idx_f) ) ) = evc2 (ig, bnd_idx_f)
  ENDDO
  CALL invfft ('Wave', psic2, dffts)
  CALL invfft ('Wave', psic1, dffts)
  arg=0
  inr=0

  do irz =0, V_d%nr3-1
    ir3mod=irz-(irz/(dffts%nr3))*dffts%nr3
    do iry =0, V_d%nr2-1
      ir2mod=iry-(iry/(dffts%nr2))*dffts%nr2
      do irx =0, V_d%nr1-1
          ir1mod=irx-(irx/(dffts%nr1))*dffts%nr1
          !!arg=tpi*(real(irx)/V_d%nr1*at_perturb(1,1)+real(iry)/V_d%nr2*at_perturb(1,2)&
          !!                                              +real(irz)/V_d%nr3*at_perturb(1,3))*(xk(1,kp_idx_i)-xk(1,kp_idx_f)) +&
          !!    tpi*(real(irx)/V_d%nr1*at_perturb(2,1)+real(iry)/V_d%nr2*at_perturb(2,2)&
          !!                                              +real(irz)/V_d%nr3*at_perturb(2,3))*(xk(2,kp_idx_i)-xk(2,kp_idx_f)) +&
          !!    tpi*(real(irx)/V_d%nr1*at_perturb(3,1)+real(iry)/V_d%nr2*at_perturb(3,2)&
          !!                                              +real(irz)/V_d%nr3*at_perturb(3,3))*(xk(3,kp_idx_i)-xk(3,kp_idx_f))   
          !
          !arg=irz*d3+iry*d2+irx*d1
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!move vloc center 
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!arg=ir3mod*d3+ir2mod*d2+ir1mod*d1
          !
          !arg=tpi*(real(irx)/dffts%nr1*at(1,1)+real(iry)/dffts%nr2*at(1,2)+&
          !                                     real(irz)/dffts%nr3*at(1,3))*(xk(1,kp_idx_f)-xk(1,kp_idx_i)) +&
          !    tpi*(real(irx)/dffts%nr1*at(2,1)+real(iry)/dffts%nr2*at(2,2)+&
          !                                     real(irz)/dffts%nr3*at(2,3))*(xk(2,kp_idx_f)-xk(2,kp_idx_i)) +&
          !    tpi*(real(irx)/dffts%nr1*at(3,1)+real(iry)/dffts%nr2*at(3,2)+&
          !                                     real(irz)/dffts%nr3*at(3,3))*(xk(3,kp_idx_f)-xk(3,kp_idx_i))   
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !! shift arg center
          !arg=irz*d3+(iry-iry/(V_d%nr2/2+1)*V_d%nr2)*d2+(irx-irx/(V_d%nr1/2+1)*V_d%nr1)*d1
          !!arg=irz*d3+(iry-iry/(dffts%nr2/2+1)*dffts%nr1)*d2+(irx-irx/(dffts%nr1/2+1)*dffts%nr1)*d1
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          arg=irz*d3+iry*d2+irx*d1
          
          phase=CMPLX(COS(arg),SIN(arg),kind=dp)
          inr=inr+1
          irnmod=(ir3mod)*dffts%nr1*dffts%nr2+(ir2mod)*dffts%nr1+ir1mod+1
          ml=ml+CONJG(psic1(irnmod))*psic2(irnmod)*v_colin(inr)*phase
          
          if ( irnmod<0 .or. irnmod>dffts%nnr ) then
             write (*,*) 'grid mismatch', irnmod, dffts%nnr 
          endif
      enddo
    enddo
  enddo
  
  ml=ml/dffts%nnr
1001 format(A,I9,I9,3F14.9,3F14.9," ( ",e17.9," , ",e17.9," ) ",e17.9)
  mlocal=ml
END SUBROUTINE calcmdefect_ml_rs
