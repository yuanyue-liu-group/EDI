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
  !INTEGER:: irx2,iry2,irz2
  !INTEGER:: irx1,iry1,irz1
  
  !INTEGER :: ix0,ix1,ix2
  !INTEGER :: iy0,iy1,iy2
  !INTEGER :: iz0,iz1,iz2, ikpsi0, ikpsi1, ikpsi2
  !COMPLEX(DP)::vlfft
  COMPLEX(DP) ::  ml0,ml1,ml2, ml3,ml4,ml5,ml6,ml7
  
  COMPLEX(DP) ,intent(inout)::  mlocal

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ALLOCATE(psiprod(dfftp%nnr))
  ALLOCATE(vgk(dfftp%nnr))
  ALLOCATE(vgk_perturb(dfftp%nnr))

  write(*,*)'ML0'
  write(*,*)'evc2 sizej',size(evc2)
  write(*,*)'evc2 (1,1)',evc2(1,1)
  write(*,*)'evc2 (1,2)',evc2(1,2)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!! vl in real super2prim, module
  !  write(*,*)'vrs',auxr(1)
  !auxr(:) =  vrs(:,1)
  write(*,*)'ML1'
  psiprod(:)=0.00
  vgk_perturb(:)=0.00
  ml=0
  psicprod=0
  psicprod1=0
  !mltot=0

  !mltot1=0

    !!!!!!!!!!!!!!!!!psi=u*exp(i*k*x)
    !!!!!!!!!!!!!!!!!V_g=int V(r)*exp(-i*g*r)

  d1=((1.0/dffts%nr1*at(1,1))*(xk(1,kp_idx_f)-xk(1,kp_idx_i)) +&
      (1.0/dffts%nr1*at(2,1))*(xk(2,kp_idx_f)-xk(2,kp_idx_i)) +&
      (1.0/dffts%nr1*at(3,1))*(xk(3,kp_idx_f)-xk(3,kp_idx_i)) )*tpi 
  d2=((1.0/dffts%nr2*at(1,2))*(xk(1,kp_idx_f)-xk(1,kp_idx_i)) +&
      (1.0/dffts%nr2*at(2,2))*(xk(2,kp_idx_f)-xk(2,kp_idx_i)) +&
      (1.0/dffts%nr2*at(3,2))*(xk(3,kp_idx_f)-xk(3,kp_idx_i)) )*tpi 
  d3=((1.0/dffts%nr3*at(1,3))*(xk(1,kp_idx_f)-xk(1,kp_idx_i)) +&
      (1.0/dffts%nr3*at(2,3))*(xk(2,kp_idx_f)-xk(2,kp_idx_i)) +&
      (1.0/dffts%nr3*at(3,3))*(xk(3,kp_idx_f)-xk(3,kp_idx_i)) )*tpi 
  
  !write(*,*) 'v_colin',shape(v_colin),V_d%nr1,V_d%nr2,V_d%nr3
  !write(*,*) 'psi',shape(psic1),dffts%nr1,dffts%nr2,dffts%nr3
  !write(*,*)'ML2',size(psic2),shape(psic2),shape(psic1),dffts%nnr

  psic2(1:dffts%nnr) = (0.d0,0.d0)
  psic1(1:dffts%nnr) = (0.d0,0.d0)
  ig=1
  !write(*,*)'ML2.0',shape(igk_k)
  !write(*,*)'ML2.1',igk_k(1:4,kp_idx_i)
  !write(*,*)'ML2.2',ngk(kp_idx_i)
  !write(*,*)'ML2.3',dffts%nl (igk_k(1:4,kp_idx_i) ) ,kp_idx_i, bnd_idx_i
  DO ig = 1, ngk(kp_idx_i)
     psic1 (dffts%nl (igk_k(ig,kp_idx_i) ) ) = evc1 (ig, bnd_idx_i)
  ENDDO
  !write(*,*)'ML3' , sum(psic1),sum(Evc1)
  !write(*,*)'ML3.1'
  DO ig = 1, ngk(kp_idx_f)
     psic2 (dffts%nl (igk_k(ig,kp_idx_f) ) ) = evc2 (ig, bnd_idx_f)
  ENDDO
  !write(*,*)'ML3.1'
  CALL invfft ('Wave', psic2, dffts)
  !write(*,*)'ML4' , sum(psic2),sum(Evc2)

  !write(*,*)'ML4.1' , psic1(1),Evc1(1,bnd_idx_i)
  !write(*,*)'ML4.2' , psic1(2),Evc1(2,bnd_idx_i)
  !write(*,*)'ML4.3' , psic2(1),Evc2(1,bnd_idx_f)
  !write(*,*)'ML4.4' , psic2(2),Evc2(2,bnd_idx_f)
  CALL invfft ('Wave', psic1, dffts)
  
  
  
  !d1=((1.0/V_d%nr1*at_perturb(1,1))*(xk(1,kp_idx_i)-xk(1,kp_idx_f)) +&
  !    (1.0/V_d%nr1*at_perturb(2,1))*(xk(2,kp_idx_i)-xk(2,kp_idx_f)) +&
  !    (1.0/V_d%nr1*at_perturb(3,1))*(xk(3,kp_idx_i)-xk(3,kp_idx_f)) )*tpi 
  !d2=((1.0/V_d%nr2*at_perturb(1,2))*(xk(1,kp_idx_i)-xk(1,kp_idx_f)) +&
  !    (1.0/V_d%nr2*at_perturb(2,2))*(xk(2,kp_idx_i)-xk(2,kp_idx_f)) +&
  !    (1.0/V_d%nr2*at_perturb(3,2))*(xk(3,kp_idx_i)-xk(3,kp_idx_f)) )*tpi 
  !d3=((1.0/V_d%nr3*at_perturb(1,3))*(xk(1,kp_idx_i)-xk(1,kp_idx_f)) +&
  !    (1.0/V_d%nr3*at_perturb(2,3))*(xk(2,kp_idx_i)-xk(2,kp_idx_f)) +&
  !    (1.0/V_d%nr3*at_perturb(3,3))*(xk(3,kp_idx_i)-xk(3,kp_idx_f)) )*tpi 
  !
  arg=0
  inr=0
  write(*,*) 'xk-xk01',xk(1,kp_idx_i)-xk(1,kp_idx_f)
  write(*,*) 'xk-xk02',xk(2,kp_idx_i)-xk(2,kp_idx_f)
  write(*,*) 'xk-xk03',xk(3,kp_idx_i)-xk(3,kp_idx_f)

  do irz =0, V_d%nr3-1
    ir3mod=irz-(irz/(dffts%nr3))*dffts%nr3
    do iry =0, V_d%nr2-1
      ir2mod=iry-(iry/(dffts%nr2))*dffts%nr2
      do irx =0, V_d%nr1-1
          ir1mod=irx-(irx/(dffts%nr1))*dffts%nr1
          !arg=tpi*(real(irx)/V_d%nr1*at_perturb(1,1)+real(iry)/V_d%nr2*at_perturb(1,2)&
          !                                              +real(irz)/V_d%nr3*at_perturb(1,3))*(xk(1,kp_idx_i)-xk(1,kp_idx_f)) +&
          !    tpi*(real(irx)/V_d%nr1*at_perturb(2,1)+real(iry)/V_d%nr2*at_perturb(2,2)&
          !                                              +real(irz)/V_d%nr3*at_perturb(2,3))*(xk(2,kp_idx_i)-xk(2,kp_idx_f)) +&
          !    tpi*(real(irx)/V_d%nr1*at_perturb(3,1)+real(iry)/V_d%nr2*at_perturb(3,2)&
          !                                              +real(irz)/V_d%nr3*at_perturb(3,3))*(xk(3,kp_idx_i)-xk(3,kp_idx_f))   
          
          arg=irz*d3+iry*d2+irx*d1
          !!!!!!!!!!!!!!!!!!!!!!!!!!!
          !move vloc center 
          !!!!!!!!!!!!!!!!!!!!!!!!!!!
          !arg=ir3mod*d3+ir2mod*d2+ir1mod*d1
          
          arg=tpi*(real(irx)/dffts%nr1*at(1,1)+real(iry)/dffts%nr2*at(1,2)+&
                                               real(irz)/dffts%nr3*at(1,3))*(xk(1,kp_idx_f)-xk(1,kp_idx_i)) +&
              tpi*(real(irx)/dffts%nr1*at(2,1)+real(iry)/dffts%nr2*at(2,2)+&
                                               real(irz)/dffts%nr3*at(2,3))*(xk(2,kp_idx_f)-xk(2,kp_idx_i)) +&
              tpi*(real(irx)/dffts%nr1*at(3,1)+real(iry)/dffts%nr2*at(3,2)+&
                                               real(irz)/dffts%nr3*at(3,3))*(xk(3,kp_idx_f)-xk(3,kp_idx_i))   
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! shift arg center
          arg=irz*d3+(iry-iry/(V_d%nr2/2+1)*V_d%nr2)*d2+(irx-irx/(V_d%nr1/2+1)*V_d%nr1)*d1
          !arg=irz*d3+(iry-iry/(dffts%nr2/2+1)*dffts%nr1)*d2+(irx-irx/(dffts%nr1/2+1)*dffts%nr1)*d1
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          arg=irz*d3+iry*d2+irx*d1
          
          phase=CMPLX(COS(arg),SIN(arg),kind=dp)
          inr=inr+1
          irnmod=(ir3mod)*dffts%nr1*dffts%nr2+(ir2mod)*dffts%nr1+ir1mod+1
          !write(*,*)'irnmod,inr',ir1mod,ir2mod,ir3mod,irx,iry,irz, irnmod,inr
          ml=ml+CONJG(psic1(irnmod))*psic2(irnmod)*v_colin(inr)*phase
          psicprod=psicprod+CONJG(psic1(irnmod))*psic2(irnmod)*phase
          psicprod1=psicprod1+CONJG(psic1(irnmod))*psic2(irnmod)
          !ml2=ml2+CONJG(psic1(irnmod))*psic2(irnmod)*v_colin(inr)*phase
          !mltot=mltot+CONJG(psic1(irnmod))*psic2(irnmod)*phase
          !mltot1=mltot1+CONJG(psic1(irnmod))*psic2(irnmod)
          !write (*,*) 'iri',ir1mod,ir2mod,ir3mod
          !write (*,*) 'grid ', irnmod
          !write (*,*) 'psic1 ', psic1(irnmod)
          !write (*,*) 'psic2 ', psic2(irnmod)
          !write (*,*) 'arg', arg
          
          if ( irnmod<0 .or. irnmod>dffts%nnr ) then
             write (*,*) 'grid mismatch', irnmod, dffts%nnr 
          endif
          
          
          if (irz==dffts%nr3/2) then
              argt= atan2(real(CONJG(psic1(irnmod))*psic2(irnmod)*phase),aimag(CONJG(psic1(irnmod))*psic2(irnmod)*phase))
              argt2= atan2(real(CONJG(psic1(irnmod))*psic2(irnmod)),aimag(CONJG(psic1(irnmod))*psic2(irnmod)))
              if (argt<0) argt=argt+tpi
          !            if (argt2<0) argt2=argt2+tpi
          !    write(*,*) 'psiplts kp_idx_i',kp_idx_i, 'xyz', irx,iry,irz,  'psi1', psic1(irnmod),abs( psic1(irnmod)),  'psi2', &
          !                        psic2(irnmod),abs(psic2(irnmod)),&
          !   'arg', arg,'prod', CONJG(psic1(irnmod))*psic2(irnmod)*phase, abs(CONJG(psic1(irnmod))*psic2(irnmod)*phase),argt,argt2,&
          !            real(CONJG(psic1(irnmod))*psic2(irnmod)),aimag(CONJG(psic1(irnmod))*psic2(irnmod)), psicprod1
          endif

             
      enddo
    enddo
  enddo
  
  !ml=ml/V_d%nr1/V_d%nr2/V_d%nr3
  ml=ml/dffts%nnr
  psicprod=psicprod/V_d%nr1/V_d%nr2/V_d%nr3
  psicprod1=psicprod1/V_d%nr1/V_d%nr2/V_d%nr3
  !psicprod=psicprod/dffts%nnr

  write(*,*) 'psicprods', psicprod , abs(psicprod), abs(psicprod1)

  !write (*,*) 'ml super to primitive ki->kf',kp_idx_f,kp_idx_i, ml, abs(ml), log(ml)
  !write (*,*) 'mlpsi*psi0 to primitive ki->kf',kp_idx_f,kp_idx_i, mltot1, abs(mltot1), log(mltot1)
  !write (*,*) 'mlpsi*psi0*phase to primitive ki->kf',kp_idx_f,kp_idx_i, mltot, abs(mltot), log(mltot)
  !write (*,*) 'modml super  ki->kf',kp_idx_f,kp_idx_i, abs(ml)
  !write (*,*) 'Ml ki->kf ',kp_idx_f,kp_idx_i, ml, abs(ml)
  write (*,*) 'Ml ibnd,ki->fbnd,kf ',bnd_idx_i,kp_idx_i,bnd_idx_f,kp_idx_f, xk(:,kp_idx_i),xk(:,kp_idx_f), ml, abs(ml)
  !write(*,*) 'nrx_perturb',V_d%nr1,V_d%nr2,V_d%nr3
  !write(*,*) 'nrx_perturb',at,at_perturb, alat
  !write (*,*) 'dvgk', vgk(:)-vgk_perturb(:)
  !write (*,*) 'vgk', vgk(:)
  !write (*,*) 'vgk_perturb', vgk_perturb(:)
  !write (*,*) 'sum dvgk', sum(vgk(:)-vgk_perturb(:))
  !write(*,*) 'xk(kp_idx_i)', xk(:,kp_idx_i),kp_idx_i
  !write(*,*) 'xk(kp_idx_f)', xk(:,kp_idx_f),kp_idx_f
1001 format(A,I9,I9,3F14.9,3F14.9," ( ",e17.9," , ",e17.9," ) ",e17.9)
  mlocal=ml
END SUBROUTINE calcmdefect_ml_rs
