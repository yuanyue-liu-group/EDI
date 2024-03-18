subroutine calcmdefect_ml_rs_noncolin(ibnd,ibnd0,ik,ik0,mlocal)
  Use kinds,          Only : dp
  USE cell_base,       ONLY : alat, ibrav, omega, at, bg, celldm, wmass
  USE constants, ONLY: tpi, e2, eps6,pi
  Use fft_base,       Only : dfftp, dffts
  USE fft_interfaces, ONLY : fwfft, invfft
  USE klist , ONLY: nks, nelec, xk, wk, degauss, ngauss, igk_k, ngk
  Use edic_mod, Only : evc1,evc2,evc3,evc4,&
                             psic1, psic2, psic3, psic4
  USE wvfct, ONLY: npwx, nbnd, wg, et, g2kin
  USE edic_mod, Only :  V_d, V_nc
  Use edic_mod, Only: m_loc
  USE io_files,  ONLY : prefix, tmp_dir, nwordwfc, iunwfc, restart_dir
  USE pw_restart_new,   ONLY : read_collected_wfc
  use edic_mod, only: bndkp_pair
  Implicit None 

  COMPLEX(DP) ::  mlocal
  complex(dp) :: ml_up, ml_down,ml
  real(dp) :: d1, d2, d3

  INTEGER :: ibnd, ik, ik0,ibnd0, ig!,ipair

  REAL(dp)::arg,argt,argt2
  COMPLEX(DP)::phase
  integer::inr, ig1, ig2
  INTEGER:: irx,iry,irz
  integer :: ir1mod,ir2mod,ir3mod,irnmod

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!! vl in real super2prim, module
  ml_up=0
  ml_down=0
  ml=0.0
  d1=((1.0/dffts%nr1*at(1,1))*(xk(1,ik)-xk(1,ik0)) +&
      (1.0/dffts%nr1*at(2,1))*(xk(2,ik)-xk(2,ik0)) +&
      (1.0/dffts%nr1*at(3,1))*(xk(3,ik)-xk(3,ik0)) )*tpi 
  d2=((1.0/dffts%nr2*at(1,2))*(xk(1,ik)-xk(1,ik0)) +&
      (1.0/dffts%nr2*at(2,2))*(xk(2,ik)-xk(2,ik0)) +&
      (1.0/dffts%nr2*at(3,2))*(xk(3,ik)-xk(3,ik0)) )*tpi 
  d3=((1.0/dffts%nr3*at(1,3))*(xk(1,ik)-xk(1,ik0)) +&
      (1.0/dffts%nr3*at(2,3))*(xk(2,ik)-xk(2,ik0)) +&
      (1.0/dffts%nr3*at(3,3))*(xk(3,ik)-xk(3,ik0)) )*tpi 
  
  psic2(1:dffts%nnr) = (0.d0,0.d0)
  psic1(1:dffts%nnr) = (0.d0,0.d0)
  psic4(1:dffts%nnr) = (0.d0,0.d0)
  psic3(1:dffts%nnr) = (0.d0,0.d0)
  ig=1

  DO ig = 1, ngk(ik)
    psic2 (dffts%nl (igk_k(ig,ik) ) ) = evc2 (ig, ibnd)
  ENDDO
  DO ig = 1, ngk(ik0)
    psic1 (dffts%nl (igk_k(ig,ik0) ) ) = evc1 (ig, ibnd0)
  ENDDO

  DO ig = 1, ngk(ik)
    psic4 (dffts%nl (igk_k(ig,ik) ) ) = evc2 (ig+npwx, ibnd)
  ENDDO
  DO ig = 1, ngk(ik0)
    psic3 (dffts%nl (igk_k(ig,ik0) ) ) = evc1 (ig+npwx, ibnd0)
  ENDDO
  CALL invfft ('Wave', psic2, dffts)
  CALL invfft ('Wave', psic1, dffts)
  CALL invfft ('Wave', psic4, dffts)
  CALL invfft ('Wave', psic3, dffts)
  arg=0
  inr=0
  do irz =0, V_d%nr3-1
    ir3mod=irz-(irz/(dffts%nr3))*dffts%nr3
    do iry =0, V_d%nr2-1
      ir2mod=iry-(iry/(dffts%nr2))*dffts%nr2
      do irx =0, V_d%nr1-1
          ir1mod=irx-(irx/(dffts%nr1))*dffts%nr1
          
          !arg=irz*d3+iry*d2+irx*d1
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!move vloc center 
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !
          !arg=tpi*(real(irx)/dffts%nr1*at(1,1)+real(iry)/dffts%nr2*at(1,2)+real(irz)/dffts%nr3*at(1,3))*(xk(1,ik)-xk(1,ik0)) +&
          !    tpi*(real(irx)/dffts%nr1*at(2,1)+real(iry)/dffts%nr2*at(2,2)+real(irz)/dffts%nr3*at(2,3))*(xk(2,ik)-xk(2,ik0)) +&
          !    tpi*(real(irx)/dffts%nr1*at(3,1)+real(iry)/dffts%nr2*at(3,2)+real(irz)/dffts%nr3*at(3,3))*(xk(3,ik)-xk(3,ik0))   
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !! shift arg center
          !arg=irz*d3+(iry-iry/(V_d%nr2/2+1)*V_d%nr2)*d2+(irx-irx/(V_d%nr1/2+1)*V_d%nr1)*d1
          !!arg=irz*d3+(iry-iry/(dffts%nr2/2+1)*dffts%nr1)*d2+(irx-irx/(dffts%nr1/2+1)*dffts%nr1)*d1
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          arg=irz*d3+iry*d2+irx*d1
          
          phase=CMPLX(COS(arg),SIN(arg),kind=dp)
          inr=inr+1
          irnmod=(ir3mod)*dffts%nr1*dffts%nr2+(ir2mod)*dffts%nr1+ir1mod+1
          ml=ml+CONJG(psic1(irnmod))*psic2(irnmod)*phase*(v_nc(inr, 1)+v_nc(inr,4))&
               +CONJG(psic3(irnmod))*psic4(irnmod)*phase*(v_nc(inr, 1)-v_nc(inr,4))&
               +CONJG(psic1(irnmod))*psic2(irnmod)*phase*(v_nc(inr, 2)+v_nc(inr,3))&
               +CONJG(psic3(irnmod))*psic4(irnmod)*phase*(v_nc(inr, 2)-v_nc(inr,3)) 

          if ( irnmod<0 .or. irnmod>dffts%nnr ) then
             write (*,*) 'grid mismatch', irnmod, dffts%nnr 
          endif
          
          if (irz==dffts%nr3/2) then
             argt= atan2(real(CONJG(psic1(irnmod))*psic2(irnmod)*phase),aimag(CONJG(psic1(irnmod))*psic2(irnmod)*phase))
             argt2= atan2(real(CONJG(psic1(irnmod))*psic2(irnmod)),aimag(CONJG(psic1(irnmod))*psic2(irnmod)))
             if (argt<0) argt=argt+tpi
          endif
      enddo
    enddo
  enddo
  ml=ml/dffts%nnr
1001 format(A16,I9,I9," ( ",e17.9," , ",e17.9," ) ",e17.9)
1002 format(A16,I9,I9," ( ",e17.9," , ",e17.9," ) ",e17.9/)
  mlocal=ml
   
    
    
end subroutine calcmdefect_ml_rs_noncolin



