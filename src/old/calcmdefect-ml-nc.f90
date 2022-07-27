subroutine calcmdefect_ml_rs_noncolin(ibnd0,ibnd,ik0,ik, V_0, V_loc)

      USE cell_base,       ONLY : alat, ibrav, omega, at, bg, celldm, wmass
      INTEGER :: ibnd, ik, ik0,ibnd0
      type(V_file) :: V_0
      real(DP) :: V_loc(:,:)


























      !write(*,*) V_loc(20:40,1)
      !write(*,*) V_loc(20:40,2)
      !write(*,*) V_0%plot(20:40)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!! vl in real super2prim, module
      auxr(:) =  vrs(:,1)
      psiprod(:)=0.00
      vgk_perturb(:)=0.00
      ml_up=0
      ml_down=0
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

      DO ig = 1, ngk(ikk)
         psic2 (dffts%nl (igk_k(ig,ikk) ) ) = evc2 (ig, ibnd)
      ENDDO
      DO ig = 1, ngk(ik0)
         psic1 (dffts%nl (igk_k(ig,ik0) ) ) = evc1 (ig, ibnd0)
      ENDDO

      DO ig = 1, ngk(ikk)
         psic4 (dffts%nl (igk_k(ig,ikk) ) ) = evc2 (ig+npwx, ibnd)
      ENDDO
      DO ig = 1, ngk(ik0)
         psic3 (dffts%nl (igk_k(ig,ik0) ) ) = evc1 (ig+npwx, ibnd0)
      ENDDO
      CALL invfft ('Wave', psic2, dffts)
      CALL invfft ('Wave', psic1, dffts)
      CALL invfft ('Wave', psic4, dffts)
      CALL invfft ('Wave', psic3, dffts)
      
      !write(*,*) psic1(1:10)
      !write(*,*) psic2(1:10)
      !write(*,*) psic3(1:10)
      !write(*,*) psic4(1:10)
    
      arg=0
      inr=0
      do irz =0, V_0%nr3-1
      ir3mod=irz-(irz/(dffts%nr3))*dffts%nr3
      do iry =0, V_0%nr2-1
      ir2mod=iry-(iry/(dffts%nr2))*dffts%nr2
      do irx =0, V_0%nr1-1
      ir1mod=irx-(irx/(dffts%nr1))*dffts%nr1
      
      arg=irz*d3+iry*d2+irx*d1
      !!!!!!!!!!!!!!!!!!!!!!!!!!!
      !move vloc center 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      arg=tpi*(real(irx)/dffts%nr1*at(1,1)+real(iry)/dffts%nr2*at(1,2)+real(irz)/dffts%nr3*at(1,3))*(xk(1,ik)-xk(1,ik0)) +&
          tpi*(real(irx)/dffts%nr1*at(2,1)+real(iry)/dffts%nr2*at(2,2)+real(irz)/dffts%nr3*at(2,3))*(xk(2,ik)-xk(2,ik0)) +&
          tpi*(real(irx)/dffts%nr1*at(3,1)+real(iry)/dffts%nr2*at(3,2)+real(irz)/dffts%nr3*at(3,3))*(xk(3,ik)-xk(3,ik0))   
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! shift arg center
      arg=irz*d3+(iry-iry/(V_0%nr2/2+1)*V_0%nr2)*d2+(irx-irx/(V_0%nr1/2+1)*V_0%nr1)*d1
      !arg=irz*d3+(iry-iry/(dffts%nr2/2+1)*dffts%nr1)*d2+(irx-irx/(dffts%nr1/2+1)*dffts%nr1)*d1
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !arg=irz*d3+iry*d2+irx*d1
      
      phase=CMPLX(COS(arg),SIN(arg),kind=dp)
      inr=inr+1
      irnmod=(ir3mod)*dffts%nr1*dffts%nr2+(ir2mod)*dffts%nr1+ir1mod+1
      ml_up=ml_up+CONJG(psic1(irnmod))*psic2(irnmod)*phase!*V_loc(inr, 1)
      ml_down=ml_down+CONJG(psic3(irnmod))*psic4(irnmod)*phase!*V_loc(inr, 2)

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
      ml_up=ml_up/dffts%nnr
      ml_down=ml_down/dffts%nnr
      write (*,1001) 'Ml_up ki->kf ',ik0,ik, ml_up, abs(ml_up)
      write (*,1001) 'Ml_down ki->kf ',ik0,ik,  ml_down, abs(ml_down)
      write (*,1002) 'Ml ki->kf ',ik0,ik,  ml_up+ml_down, abs(ml_up+ml_down)
  1001 format(A16,I9,I9," ( ",e17.9," , ",e17.9," ) ",e17.9)
  1002 format(A16,I9,I9," ( ",e17.9," , ",e17.9," ) ",e17.9/)
     
      
      
   end subroutine calcmdefect_ml_rs_noncolin



