!subroutine calcmdefect_ml_rs_noncolin(ibnd0,ibnd,ik0,ik, V_d, v_nc)
subroutine calcmdefect_ml_rs_noncolin(ibnd,ibnd0,ik,ik0)
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
    Implicit None 

    complex(dp) :: ml_up, ml_down
    real(dp) :: d1, d2, d3

    INTEGER :: ibnd, ik, ik0,ibnd0, ig

    REAL(dp)::arg,argt,argt2
    COMPLEX(DP)::phase
    INTEGER :: npw,  ispin, nbndup, nbnddown, &
                nk , ikk,ikk0,  inr, ig1, ig2
    INTEGER:: irx,iry,irz
    INTEGER:: irx2,iry2,irz2
    INTEGER:: irx1,iry1,irz1
    
    INTEGER :: ix0,ix1,ix2
    INTEGER :: iy0,iy1,iy2
    INTEGER :: iz0,iz1,iz2, ikpsi0, ikpsi1, ikpsi2
    integer :: ir1mod,ir2mod,ir3mod,irnmod


!ikk=ik
!                  CALL read_collected_wfc ( restart_dir(), ikk, evc2 )
!    type(V_file) :: V_d
!    real(DP) :: v_nc(:,:)


    !write(*,*) v_nc(20:40,1)
    !write(*,*) v_nc(20:40,2)
    !write(*,*) V_d%plot(20:40)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!! vl in real super2prim, module
    !auxr(:) =  vrs(:,1)
    !psiprod(:)=0.00
    !vgk_perturb(:)=0.00
      write(*,*)'ML1'
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
    
    write(*,*) 'v_colin',shape(v_nc),V_d%nr1,V_d%nr2,V_d%nr3
    write(*,*) 'psi',shape(psic1),dffts%nr1,dffts%nr2,dffts%nr3
      write(*,*)'ML2',size(psic2),shape(psic2),shape(psic1),dffts%nnr
    psic2(1:dffts%nnr) = (0.d0,0.d0)
    psic1(1:dffts%nnr) = (0.d0,0.d0)
    psic4(1:dffts%nnr) = (0.d0,0.d0)
    psic3(1:dffts%nnr) = (0.d0,0.d0)
      ig=1
     write(*,*)'ML2.0',shape(igk_k)
     write(*,*)'ML2.1',igk_k(1:4,ik)
     write(*,*)'ML2.2',ngk(ik)
     write(*,*)'ML2.3',dffts%nl (igk_k(1:4,ik) ) ,ik, ibnd

    DO ig = 1, ngk(ik)
       psic2 (dffts%nl (igk_k(ig,ik) ) ) = evc2 (ig, ibnd)
    ENDDO
      write(*,*)'ML3' , sum(psic2),sum(Evc1)
      write(*,*)'ML3.1'
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
      write(*,*)'ML4' , sum(psic2),sum(Evc1)
      write(*,*)'ML4.1' , psic1(1),Evc1(1,ibnd0)
      write(*,*)'ML4.2' , psic1(2),Evc1(2,ibnd0)
      write(*,*)'ML4.3' , psic2(1),Evc2(1,ibnd0)
      write(*,*)'ML4.4' , psic2(2),Evc2(2,ibnd0)
    CALL invfft ('Wave', psic1, dffts)
    CALL invfft ('Wave', psic4, dffts)
    CALL invfft ('Wave', psic3, dffts)
    
    write(*,*) psic1(1:10)
    write(*,*) psic2(1:10)
    write(*,*) psic3(1:10)
    write(*,*) psic4(1:10)

    arg=0
    inr=0
    write(*,*) 'xk-xk01',xk(1,ik)-xk(1,ik0)
    write(*,*) 'xk-xk02',xk(2,ik)-xk(2,ik0)
    write(*,*) 'xk-xk03',xk(3,ik)-xk(3,ik0)
    do irz =0, V_d%nr3-1
    ir3mod=irz-(irz/(dffts%nr3))*dffts%nr3
    do iry =0, V_d%nr2-1
    ir2mod=iry-(iry/(dffts%nr2))*dffts%nr2
    do irx =0, V_d%nr1-1
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
    arg=irz*d3+(iry-iry/(V_d%nr2/2+1)*V_d%nr2)*d2+(irx-irx/(V_d%nr1/2+1)*V_d%nr1)*d1
    !arg=irz*d3+(iry-iry/(dffts%nr2/2+1)*dffts%nr1)*d2+(irx-irx/(dffts%nr1/2+1)*dffts%nr1)*d1
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    arg=irz*d3+iry*d2+irx*d1
    
    phase=CMPLX(COS(arg),SIN(arg),kind=dp)
    inr=inr+1
    irnmod=(ir3mod)*dffts%nr1*dffts%nr2+(ir2mod)*dffts%nr1+ir1mod+1
    ml_up=ml_up+CONJG(psic1(irnmod))*psic2(irnmod)*phase*v_nc(inr, 1)
    ml_down=ml_down+CONJG(psic3(irnmod))*psic4(irnmod)*phase*v_nc(inr, 2)

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



