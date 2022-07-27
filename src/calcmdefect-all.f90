

    SUBROUTINE calcmdefect_all()! initialization and call M subroutines 

    write (*,*) 'enter calcmdefect module'

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !extra function: not fully implemented
    ALLOCATE (idx (ngm) )
    ALLOCATE (igtog (ngm) )
    ALLOCATE (gtoig (ngm) )
    idx(:) = 0
    igtog(:) = 0
    IF( lsda )THEN
       nbndup = nbnd
       nbnddown = nbnd
       nk = nks/2
       !     nspin = 2
    ELSE
       nbndup = nbnd
       nbnddown = 0
       nk = nks
       !     nspin = 1
    ENDIF


    DO ispin = 1, nspin
       DO ik = 1, nk
          ikk = ik + nk*(ispin-1)
          idx( igk_k(1:ngk(ikk),ikk) ) = 1
       ENDDO
    ENDDO

    ngtot_l = 0
    DO ig = 1, ngm
       IF( idx(ig) >= 1 )THEN
          ngtot_l = ngtot_l + 1
          igtog(ngtot_l) = ig
          gtoig(ig) = ngtot_l
       ENDIF
    ENDDO
    !extra function: not fully implemented
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ALLOCATE (aux(dfftp%nnr))
    ALLOCATE(auxr(dfftp%nnr))
    ALLOCATE(psiprod(dfftp%nnr))
    ALLOCATE(vgk(dfftp%nnr))
    ALLOCATE(vgk_perturb(dfftp%nnr))
    ALLOCATE( auxg( dfftp%ngm ) )
    !mltot=0
    !mnltot=0
    !mltot1=0
    !mltot2=0
    
    !write(*,*) 'use_calcmdefect', use_calcmdefect
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! extra data read in, not used
    IF ( npool > 1 .or. nimage > 1 .or. nbgrp > 1 ) &
      CALL errore('calcmdefect', 'pool/band/image parallelization not (yet) implemented',1)
    IF ( noncolin .OR. lspinorb ) &
      CALL errore('calcmdefect', 'noncollinear/spinorbit magnetism not (yet) implemented',2)
    tmp_unit = find_free_unit()
    OPEN(unit=tmp_unit,file = 'calcmdefect.dat',status='old',err=20)
    !OPEN(unit=tmp_unit,file = trim(tmp_dir)//'calcmdefect.dat',status='old',err=20)
    20 continue
        READ(tmp_unit,calcmcontrol,iostat=ios)
    CLOSE(tmp_unit)
    ! extra data read in, not used
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!! vloc.dat read
    iunpot_perturb=99 
    filpot_perturb=vperturb_filename
    !write(*,*) vperturb_filename
    open (unit = iunpot_perturb, file = filpot_perturb, form = 'formatted', &
         status = 'old', err = 99, iostat = ios_perturb)
    99 call errore ('mloc', 'opening file '//TRIM(filpot_perturb), abs (ios_perturb) )
    
    read (iunpot_perturb, '(a)') title_perturb
    read (iunpot_perturb, * ) nr1x_perturb, nr2x_perturb, nr3x_perturb,&
            nr1_perturb, nr2_perturb, nr3_perturb, nat_perturb, ntyp_perturb
    
    allocate(plot_perturb( nr1_perturb*nr2_perturb*nr3_perturb))
    allocate(ityp_perturb(nat_perturb))
    allocate(zv_perturb(ntyp_perturb))
    allocate(atm_perturb(ntyp_perturb))
    allocate(tau_perturb(3,nat_perturb))
    
    read (iunpot_perturb, * ) ibrav_perturb, celldm_perturb
    if (ibrav_perturb == 0) then
       do i_perturb = 1,3
          read ( iunpot_perturb, * ) ( at_perturb(ipol_perturb,i_perturb),ipol_perturb=1,3 )
       enddo
       alat_perturb=celldm_perturb(1)
    else
       call latgen(ibrav_perturb,celldm_perturb,at_perturb(1,1),at_perturb(1,2),at_perturb(1,3),omega_perturb)
       at_perturb(:,:)=at_perturb(:,:)/alat
    endif
    read (iunpot_perturb, * ) gcutm_perturb, dual_perturb, ecut_perturb, plot_num_perturb
    !read (iunpot_perturb, *) &
    read (iunpot_perturb, '(i4,3x,a2,3x,f5.2)') &
            (nt_perturb, atm_perturb(nt_perturb), zv_perturb(nt_perturb), nt_perturb=1, ntyp_perturb)
    read (iunpot_perturb, *) (ndum_perturb,  (tau_perturb (ipol_perturb, na_perturb), ipol_perturb = 1, 3), &
            ityp_perturb(na_perturb), na_perturb = 1, nat_perturb)
    read (iunpot_perturb, * ) (plot_perturb (ir_perturb), ir_perturb = 1, nr1_perturb * nr2_perturb * nr3_perturb)
    tau_perturb(:,:)=tau_perturb(:,:)*alat_perturb/alat

    !debug output
    !write (*,*) 'dv readin-vrs', sum(plot_perturb(:)-vrs(:,1))
    write (*,*) 'dv readin-vrs: , sum(plot_perturb(:)),sum(vrs(:,1)),sum(plot_perturb(:))-sum(vrs(:,1))'
    write (*,*)  sum(plot_perturb(:)),sum(vrs(:,1)),sum(plot_perturb(:))-sum(vrs(:,1))
    write (*,*) 'at-perturb', at_perturb
    write (*,*) 'alat-perturb', alat_perturb
    write (*,*) 'nr1_perturb ', nr1_perturb
    write (*,*) 'nr2_perturb ', nr2_perturb
    write (*,*) 'nr3_perturb ', nr3_perturb
    write (*,*) 'at', at(:,1)
    write (*,*) 'at', at(:,2)
    write (*,*) 'at', at(:,3)
    write (*,*) 'dfftp%nr1 ', dfftp%nr1
    write (*,*) 'dfftp%nr2 ', dfftp%nr2
    write (*,*) 'dfftp%nr3 ', dfftp%nr3
    write (*,*) 'dffts%nr1 ', dffts%nr1
    write (*,*) 'dffts%nr2 ', dffts%nr2
    write (*,*) 'dffts%nr3 ', dffts%nr3
    
     CALL md5_from_file(vperturb_filename, vf_md5_cksum)
    write (*,*) 'potential files:',TRIM(vperturb_filename),'  MD5 sum:',vf_md5_cksum
    !write (*,*) 'dv readin-vrs', plot_perturb(:)-vrs(:,1)
    !write (*,*) 'dv readin-vrs', vrs(:,1)
    !write (*,*) 'dv readin-vrs', plot_perturb(:)
    !
    !!!!!! vloc.dat read
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

    
    allocate(evc3(npwx,nbnd))
    allocate(evc4(npwx,nbnd))
    allocate(mlat2(dfftp%nr3))
    allocate(mlat1(dfftp%nr3))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!! eps read 
    iunpot_perturb=99 
    filpot_perturb=eps_filename
    open (unit = iunpot_perturb, file = filpot_perturb, form = 'formatted', &
         status = 'old', err = 99, iostat = ios_perturb)
    
    
    read (iunpot_perturb, '(a)') title_perturb
!    read (iunpot_perturb, * ) k0screen_read
    read (iunpot_perturb, * ) nepslines
    
    allocate(eps_data(2,nepslines))
    do ig= 1, nepslines
         read (iunpot_perturb, * ) eps_data(1,ig),eps_data(2,ig)
    enddo
    write (*,*) 'eps lines  ', nepslines
    write (*,*) 'eps data  ', eps_data(1,1),eps_data(2,1)
    write (*,*) 'eps data  ', eps_data(1,2),eps_data(2,2)
    write (*,*) 'eps data  ', eps_data(1,3),eps_data(2,3)
    write (*,*) 'eps data  ', eps_data(1,7),eps_data(2,7)
    k0screen=k0screen_read
    
     CALL md5_from_file(eps_filename, epsf_md5_cksum)
    write (*,*) 'eps files:',trim(eps_filename),'  MD5 sum:',epsf_md5_cksum
    !write (*,*) 'dv readin-vrs', plot_perturb(:)-vrs(:,1)
    !write (*,*) 'dv readin-vrs', vrs(:,1)
    !write (*,*) 'dv readin-vrs', plot_perturb(:)
    !
    !!!!!! eps read 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
 call gw_eps_init(epsmat_q0_filename,gw_ng_data_q0 ,gw_nmtx_max_data_q0 ,gw_nmtx_data_q0 ,&
gw_gind_eps2rho_data_q0 ,gw_gind_rho2eps_data_q0 ,&
                                 gw_g_components_data_q0 ,gw_bvec_data_q0 ,gw_blat_data_q0 ,gw_qpts_data_q0 ,gw_nq_data_q0 ,&
              gw_epsmat_diag_data_q0 ,gw_epsmat_full_data_q0 ,gw_q_g_commonsubset_indinrho_q0 ,&
gw_q_g_commonsubset_size_q0,gw_qabs_q0)

 call gw_eps_init(epsmat_q1_filename,gw_ng_data_q1 ,gw_nmtx_max_data_q1 ,gw_nmtx_data_q1 ,&
gw_gind_eps2rho_data_q1 ,gw_gind_rho2eps_data_q1 ,&
                                 gw_g_components_data_q1 ,gw_bvec_data_q1 ,gw_blat_data_q1 ,gw_qpts_data_q1 ,gw_nq_data_q1 ,&
              gw_epsmat_diag_data_q1 ,gw_epsmat_full_data_q1 ,gw_q_g_commonsubset_indinrho_q1 ,&
gw_q_g_commonsubset_size_q1,gw_qabs_q1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


call get_gind_rhoandpsi_gw(gind_rho2psi_gw_q1,gind_psi2rho_gw_q1,gw_ng_data_q1,&
gw_q_g_commonsubset_size_q1,gvec_gw_q1,gw_bvec_data_q1,gw_q_g_commonsubset_indinrho_q1)
call get_gind_rhoandpsi_gw(gind_rho2psi_gw_q0,gind_psi2rho_gw_q0,gw_ng_data_q0,&
gw_q_g_commonsubset_size_q0,gvec_gw_q0,gw_bvec_data_q0,gw_q_g_commonsubset_indinrho_q0)

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!eloc from rho*vloc
    mnl=0
    DO ig = 1, dffts%nnr
       mnl=mnl+rho%of_r(ig,1)
    ENDDO
    write(*,*) 'rhotot',mnl, ml/mnl*8
    
    ml=0
    auxr(:) =  vltot(:)
    DO ig = 1, dffts%nnr
       ml=ml+rho%of_r(ig,1)*auxr(ig)
    ENDDO
    write(*,*) 'el=rho*vltot', ml
    
    ml=0
    auxr(:) = v%of_r(:,1) 
    DO ig = 1, dffts%nnr
       ml=ml+rho%of_r(ig,1)*auxr(ig)
    ENDDO
    write(*,*) 'el=rho*v%of_r', ml
    
    ml=0
    auxr(:) = vrs(:,1)
    DO ig = 1, dffts%nnr
       ml=ml+rho%of_r(ig,1)*auxr(ig)
    ENDDO
    write(*,*) 'el=rho*vrs', ml
    !!!!!!!!!!!eloc from rho*vloc
    
    
    
    !allocate(evc1(size(evc)/nbnd,nbnd))
    !allocate(evc2(size(evc)/nbnd,nbnd))
    !write (*,*) 'npwx,npw',npwx,npw
    allocate(evc1(npwx,nbnd))
    allocate(evc2(npwx,nbnd))
    allocate(psic1(dfftp%nnr))
    allocate(psic2(dfftp%nnr))
    
       
    
    
!    tau_perturb(1,:)=tau_perturb(1,:)-(at(1,1)+at(2,1)+at(3,1))*1.4
!    tau_perturb(2,:)=tau_perturb(2,:)-(at(1,2)+at(2,2)+at(3,2))*1.4
!    tau(1,:)=tau(1,:)-(at(1,1)+at(2,1)+at(3,1))*1.0
!    tau(2,:)=tau(2,:)-(at(1,2)+at(2,2)+at(3,2))*1.0
!    do ig=1,nat_perturb
!    write (*,*) 'tau_perturb, ', shape(tau_perturb),tau_perturb(:,ig)
!    enddo

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! loop through k points
    write (*,*) 'start M calculation k loop'
!    write (*,*) 'xk',xk,nk
    ibnd0=bnd_initial
    ibnd=bnd_final
    write (*,*) 'ibnd0->ibnd:',ibnd0,ibnd
    DO ik0=kpoint_initial,kpoint_final
     DO ik = 1, nk
      DO ispin = 1, nspin
        ikk = ik + nk*(ispin-1)
        ikk0 = ik0 + nk*(ispin-1)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !ml=0
        !write (*,*) 'evc read',ik
        
        !IF( nks > 1 ) CALL get_buffer (evc, nwordwfc, iunwfc, ik )
        !         write(*,*) 'npw,npwx,ngk(ik0),ngk(ikk)',npw,npwx,ngk(ik0),ngk(ikk)
        !npw = ngk(ik)
        !         write(*,*) 'npw,npwx,ngk(ik0),ngk(ikk)',npw,npwx,ngk(ik0),ngk(ikk)
        !            CALL init_us_2 (npw, igk_k(1,ik), xk (1, ik), vkb)
        !            CALL calbec ( npw, vkb, evc, becp )
        
        !write (*,*) 'evc2 read',ik
        !CALL get_buffer ( evc2, nwordwfc, iunwfc, ik )
        !write (*,*) 'evc1 read',ik0
        !CALL get_buffer ( evc1, nwordwfc, iunwfc, ik0 )
        !!!!!!!!!!!!write (*,*) "size evc evc1:" , size(evc),size(evc1)
        !!!!!!!!!!!!!!! evc
        
        
        !write(*,*) 'evc1', evc1(1:10,10)
        !write(*,*) 'evc2', evc2(1:10,10)
!    write (*,*) 'ik',ikk0,ikk
!    write(*,*) 'xk,xk0,xk-xk01',xk(1,ik),xk(1,ik0),xk(1,ik)-xk(1,ik0)
!    write(*,*) 'xk,xk0,xk-xk02',xk(2,ik),xk(2,ik0),xk(2,ik)-xk(2,ik0)
!    write(*,*) 'xk,xk0,xk-xk03',xk(3,ik),xk(3,ik0),xk(3,ik)-xk(3,ik0)
        CALL get_buffer ( evc2, nwordwfc, iunwfc, ikk )
        CALL get_buffer ( evc1, nwordwfc, iunwfc, ikk0 )
    
        if (calcmlocal) then
         call calcmdefect_ml_rs(ibnd0,ibnd,ikk0,ikk)
         !call calcmdefect_ml_rd(ibnd0,ibnd,ikk0,ikk)
         !call calcmdefect_ml_ks(ibnd0,ibnd,ikk0,ikk)
         !call calcmdefect_ml_kd(ibnd0,ibnd,ikk0,ikk)
        endif
        if (calcmnonlocal) then
         call calcmdefect_mnl_ks(ibnd0,ibnd,ikk0,ikk)
         !call calcmdefect_mnl_kd(ibnd0,ibnd,ikk0,ikk)
         !call calcmdefect_mnl_rs(ibnd0,ibnd,ikk0,ikk)
         !call calcmdefect_mnl_rd(ibnd0,ibnd,ikk0,ikk)
        endif
        if (calcmcharge) then
         !call calcmdefect_charge(ibnd0,ibnd,ikk0,ikk)

         if (mcharge_dolfa) then
         call calcmdefect_charge_lfa(ibnd0,ibnd,ik0,ik)
         !!call calcmdefect_charge_2dlfa(ibnd0,ibnd,ikk0,ikk)
         !!call calcmdefect_charge_3dlfa(ibnd0,ibnd,ikk0,ikk)
         !!call calcmdefect_charge_qehlfa(ibnd0,ibnd,ikk0,ikk)
         else
         call calcmdefect_charge_nolfa(ibnd0,ibnd,ik0,ik)
         !!call calcmdefect_charge_2dnolfa(ibnd0,ibnd,ikk0,ikk)
         !!call calcmdefect_charge_3dnolfa(ibnd0,ibnd,ikk0,ikk)
         !!call calcmdefect_charge_qehnolfa(ibnd0,ibnd,ikk0,ikk)
         endif
        endif
    
    
      enddo
     enddo
    enddo
    
    END SUBROUTINE calcmdefect_all


