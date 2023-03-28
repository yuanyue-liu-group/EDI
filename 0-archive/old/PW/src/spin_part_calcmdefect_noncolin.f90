
    SUBROUTINE calcmdefect_noncolin()

      type(V_file) :: V_0, Bxc_1, Bxc_2, Bxc_3, V_p
      real(DP),allocatable ::  V_loc(:,:)

      write(*,"(//A/)") 'Enter calcmdefect_noncolin module '

      write(*,*) eband
      write(*,*) deband
      write(*,*) demet
      write(*,*) ewld
      write(*,*) etxcc
      write(*,*) etxc
      write(*,*) ehart

      ALLOCATE (idx (ngm) )
      ALLOCATE (igtog (ngm) )
      ALLOCATE (gtoig (ngm) )
      idx(:) = 0
      igtog(:) = 0
      !IF( lsda )THEN
       !  nbndup = nbnd
       !  nbnddown = nbnd
       !  nk = nks/2
       !     nspin = 2
      !ELSE
       !  nbndup = nbnd
       !  nbnddown = 0
      nk = nks
       !     nspin = 1
      !ENDIF


    
      DO ik = 1, nk
         ikk = ik
         idx( igk_k(1:ngk(ikk),ikk) ) = 1
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

      
      OPEN(unit=tmp_unit,file = 'calcmdefect.dat',status='old',err=20)
      20 continue
         READ(tmp_unit,calcmcontrol,iostat=ios)
      CLOSE(tmp_unit)
      
      

      V_0%filename = V_0_filename
      Bxc_1%filename = Bxc_1_filename
      Bxc_2%filename = Bxc_2_filename
      Bxc_3%filename = Bxc_3_filename
      V_p%filename = V_p_filename

      call read_perturb_file(V_0)
     ! call read_perturb_file(Bxc_1)
     ! call read_perturb_file(Bxc_2)
      call read_perturb_file(Bxc_3)
      call read_perturb_file(V_p)
      
      allocate(V_loc ( V_0%nr1*V_0%nr2*V_0%nr3, 2))
      
      call get_vloc_colin(V_0, Bxc_3, V_loc)

      allocate(evc1(2*npwx,nbnd))
      allocate(evc2(2*npwx,nbnd))
      !allocate(evc3(2*npwx,nbnd))
      !allocate(evc4(2*npwx,nbnd))
      allocate(psic1(dfftp%nnr))
      allocate(psic2(dfftp%nnr))
      allocate(psic3(dfftp%nnr))
      allocate(psic4(dfftp%nnr))


      ! loop through k points
      write (*,"(/A/)") ' start M calculation k loop'
      ibnd0=bnd_initial
      ibnd=bnd_final
  !    do ibnd0 = bnd_initial, bnd_final
   !      ibnd=ibnd0
         DO ik0=kpoint_initial,kpoint_final
            DO ik = 1, nk
               
               ikk = ik 
               ikk0 = ik0 
         
               CALL get_buffer ( evc2, nwordwfc, iunwfc, ikk )
               CALL get_buffer ( evc1, nwordwfc, iunwfc, ikk0 )

            ! write(*,*) evc1(11:15, 27)/evc1(11+npwx:15+npwx, 27)
            ! write(*,*) evc1(11+npwx:15+npwx, 27)
            ! write(*,*) evc1(2*npwx-3:2*npwx+3, 27)
      
               if (calcmlocal) then
                  call calcmdefect_ml_rs_noncolin(ibnd0,ibnd,ikk0,ikk, V_0, V_loc)
               endif
               if (calcmnonlocal) then
                  call calcmdefect_mnl_ks_noncolin(ibnd0,ibnd,ikk0,ikk, V_0, V_p)
               endif

      1003 format(A24,I6,I6,A6,I6,I6 " ( ",e17.9," , ",e17.9," ) ",e17.9//)
   !  write (stdout,1003) 'M_tot ni ki --> nf kf ', ibnd0,ikk0, '-->', ibnd,ikk, mnl_d-mnl_p+ml_up+ml_down, abs(mnl_d-mnl_p+ml_up+ml_down)
      write (stdout,1003) 'M_tot ni ki --> nf kf ', ibnd0,ikk0, '-->', ibnd,ikk, mnl_d+ml_up+ml_down, abs(mnl_d+ml_up+ml_down)
            enddo
         enddo 
    ! enddo
    END SUBROUTINE calcmdefect_noncolin




