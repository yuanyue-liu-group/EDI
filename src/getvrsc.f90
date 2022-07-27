
subroutine getvrsc()
    use kinds,    only: dp
USE clib_wrappers,     ONLY: md5_from_file
    use edic_mod,   only: V_file   
    USE cell_base, ONLY: omega, alat, tpiba2, at, bg, tpiba
      Use fft_base,  ONLY: dfftp, dffts
USE scf, ONLY: rho, rho_core, rhog_core, v, vltot, vrs

USE wavefunctions, ONLY : evc,evc1,evc2,evc3,evc4, psic, psic1, psic2, psic3, psic4
COMPLEX(DP), ALLOCATABLE ::  mlat1(:),mlat2(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! vl  supercell
 CHARACTER(len=32)::vf_md5_cksum="NA"
integer ::  unit_pert !rg_spin
character (len=75) ::  perturb_file_name!rg_spin
integer :: iunpot_perturb
character (len=75) :: filpot_perturb
character (len=75) :: title_perturb
character (len=3) ,allocatable :: atm_perturb(:)
integer :: nr1x_perturb, nr2x_perturb, nr3x_perturb, nr1_perturb, nr2_perturb, nr3_perturb, &
nat_perturb, ntyp_perturb, ibrav_perturb, plot_num_perturb,  i_perturb,nkb_perturb
integer :: iunplot_perturb, ios_perturb, ipol_perturb, na_perturb, nt_perturb, &
ir_perturb, ndum_perturb
real(DP) :: celldm_perturb (6), gcutm_perturb, dual_perturb, ecut_perturb,  at_perturb(3,3), omega_perturb, alat_perturb
integer, allocatable:: ityp_perturb (:)
real(DP),allocatable:: zv_perturb (:), tau_perturb (:, :)  , plot_perturb (:)
integer :: ir1mod,ir2mod,ir3mod,irnmod
real(DP):: d1,d2,d3

!call getvloc0(vperturb_filename_p,)
!contains 
!subroutine getvloc0(vperturb_filename)

!integer ::  unit_pert !rg_spin
!character (len=75) ::  perturb_file_name!rg_spin
!integer :: iunpot_perturb
!character (len=75) :: filpot_perturb
!character (len=75) :: title_perturb
!character (len=3) ,allocatable :: atm_perturb(:)
!integer :: nr1x_perturb, nr2x_perturb, nr3x_perturb, nr1_perturb, nr2_perturb, nr3_perturb, &
!nat_perturb, ntyp_perturb, ibrav_perturb, plot_num_perturb,  i_perturb,nkb_perturb
!integer :: iunplot_perturb, ios_perturb, ipol_perturb, na_perturb, nt_perturb, &
!ir_perturb, ndum_perturb
!real(DP) :: celldm_perturb (6), gcutm_perturb, dual_perturb, ecut_perturb,  at_perturb(3,3), omega_perturb, alat_perturb
!integer, allocatable:: ityp_perturb (:)
!real(DP),allocatable:: zv_perturb (:), tau_perturb (:, :)  , plot_perturb (:)
!integer :: ir1mod,ir2mod,ir3mod,irnmod
!real(DP):: d1,d2,d3


          CHARACTER(LEN=256) :: vperturb_filename='vloc.dat'
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
    
write(*,*)'nr1,nat_perturb)',-dfftp%nr1,dfftp%nr1,nat_perturb
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
!end subroutine getvloc0

end subroutine getvrsc

