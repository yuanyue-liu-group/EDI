

subroutine getepsdata()
    use kinds,    only: dp
    use edic_mod,   only: V_file,eps_data   
    USE cell_base, ONLY: omega, alat, tpiba2, at, bg, tpiba
USE clib_wrappers,     ONLY: md5_from_file
 CHARACTER(len=32)::epsf_md5_cksum="NA"
!real(DP),allocatable:: eps_data (:,:)
          CHARACTER(LEN=256) :: eps_filename='eps.dat'
character (len=75) :: filpot_perturb
    iunpot_perturb=99 
    filpot_perturb=eps_filename
    open (unit = iunpot_perturb, file = filpot_perturb, form = 'formatted', &
         status = 'old', err = 99, iostat = ios_perturb)
    99 call errore ('mloc', 'opening file '//TRIM(filpot_perturb), abs (ios_perturb) )
    
    
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
 
end subroutine getepsdata
