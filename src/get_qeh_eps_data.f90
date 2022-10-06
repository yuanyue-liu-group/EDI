subroutine get_qeh_eps_data()
    use kinds,    only: dp
    use edic_mod,   only: V_file,qeh_eps_data   
    USE cell_base, ONLY: omega, alat, tpiba2, at, bg, tpiba
    USE clib_wrappers,     ONLY: md5_from_file
    CHARACTER(len=32)::epsf_md5_cksum="NA"
    CHARACTER(LEN=256) :: qeh_eps_filename='eps.dat'
    character (len=75) :: filpot_perturb
    iunpot_perturb=99 
    filpot_perturb=qeh_eps_filename
    open (unit = iunpot_perturb, file = filpot_perturb, form = 'formatted', &
         status = 'old', err = 99, iostat = ios_perturb)
    99 call errore ('mloc', 'opening file '//TRIM(filpot_perturb), abs (ios_perturb) )
    
    
    read (iunpot_perturb, '(a)') title_perturb
    read (iunpot_perturb, * ) nepslines
    
    allocate(qeh_eps_data(2,nepslines))
    do ig= 1, nepslines
         read (iunpot_perturb, * ) qeh_eps_data(1,ig),qeh_eps_data(2,ig)
    enddo
    write (*,*) 'eps lines  ', nepslines
    write (*,*) 'eps data  ', qeh_eps_data(1,1),qeh_eps_data(2,1)
    write (*,*) 'eps data  ', qeh_eps_data(1,2),qeh_eps_data(2,2)
    write (*,*) 'eps data  ', qeh_eps_data(1,3),qeh_eps_data(2,3)
    write (*,*) 'eps data  ', qeh_eps_data(1,7),qeh_eps_data(2,7)
    !k0screen=k0screen_read
    
    CALL md5_from_file(qeh_eps_filename, epsf_md5_cksum)
    write (*,*) 'eps files:',trim(qeh_eps_filename),'  MD5 sum:',epsf_md5_cksum
    !write (*,*) 'dv readin-vrs', plot_perturb(:)-vrs(:,1)
    !write (*,*) 'dv readin-vrs', vrs(:,1)
    !write (*,*) 'dv readin-vrs', plot_perturb(:)
    !
 
end subroutine get_qeh_eps_data
