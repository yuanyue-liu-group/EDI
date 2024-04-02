subroutine get_wt_data()
    use kinds,    only: dp
    use edic_mod, only: wfc_pair
    use edic_mod, only: bndkp_pair
    use edic_mod,   only: V_file,qeh_eps_data   
    use edic_mod,   only: wt_filename,klist_filename,ev_filename
    USE cell_base, ONLY: omega, alat, tpiba2, at, bg, tpiba
    USE clib_wrappers,     ONLY: md5_from_file
    CHARACTER(len=32)::wt_md5_cksum="NA"
    character (len=75) :: filpot_perturb,buf
    integer::nlines
    iunpot_perturb=99 
    filpot_perturb=wt_filename
    open (unit = iunpot_perturb, file = filpot_perturb, form = 'formatted', &
         status = 'old', err = 99, iostat = ios_perturb)
    99 call errore ('mloc', 'opening file '//TRIM(filpot_perturb), abs (ios_perturb) )
    
    
    nlines=0
    DO
      READ (iunpot_perturb,*,end=10)
      nlines = nlines + 1
    END DO
    10 CLOSE (1)


    !write(*,*) wt_filename, 'lines',nlines
    bndkp_pair%npairs=nlines-1

    REwind (iunpot_perturb)
    read (iunpot_perturb, '(a)') title_perturb
    allocate(bndkp_pair%kp_idx(bndkp_pair%npairs,2))
    allocate(bndkp_pair%bnd_idx(bndkp_pair%npairs,2))
    allocate(bndkp_pair%k_coord(bndkp_pair%npairs,3,2))
    allocate(bndkp_pair%e_pair(bndkp_pair%npairs,2))
    allocate(bndkp_pair%v_pair(bndkp_pair%npairs,3,2))
    allocate(bndkp_pair%wt(bndkp_pair%npairs))
    allocate(bndkp_pair%m(bndkp_pair%npairs))
    allocate(bndkp_pair%mc(bndkp_pair%npairs))
    do ig= 1, bndkp_pair%npairs
      read (iunpot_perturb, * ) &
          bndkp_pair%bnd_idx(ig,1),bndkp_pair%kp_idx(ig,1),buf,&
          bndkp_pair%k_coord(ig,1,1),bndkp_pair%k_coord(ig,2,1),bndkp_pair%k_coord(ig,3,1),buf,&
          bndkp_pair%e_pair(ig,1),buf,&
          bndkp_pair%v_pair(ig,1,1),bndkp_pair%v_pair(ig,2,1),bndkp_pair%v_pair(ig,3,1),buf,buf,&
          bndkp_pair%bnd_idx(ig,2),bndkp_pair%kp_idx(ig,2),buf,&
          bndkp_pair%k_coord(ig,1,2),bndkp_pair%k_coord(ig,2,2),bndkp_pair%k_coord(ig,3,2),buf,&
          bndkp_pair%e_pair(ig,2),buf,&
          bndkp_pair%v_pair(ig,1,2),bndkp_pair%v_pair(ig,2,2),bndkp_pair%v_pair(ig,3,2),buf,buf, &
          bndkp_pair%wt(ig)
    enddo
    CALL md5_from_file(wt_filename, wt_md5_cksum)
    write (*,*) 'WT files:',trim(wt_filename),'  MD5 sum:',wt_md5_cksum
 
end subroutine get_wt_data
