Subroutine get_vloc_dat(v_file_)
    use kinds,    only: dp
    use edic_mod,   only: V_file   
    USE cell_base, ONLY: omega, alat, tpiba2, at, bg, tpiba
    USE clib_wrappers,     ONLY: md5_from_file
    Implicit none

    type(V_file) :: v_file_
    integer :: v_file_i_, v_file_ipol_, v_file_nt_, v_file_ir_, v_file_na_
    integer, external :: find_free_unit
    CHARACTER(len=32)::vf_md5_cksum="NA"


    v_file_%unit = find_free_unit()
    open (unit = v_file_%unit, file = v_file_%filename, form = 'formatted', &
      status = 'old', err = 99, iostat = v_file_%ios)
    99 call errore ('mloc', 'opening file '//TRIM(v_file_%filename), abs (v_file_%ios) )

    read (v_file_%unit, '(a)') v_file_%title
    read (v_file_%unit, * ) v_file_%nr1x, v_file_%nr2x, v_file_%nr3x,&
      v_file_%nr1, v_file_%nr2, v_file_%nr3, v_file_%nat, v_file_%ntyp
    
    allocate(v_file_%pot ( v_file_%nr1*v_file_%nr2*v_file_%nr3))
    allocate(v_file_%ityp (v_file_%nat))
    allocate(v_file_%zv (v_file_%ntyp))
    allocate(v_file_%atm (v_file_%ntyp))
    allocate(v_file_%intyp (v_file_%ntyp))
    allocate(v_file_%tau (3,v_file_%nat))
    
    read (v_file_%unit, * ) v_file_%ibrav, v_file_%celldm
    if (v_file_%ibrav == 0) then
        do v_file_i_ = 1,3
            read ( v_file_%unit, * ) ( v_file_%at(v_file_ipol_,v_file_i_),v_file_ipol_=1,3 )
        enddo
        v_file_%alat=v_file_%celldm(1)
    else
        call latgen(v_file_%ibrav,v_file_%celldm,v_file_%at(1,1),v_file_%at(1,2),v_file_%at(1,3),v_file_%omega)
        v_file_%at(:,:)=v_file_%at(:,:)/alat
    endif

    read (v_file_%unit, * ) v_file_%gcutm, v_file_%dual, v_file_%ecut, v_file_%plot_num
    read (v_file_%unit, '(i4,3x,a2,3x,f5.2)') &
            (v_file_%intyp(v_file_nt_), v_file_%atm(v_file_nt_), v_file_%zv(v_file_nt_), v_file_nt_=1, v_file_%ntyp)
    read (v_file_%unit, *) (v_file_%ndum,  (v_file_%tau (v_file_ipol_, v_file_na_), v_file_ipol_ = 1, 3 ), &
         v_file_%ityp(v_file_na_), v_file_na_ = 1, v_file_%nat)
    read (v_file_%unit, * ) (v_file_%pot (v_file_ir_), v_file_ir_ = 1, v_file_%nr1 * v_file_%nr2 * v_file_%nr3)
    v_file_%tau(:,:)=v_file_%tau(:,:)*v_file_%alat/alat

    CALL md5_from_file(v_file_%filename, vf_md5_cksum)
    write (*,*) 'Potential files:',TRIM(v_file_%filename),'  MD5 sum:',vf_md5_cksum


    
End Subroutine get_vloc_dat


