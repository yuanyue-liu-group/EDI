subroutine getwtdata()
    use kinds,    only: dp
use edic_mod, only: kpair
    use edic_mod,   only: V_file,eps_data   
    USE cell_base, ONLY: omega, alat, tpiba2, at, bg, tpiba
USE clib_wrappers,     ONLY: md5_from_file
 CHARACTER(len=32)::wt_md5_cksum="NA"
!real(DP),allocatable:: eps_data (:,:)
          CHARACTER(LEN=256) :: wt_filename='kq180.dat'
character (len=75) :: filpot_perturb,buf
integer::nlines
    iunpot_perturb=99 
    filpot_perturb=wt_filename
    open (unit = iunpot_perturb, file = filpot_perturb, form = 'formatted', &
         status = 'old', err = 99, iostat = ios_perturb)
    99 call errore ('mloc', 'opening file '//TRIM(filpot_perturb), abs (ios_perturb) )
    
    
nlines=0
DO
    READ (iunpot_perturb,*, END=10)
    nlines = nlines + 1
END DO


    read (iunpot_perturb, '(a)') title_perturb
    read (iunpot_perturb, '(a)') title_perturb
    
    allocate(kpair%idx(nlines-2,2))
    allocate(kpair%coord(nlines-2,3,2))
    allocate(kpair%wt(nlines-2))
    do ig= 1, nlines-2
         read (iunpot_perturb, * ) kpair%idx(ig,1),buf,kpair%idx(ig,2)
    enddo
    CALL md5_from_file(wt_filename, wt_md5_cksum)
    write (*,*) 'wt files:',trim(wt_filename),'  MD5 sum:',wt_md5_cksum
 
end subroutine getwtdata
