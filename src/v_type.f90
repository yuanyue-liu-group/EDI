Module v_type

        Use kinds,    only : dp
        Implicit None

        save

        type :: V_file
                integer :: unit
                character (len=75) :: filename
                character (len=75) :: title
                character (len=3) ,allocatable :: atm(:)
                integer :: nr1x, nr2x, nr3x, nr1, nr2, nr3, &
                        nat, ntyp, ibrav, plot_num,  i,nkb
                integer :: iunplot, ios, ipol, na, nt, &
                        ir, ndum
                real(dp) :: celldm(6), gcutm, dual, ecut,  at(3,3), omega, alat
                integer, allocatable:: ityp(:), intyp(:)
                real(dp),allocatable:: zv(:), tau(:, :)  , plot(:)
        end type V_file

        type(V_file), target :: V_0, Bxc_1, Bxc_2, Bxc_3, V_p

        real(dp), allocatable, target :: &
                V_loc(:, :)



End Module v_type
