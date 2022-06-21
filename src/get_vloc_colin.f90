Subroutine get_vloc_colin()
    use kinds,    only: dp
    use v_type,   only: V_file, V_loc, V_0, Bxc_3
    Implicit none




    V_loc(:, 1) = V_0%plot(:) + Bxc_3%plot(:)
    V_loc(:, 2) = V_0%plot(:) - Bxc_3%plot(:)

    
    !write(*,*) 'hello'
    !write(*,*) V_loc(1:100, 2)
End Subroutine get_vloc_colin