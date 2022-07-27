

   subroutine get_vloc_colin(V_0, Bxc_3, V_loc)
      type(V_file) :: V_0, Bxc_3
      real(DP) :: V_loc(:,:)

     ! allocate(V_loc( V_0%nr1 * V_0%nr2 * V_0%nr3, 2))

      V_loc(:, 1) = V_0%plot(:) + Bxc_3%plot(:)
      V_loc(:, 2) = V_0%plot(:) - Bxc_3%plot(:)

      
   end subroutine get_vloc_colin

!   subroutine get_vloc_noncolin(V_0, Bxc_3)
!      type(V_file) :: V_0, Bxc_1, Bxc_2, Bxc_3
!      complex(DP), allocatable:: get_vloc_noncolin(:, 4)
!
!      allocate(get_vloc_noncolin( V_0%nr1 * Bxc_1%nr2 * Bxc_3%nr3))

      
!   end subroutine get_vloc_noncolin



