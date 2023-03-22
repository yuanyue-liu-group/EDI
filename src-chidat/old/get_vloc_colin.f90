

   subroutine get_vloc_colin(V_d, Bxc_3, V_loc)
    use edic_mod,   only: V_file   
USE kinds, ONLY: DP,sgl
      type(V_file) :: V_d, Bxc_3
      real(DP) :: V_loc(:,:)

     ! allocate(V_loc( V_d%nr1 * V_d%nr2 * V_d%nr3, 2))

      V_loc(:, 1) = V_d%plot(:) + Bxc_3%plot(:)
      V_loc(:, 2) = V_d%plot(:) - Bxc_3%plot(:)

      
   end subroutine get_vloc_colin

!   subroutine get_vloc_noncolin(V_d, Bxc_3)
!      type(V_file) :: V_d, Bxc_1, Bxc_2, Bxc_3
!      complex(DP), allocatable:: get_vloc_noncolin(:, 4)
!
!      allocate(get_vloc_noncolin( V_d%nr1 * Bxc_1%nr2 * Bxc_3%nr3))

      
!   end subroutine get_vloc_noncolin



