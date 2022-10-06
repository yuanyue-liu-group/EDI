subroutine postprocess()
  
  Use kinds,    only: dp
  use edic_mod, only: bndkp_pair
  real(dp),allocatable::gamma_ki(:) 
  real(dp)::mu,nc
  integer :: nki
  integer, allocatable:: kilist(:)

  nki=0
  do ig = 1,bndkp_pair%npairs
    !if bndkp_pair%kp_idx(ig,1)
    !  bndkp_pair%mc(ig)
    write (*,*)    bndkp_pair%bnd_idx(ig,1),bndkp_pair%kp_idx(ig,1),&
                   bndkp_pair%k_coord(ig,1,1),bndkp_pair%k_coord(ig,2,1),bndkp_pair%k_coord(ig,3,1),&
                   bndkp_pair%e_pair(ig,1),&
                   bndkp_pair%v_pair(ig,1,1),bndkp_pair%v_pair(ig,2,1),bndkp_pair%v_pair(ig,3,1),&
                   bndkp_pair%bnd_idx(ig,2),bndkp_pair%kp_idx(ig,2),&
                   bndkp_pair%k_coord(ig,1,2),bndkp_pair%k_coord(ig,2,2),bndkp_pair%k_coord(ig,3,2),&
                   bndkp_pair%e_pair(ig,2),&
                   bndkp_pair%v_pair(ig,1,2),bndkp_pair%v_pair(ig,2,2),bndkp_pair%v_pair(ig,3,2), &
                   bndkp_pair%wt(ig),bndkp_pair%m(ig),bndkp_pair%mc(ig)
  enddo
end subroutine postprocess
