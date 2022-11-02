subroutine postprocess()
  
  Use kinds,    only: dp
  use edic_mod, only: bndkp_pair
  real(dp),allocatable::gamma_ki(:) 
  real(dp)::mu,nc
  integer :: nki
  integer, external :: find_free_unit
  integer, allocatable:: kilist(:)
  integer :: tmp_unit

  nki=0
  tmp_unit = find_free_unit()
  OPEN(unit=1111,file = 'pp.dat',status='old',err=20)
  20 continue
  write (*,*)    'ibnd, ik,ki(xyz),ei,vi(xyz);fbnd, fk,kf(xyz),ef,vf(xyz);  wt, m,mc'
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
  close(tmp_unit)
end subroutine postprocess
