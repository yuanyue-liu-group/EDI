
subroutine interp_eps_2d(epsmat_inted,gw_q_g_commonsubset_size,gind_psi2rho_gw,ik0,ik)
  USE kinds, ONLY: DP
  Use edic_mod,   only: gw_epsq1_data,gw_epsq0_data
  USE klist , ONLY:  xk
  real(DP) ::epsinttmp1s
  real(DP) ::epsinttmp2s
  real(DP) ::epsinttmp3s
  real(DP) ::epsinttmp4s
  real(DP),allocatable ::w1(:)
  INTEGER :: gind_gw_eps1,gind_gw_eps2
  !real(DP) ::  deltak_para

  complex(DP),intent(inout),allocatable ::epsmat_inted(:,:) ! interpolated eps matrix(epsilon^-1)
  INTEGER ,intent(inout):: gw_q_g_commonsubset_size
  integer(DP),intent(inout),allocatable ::gind_psi2rho_gw(:)
  integer ,intent(in):: ik0,ik
  real(dp)::q1(3)


    if (allocated(gind_psi2rho_gw)) deallocate(gind_psi2rho_gw)
    allocate(gind_psi2rho_gw(size(gw_epsq1_data%gind_psi2rho)))
    gind_psi2rho_gw(:)=gw_epsq1_data%gind_psi2rho(:)
    gw_q_g_commonsubset_size=gw_epsq1_data%q_g_commonsubset_size


    if (allocated(epsmat_inted)) deallocate(epsmat_inted)
    allocate(epsmat_inted(gw_q_g_commonsubset_size,gw_q_g_commonsubset_size))
    epsmat_inted(:,:)=(0.0,0.0)

    !deltak_para=norm2(xk(1:3,ik0)-xk(1:3,ik))*tpiba ! fixme

    write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(1,1)
    write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(1,2)
    write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(1,3)
    write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(1,4)
    write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(1,1)
    write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(2,2)
    write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(3,3)
    write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(4,4)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 2d simple interpolate prepare fixme
    !if(interpolate_2d) then
        !deltakG=norm2(xk(1:3,ik0)-xk(1:3,ik))*tpiba
        !deltakG_para=deltakG
        allocate(w1(gw_epsq1_data%nq_data(1)))
        w1(:)=0.0
        do iq1 = 1, gw_epsq1_data%nq_data(1)
           q1(:)= gw_epsq1_data%qpts_data(1,iq1)*gw_epsq1_data%bvec_data(:,1)+ &
                  gw_epsq1_data%qpts_data(2,iq1)*gw_epsq1_data%bvec_data(:,2)+ &
                  gw_epsq1_data%qpts_data(3,iq1)*gw_epsq1_data%bvec_data(:,3)
           if(abs(norm2((xk(1:3,ik0)-xk(1:3,ik))*tpiba)-norm2(q1(:)))<tpiba*(2*3**.5/3.0)*1.0/nqgrid_gw) then
             w1(iq1)=1/abs(norm2((xk(1:3,ik0)-xk(1:3,ik))*tpiba)-norm2(q1(:)))
           endif
        enddo

        write(*,*) 'gw_debug w1',w1(:)

        if(sum(w1(:))<machine_eps) then
            write(*,*) 'eps 2d interpolation error'
            stop -1
        endif
        !do iq1 = 1, gw_nq_data(1)
        !   if(abs(norm2((xk(1:3,ik0)-xk(1:3,ik))*tpiba)-norm2(q1))<machine_eps*1e-6) 
        !enddo
    !endif
    ! 2d simple interpolate prepare fixme
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do ig1 = 1, gw_q_g_commonsubset_size
      do ig2 = 1, gw_q_g_commonsubset_size
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! 2d simple interpolate  fixme
        !if (interpolate_2d) then
            !write(*,*) 'interp 2d'
            do iq1 = 1, gw_epsq1_data%nq_data(1)
                gind_gw_eps1=gw_epsq1_data%gind_rho2eps_data(gw_epsq1_data%q_g_commonsubset_indinrho(ig1),iq1)
                gind_gw_eps2=gw_epsq1_data%gind_rho2eps_data(gw_epsq1_data%q_g_commonsubset_indinrho(ig2),iq1)
                epsinttmp1s=gw_epsq1_data%epsmat_full_data(1,gind_gw_eps1,gind_gw_eps2,1,1,iq1)
                epsinttmp2s=gw_epsq1_data%epsmat_full_data(2,gind_gw_eps1,gind_gw_eps2,1,1,iq1)
                epsmat_inted(ig1,ig2)=epsmat_inted(ig1,ig2)+complex(epsinttmp1s,epsinttmp2s)*w1(iq1)
            enddo
            epsmat_inted(ig1,ig2)=epsmat_inted(ig1,ig2)/sum(w1(:))
        !endif
        ! 2d simple interpolate  fixme
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
        !write(*,*) 'gw_debug epsmat_inted ig1,ig2,q',epsmat_inted(ig1,ig2),'ig1',ig1,'ig2',ig2,deltakG_para
       
      enddo
    enddo
! get interpolated eps matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end subroutine  interp_eps_2d
