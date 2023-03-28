
subroutine interp_eps_1d0(epsmat_inted,gw_q_g_commonsubset_size,ik0,ik)
  USE kinds, ONLY: DP
  Use edic_mod,   only: gw_epsq1_data,gw_epsq0_data
  USE klist , ONLY:  xk
  use splinelib, only: spline,splint
  real(DP),allocatable ::epsint_q0_tmp1(:)
  real(DP),allocatable ::epsint_q0_tmp2(:)
  real(DP),allocatable ::epsint_q0_tmp3(:)
  real(DP),allocatable ::epsint_q0_tmp4(:)
  real(DP) ::epsinttmp1s
  real(DP) ::epsinttmp2s
  real(DP) ::epsinttmp3s
  real(DP) ::epsinttmp4s
  COMPLEX(DP) :: epstmp1,epstmp2
  INTEGER :: gind_gw_eps1,gind_gw_eps2
  real(DP) ::  deltak_para

  complex(DP),intent(inout),allocatable ::epsmat_inted(:,:) ! interpolated eps matrix(epsilon^-1)
  INTEGER ,intent(in):: gw_q_g_commonsubset_size
  !integer(DP),intent(inout),allocatable ::gind_psi2rho_gw(:)
  integer ,intent(in):: ik0,ik


  write(*,*) 'enter interp1d'
  !write(*,*) allocated(gind_psi2rho_gw)
  write(*,*) allocated(epsmat_inted)
!  if (allocated(gind_psi2rho_gw)) then
!      write(*,*) 'enter interp1d'
!  else
!      write(*,*) 'enter interp1d'
!  endif
!
!  deallocate(gind_psi2rho_gw)
  !if (allocated(gind_psi2rho_gw))   deallocate(gind_psi2rho_gw)
  !allocate(gind_psi2rho_gw(size(gw_epsq0_data%gind_psi2rho)))
  !gind_psi2rho_gw(:)=gw_epsq0_data%gind_psi2rho(:)
  !write(*,*) 'gind_psi2rho_gw',gind_psi2rho_gw(1:10)


  if (allocated(epsmat_inted)) deallocate(epsmat_inted)
  allocate(epsmat_inted(gw_q_g_commonsubset_size,gw_q_g_commonsubset_size))
  epsmat_inted(:,:)=(0.0,0.0)

  deltak_para=norm2(xk(1:3,ik0)-xk(1:3,ik))*tpiba ! fixme



  if (allocated(epsint_q0_tmp1)) deallocate(epsint_q0_tmp1)
  if (allocated(epsint_q0_tmp2)) deallocate(epsint_q0_tmp2)
  if (allocated(epsint_q0_tmp3)) deallocate(epsint_q0_tmp3)
  if (allocated(epsint_q0_tmp4)) deallocate(epsint_q0_tmp4)
  allocate(epsint_q0_tmp1(gw_epsq0_data%nq_data(1)))
  allocate(epsint_q0_tmp2(gw_epsq0_data%nq_data(1)))
  allocate(epsint_q0_tmp3(gw_epsq0_data%nq_data(1)))
  allocate(epsint_q0_tmp4(gw_epsq0_data%nq_data(1)))

  write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(1,1)
  write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(1,2)
  write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(1,3)
  write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(1,4)
  write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(1,1)
  write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(2,2)
  write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(3,3)
  write(*,*) 'gw-lin0',shape(epsmat_inted),epsmat_inted(4,4)
 
  do ig1 = 1, gw_q_g_commonsubset_size
    do ig2 = 1, gw_q_g_commonsubset_size

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! simple interpolate
      !epstmp1=(gw_epsmat_full_data(1,ig1,ig2,1,1,iq1),&
      !         gw_epsmat_full_data(2,ig1,ig2,1,1,iq1))
      !epstmp2=(gw_epsmat_full_data(1,ig1,ig2,1,1,iq2),&
      !         gw_epsmat_full_data(2,ig1,ig2,1,1,iq2))
      !      eps_gw=gw_epsmat_full_data(:,ig1,ig2,1,1,iq1)*wq1+&
      !       gw_epsmat_full_data(:,ig1,ig2,1,1,iq2)*wq2
      ! simple interpolate
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!spline interpolate fixme 
      !! change to 2d space matrix interpolate
      !! change boundary condition
      !if( interpolate_smallq1d) then

          !icount=0


          !write(*,*) 'interp 1d'
      
          !*gw_blat_data(1)
          !write(*,*) 'gw3.1 gind',ig1,ig2,gw_gind_rho2eps_data(1,1:5)
          !write(*,*) 'gw3.1 gind',ig1,ig2,gw_gind_rho2eps_data(2,1:5)
          !write(*,*) 'gw3.1 gind',ig1,ig2,gw_gind_rho2eps_data(3,1:5)
          !write(*,*) 'gw3.1',gw_q_g_commonsubset_indinrho(1:5)
          !write(*,*) 'gw3.1',gw_gind_rho2eps_data(1,1:5)
          !write(*,*) 'gw3.1',ig1,ig2
          do iq1=1,gw_epsq0_data%nq_data(1)
             !write(*,*) 'gw3.1.1',iq1, gind_gw_eps2,gw_q_g_commonsubset_indinrho(ig2),gw_gind_rho2eps_data(gw_q_g_commonsubset_indinrho(ig2),iq1)
             gind_gw_eps1=gw_epsq0_data%gind_rho2eps_data(gw_epsq0_data%q_g_commonsubset_indinrho(ig1),iq1)
             gind_gw_eps2=gw_epsq0_data%gind_rho2eps_data(gw_epsq0_data%q_g_commonsubset_indinrho(ig2),iq1)
             !write(*,*) 'gw3.1.2',iq1
             if  (gind_gw_eps1>gw_epsq0_data%nmtx_data(iq1).or. gind_gw_eps2>gw_epsq0_data%nmtx_data(iq1)  )  &
                  write(*,*) 'gindex of eps qpts messedup'
             !write(*,*) 'gw3.1.3',1,gind_gw_eps1,gind_gw_eps2,1,1,iq1
             !write(*,*) 'gw3.1.3',gw_epsmat_full_data(1,gind_gw_eps1,gind_gw_eps2,1,1,iq1)
             epsint_q0_tmp1(iq1)=gw_epsq0_data%epsmat_full_data(1,gind_gw_eps1,gind_gw_eps2,1,1,iq1)
             !write(*,*) 'gw3.1.3',iq1
             epsint_q0_tmp2(iq1)=gw_epsq0_data%epsmat_full_data(2,gind_gw_eps1,gind_gw_eps2,1,1,iq1)
             !write(*,*) 'gw3.1.4',iq1
          enddo
          !write(*,*) 'gw3.2',ig1,ig2
          call  spline(gw_epsq0_data%qabs(:),epsint_q0_tmp1(:),0.0_DP,0.0_DP,epsint_q0_tmp3(:))
          epsinttmp1s= splint(gw_epsq0_data%qabs(:),epsint_q0_tmp1(:),epsint_q0_tmp3(:),deltak_para)
          !write(*,*) 'gw3.2',epsinttmp1s
          if (deltak_para>maxval(gw_epsq0_data%qabs(:)))  epsinttmp1s=minval(epsint_q0_tmp1(:))
          if (deltak_para<minval(gw_epsq0_data%qabs(:)))  epsinttmp1s=1.0
              
          call  spline(gw_epsq0_data%qabs(:),epsint_q0_tmp2(:),0.0_DP,0.0_DP,epsint_q0_tmp4(:))
          epsinttmp2s= splint(gw_epsq0_data%qabs(:),epsint_q0_tmp2(:),epsint_q0_tmp4(:),deltak_para)
          !write(*,*) 'gw3.2',epsinttmp2s
          if (deltak_para>maxval(gw_epsq0_data%qabs(:)))  epsinttmp2s=minval(epsint_q0_tmp2(:))
          if (deltak_para<minval(gw_epsq0_data%qabs(:)))  epsinttmp2s=0.0
          !epsmat_inted(gw_gind_rho2eps_data(ig1,iq1),gw_gind_rho2eps_data(ig2,iq1))=complex(epsinttmp1s,epsinttmp2s)
          epsmat_inted(ig1,ig2)=complex(epsinttmp1s,epsinttmp2s)
          !write(*,*) 'gw3.2',ig1,ig2,epsinttmp1s,epsinttmp2s,epsmat_inted(ig1,ig2)
              
      !endif
      !!spline interpolate fixme 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    enddo
  enddo

  write(*,*) 'gw-lin1 inted'
  write(*,*) 'gw-lin1',shape(epsmat_inted),epsmat_inted(1,1)
  write(*,*) 'gw-lin1',shape(epsmat_inted),epsmat_inted(1,2)
  write(*,*) 'gw-lin1',shape(epsmat_inted),epsmat_inted(1,3)
  write(*,*) 'gw-lin1',shape(epsmat_inted),epsmat_inted(1,4)
  write(*,*) 'gw-lin1',shape(epsmat_inted),epsmat_inted(1,1)
  write(*,*) 'gw-lin1',shape(epsmat_inted),epsmat_inted(2,2)
  write(*,*) 'gw-lin1',shape(epsmat_inted),epsmat_inted(3,3)
  write(*,*) 'gw-lin1',shape(epsmat_inted),epsmat_inted(4,4)

end subroutine  interp_eps_1d0
