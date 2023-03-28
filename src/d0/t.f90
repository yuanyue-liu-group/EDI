 call gw_eps_init(epsmat_q0_filename,gw_ng_data_q0 ,gw_nmtx_max_data_q0 ,gw_nmtx_data_q0 ,&
gw_gind_eps2rho_data_q0 ,gw_gind_rho2eps_data_q0 ,&
                                 gw_g_components_data_q0 ,gw_bvec_data_q0 ,gw_blat_data_q0 ,gw_qpts_data_q0 ,gw_nq_data_q0 ,&
              gw_epsmat_diag_data_q0 ,gw_epsmat_full_data_q0 ,gw_q_g_commonsubset_indinrho_q0 ,&
gw_q_g_commonsubset_size_q0,gw_qabs_q0)

 call gw_eps_init(epsmat_q1_filename,gw_ng_data_q1 ,gw_nmtx_max_data_q1 ,gw_nmtx_data_q1 ,&
gw_gind_eps2rho_data_q1 ,gw_gind_rho2eps_data_q1 ,&
                                 gw_g_components_data_q1 ,gw_bvec_data_q1 ,gw_blat_data_q1 ,gw_qpts_data_q1 ,gw_nq_data_q1 ,&
              gw_epsmat_diag_data_q1 ,gw_epsmat_full_data_q1 ,gw_q_g_commonsubset_indinrho_q1 ,&
gw_q_g_commonsubset_size_q1,gw_qabs_q1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


call get_gind_rhoandpsi_gw(gind_rho2psi_gw_q1,gind_psi2rho_gw_q1,gw_ng_data_q1,&
gw_q_g_commonsubset_size_q1,gvec_gw_q1,gw_bvec_data_q1,gw_q_g_commonsubset_indinrho_q1)
call get_gind_rhoandpsi_gw(gind_rho2psi_gw_q0,gind_psi2rho_gw_q0,gw_ng_data_q0,&
gw_q_g_commonsubset_size_q0,gvec_gw_q0,gw_bvec_data_q0,gw_q_g_commonsubset_indinrho_q0)

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
