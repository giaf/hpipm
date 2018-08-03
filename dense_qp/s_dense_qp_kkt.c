/**************************************************************************************************
*                                                                                                 *
* This file is part of HPIPM.                                                                     *
*                                                                                                 *
* HPIPM -- High Performance Interior Point Method.                                                *
* Copyright (C) 2017 by Gianluca Frison.                                                          *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
*                                                                                                 *
* HPMPC is free software; you can redistribute it and/or                                          *
* modify it under the terms of the GNU Lesser General Public                                      *
* License as published by the Free Software Foundation; either                                    *
* version 2.1 of the License, or (at your option) any later version.                              *
*                                                                                                 *
* HPMPC is distributed in the hope that it will be useful,                                        *
* but WITHOUT ANY WARRANTY; without even the implied warranty of                                  *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                                            *
* See the GNU Lesser General Public License for more details.                                     *
*                                                                                                 *
* You should have received a copy of the GNU Lesser General Public                                *
* License along with HPMPC; if not, write to the Free Software                                    *
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA                  *
*                                                                                                 *
* Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             *
*                                                                                                 *
**************************************************************************************************/



#include <math.h>

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_s_aux.h>
#include <blasfeo_s_blas.h>

#include <hpipm_s_dense_qp_dim.h>
#include <hpipm_s_dense_qp.h>
#include <hpipm_s_dense_qp_sol.h>
#include <hpipm_s_dense_qp_res.h>
#include <hpipm_s_dense_qp_ipm.h>
#include <hpipm_s_core_qp_ipm.h>
#include <hpipm_s_core_qp_ipm_aux.h>



#define AXPY blasfeo_saxpy
#define COLPE blasfeo_scolpe
#define COLPEI blasfeo_scolpei
#define COLSC blasfeo_scolsc
#define COMPUTE_LAM_T_QP s_compute_lam_t_qp
#define COMPUTE_GAMMA_GAMMA_QP s_compute_Gamma_gamma_qp
#define COMPUTE_GAMMA_QP s_compute_gamma_qp
#define CORE_QP_IPM_WORKSPACE s_core_qp_ipm_workspace
#define DENSE_QP_RES s_dense_qp_res
#define DENSE_QP_RES_WORKSPACE s_dense_qp_res_workspace
#define DENSE_QP s_dense_qp
#define DENSE_QP_IPM_ARG s_dense_qp_ipm_arg
#define DENSE_QP_IPM_WORKSPACE s_dense_qp_ipm_workspace
#define DENSE_QP_SOL s_dense_qp_sol
#define DIAAD_SP blasfeo_sdiaad_sp
#define DIAEX blasfeo_sdiaex
#define DIARE blasfeo_sdiare
#define GECP blasfeo_sgecp
#define GELQF blasfeo_sgelqf
#define GELQF_PD_LA blasfeo_sgelqf_pd_la
#define GELQF_PD_LLA blasfeo_sgelqf_pd_lla
#define GELQF_PD blasfeo_sgelqf_pd
#define GELQF_WORKSIZE blasfeo_sgelqf_worksize
#define GEMM_L_DIAG blasfeo_sgemm_dn
#define GEMM_R_DIAG blasfeo_sgemm_nd
#define GEMV_DIAG blasfeo_sgemv_d
#define GEMV_N blasfeo_sgemv_n
#define GEMV_NT blasfeo_sgemv_nt
#define GEMV_T blasfeo_sgemv_t
#define GESE blasfeo_sgese
//#define GETRF blasfeo_sgetrf_rowpivot
#define POTRF_L blasfeo_spotrf_l
#define POTRF_L_MN blasfeo_spotrf_l_mn
#define PSTRF_L spstrf_l_libstr
#define REAL float
#define ROWAD_SP blasfeo_srowad_sp
#define ROWEX blasfeo_srowex
#define ROWIN blasfeo_srowin
#define ROWPE blasfeo_srowpe
#define ROWPEI blasfeo_srowpei
#define STRMAT blasfeo_smat
#define STRVEC blasfeo_svec
#define SYMV_L blasfeo_ssymv_l
#define SYRK_LN blasfeo_ssyrk_ln
#define SYRK_LN_MN blasfeo_ssyrk_ln_mn
#define SYRK_POTRF_LN blasfeo_ssyrk_spotrf_ln
#define SYRK_POTRF_LN_MN blasfeo_ssyrk_spotrf_ln_mn
#define TRCP_L blasfeo_strcp_l
#define TRSM_RLTN blasfeo_strsm_rltn
#define TRSM_RLTU blasfeo_strsm_rltu
#define TRSM_RUNN blasfeo_strsm_runn
#define TRSV_LNN blasfeo_strsv_lnn
#define TRSV_LNU blasfeo_strsv_lnu
#define TRSV_LTN blasfeo_strsv_ltn
#define TRSV_UNN blasfeo_strsv_unn
#define TRTR_L blasfeo_strtr_l
#define TRTR_U blasfeo_strtr_u
#define VECAD_SP blasfeo_svecad_sp
#define VECCP blasfeo_sveccp
#define VECEX_SP blasfeo_svecex_sp
#define VECMULACC blasfeo_svecmulacc
#define VECMULDOT blasfeo_svecmuldot
#define VECSC blasfeo_svecsc
#define VECSE blasfeo_svecse
#define VECCPSC blasfeo_sveccpsc
#define VECPE blasfeo_svecpe
#define VECPEI blasfeo_svecpei

#define INIT_VAR_DENSE_QP s_init_var_dense_qp
#define COMPUTE_RES_DENSE_QP s_compute_res_dense_qp
#define COMPUTE_LIN_RES_DENSE_QP s_compute_lin_res_dense_qp
#define FACT_SOLVE_KKT_UNCONSTR_DENSE_QP s_fact_solve_kkt_unconstr_dense_qp
#define FACT_LQ_SOLVE_KKT_STEP_DENSE_QP s_fact_lq_solve_kkt_step_dense_qp
//#define FACT_SOLVE_LU_KKT_STEP_DENSE_QP s_fact_solve_lu_kkt_step_dense_qp
#define FACT_SOLVE_KKT_STEP_DENSE_QP s_fact_solve_kkt_step_dense_qp
#define SOLVE_KKT_STEP_DENSE_QP s_solve_kkt_step_dense_qp



#include "x_dense_qp_kkt.c"

