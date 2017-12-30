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
#include <blasfeo_d_aux.h>
#include <blasfeo_d_blas.h>

#include "../include/hpipm_d_dense_qp_dim.h"
#include "../include/hpipm_d_dense_qp.h"
#include "../include/hpipm_d_dense_qp_sol.h"
#include "../include/hpipm_d_dense_qp_res.h"
#include "../include/hpipm_d_dense_qp_ipm.h"
#include "../include/hpipm_d_core_qp_ipm.h"
#include "../include/hpipm_d_core_qp_ipm_aux.h"



#define AXPY_LIBSTR daxpy_libstr
#define COLPE_LIBSTR dcolpe_libstr
#define COLPEI_LIBSTR dcolpei_libstr
#define COMPUTE_LAM_T_QP d_compute_lam_t_qp
#define COMPUTE_GAMMA_GAMMA_QP d_compute_Gamma_gamma_qp
#define COMPUTE_GAMMA_QP d_compute_gamma_qp
#define CORE_QP_IPM_WORKSPACE d_core_qp_ipm_workspace
#define DENSE_QP d_dense_qp
#define DENSE_QP_IPM_ARG d_dense_qp_ipm_arg
#define DENSE_QP_IPM_WORKSPACE d_dense_qp_ipm_workspace
#define DENSE_QP_RES d_dense_qp_res
#define DENSE_QP_RES_WORKSPACE d_dense_qp_res_workspace
#define DENSE_QP_SOL d_dense_qp_sol
#define DIAAD_SP_LIBSTR ddiaad_sp_libstr
#define DIAEX_LIBSTR ddiaex_libstr
#define DIARE_LIBSTR ddiare_libstr
#define GECP_LIBSTR dgecp_libstr
#define GEMM_L_DIAG_LIBSTR dgemm_l_diag_libstr
#define GEMM_R_DIAG_LIBSTR dgemm_r_diag_libstr
#define GEMV_DIAG_LIBSTR dgemv_diag_libstr
#define GEMV_N_LIBSTR dgemv_n_libstr
#define GEMV_NT_LIBSTR dgemv_nt_libstr
#define GEMV_T_LIBSTR dgemv_t_libstr
#define GESE_LIBSTR dgese_libstr
#define POTRF_L_LIBSTR dpotrf_l_libstr
#define POTRF_L_MN_LIBSTR dpotrf_l_mn_libstr
#define PSTRF_L_LIBSTR dpstrf_l_libstr
#define REAL double
#define ROWAD_SP_LIBSTR drowad_sp_libstr
#define ROWEX_LIBSTR drowex_libstr
#define ROWIN_LIBSTR drowin_libstr
#define ROWPE_LIBSTR drowpe_libstr
#define ROWPEI_LIBSTR drowpei_libstr
#define STRMAT blasfeo_dmat
#define STRVEC blasfeo_dvec
#define SYMV_L_LIBSTR dsymv_l_libstr
#define SYRK_LN_LIBSTR dsyrk_ln_libstr
#define SYRK_POTRF_LN_LIBSTR dsyrk_dpotrf_ln_libstr
#define TRCP_L_LIBSTR dtrcp_l_libstr
#define TRSM_RLTN_LIBSTR dtrsm_rltn_libstr
#define TRSV_LNN_LIBSTR dtrsv_lnn_libstr
#define TRSV_LTN_LIBSTR dtrsv_ltn_libstr
#define VECAD_SP_LIBSTR dvecad_sp_libstr
#define VECCP_LIBSTR dveccp_libstr
#define VECEX_SP_LIBSTR dvecex_sp_libstr
#define VECMULACC_LIBSTR dvecmulacc_libstr
#define VECMULDOT_LIBSTR dvecmuldot_libstr
#define VECSC_LIBSTR dvecsc_libstr
#define VECCPSC_LIBSTR dveccpsc_libstr
#define VECPE_LIBSTR dvecpe_libstr
#define VECPEI_LIBSTR dvecpei_libstr

#define INIT_VAR_DENSE_QP d_init_var_dense_qp
#define COMPUTE_RES_DENSE_QP d_compute_res_dense_qp
#define COMPUTE_LIN_RES_DENSE_QP d_compute_lin_res_dense_qp
#define FACT_SOLVE_KKT_UNCONSTR_DENSE_QP d_fact_solve_kkt_unconstr_dense_qp
#define FACT_SOLVE_KKT_STEP_DENSE_QP d_fact_solve_kkt_step_dense_qp
#define SOLVE_KKT_STEP_DENSE_QP d_solve_kkt_step_dense_qp



#include "x_dense_qp_kkt.c"
