/**************************************************************************************************
*                                                                                                 *
* This file is part of HPIPM.                                                                     *
*                                                                                                 *
* HPIPM -- High Performance Interior Point Method.                                                *
* Copyright (C) 2017 by Gianluca Frison.                                                          *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
*                                                                                                 *
* HPIPM is free software; you can redistribute it and/or                                          *
* modify it under the terms of the GNU Lesser General Public                                      *
* License as published by the Free Software Foundation; either                                    *
* version 2.1 of the License, or (at your option) any later version.                              *
*                                                                                                 *
* HPIPM is distributed in the hope that it will be useful,                                        *
* but WITHOUT ANY WARRANTY; without even the implied warranty of                                  *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                                            *
* See the GNU Lesser General Public License for more details.                                     *
*                                                                                                 *
* You should have received a copy of the GNU Lesser General Public                                *
* License along with HPIPM; if not, write to the Free Software                                    *
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA                  *
*                                                                                                 *
* Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             *
*                                                                                                 *
**************************************************************************************************/

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_d_aux.h>
#include <blasfeo_d_blas.h>

#include "../include/hpipm_tree.h"
#include "../include/hpipm_d_tree_ocp_qp.h"
#include "../include/hpipm_d_tree_ocp_qp_sol.h"
#include "../include/hpipm_d_tree_ocp_qp_ipm.h"
#include "../include/hpipm_d_core_qp_ipm.h"
#include "../include/hpipm_d_core_qp_ipm_aux.h"

#define AXPY_LIBSTR daxpy_libstr
#define COMPUTE_LAM_T_QP d_compute_lam_t_qp
#define COMPUTE_GAMMA_GAMMA_QP d_compute_Gamma_gamma_qp
#define COMPUTE_GAMMA_QP d_compute_gamma_qp
#define CORE_QP_IPM_WORKSPACE d_core_qp_ipm_workspace
#define DIAAD_SP_LIBSTR ddiaad_sp_libstr
#define GEAD_LIBSTR blasfeo_dgead
#define GECP_LIBSTR blasfeo_dgecp
#define GEMM_R_DIAG_LIBSTR blasfeo_dgemm_nd
#define GEMV_DIAG_LIBSTR blasfeo_dgemv_d
#define GEMV_N_LIBSTR blasfeo_dgemv_n
#define GEMV_NT_LIBSTR blasfeo_dgemv_nt
#define GEMV_T_LIBSTR blasfeo_dgemv_t
#define POTRF_L_MN_LIBSTR blasfeo_dpotrf_l_mn
#define REAL double
#define ROWAD_SP_LIBSTR drowad_sp_libstr
#define ROWIN_LIBSTR drowin_libstr
#define ROWEX_LIBSTR drowex_libstr
#define STRMAT blasfeo_dmat
#define STRVEC blasfeo_dvec
#define SYMV_L_LIBSTR blasfeo_dsymv_l
#define SYRK_LN_MN_LIBSTR blasfeo_dsyrk_ln_mn
#define SYRK_POTRF_LN_LIBSTR blasfeo_dsyrk_dpotrf_ln
#define TRCP_L_LIBSTR blasfeo_dtrcp_l
#define TREE_OCP_QP d_tree_ocp_qp
#define TREE_OCP_QP_IPM_WORKSPACE d_tree_ocp_qp_ipm_workspace
#define TREE_OCP_QP_SOL d_tree_ocp_qp_sol
#define TRMM_RLNN_LIBSTR blasfeo_dtrmm_rlnn
#define TRMV_LNN_LIBSTR blasfeo_dtrmv_lnn
#define TRMV_LTN_LIBSTR blasfeo_dtrmv_ltn
#define TRSV_LNN_LIBSTR blasfeo_dtrsv_lnn
#define TRSV_LNN_MN_LIBSTR blasfeo_dtrsv_lnn_mn
#define TRSV_LTN_LIBSTR blasfeo_dtrsv_ltn
#define TRSV_LTN_MN_LIBSTR blasfeo_dtrsv_ltn_mn
#define VECAD_SP_LIBSTR dvecad_sp_libstr
#define VECCP_LIBSTR blasfeo_dveccp
#define VECEX_SP_LIBSTR dvecex_sp_libstr
#define VECMULDOT_LIBSTR dvecmuldot_libstr
#define VECSC_LIBSTR blasfeo_dvecsc

#define INIT_VAR_TREE_OCP_QP d_init_var_tree_ocp_qp
#define COMPUTE_RES_TREE_OCP_QP d_compute_res_tree_ocp_qp
#define FACT_SOLVE_KKT_UNCONSTR_TREE_OCP_QP d_fact_solve_kkt_unconstr_tree_ocp_qp
#define FACT_SOLVE_KKT_STEP_TREE_OCP_QP d_fact_solve_kkt_step_tree_ocp_qp
#define SOLVE_KKT_STEP_TREE_OCP_QP d_solve_kkt_step_tree_ocp_qp

#define DOUBLE_PRECISION



#include "x_tree_ocp_qp_kkt.c"
