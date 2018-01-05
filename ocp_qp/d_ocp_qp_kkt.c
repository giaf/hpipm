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

#include "../include/hpipm_d_ocp_qp_dim.h"
#include "../include/hpipm_d_ocp_qp.h"
#include "../include/hpipm_d_ocp_qp_sol.h"
#include "../include/hpipm_d_ocp_qp_ipm.h"
#include "../include/hpipm_d_core_qp_ipm.h"
#include "../include/hpipm_d_core_qp_ipm_aux.h"



#define DOUBLE_PRECISION

#define AXPY_LIBSTR blasfeo_daxpy
#define COMPUTE_LAM_T_QP d_compute_lam_t_qp
#define COMPUTE_GAMMA_GAMMA_QP d_compute_Gamma_gamma_qp
#define COMPUTE_GAMMA_QP d_compute_gamma_qp
#define CORE_QP_IPM_WORKSPACE d_core_qp_ipm_workspace
#define DIAAD_SP_LIBSTR blasfeo_ddiaad_sp
#define GEAD_LIBSTR blasfeo_dgead
#define GECP_LIBSTR blasfeo_dgecp
#define GEMM_R_DIAG_LIBSTR blasfeo_dgemm_nd
#define GEMV_DIAG_LIBSTR blasfeo_dgemv_d
#define GEMV_N_LIBSTR blasfeo_dgemv_n
#define GEMV_NT_LIBSTR blasfeo_dgemv_nt
#define GEMV_T_LIBSTR blasfeo_dgemv_t
#define OCP_QP d_ocp_qp
#define OCP_QP_IPM_WORKSPACE d_ocp_qp_ipm_workspace
#define OCP_QP_RES d_ocp_qp_res
#define OCP_QP_RES_WORKSPACE d_ocp_qp_res_workspace
#define OCP_QP_DIM d_ocp_qp_dim
#define OCP_QP_SOL d_ocp_qp_sol
#define POTRF_L_MN_LIBSTR blasfeo_dpotrf_l_mn
#define PRINT_E_MAT d_print_e_mat
#define PRINT_E_STRVEC blasfeo_print_exp_dvec
#define PRINT_E_TRAN_STRVEC blasfeo_print_exp_tran_dvec
#define PRINT_STRMAT d_print_strmat
#define PRINT_STRVEC blasfeo_print_dvec
#define PRINT_TRAN_STRVEC blasfeo_print_tran_dvec
#define REAL double
#define ROWAD_SP_LIBSTR blasfeo_drowad_sp
#define ROWEX_LIBSTR blasfeo_drowex
#define ROWIN_LIBSTR blasfeo_drowin
#define STRMAT blasfeo_dmat
#define STRVEC blasfeo_dvec
#define SYMV_L_LIBSTR blasfeo_dsymv_l
#define SYRK_POTRF_LN_LIBSTR blasfeo_dsyrk_dpotrf_ln
#define TRCP_L_LIBSTR blasfeo_dtrcp_l
#define TRMM_RLNN_LIBSTR blasfeo_dtrmm_rlnn
#define TRMV_LNN_LIBSTR blasfeo_dtrmv_lnn
#define TRMV_LTN_LIBSTR blasfeo_dtrmv_ltn
#define TRSV_LNN_LIBSTR blasfeo_dtrsv_lnn
#define TRSV_LNN_MN_LIBSTR blasfeo_dtrsv_lnn_mn
#define TRSV_LTN_LIBSTR blasfeo_dtrsv_ltn
#define TRSV_LTN_MN_LIBSTR blasfeo_dtrsv_ltn_mn
#define VECAD_SP_LIBSTR blasfeo_dvecad_sp
#define VECCP_LIBSTR blasfeo_dveccp
#define VECEX_SP_LIBSTR blasfeo_dvecex_sp
#define VECMULDOT_LIBSTR blasfeo_dvecmuldot
#define VECSC_LIBSTR blasfeo_dvecsc



#define INIT_VAR_OCP_QP d_init_var_ocp_qp
#define COMPUTE_RES_OCP_QP d_compute_res_ocp_qp
#define FACT_SOLVE_KKT_UNCONSTR_OCP_QP d_fact_solve_kkt_unconstr_ocp_qp
#define FACT_SOLVE_KKT_STEP_OCP_QP d_fact_solve_kkt_step_ocp_qp
#define SOLVE_KKT_STEP_OCP_QP d_solve_kkt_step_ocp_qp



#include "x_ocp_qp_kkt.c"
