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

#include "../include/hpipm_s_ocp_qp_dim.h"
#include "../include/hpipm_s_ocp_qp.h"
#include "../include/hpipm_s_ocp_qp_sol.h"
#include "../include/hpipm_s_ocp_qp_ipm.h"
#include "../include/hpipm_s_core_qp_ipm.h"
#include "../include/hpipm_s_core_qp_ipm_aux.h"



#define SINGLE_PRECISION

#define AXPY_LIBSTR saxpy_libstr
#define COMPUTE_LAM_T_QP s_compute_lam_t_qp
#define COMPUTE_GAMMA_GAMMA_QP s_compute_Gamma_gamma_qp
#define COMPUTE_GAMMA_QP s_compute_gamma_qp
#define CORE_QP_IPM_WORKSPACE s_core_qp_ipm_workspace
#define DIAAD_SP_LIBSTR sdiaad_sp_libstr
#define GEAD_LIBSTR blasfeo_sgead
#define GECP_LIBSTR blasfeo_sgecp
#define GEMM_R_DIAG_LIBSTR blasfeo_sgemm_nd
#define GEMV_DIAG_LIBSTR blasfeo_sgemv_d
#define GEMV_N_LIBSTR blasfeo_sgemv_n
#define GEMV_NT_LIBSTR blasfeo_sgemv_nt
#define GEMV_T_LIBSTR blasfeo_sgemv_t
#define OCP_QP s_ocp_qp
#define OCP_QP_IPM_WORKSPACE s_ocp_qp_ipm_workspace
#define OCP_QP_RES s_ocp_qp_res
#define OCP_QP_RES_WORKSPACE s_ocp_qp_res_workspace
#define OCP_QP_DIM s_ocp_qp_dim
#define OCP_QP_SOL s_ocp_qp_sol
#define POTRF_L_MN_LIBSTR blasfeo_spotrf_l_mn
#define PRINT_E_MAT s_print_e_mat
#define PRINT_E_STRVEC s_print_e_strvec
#define PRINT_E_TRAN_STRVEC s_print_e_tran_strvec
#define PRINT_STRMAT s_print_strmat
#define PRINT_STRVEC s_print_strvec
#define PRINT_TRAN_STRVEC s_print_tran_strvec
#define REAL float
#define ROWAD_SP_LIBSTR srowad_sp_libstr
#define ROWEX_LIBSTR srowex_libstr
#define ROWIN_LIBSTR srowin_libstr
#define STRMAT blasfeo_smat
#define STRVEC blasfeo_svec
#define SYMV_L_LIBSTR blasfeo_ssymv_l
#define SYRK_POTRF_LN_LIBSTR ssyrk_spotrf_ln_libstr
#define TRCP_L_LIBSTR blasfeo_strcp_l
#define TRMM_RLNN_LIBSTR blasfeo_strmm_rlnn
#define TRMV_LNN_LIBSTR blasfeo_strmv_lnn
#define TRMV_LTN_LIBSTR blasfeo_strmv_ltn
#define TRSV_LNN_LIBSTR blasfeo_strsv_lnn
#define TRSV_LTN_LIBSTR blasfeo_strsv_ltn
#define TRSV_LNN_MN_LIBSTR blasfeo_strsv_lnn_mn
#define TRSV_LTN_MN_LIBSTR blasfeo_strsv_ltn_mn
#define VECAD_SP_LIBSTR svecad_sp_libstr
#define VECCP_LIBSTR blasfeo_sveccp
#define VECEX_SP_LIBSTR svecex_sp_libstr
#define VECMULDOT_LIBSTR svecmuldot_libstr
#define VECSC_LIBSTR blasfeo_svecsc



#define INIT_VAR_OCP_QP s_init_var_ocp_qp
#define COMPUTE_RES_OCP_QP s_compute_res_ocp_qp
#define FACT_SOLVE_KKT_UNCONSTR_OCP_QP s_fact_solve_kkt_unconstr_ocp_qp
#define COND_SLACKS s_cond_slacks
#define FACT_SOLVE_KKT_STEP_OCP_QP s_fact_solve_kkt_step_ocp_qp
#define SOLVE_KKT_STEP_OCP_QP s_solve_kkt_step_ocp_qp



#include "x_ocp_qp_kkt.c"

