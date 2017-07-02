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

#include "../include/hpipm_s_dense_qp.h"
#include "../include/hpipm_s_dense_qp_sol.h"
#include "../include/hpipm_s_dense_qp_ipm_hard.h"
#include "../include/hpipm_s_core_qp_ipm_hard.h"
#include "../include/hpipm_s_core_qp_ipm_hard_aux.h"



#define AXPY_LIBSTR saxpy_libstr
#define COMPUTE_LAM_T_HARD_QP s_compute_lam_t_hard_qp
#define COMPUTE_QX_HARD_QP s_compute_qx_hard_qp
#define COMPUTE_QX_QX_HARD_QP s_compute_Qx_qx_hard_qp
#define DENSE_QP s_dense_qp
#define DENSE_QP_SOL s_dense_qp_sol
#define DIAAD_SP_LIBSTR sdiaad_sp_libstr
#define GECP_LIBSTR sgecp_libstr
#define GEMM_R_DIAG_LIBSTR sgemm_r_diag_libstr
#define GEMV_N_LIBSTR sgemv_n_libstr
#define GEMV_NT_LIBSTR sgemv_nt_libstr
#define GEMV_T_LIBSTR sgemv_t_libstr
#define GESE_LIBSTR sgese_libstr
#define IPM_HARD_CORE_QP_WORKSPACE s_ipm_hard_core_qp_workspace
#define IPM_HARD_DENSE_QP_WORKSPACE s_ipm_hard_dense_qp_workspace
#define POTRF_L_LIBSTR spotrf_l_libstr
#define POTRF_L_MN_LIBSTR spotrf_l_mn_libstr
#define REAL float
#define ROWAD_SP_LIBSTR srowad_sp_libstr
#define ROWEX_LIBSTR srowex_libstr
#define ROWIN_LIBSTR srowin_libstr
#define STRMAT s_strmat
#define STRVEC s_strvec
#define SYMV_L_LIBSTR ssymv_l_libstr
#define SYRK_POTRF_LN_LIBSTR ssyrk_spotrf_ln_libstr
#define TRCP_L_LIBSTR strcp_l_libstr
#define TRSM_RLTN_LIBSTR strsm_rltn_libstr
#define TRSV_LNN_LIBSTR strsv_lnn_libstr
#define TRSV_LTN_LIBSTR strsv_ltn_libstr
#define VECAD_SP_LIBSTR svecad_sp_libstr
#define VECCP_LIBSTR sveccp_libstr
#define VECEX_SP_LIBSTR svecex_sp_libstr
#define VECMULDOT_LIBSTR svecmuldot_libstr
#define VECSC_LIBSTR svecsc_libstr

#define INIT_VAR_HARD_DENSE_QP s_init_var_hard_dense_qp
#define COMPUTE_RES_HARD_DENSE_QP s_compute_res_hard_dense_qp
#define FACT_SOLVE_KKT_UNCONSTR_DENSE_QP s_fact_solve_kkt_unconstr_dense_qp
#define FACT_SOLVE_KKT_STEP_HARD_DENSE_QP s_fact_solve_kkt_step_hard_dense_qp
#define SOLVE_KKT_STEP_HARD_DENSE_QP s_solve_kkt_step_hard_dense_qp



#include "x_dense_qp_kkt.c"

