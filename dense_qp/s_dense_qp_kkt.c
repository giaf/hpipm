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

#include "../include/hpipm_s_dense_qp_dim.h"
#include "../include/hpipm_s_dense_qp.h"
#include "../include/hpipm_s_dense_qp_sol.h"
#include "../include/hpipm_s_dense_qp_res.h"
#include "../include/hpipm_s_dense_qp_ipm.h"
#include "../include/hpipm_s_core_qp_ipm.h"
#include "../include/hpipm_s_core_qp_ipm_aux.h"



#define AXPY_LIBSTR saxpy_libstr
#define COLPE_LIBSTR scolpe_libstr
#define COLPEI_LIBSTR scolpei_libstr
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
#define DIAAD_SP_LIBSTR sdiaad_sp_libstr
#define DIAEX_LIBSTR sdiaex_libstr
#define DIARE_LIBSTR sdiare_libstr
#define GECP_LIBSTR blasfeo_sgecp
#define GEMM_L_DIAG_LIBSTR blasfeo_sgemm_dn
#define GEMM_R_DIAG_LIBSTR blasfeo_sgemm_nd
#define GEMV_DIAG_LIBSTR blasfeo_sgemv_d
#define GEMV_N_LIBSTR blasfeo_sgemv_n
#define GEMV_NT_LIBSTR blasfeo_sgemv_nt
#define GEMV_T_LIBSTR blasfeo_sgemv_t
#define GESE_LIBSTR blasfeo_sgese
#define POTRF_L_LIBSTR blasfeo_spotrf_l
#define POTRF_L_MN_LIBSTR blasfeo_spotrf_l_mn
#define PSTRF_L_LIBSTR spstrf_l_libstr
#define REAL float
#define ROWAD_SP_LIBSTR srowad_sp_libstr
#define ROWEX_LIBSTR srowex_libstr
#define ROWIN_LIBSTR srowin_libstr
#define ROWPE_LIBSTR srowpe_libstr
#define ROWPEI_LIBSTR srowpei_libstr
#define STRMAT blasfeo_smat
#define STRVEC blasfeo_svec
#define SYMV_L_LIBSTR blasfeo_ssymv_l
#define SYRK_LN_LIBSTR blasfeo_ssyrk_ln
#define SYRK_POTRF_LN_LIBSTR ssyrk_spotrf_ln_libstr
#define TRCP_L_LIBSTR blasfeo_strcp_l
#define TRSM_RLTN_LIBSTR blasfeo_strsm_rltn
#define TRSV_LNN_LIBSTR blasfeo_strsv_lnn
#define TRSV_LTN_LIBSTR blasfeo_strsv_ltn
#define VECAD_SP_LIBSTR svecad_sp_libstr
#define VECCP_LIBSTR blasfeo_sveccp
#define VECEX_SP_LIBSTR svecex_sp_libstr
#define VECMULACC_LIBSTR svecmulacc_libstr
#define VECMULDOT_LIBSTR svecmuldot_libstr
#define VECSC_LIBSTR blasfeo_svecsc
#define VECCPSC_LIBSTR blasfeo_sveccpsc
#define VECPE_LIBSTR svecpe_libstr
#define VECPEI_LIBSTR svecpei_libstr

#define INIT_VAR_DENSE_QP s_init_var_dense_qp
#define COMPUTE_RES_DENSE_QP s_compute_res_dense_qp
#define COMPUTE_LIN_RES_DENSE_QP s_compute_lin_res_dense_qp
#define FACT_SOLVE_KKT_UNCONSTR_DENSE_QP s_fact_solve_kkt_unconstr_dense_qp
#define FACT_SOLVE_KKT_STEP_DENSE_QP s_fact_solve_kkt_step_dense_qp
#define SOLVE_KKT_STEP_DENSE_QP s_solve_kkt_step_dense_qp



#include "x_dense_qp_kkt.c"

