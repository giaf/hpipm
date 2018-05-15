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



#define AXPY blasfeo_daxpy
#define COLPE blasfeo_dcolpe
#define COLPEI blasfeo_dcolpei
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
#define DIAAD_SP blasfeo_ddiaad_sp
#define DIAEX blasfeo_ddiaex
#define DIARE blasfeo_ddiare
#define GECP blasfeo_dgecp
#define GEMM_L_DIAG blasfeo_dgemm_dn
#define GEMM_R_DIAG blasfeo_dgemm_nd
#define GEMV_DIAG blasfeo_dgemv_d
#define GEMV_N blasfeo_dgemv_n
#define GEMV_NT blasfeo_dgemv_nt
#define GEMV_T blasfeo_dgemv_t
#define GESE blasfeo_dgese
#define GETRF blasfeo_dgetrf_rowpivot
#define POTRF_L blasfeo_dpotrf_l
#define POTRF_L_MN blasfeo_dpotrf_l_mn
#define PSTRF_L dpstrf_l_libstr
#define REAL double
#define ROWAD_SP blasfeo_drowad_sp
#define ROWEX blasfeo_drowex
#define ROWIN blasfeo_drowin
#define ROWPE blasfeo_drowpe
#define ROWPEI blasfeo_drowpei
#define STRMAT blasfeo_dmat
#define STRVEC blasfeo_dvec
#define SYMV_L blasfeo_dsymv_l
#define SYRK_LN blasfeo_dsyrk_ln
#define SYRK_LN_MN blasfeo_dsyrk_ln_mn
#define SYRK_POTRF_LN blasfeo_dsyrk_dpotrf_ln
#define TRCP_L blasfeo_dtrcp_l
#define TRSM_RLTN blasfeo_dtrsm_rltn
#define TRSV_LNN blasfeo_dtrsv_lnn
#define TRSV_LNU blasfeo_dtrsv_lnu
#define TRSV_LTN blasfeo_dtrsv_ltn
#define TRSV_UNN blasfeo_dtrsv_unn
#define TRTR_L blasfeo_dtrtr_l
#define VECAD_SP blasfeo_dvecad_sp
#define VECCP blasfeo_dveccp
#define VECEX_SP blasfeo_dvecex_sp
#define VECMULACC blasfeo_dvecmulacc
#define VECMULDOT blasfeo_dvecmuldot
#define VECSC blasfeo_dvecsc
#define VECCPSC blasfeo_dveccpsc
#define VECPE blasfeo_dvecpe
#define VECPEI blasfeo_dvecpei

#define INIT_VAR_DENSE_QP d_init_var_dense_qp
#define COMPUTE_RES_DENSE_QP d_compute_res_dense_qp
#define COMPUTE_LIN_RES_DENSE_QP d_compute_lin_res_dense_qp
#define FACT_SOLVE_KKT_UNCONSTR_DENSE_QP d_fact_solve_kkt_unconstr_dense_qp
#define FACT_SOLVE_KKT_STEP_DENSE_QP d_fact_solve_kkt_step_dense_qp
#define FACT_SOLVE_HA_KKT_STEP_DENSE_QP d_fact_solve_ha_kkt_step_dense_qp
#define SOLVE_KKT_STEP_DENSE_QP d_solve_kkt_step_dense_qp



#include "x_dense_qp_kkt.c"
