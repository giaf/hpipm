/**************************************************************************************************
*                                                                                                 *
* This file is part of HPIPM.                                                                     *
*                                                                                                 *
* HPIPM -- High-Performance Interior Point Method.                                                *
* Copyright (C) 2017-2018 by Gianluca Frison.                                                     *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
*                                                                                                 *
* This program is free software: you can redistribute it and/or modify                            *
* it under the terms of the GNU General Public License as published by                            *
* the Free Software Foundation, either version 3 of the License, or                               *
* (at your option) any later version                                                              *.
*                                                                                                 *
* This program is distributed in the hope that it will be useful,                                 *
* but WITHOUT ANY WARRANTY; without even the implied warranty of                                  *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                   *
* GNU General Public License for more details.                                                    *
*                                                                                                 *
* You should have received a copy of the GNU General Public License                               *
* along with this program.  If not, see <https://www.gnu.org/licenses/>.                          *
*                                                                                                 *
* The authors designate this particular file as subject to the "Classpath" exception              *
* as provided by the authors in the LICENSE file that accompained this code.                      *
*                                                                                                 *
* Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             *
*                                                                                                 *
**************************************************************************************************/



#include <math.h>

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_s_aux.h>
#include <blasfeo_s_blas.h>

#include <hpipm_s_ocp_qp_dim.h>
#include <hpipm_s_ocp_qp.h>
#include <hpipm_s_ocp_qp_sol.h>
#include <hpipm_s_ocp_qp_ipm.h>
#include <hpipm_s_core_qp_ipm.h>
#include <hpipm_s_core_qp_ipm_aux.h>



#define SINGLE_PRECISION

#define AXPY blasfeo_saxpy
#define COLSC blasfeo_scolsc
#define COMPUTE_LAM_T_QP s_compute_lam_t_qp
#define COMPUTE_GAMMA_GAMMA_QP s_compute_Gamma_gamma_qp
#define COMPUTE_GAMMA_QP s_compute_gamma_qp
#define CORE_QP_IPM_WORKSPACE s_core_qp_ipm_workspace
#define DIAAD_SP blasfeo_sdiaad_sp
#define DIARE blasfeo_sdiare
#define GEAD blasfeo_sgead
#define GECP blasfeo_sgecp
#define GELQF blasfeo_sgelqf
#define GELQF_PD blasfeo_sgelqf_pd
#define GELQF_PD_LA blasfeo_sgelqf_pd_la
#define GELQF_PD_LLA blasfeo_sgelqf_pd_lla
#define GEMM_NT blasfeo_sgemm_nt
#define GEMM_R_DIAG blasfeo_sgemm_nd
#define GEMV_DIAG blasfeo_sgemv_d
#define GEMV_N blasfeo_sgemv_n
#define GEMV_NT blasfeo_sgemv_nt
#define GEMV_T blasfeo_sgemv_t
#define GESE blasfeo_sgese
#define OCP_QP s_ocp_qp
#define OCP_QP_IPM_ARG s_ocp_qp_ipm_arg
#define OCP_QP_IPM_WORKSPACE s_ocp_qp_ipm_workspace
#define OCP_QP_RES s_ocp_qp_res
#define OCP_QP_RES_WORKSPACE s_ocp_qp_res_workspace
#define OCP_QP_DIM s_ocp_qp_dim
#define OCP_QP_SOL s_ocp_qp_sol
#define POTRF_L blasfeo_spotrf_l
#define POTRF_L_MN blasfeo_spotrf_l_mn
#define PRINT_E_MAT s_print_exp_mat
#define PRINT_E_STRVEC blasfeo_print_exp_svec
#define PRINT_E_TRAN_STRVEC blasfeo_print_exp_tran_svec
#define PRINT_STRMAT blasfeo_print_smat
#define PRINT_STRVEC blasfeo_print_svec
#define PRINT_TRAN_STRVEC blasfeo_print_tran_svec
#define REAL float
#define ROWAD_SP blasfeo_srowad_sp
#define ROWEX blasfeo_srowex
#define ROWIN blasfeo_srowin
#define STRMAT blasfeo_smat
#define STRVEC blasfeo_svec
#define SYMV_L blasfeo_ssymv_l
#define SYRK_LN blasfeo_ssyrk_ln
#define SYRK_LN_MN blasfeo_ssyrk_ln_mn
#define SYRK_POTRF_LN_MN blasfeo_ssyrk_spotrf_ln_mn
#define TRCP_L blasfeo_strcp_l
#define TRMM_RLNN blasfeo_strmm_rlnn
#define TRMV_LNN blasfeo_strmv_lnn
#define TRMV_LTN blasfeo_strmv_ltn
#define TRSV_LNN blasfeo_strsv_lnn
#define TRSV_LTN blasfeo_strsv_ltn
#define TRSV_LNN_MN blasfeo_strsv_lnn_mn
#define TRSV_LTN_MN blasfeo_strsv_ltn_mn
#define TRTR_L blasfeo_strtr_l
#define VECAD_SP blasfeo_svecad_sp
#define VECCP blasfeo_sveccp
#define VECCPSC blasfeo_sveccpsc
#define VECEX_SP blasfeo_svecex_sp
#define VECMULACC blasfeo_svecmulacc
#define VECMULDOT blasfeo_svecmuldot
#define VECSC blasfeo_svecsc



#define INIT_VAR_OCP_QP s_init_var_ocp_qp
#define COMPUTE_LIN_RES_OCP_QP s_compute_lin_res_ocp_qp
#define COMPUTE_RES_OCP_QP s_compute_res_ocp_qp
#define FACT_SOLVE_KKT_UNCONSTR_OCP_QP s_fact_solve_kkt_unconstr_ocp_qp
#define COND_SLACKS s_cond_slacks
#define FACT_SOLVE_KKT_STEP_OCP_QP s_fact_solve_kkt_step_ocp_qp
#define FACT_LQ_SOLVE_KKT_STEP_OCP_QP s_fact_lq_solve_kkt_step_ocp_qp
#define SOLVE_KKT_STEP_OCP_QP s_solve_kkt_step_ocp_qp



#include "x_ocp_qp_kkt.c"

