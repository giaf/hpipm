/**************************************************************************************************
*                                                                                                 *
* This file is part of HPIPM.                                                                     *
*                                                                                                 *
* HPIPM -- High-Performance Interior Point Method.                                                *
* Copyright (C) 2019 by Gianluca Frison.                                                          *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
*                                                                                                 *
* The 2-Clause BSD License                                                                        *
*                                                                                                 *
* Redistribution and use in source and binary forms, with or without                              *
* modification, are permitted provided that the following conditions are met:                     *
*                                                                                                 *
* 1. Redistributions of source code must retain the above copyright notice, this                  *
*    list of conditions and the following disclaimer.                                             *
* 2. Redistributions in binary form must reproduce the above copyright notice,                    *
*    this list of conditions and the following disclaimer in the documentation                    *
*    and/or other materials provided with the distribution.                                       *
*                                                                                                 *
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND                 *
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED                   *
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE                          *
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR                 *
* ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES                  *
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;                    *
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND                     *
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT                      *
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS                   *
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                                    *
*                                                                                                 *
* Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             *
*                                                                                                 *
**************************************************************************************************/



#include <stdlib.h>
#include <stdio.h>
#ifdef USE_C99_MATH
#include <math.h>
#endif

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_s_aux.h>
#include <blasfeo_s_blas.h>

#include <hpipm_tree.h>
#include <hpipm_aux_string.h>
#include <hpipm_aux_mem.h>
#include <hpipm_s_core_qp_ipm.h>
#include <hpipm_s_core_qp_ipm_aux.h>
#include <hpipm_s_tree_ocp_qp_dim.h>
#include <hpipm_s_tree_ocp_qp.h>
#include <hpipm_s_tree_ocp_qp_sol.h>
#include <hpipm_s_tree_ocp_qp_res.h>
#include <hpipm_s_tree_ocp_qp_ipm.h>
#include <hpipm_s_tree_ocp_qp_kkt.h>
//#include <hpipm_s_tree_ocp_qp_utils.h>
#include <hpipm_s_tree_ocp_qcqp_dim.h>
#include <hpipm_s_tree_ocp_qcqp.h>
#include <hpipm_s_tree_ocp_qcqp_sol.h>
#include <hpipm_s_tree_ocp_qcqp_res.h>
#include <hpipm_s_tree_ocp_qcqp_ipm.h>
//#include <hpipm_s_tree_ocp_qcqp_utils.h>



#define SINGLE_PRECISION
#define BLASFEO_VECEL BLASFEO_SVECEL
#define BLASFEO_MATEL BLASFEO_SMATEL



#define AXPY blasfeo_saxpy
#define AXPBY blasfeo_saxpby
#define BACKUP_RES_M s_backup_res_m
#define COLEX blasfeo_scolex
#define COMPUTE_ALPHA_QCQP s_compute_alpha_qcqp
#define COMPUTE_CENTERING_CORRECTION_QCQP s_compute_centering_correction_qcqp
#define COMPUTE_CENTERING_QCQP s_compute_centering_qcqp
#define COMPUTE_MU_AFF_QCQP s_compute_mu_aff_qcqp
#define COLIN blasfeo_scolin
#define CORE_QP_IPM_WORKSPACE s_core_qp_ipm_workspace
//#define CREATE_CORE_QP_IPM s_create_core_qp_ipm
#define CREATE_STRMAT blasfeo_create_smat
#define CREATE_STRVEC blasfeo_create_svec
#define DOT blasfeo_sdot
#define TREE_OCP_QP_FACT_LQ_SOLVE_KKT_STEP s_tree_ocp_qp_fact_lq_solve_kkt_step
#define TREE_OCP_QP_FACT_SOLVE_LU_KKT_STEP s_tree_ocp_qp_fact_solve_lu_kkt_step
#define TREE_OCP_QP_FACT_SOLVE_KKT_STEP s_tree_ocp_qp_fact_solve_kkt_step
#define TREE_OCP_QP_FACT_SOLVE_KKT_UNCONSTR s_tree_ocp_qp_fact_solve_kkt_unconstr
#define GEAD blasfeo_sgead
#define GECP blasfeo_sgecp
#define GELQF_WORKSIZE blasfeo_dgelqf_worksize
#define GEMV_T blasfeo_sgemv_t
#define HPIPM_MODE hpipm_mode
#define INIT_VAR_OCP_QCQP s_init_var_ocp_qcqp
//#define MEMSIZE_CORE_QP_IPM s_memsize_core_qp_ipm
#define TREE_OCP_QCQP s_tree_ocp_qcqp
#define TREE_OCP_QCQP_IPM_ARG s_tree_ocp_qcqp_ipm_arg
#define TREE_OCP_QCQP_IPM_WS s_tree_ocp_qcqp_ipm_ws
#define TREE_OCP_QCQP_DIM s_tree_ocp_qcqp_dim
#define TREE_OCP_QCQP_RES s_tree_ocp_qcqp_res
#define TREE_OCP_QCQP_RES_COMPUTE s_tree_ocp_qcqp_res_compute
#define TREE_OCP_QCQP_RES_COMPUTE_INF_NORM s_tree_ocp_qcqp_res_compute_inf_norm
#define TREE_OCP_QCQP_RES_COMPUTE_LIN s_tree_ocp_qcqp_res_compute_lin
#define TREE_OCP_QCQP_RES_CREATE s_tree_ocp_qcqp_res_create
#define TREE_OCP_QCQP_RES_MEMSIZE s_tree_ocp_qcqp_res_memsize
#define TREE_OCP_QCQP_RES_WS s_tree_ocp_qcqp_res_ws
#define TREE_OCP_QCQP_RES_WS_CREATE s_tree_ocp_qcqp_res_ws_create
#define TREE_OCP_QCQP_RES_WS_MEMSIZE s_tree_ocp_qcqp_res_ws_memsize
#define TREE_OCP_QCQP_SOL s_tree_ocp_qcqp_sol
#define TREE_OCP_QCQP_SOL_CREATE s_tree_ocp_qcqp_sol_create
#define TREE_OCP_QCQP_SOL_MEMSIZE s_tree_ocp_qcqp_sol_memsize
#define TREE_OCP_QP s_tree_ocp_qp
#define TREE_OCP_QP_CREATE s_tree_ocp_qp_create
#define TREE_OCP_QP_MEMSIZE s_tree_ocp_qp_memsize
#define TREE_OCP_QP_IPM_ARG s_tree_ocp_qp_ipm_arg
#define TREE_OCP_QP_IPM_ARG_CREATE s_tree_ocp_qp_ipm_arg_create
#define TREE_OCP_QP_IPM_ARG_MEMSIZE s_tree_ocp_qp_ipm_arg_memsize
#define TREE_OCP_QP_IPM_ARG_SET_DEFAULT s_tree_ocp_qp_ipm_arg_set_default
#define TREE_OCP_QP_IPM_ARG_SET s_tree_ocp_qp_ipm_arg_set
#define TREE_OCP_QP_IPM_ARG_SET_ITER_MAX s_tree_ocp_qp_ipm_arg_set_iter_max
#define TREE_OCP_QP_IPM_ARG_SET_ALPHA_MIN s_tree_ocp_qp_ipm_arg_set_alpha_min
#define TREE_OCP_QP_IPM_ARG_SET_MU0 s_tree_ocp_qp_ipm_arg_set_mu0
#define TREE_OCP_QP_IPM_ARG_SET_TOL_STAT s_tree_ocp_qp_ipm_arg_set_tol_stat
#define TREE_OCP_QP_IPM_ARG_SET_TOL_EQ s_tree_ocp_qp_ipm_arg_set_tol_eq
#define TREE_OCP_QP_IPM_ARG_SET_TOL_INEQ s_tree_ocp_qp_ipm_arg_set_tol_ineq
#define TREE_OCP_QP_IPM_ARG_SET_TOL_COMP s_tree_ocp_qp_ipm_arg_set_tol_comp
#define TREE_OCP_QP_IPM_ARG_SET_TOL_DUAL_GAP s_tree_ocp_qp_ipm_arg_set_tol_dual_gap
#define TREE_OCP_QP_IPM_ARG_SET_REG_PRIM s_tree_ocp_qp_ipm_arg_set_reg_prim
#define TREE_OCP_QP_IPM_ARG_SET_REG_DUAL s_tree_ocp_qp_ipm_arg_set_reg_dual
#define TREE_OCP_QP_IPM_ARG_SET_WARM_START s_tree_ocp_qp_ipm_arg_set_warm_start
#define TREE_OCP_QP_IPM_ARG_SET_PRED_CORR s_tree_ocp_qp_ipm_arg_set_pred_corr
#define TREE_OCP_QP_IPM_ARG_SET_COND_PRED_CORR s_tree_ocp_qp_ipm_arg_set_cond_pred_corr
#define TREE_OCP_QP_IPM_ARG_SET_COMP_RES_PRED s_tree_ocp_qp_ipm_arg_set_comp_res_pred
#define TREE_OCP_QP_IPM_ARG_SET_COMP_RES_EXIT s_tree_ocp_qp_ipm_arg_set_comp_res_exit
#define TREE_OCP_QP_IPM_ARG_SET_RIC_ALG s_tree_ocp_qp_ipm_arg_set_ric_alg
#define TREE_OCP_QP_IPM_ARG_SET_LAM_MIN s_tree_ocp_qp_ipm_arg_set_lam_min
#define TREE_OCP_QP_IPM_ARG_SET_T_MIN s_tree_ocp_qp_ipm_arg_set_t_min
#define TREE_OCP_QP_IPM_ARG_SET_TAU_MIN s_tree_ocp_qp_ipm_arg_set_tau_min
#define TREE_OCP_QP_IPM_ARG_SET_SPLIT_STEP s_tree_ocp_qp_ipm_arg_set_split_step
#define TREE_OCP_QP_IPM_ARG_SET_T_LAM_MIN s_tree_ocp_qp_ipm_arg_set_t_lam_min
#define TREE_OCP_QP_IPM_ABS_STEP s_tree_ocp_qp_ipm_abs_step
#define TREE_OCP_QP_IPM_DELTA_STEP s_tree_ocp_qp_ipm_delta_step
#define TREE_OCP_QP_IPM_GET_STAT s_tree_ocp_qp_ipm_get_stat
#define TREE_OCP_QP_IPM_GET_STAT_M s_tree_ocp_qp_ipm_get_stat_m
#define TREE_OCP_QP_IPM_WS s_tree_ocp_qp_ipm_ws
#define TREE_OCP_QP_IPM_WS_CREATE s_tree_ocp_qp_ipm_ws_create
#define TREE_OCP_QP_IPM_WS_MEMSIZE s_tree_ocp_qp_ipm_ws_memsize
#define TREE_OCP_QP_RES s_tree_ocp_qp_res
#define TREE_OCP_QP_SOL s_tree_ocp_qp_sol
#define TREE_OCP_QP_SOL_CREATE s_tree_ocp_qp_sol_create
#define TREE_OCP_QP_SOL_MEMSIZE s_tree_ocp_qp_sol_memsize
#define REAL float
#define SIZE_STRMAT blasfeo_memsize_smat
#define SIZE_STRVEC blasfeo_memsize_svec
#define SOLVE_KKT_STEP_OCP_QCQP s_solve_kkt_step_ocp_qcqp
#define STRMAT blasfeo_smat
#define STRVEC blasfeo_svec
#define SYMV_L blasfeo_ssymv_l
#define UPDATE_VAR_QCQP s_update_var_qcqp
#define VECCP blasfeo_sveccp
#define VECMUL blasfeo_svecmul
#define VECMULDOT blasfeo_svecmuldot
#define VECNRM_INF blasfeo_svecnrm_inf
#define VECSC blasfeo_svecsc
#define VECSE blasfeo_svecse



// arg
#define TREE_OCP_QCQP_IPM_ARG_STRSIZE s_tree_ocp_qcqp_ipm_arg_strsize
#define TREE_OCP_QCQP_IPM_ARG_MEMSIZE s_tree_ocp_qcqp_ipm_arg_memsize
#define TREE_OCP_QCQP_IPM_ARG_CREATE s_tree_ocp_qcqp_ipm_arg_create
#define TREE_OCP_QCQP_IPM_ARG_SET_DEFAULT s_tree_ocp_qcqp_ipm_arg_set_default
#define TREE_OCP_QCQP_IPM_ARG_SET s_tree_ocp_qcqp_ipm_arg_set
#define TREE_OCP_QCQP_IPM_ARG_SET_ITER_MAX s_tree_ocp_qcqp_ipm_arg_set_iter_max
#define TREE_OCP_QCQP_IPM_ARG_SET_ALPHA_MIN s_tree_ocp_qcqp_ipm_arg_set_alpha_min
#define TREE_OCP_QCQP_IPM_ARG_SET_MU0 s_tree_ocp_qcqp_ipm_arg_set_mu0
#define TREE_OCP_QCQP_IPM_ARG_SET_TOL_STAT s_tree_ocp_qcqp_ipm_arg_set_tol_stat
#define TREE_OCP_QCQP_IPM_ARG_SET_TOL_EQ s_tree_ocp_qcqp_ipm_arg_set_tol_eq
#define TREE_OCP_QCQP_IPM_ARG_SET_TOL_INEQ s_tree_ocp_qcqp_ipm_arg_set_tol_ineq
#define TREE_OCP_QCQP_IPM_ARG_SET_TOL_COMP s_tree_ocp_qcqp_ipm_arg_set_tol_comp
#define TREE_OCP_QCQP_IPM_ARG_SET_TOL_DUAL_GAP s_tree_ocp_qcqp_ipm_arg_set_tol_dual_gap
#define TREE_OCP_QCQP_IPM_ARG_SET_REG_PRIM s_tree_ocp_qcqp_ipm_arg_set_reg_prim
#define TREE_OCP_QCQP_IPM_ARG_SET_WARM_START s_tree_ocp_qcqp_ipm_arg_set_warm_start
#define TREE_OCP_QCQP_IPM_ARG_SET_PRED_CORR s_tree_ocp_qcqp_ipm_arg_set_pred_corr
#define TREE_OCP_QCQP_IPM_ARG_SET_COND_PRED_CORR s_tree_ocp_qcqp_ipm_arg_set_cond_pred_corr
#define TREE_OCP_QCQP_IPM_ARG_SET_COMP_RES_PRED s_tree_ocp_qcqp_ipm_arg_set_comp_res_pred
#define TREE_OCP_QCQP_IPM_ARG_SET_COMP_RES_EXIT s_tree_ocp_qcqp_ipm_arg_set_comp_res_exit
#define TREE_OCP_QCQP_IPM_ARG_SET_RIC_ALG s_tree_ocp_qcqp_ipm_arg_set_ric_alg
#define TREE_OCP_QCQP_IPM_ARG_SET_LAM_MIN s_tree_ocp_qcqp_ipm_arg_set_lam_min
#define TREE_OCP_QCQP_IPM_ARG_SET_T_MIN s_tree_ocp_qcqp_ipm_arg_set_t_min
#define TREE_OCP_QCQP_IPM_ARG_SET_TAU_MIN s_tree_ocp_qcqp_ipm_arg_set_tau_min
#define TREE_OCP_QCQP_IPM_ARG_SET_SPLIT_STEP s_tree_ocp_qcqp_ipm_arg_set_split_step
#define TREE_OCP_QCQP_IPM_ARG_SET_T_LAM_MIN s_tree_ocp_qcqp_ipm_arg_set_t_lam_min
// ipm
#define TREE_OCP_QCQP_IPM_WS_STRSIZE s_tree_ocp_qcqp_ipm_ws_strsize
#define TREE_OCP_QCQP_IPM_WS_MEMSIZE s_tree_ocp_qcqp_ipm_ws_memsize
#define TREE_OCP_QCQP_IPM_WS_CREATE s_tree_ocp_qcqp_ipm_ws_create
#define TREE_OCP_QCQP_IPM_GET s_tree_ocp_qcqp_ipm_get
#define TREE_OCP_QCQP_IPM_GET_STATUS s_tree_ocp_qcqp_ipm_get_status
#define TREE_OCP_QCQP_IPM_GET_ITER s_tree_ocp_qcqp_ipm_get_iter
#define TREE_OCP_QCQP_IPM_GET_MAX_RES_STAT s_tree_ocp_qcqp_ipm_get_max_res_stat
#define TREE_OCP_QCQP_IPM_GET_MAX_RES_EQ s_tree_ocp_qcqp_ipm_get_max_res_eq
#define TREE_OCP_QCQP_IPM_GET_MAX_RES_INEQ s_tree_ocp_qcqp_ipm_get_max_res_ineq
#define TREE_OCP_QCQP_IPM_GET_MAX_RES_COMP s_tree_ocp_qcqp_ipm_get_max_res_comp
#define TREE_OCP_QCQP_IPM_GET_DUAL_GAP s_tree_ocp_qcqp_ipm_get_dual_gap
#define TREE_OCP_QCQP_IPM_GET_OBJ s_tree_ocp_qcqp_ipm_get_obj
#define TREE_OCP_QCQP_IPM_GET_STAT s_tree_ocp_qcqp_ipm_get_stat
#define TREE_OCP_QCQP_IPM_GET_STAT_M s_tree_ocp_qcqp_ipm_get_stat_m
#define TREE_OCP_QCQP_INIT_VAR s_tree_ocp_qcqp_init_var
#define TREE_OCP_QCQP_APPROX_QP s_tree_ocp_qcqp_approx_qp
#define TREE_OCP_QCQP_UPDATE_QP s_tree_ocp_qcqp_update_qp
#define TREE_OCP_QCQP_UPDATE_QP_ABS_STEP s_tree_ocp_qcqp_update_qp_abs_step
#define TREE_OCP_QCQP_SOL_CONV_QP_SOL s_tree_ocp_qcqp_sol_conv_qp_sol
#define TREE_OCP_QP_SOL_CONV_QCQP_SOL s_tree_ocp_qp_sol_conv_qcqp_sol
#define TREE_OCP_QCQP_RES_CONV_QP_RES s_tree_ocp_qcqp_res_conv_qp_res
#define TREE_OCP_QCQP_IPM_SOLVE s_tree_ocp_qcqp_ipm_solve
#define TREE_OCP_QCQP_IPM_PREDICT s_tree_ocp_qcqp_ipm_predict
#define TREE_OCP_QCQP_IPM_SENS s_tree_ocp_qcqp_ipm_sens

#include "x_tree_ocp_qcqp_ipm.c"




