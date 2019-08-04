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

#include <hpipm_s_tree_ocp_qp.h>
#include <hpipm_s_tree_ocp_qp_sol.h>
#include <hpipm_s_tree_ocp_qp_res.h>
#include <hpipm_s_tree_ocp_qp_ipm.h>
#include <hpipm_s_tree_ocp_qp_kkt.h>
#include <hpipm_s_core_qp_ipm.h>
#include <hpipm_s_core_qp_ipm_aux.h>



#define AXPY blasfeo_saxpy
#define BACKUP_RES_M s_backup_res_m
#define COMPUTE_ALPHA_QP s_compute_alpha_qp
#define COMPUTE_CENTERING_CORRECTION_QP s_compute_centering_correction_qp
#define COMPUTE_CENTERING_QP s_compute_centering_qp
#define COMPUTE_LIN_RES_TREE_OCP_QP s_compute_lin_res_tree_ocp_qp
#define COMPUTE_MU_AFF_QP s_compute_mu_aff_qp
#define COMPUTE_RES_TREE_OCP_QP s_compute_res_tree_ocp_qp
#define CORE_QP_IPM_WORKSPACE s_core_qp_ipm_workspace
#define CREATE_TREE_OCP_QP_RES s_create_tree_ocp_qp_res
#define CREATE_TREE_OCP_QP_SOL s_create_tree_ocp_qp_sol
#define CREATE_STRMAT blasfeo_create_smat
#define CREATE_STRVEC blasfeo_create_svec
#define CREATE_CORE_QP_IPM s_create_core_qp_ipm
#define FACT_LQ_SOLVE_KKT_STEP_TREE_OCP_QP s_fact_lq_solve_kkt_step_tree_ocp_qp
#define FACT_SOLVE_KKT_STEP_TREE_OCP_QP s_fact_solve_kkt_step_tree_ocp_qp
#define FACT_SOLVE_KKT_UNCONSTR_TREE_OCP_QP s_fact_solve_kkt_unconstr_tree_ocp_qp
#define GELQF_WORKSIZE blasfeo_sgelqf_worksize
#define INIT_VAR_TREE_OCP_QP s_init_var_tree_ocp_qp
#define MEMSIZE_CORE_QP_IPM s_memsize_core_qp_ipm
#define MEMSIZE_TREE_OCP_QP_RES s_memsize_tree_ocp_qp_res
#define MEMSIZE_TREE_OCP_QP_SOL s_memsize_tree_ocp_qp_sol
#define REAL float
#define SIZE_STRMAT blasfeo_memsize_smat
#define SIZE_STRVEC blasfeo_memsize_svec
#define SOLVE_KKT_STEP_TREE_OCP_QP s_solve_kkt_step_tree_ocp_qp
#define STRMAT blasfeo_smat
#define STRVEC blasfeo_svec
#define TREE_OCP_QP s_tree_ocp_qp
#define TREE_OCP_QP_DIM s_tree_ocp_qp_dim
#define TREE_OCP_QP_IPM_ARG s_tree_ocp_qp_ipm_arg
#define HPIPM_MODE hpipm_mode
#define TREE_OCP_QP_IPM_WORKSPACE s_tree_ocp_qp_ipm_workspace
#define TREE_OCP_QP_RES s_tree_ocp_qp_res
#define TREE_OCP_QP_RES_WORKSPACE s_tree_ocp_qp_res_workspace
#define TREE_OCP_QP_SOL s_tree_ocp_qp_sol
#define UPDATE_VAR_QP s_update_var_qp
#define VECMULDOT blasfeo_svecmuldot
#define VECNRM_INF blasfeo_svecnrm_inf
#define VECSC blasfeo_svecsc

#define MEMSIZE_TREE_OCP_QP_IPM_ARG s_memsize_tree_ocp_qp_ipm_arg
#define CREATE_TREE_OCP_QP_IPM_ARG s_create_tree_ocp_qp_ipm_arg
#define SET_DEFAULT_TREE_OCP_QP_IPM_ARG s_set_default_tree_ocp_qp_ipm_arg
#define SET_TREE_OCP_QP_IPM_ARG_ITER_MAX s_set_tree_ocp_qp_ipm_arg_iter_max
#define SET_TREE_OCP_QP_IPM_ARG_MU0 s_set_tree_ocp_qp_ipm_arg_mu0
#define SET_TREE_OCP_QP_IPM_ARG_TOL_STAT s_set_tree_ocp_qp_ipm_arg_tol_stat
#define SET_TREE_OCP_QP_IPM_ARG_TOL_EQ s_set_tree_ocp_qp_ipm_arg_tol_eq
#define SET_TREE_OCP_QP_IPM_ARG_TOL_INEQ s_set_tree_ocp_qp_ipm_arg_tol_ineq
#define SET_TREE_OCP_QP_IPM_ARG_TOL_COMP s_set_tree_ocp_qp_ipm_arg_tol_comp
#define SET_TREE_OCP_QP_IPM_ARG_REG_PRIM s_set_tree_ocp_qp_ipm_arg_reg_prim
// ipm
#define MEMSIZE_TREE_OCP_QP_IPM s_memsize_tree_ocp_qp_ipm
#define CREATE_TREE_OCP_QP_IPM s_create_tree_ocp_qp_ipm
#define GET_TREE_OCP_QP_IPM_ITER s_get_tree_ocp_qp_ipm_iter
#define GET_TREE_OCP_QP_IPM_RES_STAT s_get_tree_ocp_qp_ipm_res_stat
#define GET_TREE_OCP_QP_IPM_RES_EQ s_get_tree_ocp_qp_ipm_res_eq
#define GET_TREE_OCP_QP_IPM_RES_INEQ s_get_tree_ocp_qp_ipm_res_ineq
#define GET_TREE_OCP_QP_IPM_RES_COMP s_get_tree_ocp_qp_ipm_res_comp
#define GET_TREE_OCP_QP_IPM_STAT s_get_tree_ocp_qp_ipm_stat
#define SOLVE_TREE_OCP_QP_IPM s_solve_tree_ocp_qp_ipm
#define SOLVE_TREE_OCP_QP_IPM2 s_solve_tree_ocp_qp_ipm2



#include "x_tree_ocp_qp_ipm.c"

