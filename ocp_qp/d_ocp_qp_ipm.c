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



#include <stdlib.h>
#include <stdio.h>
#ifdef USE_C99_MATH
#include <math.h>
#endif

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_d_aux.h>
#include <blasfeo_d_blas.h>

#include <hpipm_d_ocp_qp_dim.h>
#include <hpipm_d_ocp_qp.h>
#include <hpipm_d_ocp_qp_sol.h>
#include <hpipm_d_ocp_qp_ipm.h>
#include <hpipm_d_core_qp_ipm.h>
#include <hpipm_d_core_qp_ipm_aux.h>
#include <hpipm_d_ocp_qp_res.h>
#include <hpipm_d_ocp_qp_kkt.h>
#include <hpipm_aux_string.h>



#define AXPY blasfeo_daxpy
#define BACKUP_RES_M d_backup_res_m
#define COMPUTE_ALPHA_QP d_compute_alpha_qp
#define COMPUTE_CENTERING_CORRECTION_QP d_compute_centering_correction_qp
#define COMPUTE_CENTERING_QP d_compute_centering_qp
#define COMPUTE_MU_AFF_QP d_compute_mu_aff_qp
#define COMPUTE_LIN_RES_OCP_QP d_compute_lin_res_ocp_qp
#define COMPUTE_RES_OCP_QP d_compute_res_ocp_qp
#define CORE_QP_IPM_WORKSPACE d_core_qp_ipm_workspace
#define CREATE_CORE_QP_IPM d_create_core_qp_ipm
#define CREATE_OCP_QP_RES d_create_ocp_qp_res
#define OCP_QP_SOL_CREATE d_ocp_qp_sol_create
#define CREATE_STRMAT blasfeo_create_dmat
#define CREATE_STRVEC blasfeo_create_dvec
#define CVT_STRVEC2VEC blasfeo_unpack_dvec
#define FACT_SOLVE_KKT_STEP_OCP_QP d_fact_solve_kkt_step_ocp_qp
#define FACT_LQ_SOLVE_KKT_STEP_OCP_QP d_fact_lq_solve_kkt_step_ocp_qp
#define FACT_SOLVE_KKT_UNCONSTR_OCP_QP d_fact_solve_kkt_unconstr_ocp_qp
#define GELQF_WORKSIZE blasfeo_dgelqf_worksize
#define INIT_VAR_OCP_QP d_init_var_ocp_qp
#define MEMSIZE_CORE_QP_IPM d_memsize_core_qp_ipm
#define MEMSIZE_OCP_QP_RES d_memsize_ocp_qp_res
#define OCP_QP_SOL_MEMSIZE d_ocp_qp_sol_memsize
#define OCP_QP d_ocp_qp
#define OCP_QP_IPM_ARG d_ocp_qp_ipm_arg
#define HPIPM_MODE hpipm_mode
#define OCP_QP_IPM_WS d_ocp_qp_ipm_ws
#define OCP_QP_RES d_ocp_qp_res
#define OCP_QP_RES_WORKSPACE d_ocp_qp_res_workspace
#define OCP_QP_DIM d_ocp_qp_dim
#define OCP_QP_SOL d_ocp_qp_sol
#define PRINT_E_MAT d_print_exp_mat
#define PRINT_E_STRVEC blasfeo_print_exp_dvec
#define PRINT_E_TRAN_STRVEC blasfeo_print_exp_tran_dvec
#define PRINT_STRMAT d_print_strmat
#define PRINT_STRVEC blasfeo_print_dvec
#define PRINT_TRAN_STRVEC blasfeo_print_tran_dvec
#define REAL double
#define SIZE_STRMAT blasfeo_memsize_dmat
#define SIZE_STRVEC blasfeo_memsize_dvec
#define SOLVE_KKT_STEP_OCP_QP d_solve_kkt_step_ocp_qp
#define STRMAT blasfeo_dmat
#define STRVEC blasfeo_dvec
#define UPDATE_VAR_QP d_update_var_qp
#define VECMULDOT blasfeo_dvecmuldot
#define VECNRM_INF blasfeo_dvecnrm_inf
#define VECSC blasfeo_dvecsc



// arg
#define OCP_QP_IPM_ARG_STRSIZE d_ocp_qp_ipm_arg_strsize
#define SIZEOF_OCP_QP_IPM_ARG d_sizeof_ocp_qp_ipm_arg
#define OCP_QP_IPM_ARG_MEMSIZE d_ocp_qp_ipm_arg_memsize
#define OCP_QP_IPM_ARG_CREATE d_ocp_qp_ipm_arg_create
#define OCP_QP_IPM_ARG_SET_DEFAULT d_ocp_qp_ipm_arg_set_default
#define OCP_QP_IPM_ARG_SET d_ocp_qp_ipm_arg_set
#define OCP_QP_IPM_ARG_SET_ITER_MAX d_ocp_qp_ipm_arg_set_iter_max
#define OCP_QP_IPM_ARG_SET_ALPHA_MIN d_ocp_qp_ipm_arg_set_alpha_min
#define OCP_QP_IPM_ARG_SET_MU0 d_ocp_qp_ipm_arg_set_mu0
#define OCP_QP_IPM_ARG_SET_TOL_STAT d_ocp_qp_ipm_arg_set_tol_stat
#define OCP_QP_IPM_ARG_SET_TOL_EQ d_ocp_qp_ipm_arg_set_tol_eq
#define OCP_QP_IPM_ARG_SET_TOL_INEQ d_ocp_qp_ipm_arg_set_tol_ineq
#define OCP_QP_IPM_ARG_SET_TOL_COMP d_ocp_qp_ipm_arg_set_tol_comp
#define OCP_QP_IPM_ARG_SET_REG_PRIM d_ocp_qp_ipm_arg_set_reg_prim
#define OCP_QP_IPM_ARG_SET_WARM_START d_ocp_qp_ipm_arg_set_warm_start
#define OCP_QP_IPM_ARG_SET_PRED_CORR d_ocp_qp_ipm_arg_set_pred_corr
#define OCP_QP_IPM_ARG_SET_RIC_ALG d_ocp_qp_ipm_arg_set_ric_alg
// ipm
#define OCP_QP_IPM_WS_STRSIZE d_ocp_qp_ipm_ws_strsize
#define OCP_QP_IPM_WS_MEMSIZE d_ocp_qp_ipm_ws_memsize
#define OCP_QP_IPM_WS_CREATE d_ocp_qp_ipm_ws_create
#define OCP_QP_IPM_GET d_ocp_qp_ipm_get
#define OCP_QP_IPM_GET_STATUS d_ocp_qp_ipm_get_status
#define OCP_QP_IPM_GET_ITER d_ocp_qp_ipm_get_iter
#define OCP_QP_IPM_GET_RES_STAT d_ocp_qp_ipm_get_res_stat
#define OCP_QP_IPM_GET_RES_EQ d_ocp_qp_ipm_get_res_eq
#define OCP_QP_IPM_GET_RES_INEQ d_ocp_qp_ipm_get_res_ineq
#define OCP_QP_IPM_GET_RES_COMP d_ocp_qp_ipm_get_res_comp
#define OCP_QP_IPM_GET_STAT d_ocp_qp_ipm_get_stat
#define OCP_QP_IPM_SOLVE d_ocp_qp_ipm_solve



#include "x_ocp_qp_ipm.c"
