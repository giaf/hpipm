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
#include <blasfeo_s_aux.h>
#include <blasfeo_s_blas.h>

#include <hpipm_s_dense_qp_dim.h>
#include <hpipm_s_dense_qp.h>
#include <hpipm_s_dense_qp_sol.h>
#include <hpipm_s_dense_qp_res.h>
#include <hpipm_s_dense_qp_ipm.h>
#include <hpipm_s_core_qp_ipm.h>
#include <hpipm_s_core_qp_ipm_aux.h>
#include <hpipm_s_dense_qp_kkt.h>



#define AXPY blasfeo_saxpy
#define BACKUP_RES_M s_backup_res_m
#define COMPUTE_ALPHA_QP s_compute_alpha_qp
#define COMPUTE_CENTERING_CORRECTION_QP s_compute_centering_correction_qp
#define COMPUTE_CENTERING_QP s_compute_centering_qp
#define COMPUTE_LIN_RES_DENSE_QP s_compute_lin_res_dense_qp
#define COMPUTE_MU_AFF_QP s_compute_mu_aff_qp
#define COMPUTE_RES_DENSE_QP s_compute_res_dense_qp
#define CORE_QP_IPM_WORKSPACE s_core_qp_ipm_workspace
#define CREATE_CORE_QP_IPM s_create_core_qp_ipm
#define CREATE_DENSE_QP_RES s_create_dense_qp_res
#define CREATE_DENSE_QP_SOL s_create_dense_qp_sol
#define CREATE_STRMAT blasfeo_create_smat
#define CREATE_STRVEC blasfeo_create_svec
#define DENSE_QP s_dense_qp
#define DENSE_QP_IPM_ARG s_dense_qp_ipm_arg
#define HPIPM_MODE hpipm_mode
#define DENSE_QP_IPM_WORKSPACE s_dense_qp_ipm_workspace
#define DENSE_QP_DIM s_dense_qp_dim
#define DENSE_QP_RES s_dense_qp_res
#define DENSE_QP_RES_WORKSPACE s_dense_qp_res_workspace
#define DENSE_QP_SOL s_dense_qp_sol
#define DOT blasfeo_sdot
#define FACT_LQ_SOLVE_KKT_STEP_DENSE_QP s_fact_lq_solve_kkt_step_dense_qp
#define FACT_SOLVE_LU_KKT_STEP_DENSE_QP s_fact_solve_lu_kkt_step_dense_qp
#define FACT_SOLVE_KKT_STEP_DENSE_QP s_fact_solve_kkt_step_dense_qp
#define FACT_SOLVE_KKT_UNCONSTR_DENSE_QP s_fact_solve_kkt_unconstr_dense_qp
#define GELQF_WORKSIZE blasfeo_sgelqf_worksize
#define INIT_VAR_DENSE_QP s_init_var_dense_qp
#define MEMSIZE_CORE_QP_IPM s_memsize_core_qp_ipm
#define MEMSIZE_DENSE_QP_RES s_memsize_dense_qp_res
#define MEMSIZE_DENSE_QP_SOL s_memsize_dense_qp_sol
#define REAL float
#define SIZE_STRMAT blasfeo_memsize_smat
#define SIZE_STRVEC blasfeo_memsize_svec
#define SOLVE_KKT_STEP_DENSE_QP s_solve_kkt_step_dense_qp
#define STRMAT blasfeo_smat
#define STRVEC blasfeo_svec
#define UPDATE_VAR_QP s_update_var_qp
#define VECMUL blasfeo_svecmul
#define VECMULDOT blasfeo_svecmuldot
#define VECNRM_INF blasfeo_svecnrm_inf
#define VECSC blasfeo_svecsc



// arg
#define MEMSIZE_DENSE_QP_IPM_ARG s_memsize_dense_qp_ipm_arg
#define CREATE_DENSE_QP_IPM_ARG s_create_dense_qp_ipm_arg
#define SET_DEFAULT_DENSE_QP_IPM_ARG s_set_default_dense_qp_ipm_arg
#define SET_DENSE_QP_IPM_ARG_ITER_MAX s_set_dense_qp_ipm_arg_iter_max
#define SET_DENSE_QP_IPM_ARG_MU0 s_set_dense_qp_ipm_arg_mu0
#define SET_DENSE_QP_IPM_ARG_TOL_STAT s_set_dense_qp_ipm_arg_tol_stat
#define SET_DENSE_QP_IPM_ARG_TOL_EQ s_set_dense_qp_ipm_arg_tol_eq
#define SET_DENSE_QP_IPM_ARG_TOL_INEQ s_set_dense_qp_ipm_arg_tol_ineq
#define SET_DENSE_QP_IPM_ARG_TOL_COMP s_set_dense_qp_ipm_arg_tol_comp
#define SET_DENSE_QP_IPM_ARG_REG_PRIM s_set_dense_qp_ipm_arg_reg_prim
#define SET_DENSE_QP_IPM_ARG_REG_DUAL s_set_dense_qp_ipm_arg_reg_dual
// ipm
#define MEMSIZE_DENSE_QP_IPM s_memsize_dense_qp_ipm
#define CREATE_DENSE_QP_IPM s_create_dense_qp_ipm
#define GET_DENSE_QP_IPM_ITER s_get_dense_qp_ipm_iter
#define GET_DENSE_QP_IPM_RES_STAT s_get_dense_qp_ipm_res_stat
#define GET_DENSE_QP_IPM_RES_EQ s_get_dense_qp_ipm_res_eq
#define GET_DENSE_QP_IPM_RES_INEQ s_get_dense_qp_ipm_res_ineq
#define GET_DENSE_QP_IPM_RES_COMP s_get_dense_qp_ipm_res_comp
#define GET_DENSE_QP_IPM_STAT s_get_dense_qp_ipm_stat
#define SOLVE_DENSE_QP_IPM s_solve_dense_qp_ipm
#define SOLVE_DENSE_QP_IPM2 s_solve_dense_qp_ipm2



#include "x_dense_qp_ipm.c"

