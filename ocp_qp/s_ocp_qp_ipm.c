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



#include <stdlib.h>
#include <stdio.h>
#ifdef USE_C99_MATH
#include <math.h>
#endif

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
#include "../include/hpipm_s_ocp_qp_res.h"
#include "../include/hpipm_s_ocp_qp_kkt.h"



#define AXPY_LIBSTR blasfeo_saxpy
#define BACKUP_RES_M s_backup_res_m
#define COMPUTE_ALPHA_QP s_compute_alpha_qp
#define COMPUTE_CENTERING_CORRECTION_QP s_compute_centering_correction_qp
#define COMPUTE_CENTERING_QP s_compute_centering_qp
#define COMPUTE_MU_AFF_QP s_compute_mu_aff_qp
#define COMPUTE_LIN_RES_OCP_QP s_compute_lin_res_ocp_qp
#define COMPUTE_RES_OCP_QP s_compute_res_ocp_qp
#define CORE_QP_IPM_WORKSPACE s_core_qp_ipm_workspace
#define CREATE_CORE_QP_IPM s_create_core_qp_ipm
#define CREATE_OCP_QP_RES s_create_ocp_qp_res
#define CREATE_OCP_QP_SOL s_create_ocp_qp_sol
#define CREATE_STRMAT blasfeo_create_smat
#define CREATE_STRVEC blasfeo_create_svec
#define CVT_STRVEC2VEC blasfeo_unpack_svec
#define FACT_SOLVE_KKT_STEP_OCP_QP s_fact_solve_kkt_step_ocp_qp
#define FACT_LQ_SOLVE_KKT_STEP_OCP_QP s_fact_lq_solve_kkt_step_ocp_qp
#define FACT_SOLVE_KKT_UNCONSTR_OCP_QP s_fact_solve_kkt_unconstr_ocp_qp
#define GELQF_WORKSIZE blasfeo_sgelqf_worksize
#define INIT_VAR_OCP_QP s_init_var_ocp_qp
#define MEMSIZE_CORE_QP_IPM s_memsize_core_qp_ipm
#define MEMSIZE_OCP_QP_RES s_memsize_ocp_qp_res
#define MEMSIZE_OCP_QP_SOL s_memsize_ocp_qp_sol
#define OCP_QP s_ocp_qp
#define OCP_QP_IPM_ARG s_ocp_qp_ipm_arg
#define OCP_QP_IPM_MODE s_ocp_qp_ipm_mode
#define OCP_QP_IPM_WORKSPACE s_ocp_qp_ipm_workspace
#define OCP_QP_RES s_ocp_qp_res
#define OCP_QP_RES_WORKSPACE s_ocp_qp_res_workspace
#define OCP_QP_DIM s_ocp_qp_dim
#define OCP_QP_SOL s_ocp_qp_sol
#define PRINT_E_MAT s_print_e_mat
#define PRINT_E_STRVEC blasfeo_print_exp_svec
#define PRINT_E_TRAN_STRVEC blasfeo_print_exp_tran_svec
#define PRINT_STRMAT blasfeo_print_smat
#define PRINT_STRVEC blasfeo_print_svec
#define PRINT_TRAN_STRVEC blasfeo_print_tran_svec
#define REAL float
#define SIZE_STRMAT blasfeo_memsize_smat
#define SIZE_STRVEC blasfeo_memsize_svec
#define SOLVE_KKT_STEP_OCP_QP s_solve_kkt_step_ocp_qp
#define STRMAT blasfeo_smat
#define STRVEC blasfeo_svec
#define UPDATE_VAR_QP s_update_var_qp
#define VECNRM_INF_LIBSTR blasfeo_svecnrm_inf



#define MEMSIZE_OCP_QP_IPM_ARG s_memsize_ocp_qp_ipm_arg
#define CREATE_OCP_QP_IPM_ARG s_create_ocp_qp_ipm_arg
#define SET_DEFAULT_OCP_QP_IPM_ARG s_set_default_ocp_qp_ipm_arg
#define MEMSIZE_OCP_QP_IPM s_memsize_ocp_qp_ipm
#define CREATE_OCP_QP_IPM s_create_ocp_qp_ipm
#define SOLVE_OCP_QP_IPM s_solve_ocp_qp_ipm
#define SOLVE_OCP_QP_IPM2 s_solve_ocp_qp_ipm2



#include "x_ocp_qp_ipm.c"

