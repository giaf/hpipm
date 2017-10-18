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



#if defined(RUNTIME_CHECKS)
#include <stdlib.h>
#include <stdio.h>
#endif

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_s_aux.h>

#include "../include/hpipm_s_ocp_qp.h"
#include "../include/hpipm_s_ocp_qp_sol.h"
#include "../include/hpipm_s_ocp_qp_ipm.h"
#include "../include/hpipm_s_core_qp_ipm.h"
#include "../include/hpipm_s_core_qp_ipm_aux.h"
#include "../include/hpipm_s_ocp_qp_kkt.h"



#define COMPUTE_ALPHA_QP s_compute_alpha_qp
#define COMPUTE_CENTERING_CORRECTION_QP s_compute_centering_correction_qp
#define COMPUTE_MU_AFF_QP s_compute_mu_aff_qp
#define COMPUTE_RES_OCP_QP s_compute_res_ocp_qp
#define CORE_QP_IPM_WORKSPACE s_core_qp_ipm_workspace
#define CREATE_CORE_QP_IPM s_create_core_qp_ipm
#define CREATE_STRMAT s_create_strmat
#define CREATE_STRVEC s_create_strvec
#define CVT_STRVEC2VEC s_cvt_strvec2vec
#define FACT_SOLVE_KKT_STEP_OCP_QP s_fact_solve_kkt_step_ocp_qp
#define FACT_SOLVE_KKT_UNCONSTR_OCP_QP s_fact_solve_kkt_unconstr_ocp_qp
#define INIT_VAR_OCP_QP s_init_var_ocp_qp
#define MEMSIZE_CORE_QP_IPM s_memsize_core_qp_ipm
#define OCP_QP s_ocp_qp
#define OCP_QP_IPM_WORKSPACE s_ocp_qp_ipm_workspace
#define OCP_QP_IPM_ARG s_ocp_qp_ipm_arg
#define OCP_QP_RES s_ocp_qp_res
#define OCP_QP_SOL s_ocp_qp_sol
#define PRINT_E_MAT s_print_e_mat
#define PRINT_E_STRVEC s_print_e_strvec
#define PRINT_E_TRAN_STRVEC s_print_e_tran_strvec
#define PRINT_STRMAT s_print_strmat
#define PRINT_STRVEC s_print_strvec
#define PRINT_TRAN_STRVEC s_print_tran_strvec
#define REAL float
#define SIZE_STRMAT s_size_strmat
#define SIZE_STRVEC s_size_strvec
#define SOLVE_KKT_STEP_OCP_QP s_solve_kkt_step_ocp_qp
#define STRMAT s_strmat
#define STRVEC s_strvec
#define UPDATE_VAR_QP s_update_var_qp
#define VECNRM_INF_LIBSTR svecnrm_inf_libstr



#define MEMSIZE_OCP_QP_IPM_ARG s_memsize_ocp_qp_ipm_arg
#define CREATE_OCP_QP_IPM_ARG s_create_ocp_qp_ipm_arg
#define SET_DEFAULT_OCP_QP_IPM_ARG s_set_default_ocp_qp_ipm_arg
#define MEMSIZE_OCP_QP_RES s_memsize_ocp_qp_res
#define CREATE_OCP_QP_RES s_create_ocp_qp_res
#define MEMSIZE_OCP_QP_IPM s_memsize_ocp_qp_ipm
#define CREATE_OCP_QP_IPM s_create_ocp_qp_ipm
#define SOLVE_OCP_QP_IPM s_solve_ocp_qp_ipm
#define SOLVE_OCP_QP_IPM2 s_solve_ocp_qp_ipm2
#define CVT_OCP_QP_RES_TO_COLMAJ s_cvt_ocp_qp_res_to_colmaj
#define CVT_OCP_QP_RES_TO_ROWMAJ s_cvt_ocp_qp_res_to_rowmaj



#include "x_ocp_qp_ipm.c"

