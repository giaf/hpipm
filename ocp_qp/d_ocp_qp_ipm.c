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
#include <blasfeo_d_aux.h>

#include "../include/hpipm_d_ocp_qp_size.h"
#include "../include/hpipm_d_ocp_qp.h"
#include "../include/hpipm_d_ocp_qp_sol.h"
#include "../include/hpipm_d_ocp_qp_ipm.h"
#include "../include/hpipm_d_core_qp_ipm.h"
#include "../include/hpipm_d_core_qp_ipm_aux.h"
#include "../include/hpipm_d_ocp_qp_kkt.h"



#define COMPUTE_ALPHA_QP d_compute_alpha_qp
#define COMPUTE_CENTERING_CORRECTION_QP d_compute_centering_correction_qp
#define COMPUTE_MU_AFF_QP d_compute_mu_aff_qp
#define COMPUTE_RES_OCP_QP d_compute_res_ocp_qp
#define CORE_QP_IPM_WORKSPACE d_core_qp_ipm_workspace
#define CREATE_CORE_QP_IPM d_create_core_qp_ipm
#define CREATE_STRMAT d_create_strmat
#define CREATE_STRVEC d_create_strvec
#define CVT_STRVEC2VEC d_cvt_strvec2vec
#define FACT_SOLVE_KKT_STEP_OCP_QP d_fact_solve_kkt_step_ocp_qp
#define FACT_SOLVE_KKT_UNCONSTR_OCP_QP d_fact_solve_kkt_unconstr_ocp_qp
#define INIT_VAR_OCP_QP d_init_var_ocp_qp
#define MEMSIZE_CORE_QP_IPM d_memsize_core_qp_ipm
#define OCP_QP d_ocp_qp
#define OCP_QP_IPM_WORKSPACE d_ocp_qp_ipm_workspace
#define OCP_QP_IPM_ARG d_ocp_qp_ipm_arg
#define OCP_QP_RES d_ocp_qp_res
#define OCP_QP_SIZE d_ocp_qp_size
#define OCP_QP_SOL d_ocp_qp_sol
#define PRINT_E_MAT d_print_e_mat
#define PRINT_E_STRVEC d_print_e_strvec
#define PRINT_E_TRAN_STRVEC d_print_e_tran_strvec
#define PRINT_STRMAT d_print_strmat
#define PRINT_STRVEC d_print_strvec
#define PRINT_TRAN_STRVEC d_print_tran_strvec
#define REAL double
#define SIZE_STRMAT d_size_strmat
#define SIZE_STRVEC d_size_strvec
#define SOLVE_KKT_STEP_OCP_QP d_solve_kkt_step_ocp_qp
#define STRMAT d_strmat
#define STRVEC d_strvec
#define UPDATE_VAR_QP d_update_var_qp
#define VECNRM_INF_LIBSTR dvecnrm_inf_libstr



#define MEMSIZE_OCP_QP_IPM_ARG d_memsize_ocp_qp_ipm_arg
#define CREATE_OCP_QP_IPM_ARG d_create_ocp_qp_ipm_arg
#define SET_DEFAULT_OCP_QP_IPM_ARG d_set_default_ocp_qp_ipm_arg
#define MEMSIZE_OCP_QP_RES d_memsize_ocp_qp_res
#define CREATE_OCP_QP_RES d_create_ocp_qp_res
#define MEMSIZE_OCP_QP_IPM d_memsize_ocp_qp_ipm
#define CREATE_OCP_QP_IPM d_create_ocp_qp_ipm
#define SOLVE_OCP_QP_IPM d_solve_ocp_qp_ipm
#define SOLVE_OCP_QP_IPM2 d_solve_ocp_qp_ipm2
#define CVT_OCP_QP_RES_TO_COLMAJ d_cvt_ocp_qp_res_to_colmaj
#define CVT_OCP_QP_RES_TO_ROWMAJ d_cvt_ocp_qp_res_to_rowmaj



#include "x_ocp_qp_ipm.c"
