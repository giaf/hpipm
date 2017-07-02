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



#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_s_aux.h>

#include "../include/hpipm_s_dense_qp.h"
#include "../include/hpipm_s_dense_qp_sol.h"
#include "../include/hpipm_s_dense_qp_ipm_hard.h"
#include "../include/hpipm_s_core_qp_ipm_hard.h"
#include "../include/hpipm_s_core_qp_ipm_hard_aux.h"
#include "../include/hpipm_s_dense_qp_kkt.h"



#define COMPUTE_ALPHA_HARD_QP s_compute_alpha_hard_qp
#define COMPUTE_CENTERING_CORRECTION_HARD_QP s_compute_centering_correction_hard_qp
#define COMPUTE_MU_AFF_HARD_QP s_compute_mu_aff_hard_qp
#define COMPUTE_RES_HARD_DENSE_QP s_compute_res_hard_dense_qp
#define CREATE_IPM_HARD_CORE_QP s_create_ipm_hard_core_qp
#define CREATE_STRMAT s_create_strmat
#define CREATE_STRVEC s_create_strvec
#define DENSE_QP s_dense_qp
#define DENSE_QP_SOL s_dense_qp_sol
#define FACT_SOLVE_KKT_STEP_HARD_DENSE_QP s_fact_solve_kkt_step_hard_dense_qp
#define FACT_SOLVE_KKT_UNCONSTR_DENSE_QP s_fact_solve_kkt_unconstr_dense_qp
#define INIT_VAR_HARD_DENSE_QP s_init_var_hard_dense_qp
#define IPM_HARD_CORE_QP_WORKSPACE s_ipm_hard_core_qp_workspace
#define IPM_HARD_DENSE_QP_ARG s_ipm_hard_dense_qp_arg
#define IPM_HARD_DENSE_QP_WORKSPACE s_ipm_hard_dense_qp_workspace
#define MEMSIZE_IPM_HARD_CORE_QP s_memsize_ipm_hard_core_qp
#define REAL float
#define SIZE_STRMAT s_size_strmat
#define SIZE_STRVEC s_size_strvec
#define SOLVE_KKT_STEP_HARD_DENSE_QP s_solve_kkt_step_hard_dense_qp
#define STRMAT s_strmat
#define STRVEC s_strvec
#define UPDATE_VAR_HARD_QP s_update_var_hard_qp



#define MEMSIZE_IPM_HARD_DENSE_QP s_memsize_ipm_hard_dense_qp
#define CREATE_IPM_HARD_DENSE_QP s_create_ipm_hard_dense_qp
#define SOLVE_IPM_HARD_DENSE_QP s_solve_ipm_hard_dense_qp
#define SOLVE_IPM2_HARD_DENSE_QP s_solve_ipm2_hard_dense_qp



#include "x_dense_qp_ipm_hard.c"

