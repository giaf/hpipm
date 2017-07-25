/**************************************************************************************************
*                                                                                                 *
* This file is part of HPIPM.                                                                     *
*                                                                                                 *
* HPIPM -- High Performance Interior Point Method.                                                *
* Copyright (C) 2017 by Gianluca Frison.                                                          *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
*                                                                                                 *
* HPIPM is free software; you can redistribute it and/or                                          *
* modify it under the terms of the GNU Lesser General Public                                      *
* License as published by the Free Software Foundation; either                                    *
* version 2.1 of the License, or (at your option) any later version.                              *
*                                                                                                 *
* HPIPM is distributed in the hope that it will be useful,                                        *
* but WITHOUT ANY WARRANTY; without even the implied warranty of                                  *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                                            *
* See the GNU Lesser General Public License for more details.                                     *
*                                                                                                 *
* You should have received a copy of the GNU Lesser General Public                                *
* License along with HPIPM; if not, write to the Free Software                                    *
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA                  *
*                                                                                                 *
* Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             *
*                                                                                                 *
**************************************************************************************************/

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_d_aux.h>

#include "../include/hpipm_d_tree_ocp_qp.h"
#include "../include/hpipm_d_tree_ocp_qp_sol.h"
#include "../include/hpipm_d_tree_ocp_qp_ipm_hard.h"
#include "../include/hpipm_d_tree_ocp_qp_kkt.h"
#include "../include/hpipm_d_core_qp_ipm_hard.h"
#include "../include/hpipm_d_core_qp_ipm_hard_aux.h"

#define COMPUTE_ALPHA_HARD_QP d_compute_alpha_hard_qp
#define COMPUTE_CENTERING_CORRECTION_HARD_QP d_compute_centering_correction_hard_qp
#define COMPUTE_MU_AFF_HARD_QP d_compute_mu_aff_hard_qp
#define COMPUTE_RES_HARD_TREE_OCP_QP d_compute_res_hard_tree_ocp_qp
#define CREATE_STRMAT d_create_strmat
#define CREATE_STRVEC d_create_strvec
#define CREATE_IPM_HARD_CORE_QP d_create_ipm_hard_core_qp
#define FACT_SOLVE_KKT_STEP_HARD_TREE_OCP_QP d_fact_solve_kkt_step_hard_tree_ocp_qp
#define FACT_SOLVE_KKT_UNCONSTR_TREE_OCP_QP d_fact_solve_kkt_unconstr_tree_ocp_qp
#define INIT_VAR_HARD_TREE_OCP_QP d_init_var_hard_tree_ocp_qp
#define IPM_HARD_CORE_QP_WORKSPACE d_ipm_hard_core_qp_workspace
#define IPM_HARD_TREE_OCP_QP_ARG d_ipm_hard_tree_ocp_qp_arg
#define IPM_HARD_TREE_OCP_QP_WORKSPACE d_ipm_hard_tree_ocp_qp_workspace
#define MEMSIZE_IPM_HARD_CORE_QP d_memsize_ipm_hard_core_qp
#define REAL double
#define SIZE_STRMAT d_size_strmat
#define SIZE_STRVEC d_size_strvec
#define SOLVE_KKT_STEP_HARD_TREE_OCP_QP d_solve_kkt_step_hard_tree_ocp_qp
#define STRMAT d_strmat
#define STRVEC d_strvec
#define TREE_OCP_QP d_tree_ocp_qp
#define TREE_OCP_QP_SOL d_tree_ocp_qp_sol
#define UPDATE_VAR_HARD_QP d_update_var_hard_qp

#define MEMSIZE_IPM_HARD_TREE_OCP_QP d_memsize_ipm_hard_tree_ocp_qp
#define CREATE_IPM_HARD_TREE_OCP_QP d_create_ipm_hard_tree_ocp_qp
#define SOLVE_IPM_HARD_TREE_OCP_QP d_solve_ipm_hard_tree_ocp_qp
#define SOLVE_IPM2_HARD_TREE_OCP_QP d_solve_ipm2_hard_tree_ocp_qp



#include "x_tree_ocp_qp_ipm_hard.c"
