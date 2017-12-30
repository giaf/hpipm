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

#include <stdlib.h>
#include <stdio.h>

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_d_aux.h>

#include "../include/hpipm_d_tree_ocp_qp.h"
#include "../include/hpipm_d_tree_ocp_qp_sol.h"
#include "../include/hpipm_d_tree_ocp_qp_ipm.h"
#include "../include/hpipm_d_tree_ocp_qp_kkt.h"
#include "../include/hpipm_d_core_qp_ipm.h"
#include "../include/hpipm_d_core_qp_ipm_aux.h"

#define BACKUP_RES_M d_backup_res_m
#define COMPUTE_ALPHA_QP d_compute_alpha_qp
#define COMPUTE_CENTERING_CORRECTION_QP d_compute_centering_correction_qp
#define COMPUTE_CENTERING_QP d_compute_centering_qp
#define COMPUTE_MU_AFF_QP d_compute_mu_aff_qp
#define COMPUTE_RES_TREE_OCP_QP d_compute_res_tree_ocp_qp
#define CORE_QP_IPM_WORKSPACE d_core_qp_ipm_workspace
#define CREATE_STRMAT d_create_strmat
#define CREATE_STRVEC d_create_strvec
#define CREATE_CORE_QP_IPM d_create_core_qp_ipm
#define FACT_SOLVE_KKT_STEP_TREE_OCP_QP d_fact_solve_kkt_step_tree_ocp_qp
#define FACT_SOLVE_KKT_UNCONSTR_TREE_OCP_QP d_fact_solve_kkt_unconstr_tree_ocp_qp
#define INIT_VAR_TREE_OCP_QP d_init_var_tree_ocp_qp
#define MEMSIZE_CORE_QP_IPM d_memsize_core_qp_ipm
#define REAL double
#define SIZE_STRMAT d_size_strmat
#define SIZE_STRVEC d_size_strvec
#define SOLVE_KKT_STEP_TREE_OCP_QP d_solve_kkt_step_tree_ocp_qp
#define STRMAT blasfeo_dmat
#define STRVEC blasfeo_dvec
#define TREE_OCP_QP d_tree_ocp_qp
#define TREE_OCP_QP_IPM_ARG d_tree_ocp_qp_ipm_arg
#define TREE_OCP_QP_IPM_WORKSPACE d_tree_ocp_qp_ipm_workspace
#define TREE_OCP_QP_SOL d_tree_ocp_qp_sol
#define UPDATE_VAR_QP d_update_var_qp
#define VECNRM_INF_LIBSTR dvecnrm_inf_libstr

#define MEMSIZE_TREE_OCP_QP_IPM_ARG d_memsize_tree_ocp_qp_ipm_arg
#define CREATE_TREE_OCP_QP_IPM_ARG d_create_tree_ocp_qp_ipm_arg
#define SET_DEFAULT_TREE_OCP_QP_IPM_ARG d_set_default_tree_ocp_qp_ipm_arg
#define MEMSIZE_TREE_OCP_QP_IPM d_memsize_tree_ocp_qp_ipm
#define CREATE_TREE_OCP_QP_IPM d_create_tree_ocp_qp_ipm
#define SOLVE_TREE_OCP_QP_IPM d_solve_tree_ocp_qp_ipm
#define SOLVE_TREE_OCP_QP_IPM2 d_solve_tree_ocp_qp_ipm2



#include "x_tree_ocp_qp_ipm.c"
