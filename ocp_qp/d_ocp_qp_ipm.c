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



#define AXPY_LIBSTR blasfeo_daxpy
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
#define CREATE_OCP_QP_SOL d_create_ocp_qp_sol
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
#define MEMSIZE_OCP_QP_SOL d_memsize_ocp_qp_sol
#define OCP_QP d_ocp_qp
#define OCP_QP_IPM_ARG d_ocp_qp_ipm_arg
#define OCP_QP_IPM_MODE ocp_qp_ipm_mode
#define OCP_QP_IPM_WORKSPACE d_ocp_qp_ipm_workspace
#define OCP_QP_RES d_ocp_qp_res
#define OCP_QP_RES_WORKSPACE d_ocp_qp_res_workspace
#define OCP_QP_DIM d_ocp_qp_dim
#define OCP_QP_SOL d_ocp_qp_sol
#define PRINT_E_MAT d_print_e_mat
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
#define VECNRM_INF_LIBSTR blasfeo_dvecnrm_inf



#define MEMSIZE_OCP_QP_IPM_ARG d_memsize_ocp_qp_ipm_arg
#define CREATE_OCP_QP_IPM_ARG d_create_ocp_qp_ipm_arg
#define SET_DEFAULT_OCP_QP_IPM_ARG d_set_default_ocp_qp_ipm_arg
#define MEMSIZE_OCP_QP_IPM d_memsize_ocp_qp_ipm
#define CREATE_OCP_QP_IPM d_create_ocp_qp_ipm
#define SOLVE_OCP_QP_IPM d_solve_ocp_qp_ipm
#define SOLVE_OCP_QP_IPM2 d_solve_ocp_qp_ipm2



#include "x_ocp_qp_ipm.c"
