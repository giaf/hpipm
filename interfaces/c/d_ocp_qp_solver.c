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

#include <hpipm_d_ocp_qp_ipm.h>
#include <hpipm_d_ocp_qp_solver.h>
#include <hpipm_aux_string.h>
#include <hpipm_aux_mem.h>



#define OCP_QP d_ocp_qp
#define OCP_QP_DIM d_ocp_qp_dim
#define OCP_QP_IPM_ARG d_ocp_qp_ipm_arg
#define OCP_QP_IPM_ARG_CREATE d_ocp_qp_ipm_arg_create
#define OCP_QP_IPM_ARG_DEEPCOPY d_ocp_qp_ipm_arg_deepcopy
#define OCP_QP_IPM_ARG_MEMSIZE d_ocp_qp_ipm_arg_memsize
#define OCP_QP_IPM_ARG_SET d_ocp_qp_ipm_arg_set
#define OCP_QP_IPM_ARG_SET_ALPHA_MIN d_ocp_qp_ipm_arg_set_alpha_min
#define OCP_QP_IPM_ARG_SET_DEFAULT d_ocp_qp_ipm_arg_set_default
#define OCP_QP_IPM_ARG_SET_ITER_MAX d_ocp_qp_ipm_arg_set_iter_max
#define OCP_QP_IPM_ARG_SET_MU0 d_ocp_qp_ipm_arg_set_mu0
#define OCP_QP_IPM_ARG_SET_TOL_COMP d_ocp_qp_ipm_arg_set_tol_comp
#define OCP_QP_IPM_ARG_SET_TOL_EQ d_ocp_qp_ipm_arg_set_tol_eq
#define OCP_QP_IPM_ARG_SET_TOL_INEQ d_ocp_qp_ipm_arg_set_tol_ineq
#define OCP_QP_IPM_ARG_SET_TOL_STAT d_ocp_qp_ipm_arg_set_tol_stat
#define OCP_QP_IPM_GET d_ocp_qp_ipm_get
#define OCP_QP_IPM_GET_STATUS d_ocp_qp_ipm_get_status
#define OCP_QP_IPM_SOLVE d_ocp_qp_ipm_solve
#define OCP_QP_IPM_WS d_ocp_qp_ipm_ws
#define OCP_QP_IPM_WS_CREATE d_ocp_qp_ipm_ws_create
#define OCP_QP_IPM_WS_MEMSIZE d_ocp_qp_ipm_ws_memsize
#define OCP_QP_SOL d_ocp_qp_sol
#define OCP_QP_SOLVER_ARG d_ocp_qp_solver_arg
#define OCP_QP_SOLVER_WS d_ocp_qp_solver_ws
#define REAL double



// arg
#define OCP_QP_SOLVER_ARG_STRSIZE d_ocp_qp_solver_arg_strsize
#define OCP_QP_SOLVER_ARG_MEMSIZE d_ocp_qp_solver_arg_memsize
#define OCP_QP_SOLVER_ARG_CREATE d_ocp_qp_solver_arg_create
#define OCP_QP_SOLVER_ARG_SET_DEFAULT d_ocp_qp_solver_arg_set_default
#define OCP_QP_SOLVER_ARG_SET d_ocp_qp_solver_arg_set
#define OCP_QP_SOLVER_ARG_DEEPCOPY d_ocp_qp_solver_arg_deepcopy
// ws
#define OCP_QP_SOLVER_WS_STRSIZE d_ocp_qp_solver_ws_strsize
#define OCP_QP_SOLVER_WS_MEMSIZE d_ocp_qp_solver_ws_memsize
#define OCP_QP_SOLVER_WS_CREATE d_ocp_qp_solver_ws_create
#define OCP_QP_SOLVER_GET d_ocp_qp_solver_get
#define OCP_QP_SOLVER_GET_STATUS d_ocp_qp_solver_get_status
#define OCP_QP_SOLVER_SET d_ocp_qp_solver_set
#define OCP_QP_SOLVER_SET_ITER_MAX d_ocp_qp_solver_set_iter_max
#define OCP_QP_SOLVER_SET_ALPHA_MIN d_ocp_qp_solver_set_alpha_min
#define OCP_QP_SOLVER_SET_MU0 d_ocp_qp_solver_set_mu0
#define OCP_QP_SOLVER_SET_TOL_STAT d_ocp_qp_solver_set_tol_stat
#define OCP_QP_SOLVER_SET_TOL_EQ d_ocp_qp_solver_set_tol_eq
#define OCP_QP_SOLVER_SET_TOL_INEQ d_ocp_qp_solver_set_tol_ineq
#define OCP_QP_SOLVER_SET_TOL_COMP d_ocp_qp_solver_set_tol_comp
#define OCP_QP_SOLVER_SOLVE d_ocp_qp_solver_solve



#include "x_ocp_qp_solver.c"
