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

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_d_blas.h>
#include <blasfeo_d_aux.h>

#include "../include/hpipm_d_ocp_qp_dim.h"
#include "../include/hpipm_d_ocp_qp.h"
#include "../include/hpipm_d_ocp_qp_sol.h"
#include "../include/hpipm_d_dense_qp.h"
#include "../include/hpipm_d_dense_qp_sol.h"
#include "../include/hpipm_d_cond.h"
#include "../include/hpipm_d_cond_aux.h"



#define COND_DCTD d_cond_DCtd
#define COND_D d_cond_d
#define COND_B d_cond_b
#define COND_BABT d_cond_BAbt
#define COND_RQ_N2NX3 d_cond_rq_N2nx3
#define COND_RSQRQ_N2NX3 d_cond_RSQrq_N2nx3
#define COND_QP_OCP2DENSE_ARG d_cond_qp_ocp2dense_arg
#define COND_QP_OCP2DENSE_WORKSPACE d_cond_qp_ocp2dense_workspace
#define CREATE_STRMAT blasfeo_create_dmat
#define CREATE_STRVEC blasfeo_create_dvec
#define DENSE_QP d_dense_qp
#define DENSE_QP_DIM d_dense_qp_dim
#define DENSE_QP_SOL d_dense_qp_sol
#define EXPAND_SOL d_expand_sol
#define EXPAND_PRIMAL_SOL d_expand_primal_sol
#define OCP_QP d_ocp_qp
#define OCP_QP_DIM d_ocp_qp_dim
#define OCP_QP_SOL d_ocp_qp_sol
#define SIZE_STRMAT blasfeo_memsize_dmat
#define SIZE_STRVEC blasfeo_memsize_dvec
#define STRMAT blasfeo_dmat
#define STRVEC blasfeo_dvec
#define UPDATE_COND_DCTD d_update_cond_DCtd
#define UPDATE_COND_BABT d_update_cond_BAbt
#define UPDATE_COND_RSQRQ_N2NX3 d_update_cond_RSQrq_N2nx3

#define COMPUTE_QP_DIM_OCP2DENSE d_compute_qp_dim_ocp2dense
#define MEMSIZE_COND_QP_OCP2DENSE_ARG d_memsize_cond_qp_ocp2dense_arg
#define CREATE_COND_QP_OCP2DENSE_ARG d_create_cond_qp_ocp2dense_arg
#define SET_DEFAULT_COND_QP_OCP2DENSE_ARG d_set_default_cond_qp_ocp2dense_arg
#define MEMSIZE_COND_QP_OCP2DENSE d_memsize_cond_qp_ocp2dense
#define CREATE_COND_QP_OCP2DENSE d_create_cond_qp_ocp2dense
#define COND_QP_OCP2DENSE d_cond_qp_ocp2dense
#define COND_RHS_QP_OCP2DENSE d_cond_rhs_qp_ocp2dense
#define EXPAND_SOL_DENSE2OCP d_expand_sol_dense2ocp
#define EXPAND_PRIMAL_SOL_DENSE2OCP d_expand_primal_sol_dense2ocp
#define UPDATE_COND_QP_OCP2DENSE d_update_cond_qp_ocp2dense



#include "x_cond.c"
