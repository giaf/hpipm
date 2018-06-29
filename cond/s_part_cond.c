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
#include <blasfeo_s_blas.h>
#include <blasfeo_s_aux.h>

#include "../include/hpipm_s_ocp_qp_dim.h"
#include "../include/hpipm_s_ocp_qp.h"
#include "../include/hpipm_s_ocp_qp_sol.h"
#include "../include/hpipm_s_dense_qp.h"
#include "../include/hpipm_s_dense_qp_sol.h"
#include "../include/hpipm_s_cond.h"
#include "../include/hpipm_s_part_cond.h"
#include "../include/hpipm_s_cond_aux.h"



#define COND_B s_cond_b
#define COND_BABT s_cond_BAbt
#define COND_D s_cond_d
#define COND_DCTD s_cond_DCtd
#define COND_RQ_N2NX3 s_cond_rq_N2nx3
#define COND_RSQRQ_N2NX3 s_cond_RSQrq_N2nx3
#define COND_QP_OCP2DENSE_ARG s_cond_qp_ocp2dense_arg
#define COND_QP_OCP2DENSE_WORKSPACE s_cond_qp_ocp2dense_workspace
#define COND_QP_OCP2OCP_ARG s_cond_qp_ocp2ocp_arg
#define COND_QP_OCP2OCP_WORKSPACE s_cond_qp_ocp2ocp_workspace
#define CREATE_COND_QP_OCP2DENSE_ARG s_create_cond_qp_ocp2dense_arg
#define CREATE_COND_QP_OCP2DENSE s_create_cond_qp_ocp2dense
#define CREATE_STRVEC blasfeo_create_svec
#define DENSE_QP s_dense_qp
#define DENSE_QP_SOL s_dense_qp_sol
#define EXPAND_SOL s_expand_sol
#define GECP_LIBSTR blasfeo_sgecp
#define MEMSIZE_COND_QP_OCP2DENSE_ARG s_memsize_cond_qp_ocp2dense_arg
#define MEMSIZE_COND_QP_OCP2DENSE s_memsize_cond_qp_ocp2dense
#define OCP_QP s_ocp_qp
#define OCP_QP_DIM s_ocp_qp_dim
#define OCP_QP_SOL s_ocp_qp_sol
#define SET_DEFAULT_COND_QP_OCP2DENSE_ARG s_set_default_cond_qp_ocp2dense_arg
#define STRVEC blasfeo_svec
#define UPDATE_COND_BABT s_update_cond_BAbt
#define UPDATE_COND_DCTD s_update_cond_DCtd
#define UPDATE_COND_RSQRQ_N2NX3 s_update_cond_RSQrq_N2nx3
#define VECCP_LIBSTR blasfeo_sveccp

#define COMPUTE_BLOCK_SIZE_COND_QP_OCP2OCP s_compute_block_size_cond_qp_ocp2ocp
#define COMPUTE_QP_DIM_OCP2OCP s_compute_qp_dim_ocp2ocp
#define MEMSIZE_COND_QP_OCP2OCP_ARG s_memsize_cond_qp_ocp2ocp_arg
#define CREATE_COND_QP_OCP2OCP_ARG s_create_cond_qp_ocp2ocp_arg
#define SET_DEFAULT_COND_QP_OCP2OCP_ARG s_set_default_cond_qp_ocp2ocp_arg
#define MEMSIZE_COND_QP_OCP2OCP s_memsize_cond_qp_ocp2ocp
#define CREATE_COND_QP_OCP2OCP s_create_cond_qp_ocp2ocp
#define COND_QP_OCP2OCP s_cond_qp_ocp2ocp
#define COND_RHS_QP_OCP2OCP s_cond_rhs_qp_ocp2ocp
#define EXPAND_SOL_OCP2OCP s_expand_sol_ocp2ocp
#define UPDATE_COND_QP_OCP2OCP s_update_cond_qp_ocp2ocp



#include "x_part_cond.c"

