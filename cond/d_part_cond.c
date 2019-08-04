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
#include "../include/hpipm_d_part_cond.h"
#include "../include/hpipm_d_cond_aux.h"



#define COND_B d_cond_b
#define COND_BABT d_cond_BAbt
#define COND_D d_cond_d
#define COND_DCTD d_cond_DCtd
#define COND_RQ_N2NX3 d_cond_rq_N2nx3
#define COND_RSQRQ_N2NX3 d_cond_RSQrq_N2nx3
#define COND_QP_OCP2DENSE_ARG d_cond_qp_ocp2dense_arg
#define COND_QP_OCP2DENSE_WORKSPACE d_cond_qp_ocp2dense_workspace
#define COND_QP_OCP2OCP_ARG d_cond_qp_ocp2ocp_arg
#define COND_QP_OCP2OCP_WORKSPACE d_cond_qp_ocp2ocp_workspace
#define CREATE_COND_QP_OCP2DENSE_ARG d_create_cond_qp_ocp2dense_arg
#define CREATE_COND_QP_OCP2DENSE d_create_cond_qp_ocp2dense
#define CREATE_STRVEC blasfeo_create_dvec
#define DENSE_QP d_dense_qp
#define DENSE_QP_SOL d_dense_qp_sol
#define EXPAND_SOL d_expand_sol
#define GECP_LIBSTR blasfeo_dgecp
#define MEMSIZE_COND_QP_OCP2DENSE_ARG d_memsize_cond_qp_ocp2dense_arg
#define MEMSIZE_COND_QP_OCP2DENSE d_memsize_cond_qp_ocp2dense
#define OCP_QP d_ocp_qp
#define OCP_QP_DIM d_ocp_qp_dim
#define OCP_QP_SOL d_ocp_qp_sol
#define SET_COND_QP_OCP2DENSE_ARG_RIC_ALG d_set_cond_qp_ocp2dense_arg_ric_alg
#define SET_DEFAULT_COND_QP_OCP2DENSE_ARG d_set_default_cond_qp_ocp2dense_arg
#define STRVEC blasfeo_dvec
#define UPDATE_COND_BABT d_update_cond_BAbt
#define UPDATE_COND_DCTD d_update_cond_DCtd
#define UPDATE_COND_RSQRQ_N2NX3 d_update_cond_RSQrq_N2nx3
#define VECCP_LIBSTR blasfeo_dveccp

#define COMPUTE_BLOCK_SIZE_COND_QP_OCP2OCP d_compute_block_size_cond_qp_ocp2ocp
#define COMPUTE_QP_DIM_OCP2OCP d_compute_qp_dim_ocp2ocp
#define MEMSIZE_COND_QP_OCP2OCP_ARG d_memsize_cond_qp_ocp2ocp_arg
#define CREATE_COND_QP_OCP2OCP_ARG d_create_cond_qp_ocp2ocp_arg
#define SET_DEFAULT_COND_QP_OCP2OCP_ARG d_set_default_cond_qp_ocp2ocp_arg
#define SET_COND_QP_OCP2OCP_ARG_RIC_ALG d_set_cond_qp_ocp2ocp_arg_ric_alg
#define MEMSIZE_COND_QP_OCP2OCP d_memsize_cond_qp_ocp2ocp
#define CREATE_COND_QP_OCP2OCP d_create_cond_qp_ocp2ocp
#define COND_QP_OCP2OCP d_cond_qp_ocp2ocp
#define COND_RHS_QP_OCP2OCP d_cond_rhs_qp_ocp2ocp
#define EXPAND_SOL_OCP2OCP d_expand_sol_ocp2ocp
#define UPDATE_COND_QP_OCP2OCP d_update_cond_qp_ocp2ocp



#include "x_part_cond.c"
