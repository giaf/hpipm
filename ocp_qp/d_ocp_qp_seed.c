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
#include <blasfeo_d_aux.h>

#include <hpipm_d_ocp_qp_dim.h>
#include <hpipm_d_ocp_qp_seed.h>
#include <hpipm_aux_string.h>
#include <hpipm_aux_mem.h>



#define BLASFEO_VECEL BLASFEO_DVECEL



#define CREATE_STRVEC blasfeo_create_dvec
#define OCP_QP_DIM d_ocp_qp_dim
#define OCP_QP_SEED d_ocp_qp_seed
#define PACK_VEC blasfeo_pack_dvec
#define REAL double
#define SIZE_STRVEC blasfeo_memsize_dvec
#define STRVEC blasfeo_dvec
#define VECSE blasfeo_dvecse
#define VECSC blasfeo_dvecsc



#define OCP_QP_SEED_STRSIZE d_ocp_qp_seed_strsize
#define OCP_QP_SEED_MEMSIZE d_ocp_qp_seed_memsize
#define OCP_QP_SEED_CREATE d_ocp_qp_seed_create
#define OCP_QP_SEED_WS_MEMSIZE d_ocp_qp_seed_ws_memsize
#define OCP_QP_SEED_WS_CREATE d_ocp_qp_seed_ws_create
#define OCP_QP_SEED_COMPUTE d_ocp_qp_seed_compute
#define OCP_QP_SEED_COMPUTE_LIN d_ocp_qp_seed_compute_lin
#define OCP_QP_SEED_COMPUTE_INF_NORM d_ocp_qp_seed_compute_inf_norm
#define OCP_QP_SEED_GET_ALL d_ocp_qp_seed_get_all
#define OCP_QP_SEED_GET_MAX_SEED_STAT d_ocp_qp_seed_get_max_seed_stat
#define OCP_QP_SEED_GET_MAX_SEED_EQ d_ocp_qp_seed_get_max_seed_eq
#define OCP_QP_SEED_GET_MAX_SEED_INEQ d_ocp_qp_seed_get_max_seed_ineq
#define OCP_QP_SEED_GET_MAX_SEED_COMP d_ocp_qp_seed_get_max_seed_comp
#define OCP_QP_SEED_SET_ZERO d_ocp_qp_seed_set_zero
#define OCP_QP_SEED_SET d_ocp_qp_seed_set
#define OCP_QP_SEED_SET_SEED_R d_ocp_qp_seed_set_seed_r
#define OCP_QP_SEED_SET_SEED_Q d_ocp_qp_seed_set_seed_q
#define OCP_QP_SEED_SET_SEED_ZL d_ocp_qp_seed_set_seed_zl
#define OCP_QP_SEED_SET_SEED_ZU d_ocp_qp_seed_set_seed_zu
#define OCP_QP_SEED_SET_SEED_B d_ocp_qp_seed_set_seed_b
#define OCP_QP_SEED_SET_SEED_LB d_ocp_qp_seed_set_seed_lb
#define OCP_QP_SEED_SET_SEED_LBU d_ocp_qp_seed_set_seed_lbu
#define OCP_QP_SEED_SET_SEED_LBX d_ocp_qp_seed_set_seed_lbx
#define OCP_QP_SEED_SET_SEED_UB d_ocp_qp_seed_set_seed_ub
#define OCP_QP_SEED_SET_SEED_UBU d_ocp_qp_seed_set_seed_ubu
#define OCP_QP_SEED_SET_SEED_UBX d_ocp_qp_seed_set_seed_ubx
#define OCP_QP_SEED_SET_SEED_LG d_ocp_qp_seed_set_seed_lg
#define OCP_QP_SEED_SET_SEED_UG d_ocp_qp_seed_set_seed_ug
#define OCP_QP_SEED_SET_SEED_LS d_ocp_qp_seed_set_seed_ls
#define OCP_QP_SEED_SET_SEED_US d_ocp_qp_seed_set_seed_us



#include "x_ocp_qp_seed.c"

