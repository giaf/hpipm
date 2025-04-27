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
#include <blasfeo_d_blas.h>

#include <hpipm_d_ocp_qp_dim.h>
#include <hpipm_d_ocp_qp_res.h>
#include <hpipm_aux_string.h>
#include <hpipm_aux_mem.h>



#define BLASFEO_VECEL BLASFEO_DVECEL



#define AXPY blasfeo_daxpy
#define CREATE_STRVEC blasfeo_create_dvec
#define DOT blasfeo_ddot
#define GEMV_DIAG blasfeo_dgemv_d
#define GEMV_NT blasfeo_dgemv_nt
#define OCP_QP d_ocp_qp
#define OCP_QP_DIM d_ocp_qp_dim
#define OCP_QP_RES d_ocp_qp_res
#define OCP_QP_RES_WS d_ocp_qp_res_ws
#define OCP_QP_SOL d_ocp_qp_sol
#define PACK_VEC blasfeo_pack_dvec
#define REAL double
#define SIZE_STRVEC blasfeo_memsize_dvec
#define STRMAT blasfeo_dmat
#define STRVEC blasfeo_dvec
#define SYMV_L blasfeo_dsymv_l
#define UNPACK_VEC blasfeo_unpack_dvec
#define VECAD_SP blasfeo_dvecad_sp
#define VECCP blasfeo_dveccp
#define VECEX_SP blasfeo_dvecex_sp
#define VECMULACC blasfeo_dvecmulacc
#define VECMULDOT blasfeo_dvecmuldot
#define VECNRM_INF blasfeo_dvecnrm_inf
#define VECSE blasfeo_dvecse



#define OCP_QP_RES_MEMSIZE d_ocp_qp_res_memsize
#define OCP_QP_RES_CREATE d_ocp_qp_res_create
#define OCP_QP_RES_WS_MEMSIZE d_ocp_qp_res_ws_memsize
#define OCP_QP_RES_WS_CREATE d_ocp_qp_res_ws_create
#define OCP_QP_RES_COMPUTE d_ocp_qp_res_compute
#define OCP_QP_RES_COMPUTE_LIN d_ocp_qp_res_compute_lin
#define OCP_QP_RES_COMPUTE_INF_NORM d_ocp_qp_res_compute_inf_norm
#define OCP_QP_RES_GET_ALL d_ocp_qp_res_get_all
#define OCP_QP_RES_GET_MAX_RES_STAT d_ocp_qp_res_get_max_res_stat
#define OCP_QP_RES_GET_MAX_RES_EQ d_ocp_qp_res_get_max_res_eq
#define OCP_QP_RES_GET_MAX_RES_INEQ d_ocp_qp_res_get_max_res_ineq
#define OCP_QP_RES_GET_MAX_RES_COMP d_ocp_qp_res_get_max_res_comp
#define OCP_QP_RES_SET_ZERO d_ocp_qp_res_set_zero
#define OCP_QP_RES_SET d_ocp_qp_res_set
#define OCP_QP_RES_SET_RES_R d_ocp_qp_res_set_res_r
#define OCP_QP_RES_SET_RES_Q d_ocp_qp_res_set_res_q
#define OCP_QP_RES_SET_RES_ZL d_ocp_qp_res_set_res_zl
#define OCP_QP_RES_SET_RES_ZU d_ocp_qp_res_set_res_zu
#define OCP_QP_RES_SET_RES_B d_ocp_qp_res_set_res_b
#define OCP_QP_RES_SET_RES_LB d_ocp_qp_res_set_res_lb
#define OCP_QP_RES_SET_RES_LBU d_ocp_qp_res_set_res_lbu
#define OCP_QP_RES_SET_RES_LBX d_ocp_qp_res_set_res_lbx
#define OCP_QP_RES_SET_RES_UB d_ocp_qp_res_set_res_ub
#define OCP_QP_RES_SET_RES_UBU d_ocp_qp_res_set_res_ubu
#define OCP_QP_RES_SET_RES_UBX d_ocp_qp_res_set_res_ubx
#define OCP_QP_RES_SET_RES_LG d_ocp_qp_res_set_res_lg
#define OCP_QP_RES_SET_RES_UG d_ocp_qp_res_set_res_ug
#define OCP_QP_RES_SET_RES_LS d_ocp_qp_res_set_res_ls
#define OCP_QP_RES_SET_RES_US d_ocp_qp_res_set_res_us



#include "x_ocp_qp_res.c"
