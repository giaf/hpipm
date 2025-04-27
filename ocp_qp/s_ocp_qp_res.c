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
#include <blasfeo_s_aux.h>
#include <blasfeo_s_blas.h>

#include <hpipm_s_ocp_qp_dim.h>
#include <hpipm_s_ocp_qp_res.h>
#include <hpipm_aux_string.h>
#include <hpipm_aux_mem.h>



#define BLASFEO_VECEL BLASFEO_SVECEL



#define AXPY blasfeo_saxpy
#define CREATE_STRVEC blasfeo_create_svec
#define DOT blasfeo_sdot
#define GEMV_DIAG blasfeo_sgemv_d
#define GEMV_NT blasfeo_sgemv_nt
#define OCP_QP s_ocp_qp
#define OCP_QP_DIM s_ocp_qp_dim
#define OCP_QP_RES s_ocp_qp_res
#define OCP_QP_RES_WS s_ocp_qp_res_ws
#define OCP_QP_SOL s_ocp_qp_sol
#define PACK_VEC blasfeo_pack_svec
#define REAL float
#define SIZE_STRVEC blasfeo_memsize_svec
#define STRMAT blasfeo_smat
#define STRVEC blasfeo_svec
#define SYMV_L blasfeo_ssymv_l
#define UNPACK_VEC blasfeo_unpack_svec
#define VECAD_SP blasfeo_svecad_sp
#define VECCP blasfeo_sveccp
#define VECEX_SP blasfeo_svecex_sp
#define VECMULACC blasfeo_svecmulacc
#define VECMULDOT blasfeo_svecmuldot
#define VECNRM_INF blasfeo_svecnrm_inf
#define VECSE blasfeo_svecse



#define OCP_QP_RES_MEMSIZE s_ocp_qp_res_memsize
#define OCP_QP_RES_CREATE s_ocp_qp_res_create
#define OCP_QP_RES_WS_MEMSIZE s_ocp_qp_res_ws_memsize
#define OCP_QP_RES_WS_CREATE s_ocp_qp_res_ws_create
#define OCP_QP_RES_COMPUTE s_ocp_qp_res_compute
#define OCP_QP_RES_COMPUTE_LIN s_ocp_qp_res_compute_lin
#define OCP_QP_RES_COMPUTE_INF_NORM s_ocp_qp_res_compute_inf_norm
#define OCP_QP_RES_GET_ALL s_ocp_qp_res_get_all
#define OCP_QP_RES_GET_MAX_RES_STAT s_ocp_qp_res_get_max_res_stat
#define OCP_QP_RES_GET_MAX_RES_EQ s_ocp_qp_res_get_max_res_eq
#define OCP_QP_RES_GET_MAX_RES_INEQ s_ocp_qp_res_get_max_res_ineq
#define OCP_QP_RES_GET_MAX_RES_COMP s_ocp_qp_res_get_max_res_comp
#define OCP_QP_RES_SET_ZERO s_ocp_qp_res_set_zero
#define OCP_QP_RES_SET s_ocp_qp_res_set
#define OCP_QP_RES_SET_RES_R s_ocp_qp_res_set_res_r
#define OCP_QP_RES_SET_RES_Q s_ocp_qp_res_set_res_q
#define OCP_QP_RES_SET_RES_ZL s_ocp_qp_res_set_res_zl
#define OCP_QP_RES_SET_RES_ZU s_ocp_qp_res_set_res_zu
#define OCP_QP_RES_SET_RES_B s_ocp_qp_res_set_res_b
#define OCP_QP_RES_SET_RES_LB s_ocp_qp_res_set_res_lb
#define OCP_QP_RES_SET_RES_LBU s_ocp_qp_res_set_res_lbu
#define OCP_QP_RES_SET_RES_LBX s_ocp_qp_res_set_res_lbx
#define OCP_QP_RES_SET_RES_UB s_ocp_qp_res_set_res_ub
#define OCP_QP_RES_SET_RES_UBU s_ocp_qp_res_set_res_ubu
#define OCP_QP_RES_SET_RES_UBX s_ocp_qp_res_set_res_ubx
#define OCP_QP_RES_SET_RES_LG s_ocp_qp_res_set_res_lg
#define OCP_QP_RES_SET_RES_UG s_ocp_qp_res_set_res_ug
#define OCP_QP_RES_SET_RES_LS s_ocp_qp_res_set_res_ls
#define OCP_QP_RES_SET_RES_US s_ocp_qp_res_set_res_us



#include "x_ocp_qp_res.c"

