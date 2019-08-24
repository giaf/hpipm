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



#include <stdio.h>
#include <stdlib.h>

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_d_aux.h>

#include <hpipm_d_ocp_qp_dim.h>
#include <hpipm_d_ocp_qp.h>
#include <hpipm_aux_string.h>



#define DOUBLE_PRECISION



#define CREATE_STRMAT blasfeo_create_dmat
#define CREATE_STRVEC blasfeo_create_dvec
#define CVT_MAT2STRMAT blasfeo_pack_dmat
#define CVT_STRMAT2MAT blasfeo_unpack_dmat
#define CVT_TRAN_MAT2STRMAT blasfeo_pack_tran_dmat
#define CVT_TRAN_STRMAT2MAT blasfeo_unpack_tran_dmat
#define CVT_VEC2STRVEC blasfeo_pack_dvec
#define CVT_STRVEC2VEC blasfeo_unpack_dvec
#define GECP blasfeo_dgecp
#define GESE blasfeo_dgese
#define OCP_QP d_ocp_qp
#define OCP_QP_DIM d_ocp_qp_dim
#define REAL double
#define STRMAT blasfeo_dmat
#define STRVEC blasfeo_dvec
#define SIZE_STRMAT blasfeo_memsize_dmat
#define SIZE_STRVEC blasfeo_memsize_dvec
#define VECCP blasfeo_dveccp
#define VECSC blasfeo_dvecsc
#define VECSE blasfeo_dvecse

#define OCP_QP_STRSIZE d_ocp_qp_strsize
#define OCP_QP_MEMSIZE d_ocp_qp_memsize
#define OCP_QP_CREATE d_ocp_qp_create
#define OCP_QP_COPY_ALL d_ocp_qp_copy_all
#define OCP_QP_SET_ALL_ZERO d_ocp_qp_set_all_zero
#define OCP_QP_SET_RHS_ZERO d_ocp_qp_set_rhs_zero
#define OCP_QP_SET_ALL d_ocp_qp_set_all
#define OCP_QP_SET d_ocp_qp_set
#define OCP_QP_SET_EL d_ocp_qp_set_el
#define OCP_QP_SET_A d_ocp_qp_set_A
#define OCP_QP_SET_B d_ocp_qp_set_B
#define OCP_QP_GET_B d_cvt_ocp_qp_B
#define OCP_QP_SET_BVEC d_ocp_qp_set_b
#define OCP_QP_GET_BVEC d_cvt_ocp_qp_b
#define OCP_QP_SET_Q d_ocp_qp_set_Q
#define OCP_QP_GET_Q d_cvt_ocp_qp_Q
#define OCP_QP_SET_S d_ocp_qp_set_S
#define OCP_QP_GET_S d_cvt_ocp_qp_S
#define OCP_QP_SET_R d_ocp_qp_set_R
#define OCP_QP_GET_R d_cvt_ocp_qp_R
#define OCP_QP_SET_QVEC d_ocp_qp_set_q
#define OCP_QP_GET_QVEC d_cvt_ocp_qp_q
#define OCP_QP_SET_RVEC d_ocp_qp_set_r
#define OCP_QP_GET_RVEC d_cvt_ocp_qp_r
#define OCP_QP_SET_LBX d_ocp_qp_set_lbx
#define OCP_QP_SET_EL_LBX d_ocp_qp_set_el_lbx
#define OCP_QP_SET_UBX d_ocp_qp_set_ubx
#define OCP_QP_SET_EL_UBX d_ocp_qp_set_el_ubx
#define OCP_QP_SET_LBU d_ocp_qp_set_lbu
#define OCP_QP_GET_LBU d_cvt_ocp_qp_lbu
#define OCP_QP_SET_UBU d_ocp_qp_set_ubu
#define OCP_QP_GET_UBU d_cvt_ocp_qp_ubu
#define OCP_QP_SET_IDXB d_ocp_qp_set_idxb
#define OCP_QP_GET_IDXB d_cvt_ocp_qp_idxb
#define OCP_QP_SET_IDXBX d_ocp_qp_set_idxbx
#define OCP_QP_GET_IDXBX d_cvt_ocp_qp_idxbx
#define OCP_QP_SET_JBX d_ocp_qp_set_Jbx
#define OCP_QP_GET_JBX d_ocp_qp_get_Jbx
#define OCP_QP_SET_IDXBU d_ocp_qp_set_idxbu
#define OCP_QP_GET_IDXBU d_cvt_ocp_qp_idxbu
#define OCP_QP_SET_JBU d_ocp_qp_set_Jbu
#define OCP_QP_GET_JBU d_ocp_qp_get_Jbu
#define OCP_QP_SET_LB d_ocp_qp_set_lb
#define OCP_QP_GET_LB d_cvt_ocp_qp_lb
#define OCP_QP_SET_UB d_ocp_qp_set_ub
#define OCP_QP_GET_UB d_cvt_ocp_qp_ub
#define OCP_QP_SET_C d_ocp_qp_set_C
#define OCP_QP_GET_C d_cvt_ocp_qp_C
#define OCP_QP_SET_D d_ocp_qp_set_D
#define OCP_QP_GET_D d_cvt_ocp_qp_D
#define OCP_QP_SET_LG d_ocp_qp_set_lg
#define OCP_QP_GET_LG d_cvt_ocp_qp_lg
#define OCP_QP_SET_UG d_ocp_qp_set_ug
#define OCP_QP_GET_UG d_cvt_ocp_qp_ug
#define OCP_QP_SET_ZL d_ocp_qp_set_Zl
#define OCP_QP_GET_ZL d_cvt_ocp_qp_Zl
#define OCP_QP_SET_ZU d_ocp_qp_set_Zu
#define OCP_QP_GET_ZU d_cvt_ocp_qp_Zu
#define OCP_QP_SET_ZLVEC d_ocp_qp_set_zl
#define OCP_QP_GET_ZLVEC d_cvt_ocp_qp_zl
#define OCP_QP_SET_ZUVEC d_ocp_qp_set_zu
#define OCP_QP_GET_ZUVEC d_cvt_ocp_qp_zu
#define OCP_QP_SET_IDXS d_ocp_qp_set_idxs
#define OCP_QP_GET_IDXS d_cvt_ocp_qp_idxs
#define OCP_QP_SET_JSBU d_ocp_qp_set_Jsbu
#define OCP_QP_GET_JSBU d_ocp_qp_get_Jsbu
#define OCP_QP_SET_JSBX d_ocp_qp_set_Jsbx
#define OCP_QP_GET_JSBX d_ocp_qp_get_Jsbx
#define OCP_QP_SET_JSG d_ocp_qp_set_Jsg
#define OCP_QP_GET_JSG d_ocp_qp_get_Jsg
#define OCP_QP_SET_LLS d_ocp_qp_set_lls
#define OCP_QP_GET_LLS d_cvt_ocp_qp_lls
#define OCP_QP_SET_LUS d_ocp_qp_set_lus
#define OCP_QP_GET_LUS d_cvt_ocp_qp_lus
// getters
#define OCP_QP_GET d_ocp_qp_get
#define OCP_QP_GET_A d_ocp_qp_get_A
#define OCP_QP_GET_LBX d_ocp_qp_get_lbx
#define OCP_QP_GET_UBX d_ocp_qp_get_ubx
// TODO remove ???
#define CHANGE_BOUNDS_DIMENSIONS_OCP_QP d_change_bounds_dimensions_ocp_qp



#include "x_ocp_qp.c"
