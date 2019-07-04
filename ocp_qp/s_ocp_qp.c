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



#include <stdio.h>
#include <stdlib.h>

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_s_aux.h>

#include <hpipm_s_ocp_qp_dim.h>
#include <hpipm_s_ocp_qp.h>
#include <hpipm_aux_string.h>



#define CREATE_STRMAT blasfeo_create_smat
#define CREATE_STRVEC blasfeo_create_svec
#define CVT_MAT2STRMAT blasfeo_pack_smat
#define CVT_STRMAT2MAT blasfeo_unpack_smat
#define CVT_TRAN_MAT2STRMAT blasfeo_pack_tran_smat
#define CVT_TRAN_STRMAT2MAT blasfeo_unpack_tran_smat
#define CVT_VEC2STRVEC blasfeo_pack_svec
#define CVT_STRVEC2VEC blasfeo_unpack_svec
#define GESE blasfeo_sgese
#define OCP_QP s_ocp_qp
#define OCP_QP_DIM s_ocp_qp_dim
#define REAL float
#define STRMAT blasfeo_smat
#define STRVEC blasfeo_svec
#define SIZE_STRMAT blasfeo_memsize_smat
#define SIZE_STRVEC blasfeo_memsize_svec
#define VECCP_LIBSTR blasfeo_sveccp
#define VECSC_LIBSTR blasfeo_svecsc
#define VECSE blasfeo_svecse

#define OCP_QP_STRSIZE s_ocp_qp_strsize
#define OCP_QP_MEMSIZE s_ocp_qp_memsize
#define OCP_QP_CREATE s_ocp_qp_create
#define OCP_QP_SET_ALL s_ocp_qp_set_all
#define OCP_QP_SET s_ocp_qp_set
#define OCP_QP s_ocp_qp
#define OCP_QP_SET_A s_ocp_qp_set_A
#define OCP_QP_GET_A s_cvt_ocp_qp_A
#define OCP_QP_SET_B s_ocp_qp_set_B
#define OCP_QP_GET_B s_cvt_ocp_qp_B
#define OCP_QP_SET_BVEC s_ocp_qp_set_b
#define OCP_QP_GET_BVEC s_cvt_ocp_qp_b
#define OCP_QP_SET_Q s_ocp_qp_set_Q
#define OCP_QP_GET_Q s_cvt_ocp_qp_Q
#define OCP_QP_SET_S s_ocp_qp_set_S
#define OCP_QP_GET_S s_cvt_ocp_qp_S
#define OCP_QP_SET_R s_ocp_qp_set_R
#define OCP_QP_GET_R s_cvt_ocp_qp_R
#define OCP_QP_SET_QVEC s_ocp_qp_set_q
#define OCP_QP_GET_QVEC s_cvt_ocp_qp_q
#define OCP_QP_SET_RVEC s_ocp_qp_set_r
#define OCP_QP_GET_RVEC s_cvt_ocp_qp_r
#define OCP_QP_SET_LBX s_ocp_qp_set_lbx
#define OCP_QP_GET_LBX s_cvt_ocp_qp_lbx
#define OCP_QP_SET_UBX s_ocp_qp_set_ubx
#define OCP_QP_GET_UBX s_cvt_ocp_qp_ubx
#define OCP_QP_SET_LBU s_ocp_qp_set_lbu
#define OCP_QP_GET_LBU s_cvt_ocp_qp_lbu
#define OCP_QP_SET_UBU s_ocp_qp_set_ubu
#define OCP_QP_GET_UBU s_cvt_ocp_qp_ubu
#define OCP_QP_SET_IDXB s_ocp_qp_set_idxb
#define OCP_QP_GET_IDXB s_cvt_ocp_qp_idxb
#define OCP_QP_SET_IDXBX s_ocp_qp_set_idxbx
#define OCP_QP_GET_IDXBX s_cvt_ocp_qp_idxbx
#define OCP_QP_SET_JBX s_ocp_qp_set_Jbx
#define OCP_QP_GET_JBX s_ocp_qp_get_Jbx
#define OCP_QP_SET_IDXBU s_ocp_qp_set_idxbu
#define OCP_QP_GET_IDXBU s_cvt_ocp_qp_idxbu
#define OCP_QP_SET_JBU s_ocp_qp_set_Jbu
#define OCP_QP_GET_JBU s_ocp_qp_get_Jbu
#define OCP_QP_SET_LB s_ocp_qp_set_lb
#define OCP_QP_GET_LB s_cvt_ocp_qp_lb
#define OCP_QP_SET_UB s_ocp_qp_set_ub
#define OCP_QP_GET_UB s_cvt_ocp_qp_ub
#define OCP_QP_SET_C s_ocp_qp_set_C
#define OCP_QP_GET_C s_cvt_ocp_qp_C
#define OCP_QP_SET_D s_ocp_qp_set_D
#define OCP_QP_GET_D s_cvt_ocp_qp_D
#define OCP_QP_SET_LG s_ocp_qp_set_lg
#define OCP_QP_GET_LG s_cvt_ocp_qp_lg
#define OCP_QP_SET_UG s_ocp_qp_set_ug
#define OCP_QP_GET_UG s_cvt_ocp_qp_ug
#define OCP_QP_SET_ZL s_ocp_qp_set_Zl
#define OCP_QP_GET_ZL s_cvt_ocp_qp_Zl
#define OCP_QP_SET_ZU s_ocp_qp_set_Zu
#define OCP_QP_GET_ZU s_cvt_ocp_qp_Zu
#define OCP_QP_SET_ZLVEC s_ocp_qp_set_zl
#define OCP_QP_GET_ZLVEC s_cvt_ocp_qp_zl
#define OCP_QP_SET_ZUVEC s_ocp_qp_set_zu
#define OCP_QP_GET_ZUVEC s_cvt_ocp_qp_zu
#define OCP_QP_SET_IDXS s_ocp_qp_set_idxs
#define OCP_QP_GET_IDXS s_cvt_ocp_qp_idxs
#define OCP_QP_SET_LLS s_ocp_qp_set_lls
#define OCP_QP_GET_LLS s_cvt_ocp_qp_lls
#define OCP_QP_SET_LUS s_ocp_qp_set_lus
#define OCP_QP_GET_LUS s_cvt_ocp_qp_lus
// TODO remove ???
#define CHANGE_BOUNDS_DIMENSIONS_OCP_QP s_change_bounds_dimensions_ocp_qp

#include "x_ocp_qp.c"

