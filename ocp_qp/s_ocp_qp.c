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



#if defined(RUNTIME_CHECKS)
#include <stdlib.h>
#include <stdio.h>
#endif

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_s_aux.h>

#include "../include/hpipm_s_ocp_qp_dim.h"
#include "../include/hpipm_s_ocp_qp.h"



#define CREATE_STRMAT blasfeo_create_smat
#define CREATE_STRVEC blasfeo_create_svec
#define CVT_MAT2STRMAT blasfeo_pack_smat
#define CVT_STRMAT2MAT blasfeo_unpack_smat
#define CVT_TRAN_MAT2STRMAT blasfeo_pack_tran_smat
#define CVT_TRAN_STRMAT2MAT blasfeo_unpack_tran_smat
#define CVT_VEC2STRVEC blasfeo_pack_svec
#define CVT_STRVEC2VEC blasfeo_unpack_svec
#define OCP_QP s_ocp_qp
#define OCP_QP_DIM s_ocp_qp_dim
#define REAL float
#define STRMAT blasfeo_smat
#define STRVEC blasfeo_svec
#define SIZE_STRMAT blasfeo_memsize_smat
#define SIZE_STRVEC blasfeo_memsize_svec
#define VECCP_LIBSTR blasfeo_sveccp
#define VECSC_LIBSTR blasfeo_svecsc

#define CAST_OCP_QP s_cast_ocp_qp
#define COPY_OCP_QP s_copy_ocp_qp
#define CREATE_OCP_QP s_create_ocp_qp
#define CVT_COLMAJ_TO_OCP_QP s_cvt_colmaj_to_ocp_qp
#define CVT_ROWMAJ_TO_OCP_QP s_cvt_rowmaj_to_ocp_qp
#define MEMSIZE_OCP_QP s_memsize_ocp_qp

#define UPDATE_Q s_update_Q
#define UPDATE_S s_update_S
#define UPDATE_R s_update_R
#define UPDATE_QVEC s_update_q
#define UPDATE_RVEC s_update_r
#define UPDATE_A s_update_A
#define UPDATE_B s_update_B
#define UPDATE_BVEC s_update_b
#define UPDATE_LBX s_update_lbx
#define UPDATE_UBX s_update_ubx
#define UPDATE_LBU s_update_lbu
#define UPDATE_UBU s_update_ubu
#define UPDATE_C s_update_C
#define UPDATE_D s_update_D
#define UPDATE_LG s_update_lg
#define UPDATE_UG s_update_ug

#define COPY_Q s_copy_Q
#define COPY_S s_copy_S
#define COPY_R s_copy_R
#define COPY_QVEC s_copy_q
#define COPY_RVEC s_copy_r
#define COPY_A s_copy_A
#define COPY_B s_copy_B
#define COPY_BVEC s_copy_b
#define COPY_LBX s_copy_lbx
#define COPY_UBX s_copy_ubx
#define COPY_LBU s_copy_lbu
#define COPY_UBU s_copy_ubu
#define COPY_C s_copy_C
#define COPY_D s_copy_D
#define COPY_LG s_copy_lg
#define COPY_UG s_copy_ug

#include "x_ocp_qp.c"

