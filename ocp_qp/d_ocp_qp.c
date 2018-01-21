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
#include <blasfeo_d_aux.h>

#include "../include/hpipm_d_ocp_qp_dim.h"
#include "../include/hpipm_d_ocp_qp.h"



#define CREATE_STRMAT blasfeo_create_dmat
#define CREATE_STRVEC blasfeo_create_dvec
#define CVT_MAT2STRMAT blasfeo_pack_dmat
#define CVT_TRAN_MAT2STRMAT blasfeo_pack_tran_dmat
#define CVT_VEC2STRVEC blasfeo_pack_dvec
#define OCP_QP d_ocp_qp
#define OCP_QP_DIM d_ocp_qp_dim
#define REAL double
#define STRMAT blasfeo_dmat
#define STRVEC blasfeo_dvec
#define SIZE_STRMAT blasfeo_memsize_dmat
#define SIZE_STRVEC blasfeo_memsize_dvec
#define VECCP_LIBSTR blasfeo_dveccp
#define VECSC_LIBSTR blasfeo_dvecsc

#define CAST_OCP_QP d_cast_ocp_qp
#define COPY_OCP_QP d_copy_ocp_qp
#define CREATE_OCP_QP d_create_ocp_qp
#define CVT_COLMAJ_TO_OCP_QP d_cvt_colmaj_to_ocp_qp
#define CVT_ROWMAJ_TO_OCP_QP d_cvt_rowmaj_to_ocp_qp
#define MEMSIZE_OCP_QP d_memsize_ocp_qp

#define UPDATE_Q d_update_Q
#define UPDATE_S d_update_S
#define UPDATE_R d_update_R
#define UPDATE_QVEC d_update_q
#define UPDATE_RVEC d_update_r
#define UPDATE_A d_update_A
#define UPDATE_B d_update_B
#define UPDATE_BVEC d_update_b
#define UPDATE_LBX d_update_lbx
#define UPDATE_UBX d_update_ubx
#define UPDATE_LBU d_update_lbu
#define UPDATE_UBU d_update_ubu
#define UPDATE_C d_update_C
#define UPDATE_D d_update_D
#define UPDATE_LG d_update_lg
#define UPDATE_UG d_update_ug

#include "x_ocp_qp.c"
