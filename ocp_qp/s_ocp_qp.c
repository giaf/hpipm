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
#define CVT_TRAN_MAT2STRMAT blasfeo_pack_tran_smat
#define CVT_VEC2STRVEC blasfeo_pack_svec
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



#include "x_ocp_qp.c"

