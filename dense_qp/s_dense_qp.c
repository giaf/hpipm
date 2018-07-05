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




#include <stdlib.h>
#include <stdio.h>

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_s_aux.h>

#include <hpipm_s_dense_qp_dim.h>
#include <hpipm_s_dense_qp.h>


#define CREATE_STRMAT blasfeo_create_smat
#define CREATE_STRVEC blasfeo_create_svec
#define CVT_MAT2STRMAT blasfeo_pack_smat
#define CVT_TRAN_MAT2STRMAT blasfeo_pack_tran_smat
#define CVT_TRAN_STRMAT2MAT blasfeo_unpack_tran_smat
#define CVT_VEC2STRVEC blasfeo_pack_svec
#define CVT_STRMAT2MAT blasfeo_unpack_smat
#define CVT_STRVEC2VEC blasfeo_unpack_svec
#define DENSE_QP s_dense_qp
#define DENSE_QP_DIM s_dense_qp_dim
#define GECP_LIBSTR blasfeo_sgecp
#define GETR_LIBSTR blasfeo_sgetr
#define REAL float
#define ROWIN_LIBSTR blasfeo_srowin
#define SIZE_STRMAT blasfeo_memsize_smat
#define SIZE_STRVEC blasfeo_memsize_svec
#define STRMAT blasfeo_smat
#define STRVEC blasfeo_svec
#define VECCP_LIBSTR blasfeo_sveccp
#define VECSC_LIBSTR blasfeo_svecsc
#define VECSE_LIBSTR blasfeo_svecse

#define MEMSIZE_DENSE_QP s_memsize_dense_qp
#define CREATE_DENSE_QP s_create_dense_qp
#define CVT_COLMAJ_TO_DENSE_QP s_cvt_colmaj_to_dense_qp
#define CVT_DENSE_QP_TO_COLMAJ s_cvt_dense_qp_to_colmaj
#define CVT_ROWMAJ_TO_DENSE_QP s_cvt_rowmaj_to_dense_qp
#define CVT_DENSE_QP_TO_ROWMAJ s_cvt_dense_qp_to_rowmaj
#define CVT_LIBSTR_TO_DENSE_QP s_cvt_libstr_to_dense_qp
#define CVT_DENSE_QP_TO_LIBSTR s_cvt_dense_qp_to_libstr
#define CAST_DENSE_QP_DIM s_cast_dense_qp_dim



#include "x_dense_qp.c"

