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
#include <blasfeo_d_aux.h>

#include "../include/hpipm_d_dense_qp_dim.h"
#include "../include/hpipm_d_dense_qp.h"


#define CREATE_STRMAT blasfeo_create_dmat
#define CREATE_STRVEC blasfeo_create_dvec
#define CVT_MAT2STRMAT blasfeo_pack_dmat
#define CVT_TRAN_MAT2STRMAT blasfeo_pack_tran_dmat
#define CVT_TRAN_STRMAT2MAT blasfeo_unpack_tran_dmat
#define CVT_VEC2STRVEC blasfeo_pack_dvec
#define CVT_STRMAT2MAT blasfeo_unpack_dmat
#define CVT_STRVEC2VEC blasfeo_unpack_dvec
#define DENSE_QP d_dense_qp
#define DENSE_QP_DIM d_dense_qp_dim
#define GECP_LIBSTR blasfeo_dgecp
#define GETR_LIBSTR blasfeo_dgetr
#define REAL double
#define ROWIN_LIBSTR blasfeo_drowin
#define SIZE_STRMAT blasfeo_memsize_dmat
#define SIZE_STRVEC blasfeo_memsize_dvec
#define STRMAT blasfeo_dmat
#define STRVEC blasfeo_dvec
#define VECCP_LIBSTR blasfeo_dveccp
#define VECSC_LIBSTR blasfeo_dvecsc
#define VECSE_LIBSTR blasfeo_dvecse

#define MEMSIZE_DENSE_QP d_memsize_dense_qp
#define CREATE_DENSE_QP d_create_dense_qp
#define CVT_COLMAJ_TO_DENSE_QP d_cvt_colmaj_to_dense_qp
#define CVT_DENSE_QP_TO_COLMAJ d_cvt_dense_qp_to_colmaj
#define CVT_ROWMAJ_TO_DENSE_QP d_cvt_rowmaj_to_dense_qp
#define CVT_DENSE_QP_TO_ROWMAJ d_cvt_dense_qp_to_rowmaj
#define CVT_LIBSTR_TO_DENSE_QP d_cvt_libstr_to_dense_qp
#define CVT_DENSE_QP_TO_LIBSTR d_cvt_dense_qp_to_libstr
#define CAST_DENSE_QP_DIM d_cast_dense_qp_dim



#include "x_dense_qp.c"
