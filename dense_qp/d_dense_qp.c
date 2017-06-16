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
#endif

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_d_aux.h>

#include "../include/hpipm_d_dense_qp.h"


#define CREATE_STRMAT d_create_strmat
#define CREATE_STRVEC d_create_strvec
#define CVT_MAT2STRMAT d_cvt_mat2strmat
#define CVT_TRAN_MAT2STRMAT d_cvt_tran_mat2strmat
#define CVT_TRAN_STRMAT2MAT d_cvt_tran_strmat2mat
#define CVT_VEC2STRVEC d_cvt_vec2strvec
#define CVT_STRMAT2MAT d_cvt_strmat2mat
#define CVT_STRVEC2VEC d_cvt_strvec2vec
#define DENSE_QP_DIM d_dense_qp_dim
#define DENSE_QP_VEC d_dense_qp_vec
#define DENSE_QP_MAT d_dense_qp_mat
#define DENSE_QP d_dense_qp
#define GECP_LIBSTR dgecp_libstr
#define GETR_LIBSTR dgetr_libstr
#define REAL double
#define ROWIN_LIBSTR drowin_libstr
#define SIZE_STRMAT d_size_strmat
#define SIZE_STRVEC d_size_strvec
#define STRMAT d_strmat
#define STRVEC d_strvec
#define VECCP_LIBSTR dveccp_libstr

#define MEMSIZE_DENSE_QP d_memsize_dense_qp
#define CREATE_DENSE_QP d_create_dense_qp
#define CVT_COLMAJ_TO_DENSE_QP d_cvt_colmaj_to_dense_qp
#define CVT_DENSE_QP_TO_COLMAJ d_cvt_dense_qp_to_colmaj
#define CVT_ROWMAJ_TO_DENSE_QP d_cvt_rowmaj_to_dense_qp
#define CVT_DENSE_QP_TO_ROWMAJ d_cvt_dense_qp_to_rowmaj
#define CVT_LIBSTR_TO_DENSE_QP d_cvt_libstr_to_dense_qp
#define CVT_DENSE_QP_TO_LIBSTR d_cvt_dense_qp_to_libstr
#define CAST_DENSE_QP_DIM d_cast_dense_qp_dim
//#define CREATE_DENSE_QP d_create_dense_qp
//#define COPY_DENSE_QP d_copy_dense_qp



#include "x_dense_qp.c"
