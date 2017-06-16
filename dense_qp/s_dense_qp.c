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
#include <blasfeo_s_aux.h>

#include "../include/hpipm_s_dense_qp.h"


#define CREATE_STRMAT s_create_strmat
#define CREATE_STRVEC s_create_strvec
#define CVT_MAT2STRMAT s_cvt_mat2strmat
#define CVT_TRAN_MAT2STRMAT s_cvt_tran_mat2strmat
#define CVT_TRAN_STRMAT2MAT s_cvt_tran_strmat2mat
#define CVT_VEC2STRVEC s_cvt_vec2strvec
#define CVT_STRMAT2MAT s_cvt_strmat2mat
#define CVT_STRVEC2VEC s_cvt_strvec2vec
#define DENSE_QP_DIM s_dense_qp_dim
#define DENSE_QP_VEC s_dense_qp_vec
#define DENSE_QP_MAT s_dense_qp_mat
#define DENSE_QP s_dense_qp
#define GECP_LIBSTR sgecp_libstr
#define GETR_LIBSTR sgetr_libstr
#define REAL float
#define ROWIN_LIBSTR srowin_libstr
#define SIZE_STRMAT s_size_strmat
#define SIZE_STRVEC s_size_strvec
#define STRMAT s_strmat
#define STRVEC s_strvec
#define VECCP_LIBSTR sveccp_libstr

#define MEMSIZE_DENSE_QP s_memsize_dense_qp
#define CREATE_DENSE_QP s_create_dense_qp
#define CVT_COLMAJ_TO_DENSE_QP s_cvt_colmaj_to_dense_qp
#define CVT_DENSE_QP_TO_COLMAJ s_cvt_dense_qp_to_colmaj
#define CVT_ROWMAJ_TO_DENSE_QP s_cvt_rowmaj_to_dense_qp
#define CVT_DENSE_QP_TO_ROWMAJ s_cvt_dense_qp_to_rowmaj
#define CVT_LIBSTR_TO_DENSE_QP s_cvt_libstr_to_dense_qp
#define CVT_DENSE_QP_TO_LIBSTR s_cvt_dense_qp_to_libstr
#define CAST_DENSE_QP_DIM s_cast_dense_qp_dim
//#define CREATE_DENSE_QP s_create_dense_qp
//#define COPY_DENSE_QP s_copy_dense_qp



#include "x_dense_qp.c"

