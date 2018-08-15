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

#define DENSE_QP_SET_H s_dense_qp_set_H
#define DENSE_QP_SET_G s_dense_qp_set_g
#define DENSE_QP_SET_A s_dense_qp_set_A
#define DENSE_QP_SET_B s_dense_qp_set_b
#define DENSE_QP_SET_IDXB s_dense_qp_set_idxb
#define DENSE_QP_SET_LB s_dense_qp_set_lb
#define DENSE_QP_SET_UB s_dense_qp_set_ub
#define DENSE_QP_SET_C s_dense_qp_set_C
#define DENSE_QP_SET_LG s_dense_qp_set_lg
#define DENSE_QP_SET_UG s_dense_qp_set_ug
#define DENSE_QP_SET_IDXS s_dense_qp_set_idxs
#define DENSE_QP_SET_ZZL s_dense_qp_set_Zl
#define DENSE_QP_SET_ZZU s_dense_qp_set_Zu
#define DENSE_QP_SET_ZL s_dense_qp_set_zl
#define DENSE_QP_SET_ZU s_dense_qp_set_zu
#define DENSE_QP_SET_LS s_dense_qp_set_ls
#define DENSE_QP_SET_US s_dense_qp_set_us

#define DENSE_QP_GET_H s_dense_qp_get_H
#define DENSE_QP_GET_G s_dense_qp_get_g
#define DENSE_QP_GET_A s_dense_qp_get_A
#define DENSE_QP_GET_B s_dense_qp_get_b
#define DENSE_QP_GET_IDXB s_dense_qp_get_idxb
#define DENSE_QP_GET_LB s_dense_qp_get_lb
#define DENSE_QP_GET_UB s_dense_qp_get_ub
#define DENSE_QP_GET_C s_dense_qp_get_C
#define DENSE_QP_GET_LG s_dense_qp_get_lg
#define DENSE_QP_GET_UG s_dense_qp_get_ug
#define DENSE_QP_GET_IDXS s_dense_qp_get_idxs
#define DENSE_QP_GET_ZZL s_dense_qp_get_Zl
#define DENSE_QP_GET_ZZU s_dense_qp_get_Zu
#define DENSE_QP_GET_ZL s_dense_qp_get_zl
#define DENSE_QP_GET_ZU s_dense_qp_get_zu
#define DENSE_QP_GET_LS s_dense_qp_get_ls
#define DENSE_QP_GET_US s_dense_qp_get_us

#include "x_dense_qp.c"

