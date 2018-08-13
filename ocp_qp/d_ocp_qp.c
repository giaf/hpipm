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

#include <hpipm_d_ocp_qp_dim.h>
#include <hpipm_d_ocp_qp.h>



#define CREATE_STRMAT blasfeo_create_dmat
#define CREATE_STRVEC blasfeo_create_dvec
#define CVT_MAT2STRMAT blasfeo_pack_dmat
#define CVT_STRMAT2MAT blasfeo_unpack_dmat
#define CVT_TRAN_MAT2STRMAT blasfeo_pack_tran_dmat
#define CVT_TRAN_STRMAT2MAT blasfeo_unpack_tran_dmat
#define CVT_VEC2STRVEC blasfeo_pack_dvec
#define CVT_STRVEC2VEC blasfeo_unpack_dvec
#define GESE blasfeo_dgese
#define OCP_QP d_ocp_qp
#define OCP_QP_DIM d_ocp_qp_dim
#define REAL double
#define STRMAT blasfeo_dmat
#define STRVEC blasfeo_dvec
#define SIZE_STRMAT blasfeo_memsize_dmat
#define SIZE_STRVEC blasfeo_memsize_dvec
#define VECCP_LIBSTR blasfeo_dveccp
#define VECSC_LIBSTR blasfeo_dvecsc
#define VECSE blasfeo_dvecse

#define CAST_OCP_QP d_cast_ocp_qp
#define COPY_OCP_QP d_copy_ocp_qp
#define CREATE_OCP_QP d_create_ocp_qp
#define CVT_COLMAJ_TO_OCP_QP d_cvt_colmaj_to_ocp_qp
#define CVT_ROWMAJ_TO_OCP_QP d_cvt_rowmaj_to_ocp_qp
#define MEMSIZE_OCP_QP d_memsize_ocp_qp
#define CHANGE_BOUNDS_DIMENSIONS_OCP_QP d_change_bounds_dimensions_ocp_qp
#define CVT_COLMAJ_TO_OCP_QP_A d_cvt_colmaj_to_ocp_qp_A
#define CVT_OCP_QP_TO_COLMAJ_A d_cvt_ocp_qp_to_colmaj_A
#define CVT_COLMAJ_TO_OCP_QP_B d_cvt_colmaj_to_ocp_qp_B
#define CVT_OCP_QP_TO_COLMAJ_B d_cvt_ocp_qp_to_colmaj_B
#define CVT_COLMAJ_TO_OCP_QP_BVEC d_cvt_colmaj_to_ocp_qp_b
#define CVT_OCP_QP_TO_COLMAJ_BVEC d_cvt_ocp_qp_to_colmaj_b
#define CVT_COLMAJ_TO_OCP_QP_Q d_cvt_colmaj_to_ocp_qp_Q
#define CVT_OCP_QP_TO_COLMAJ_Q d_cvt_ocp_qp_to_colmaj_Q
#define CVT_COLMAJ_TO_OCP_QP_S d_cvt_colmaj_to_ocp_qp_S
#define CVT_OCP_QP_TO_COLMAJ_S d_cvt_ocp_qp_to_colmaj_S
#define CVT_COLMAJ_TO_OCP_QP_R d_cvt_colmaj_to_ocp_qp_R
#define CVT_OCP_QP_TO_COLMAJ_R d_cvt_ocp_qp_to_colmaj_R
#define CVT_COLMAJ_TO_OCP_QP_QVEC d_cvt_colmaj_to_ocp_qp_q
#define CVT_OCP_QP_TO_COLMAJ_QVEC d_cvt_ocp_qp_to_colmaj_q
#define CVT_COLMAJ_TO_OCP_QP_RVEC d_cvt_colmaj_to_ocp_qp_r
#define CVT_OCP_QP_TO_COLMAJ_RVEC d_cvt_ocp_qp_to_colmaj_r
#define CVT_COLMAJ_TO_OCP_QP_LBX d_cvt_colmaj_to_ocp_qp_lbx
#define CVT_OCP_QP_TO_COLMAJ_LBX d_cvt_ocp_qp_to_colmaj_lbx
#define CVT_COLMAJ_TO_OCP_QP_UBX d_cvt_colmaj_to_ocp_qp_ubx
#define CVT_OCP_QP_TO_COLMAJ_UBX d_cvt_ocp_qp_to_colmaj_ubx
#define CVT_COLMAJ_TO_OCP_QP_LBU d_cvt_colmaj_to_ocp_qp_lbu
#define CVT_OCP_QP_TO_COLMAJ_LBU d_cvt_ocp_qp_to_colmaj_lbu
#define CVT_COLMAJ_TO_OCP_QP_UBU d_cvt_colmaj_to_ocp_qp_ubu
#define CVT_OCP_QP_TO_COLMAJ_UBU d_cvt_ocp_qp_to_colmaj_ubu
#define CVT_COLMAJ_TO_OCP_QP_IDXB d_cvt_colmaj_to_ocp_qp_idxb
#define CVT_OCP_QP_TO_COLMAJ_IDXB d_cvt_ocp_qp_to_colmaj_idxb
#define CVT_COLMAJ_TO_OCP_QP_LB d_cvt_colmaj_to_ocp_qp_lb
#define CVT_OCP_QP_TO_COLMAJ_LB d_cvt_ocp_qp_to_colmaj_lb
#define CVT_COLMAJ_TO_OCP_QP_UB d_cvt_colmaj_to_ocp_qp_ub
#define CVT_OCP_QP_TO_COLMAJ_UB d_cvt_ocp_qp_to_colmaj_ub
#define CVT_COLMAJ_TO_OCP_QP_C d_cvt_colmaj_to_ocp_qp_C
#define CVT_OCP_QP_TO_COLMAJ_C d_cvt_ocp_qp_to_colmaj_C
#define CVT_COLMAJ_TO_OCP_QP_D d_cvt_colmaj_to_ocp_qp_D
#define CVT_OCP_QP_TO_COLMAJ_D d_cvt_ocp_qp_to_colmaj_D
#define CVT_COLMAJ_TO_OCP_QP_LG d_cvt_colmaj_to_ocp_qp_lg
#define CVT_OCP_QP_TO_COLMAJ_LG d_cvt_ocp_qp_to_colmaj_lg
#define CVT_COLMAJ_TO_OCP_QP_UG d_cvt_colmaj_to_ocp_qp_ug
#define CVT_OCP_QP_TO_COLMAJ_UG d_cvt_ocp_qp_to_colmaj_ug
#define CVT_COLMAJ_TO_OCP_QP_ZL d_cvt_colmaj_to_ocp_qp_Zl
#define CVT_OCP_QP_TO_COLMAJ_ZL d_cvt_ocp_qp_to_colmaj_Zl
#define CVT_COLMAJ_TO_OCP_QP_ZU d_cvt_colmaj_to_ocp_qp_Zu
#define CVT_OCP_QP_TO_COLMAJ_ZU d_cvt_ocp_qp_to_colmaj_Zu
#define CVT_COLMAJ_TO_OCP_QP_ZLVEC d_cvt_colmaj_to_ocp_qp_zl
#define CVT_OCP_QP_TO_COLMAJ_ZLVEC d_cvt_ocp_qp_to_colmaj_zl
#define CVT_COLMAJ_TO_OCP_QP_ZUVEC d_cvt_colmaj_to_ocp_qp_zu
#define CVT_OCP_QP_TO_COLMAJ_ZUVEC d_cvt_ocp_qp_to_colmaj_zu
#define CVT_COLMAJ_TO_OCP_QP_IDXS d_cvt_colmaj_to_ocp_qp_idxs
#define CVT_OCP_QP_TO_COLMAJ_IDXS d_cvt_ocp_qp_to_colmaj_idxs
#define CVT_COLMAJ_TO_OCP_QP_LS d_cvt_colmaj_to_ocp_qp_ls
#define CVT_OCP_QP_TO_COLMAJ_LS d_cvt_ocp_qp_to_colmaj_ls
#define CVT_COLMAJ_TO_OCP_QP_US d_cvt_colmaj_to_ocp_qp_us
#define CVT_OCP_QP_TO_COLMAJ_US d_cvt_ocp_qp_to_colmaj_us

// interface functions
#define SIZEOF_OCP_QP d_sizeof_ocp_qp

#include "x_ocp_qp.c"
