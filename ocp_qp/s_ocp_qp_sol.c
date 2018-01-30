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
#include "../include/hpipm_s_ocp_qp_sol.h"



#define CREATE_STRVEC blasfeo_create_svec
#define CVT_STRVEC2VEC blasfeo_unpack_svec
#define CVT_VEC2STRVEC blasfeo_pack_svec
#define OCP_QP s_ocp_qp
#define OCP_QP_DIM s_ocp_qp_dim
#define OCP_QP_SOL s_ocp_qp_sol
#define REAL float
#define STRVEC blasfeo_svec
#define SIZE_STRVEC blasfeo_memsize_svec
#define VECCP_LIBSTR blasfeo_sveccp

#define CREATE_OCP_QP_SOL s_create_ocp_qp_sol
#define MEMSIZE_OCP_QP_SOL s_memsize_ocp_qp_sol
#define CVT_OCP_QP_SOL_TO_COLMAJ s_cvt_ocp_qp_sol_to_colmaj
#define CVT_COLMAJ_TO_OCP_QP_SOL s_cvt_colmaj_to_ocp_qp_sol
#define CVT_OCP_QP_SOL_TO_ROWMAJ s_cvt_ocp_qp_sol_to_rowmaj
#define CVT_OCP_QP_SOL_TO_LIBSTR s_cvt_ocp_qp_sol_to_libstr
#define CVT_OCP_QP_SOL_TO_COLMAJ_X s_cvt_ocp_qp_sol_to_colmaj_x
#define CVT_OCP_QP_SOL_TO_COLMAJ_U s_cvt_ocp_qp_sol_to_colmaj_u
#define CVT_OCP_QP_SOL_TO_COLMAJ_PI s_cvt_ocp_qp_sol_to_colmaj_pi
#define CVT_OCP_QP_SOL_TO_COLMAJ_LAM_LB s_cvt_ocp_qp_sol_to_colmaj_lam_lb
#define CVT_OCP_QP_SOL_TO_COLMAJ_LAM_UB s_cvt_ocp_qp_sol_to_colmaj_lam_ub
#define CVT_OCP_QP_SOL_TO_COLMAJ_LAM_LG s_cvt_ocp_qp_sol_to_colmaj_lam_lg
#define CVT_OCP_QP_SOL_TO_COLMAJ_LAM_UG s_cvt_ocp_qp_sol_to_colmaj_lam_ug


#include "x_ocp_qp_sol.c"

