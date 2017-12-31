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
#include "../include/hpipm_d_ocp_qp_sol.h"



#define CREATE_STRVEC blasfeo_create_dvec
#define CVT_STRVEC2VEC blasfeo_unpack_dvec
#define CVT_VEC2STRVEC blasfeo_pack_dvec
#define OCP_QP d_ocp_qp
#define OCP_QP_DIM d_ocp_qp_dim
#define OCP_QP_SOL d_ocp_qp_sol
#define REAL double
#define STRVEC blasfeo_dvec
#define SIZE_STRVEC blasfeo_memsize_dvec
#define VECCP_LIBSTR blasfeo_dveccp

#define CREATE_OCP_QP_SOL d_create_ocp_qp_sol
#define MEMSIZE_OCP_QP_SOL d_memsize_ocp_qp_sol
#define CVT_OCP_QP_SOL_TO_COLMAJ d_cvt_ocp_qp_sol_to_colmaj
#define CVT_COLMAJ_TO_OCP_QP_SOL d_cvt_colmaj_to_ocp_qp_sol
#define CVT_OCP_QP_SOL_TO_ROWMAJ d_cvt_ocp_qp_sol_to_rowmaj
#define CVT_OCP_QP_SOL_TO_LIBSTR d_cvt_ocp_qp_sol_to_libstr



#include "x_ocp_qp_sol.c"
