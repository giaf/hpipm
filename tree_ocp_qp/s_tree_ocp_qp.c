/**************************************************************************************************
*                                                                                                 *
* This file is part of HPIPM.                                                                     *
*                                                                                                 *
* HPIPM -- High Performance Interior Point Method.                                                *
* Copyright (C) 2017 by Gianluca Frison.                                                          *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
*                                                                                                 *
* HPIPM is free software; you can redistribute it and/or                                          *
* modify it under the terms of the GNU Lesser General Public                                      *
* License as published by the Free Software Foundation; either                                    *
* version 2.1 of the License, or (at your option) any later version.                              *
*                                                                                                 *
* HPIPM is distributed in the hope that it will be useful,                                        *
* but WITHOUT ANY WARRANTY; without even the implied warranty of                                  *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                                            *
* See the GNU Lesser General Public License for more details.                                     *
*                                                                                                 *
* You should have received a copy of the GNU Lesser General Public                                *
* License along with HPIPM; if not, write to the Free Software                                    *
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

#include "../include/hpipm_tree.h"
#include "../include/hpipm_scenario_tree.h"
#include "../include/hpipm_s_tree_ocp_qp_dim.h"
#include "../include/hpipm_s_tree_ocp_qp.h"

#define CREATE_STRMAT blasfeo_create_smat
#define CREATE_STRVEC blasfeo_create_svec
#define CVT_MAT2STRMAT blasfeo_pack_smat
#define CVT_TRAN_MAT2STRMAT blasfeo_pack_tran_smat
#define CVT_VEC2STRVEC blasfeo_pack_svec
#define REAL float
#define SIZE_STRMAT blasfeo_memsize_smat
#define SIZE_STRVEC blasfeo_memsize_svec
#define STRMAT blasfeo_smat
#define STRVEC blasfeo_svec
#define TREE_OCP_QP s_tree_ocp_qp
#define TREE_OCP_QP_DIM s_tree_ocp_qp_dim
#define VECSC_LIBSTR blasfeo_svecsc

#define MEMSIZE_TREE_OCP_QP s_memsize_tree_ocp_qp
#define CREATE_TREE_OCP_QP s_create_tree_ocp_qp
#define CVT_COLMAJ_TO_TREE_OCP_QP s_cvt_colmaj_to_tree_ocp_qp



#include "x_tree_ocp_qp.c"

