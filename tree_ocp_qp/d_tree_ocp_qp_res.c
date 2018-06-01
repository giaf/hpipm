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

#include "../include/hpipm_d_tree_ocp_qp_dim.h"
#include "../include/hpipm_d_tree_ocp_qp_res.h"



#define CREATE_STRVEC blasfeo_create_dvec
#define CVT_STRVEC2VEC blasfeo_unpack_dvec
#define TREE_OCP_QP_DIM d_tree_ocp_qp_dim
#define TREE_OCP_QP_RES d_tree_ocp_qp_res
#define TREE_OCP_QP_RES_WORKSPACE d_tree_ocp_qp_res_workspace
#define REAL double
#define SIZE_STRVEC blasfeo_memsize_dvec
#define STRVEC blasfeo_dvec



#define MEMSIZE_TREE_OCP_QP_RES d_memsize_tree_ocp_qp_res
#define CREATE_TREE_OCP_QP_RES d_create_tree_ocp_qp_res
#define MEMSIZE_TREE_OCP_QP_RES_WORKSPACE d_memsize_tree_ocp_qp_res_workspace
#define CREATE_TREE_OCP_QP_RES_WORKSPACE d_create_tree_ocp_qp_res_workspace
#define CVT_TREE_OCP_QP_RES_TO_COLMAJ d_cvt_tree_ocp_qp_res_to_colmaj
#define CVT_TREE_OCP_QP_RES_TO_ROWMAJ d_cvt_tree_ocp_qp_res_to_rowmaj



#include "x_tree_ocp_qp_res.c"

