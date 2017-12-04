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
#include "../include/hpipm_s_tree_ocp_qp.h"
#include "../include/hpipm_s_tree_ocp_qp_sol.h"

#define CREATE_STRVEC s_create_strvec
#define CVT_STRVEC2VEC s_cvt_strvec2vec
#define REAL float
#define SIZE_STRVEC s_size_strvec
#define STRVEC s_strvec
#define TREE_OCP_QP s_tree_ocp_qp
#define TREE_OCP_QP_DIM s_tree_ocp_qp_dim
#define TREE_OCP_QP_SOL s_tree_ocp_qp_sol

#define MEMSIZE_TREE_OCP_QP_SOL s_memsize_tree_ocp_qp_sol
#define CREATE_TREE_OCP_QP_SOL s_create_tree_ocp_qp_sol
#define CVT_TREE_OCP_QP_SOL_TO_COLMAJ s_cvt_tree_ocp_qp_sol_to_colmaj
#define CVT_TREE_OCP_QP_SOL_TO_ROWMAJ s_cvt_tree_ocp_qp_sol_to_rowmaj



#include "x_tree_ocp_qp_sol.c"

