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



#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_d_aux.h>

#include "../include/hpipm_tree.h"
#include "../include/hpipm_scenario_tree.h"
#include "../include/hpipm_d_tree_ocp_qp.h"

#define CREATE_STRMAT d_create_strmat
#define CREATE_STRVEC d_create_strvec
#define CVT_MAT2STRMAT d_cvt_mat2strmat
#define CVT_TRAN_MAT2STRMAT d_cvt_tran_mat2strmat
#define CVT_VEC2STRVEC d_cvt_vec2strvec
#define REAL double
#define SIZE_STRMAT d_size_strmat
#define SIZE_STRVEC d_size_strvec
#define STRMAT d_strmat
#define STRVEC d_strvec
#define TREE_OCP_QP d_tree_ocp_qp

#define MEMSIZE_TREE_OCP_QP d_memsize_tree_ocp_qp
#define CREATE_TREE_OCP_QP d_create_tree_ocp_qp
#define CVT_COLMAJ_TO_TREE_OCP_QP d_cvt_colmaj_to_tree_ocp_qp



#include "x_tree_ocp_qp.c"
