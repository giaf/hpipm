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
#include <blasfeo_d_blas.h>
#include <blasfeo_d_aux.h>

#include "../include/hpipm_d_ocp_qp.h"
#include "../include/hpipm_d_ocp_qp_sol.h"
#include "../include/hpipm_d_dense_qp.h"
#include "../include/hpipm_d_dense_qp_sol.h"
#include "../include/hpipm_d_cond.h"
#include "../include/hpipm_d_cond_aux.h"



#define COND_DCTD d_cond_DCtd
#define COND_D d_cond_d
#define COND_B d_cond_b
#define COND_BABT d_cond_BAbt
#define COND_RQ_N2NX3 d_cond_rq_N2nx3
#define COND_RSQRQ_N2NX3 d_cond_RSQrq_N2nx3
#define COND_QP_OCP2DENSE_WORKSPACE d_cond_qp_ocp2dense_workspace
#define CREATE_STRMAT d_create_strmat
#define CREATE_STRVEC d_create_strvec
#define DENSE_QP d_dense_qp
#define DENSE_QP_SOL d_dense_qp_sol
#define EXPAND_SOL d_expand_sol
#define OCP_QP d_ocp_qp
#define OCP_QP_SOL d_ocp_qp_sol
#define SIZE_STRMAT d_size_strmat
#define SIZE_STRVEC d_size_strvec
#define STRMAT d_strmat
#define STRVEC d_strvec

#define COMPUTE_QP_SIZE_OCP2DENSE d_compute_qp_size_ocp2dense
#define MEMSIZE_COND_QP_OCP2DENSE d_memsize_cond_qp_ocp2dense
#define CREATE_COND_QP_OCP2DENSE d_create_cond_qp_ocp2dense
#define COND_QP_OCP2DENSE d_cond_qp_ocp2dense
#define COND_RHS_QP_OCP2DENSE d_cond_rhs_qp_ocp2dense
#define EXPAND_SOL_DENSE2OCP d_expand_sol_dense2ocp



#include "x_cond.c"
