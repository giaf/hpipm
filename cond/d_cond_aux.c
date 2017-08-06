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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_d_blas.h>
#include <blasfeo_d_aux.h>

#include "../include/hpipm_d_ocp_qp.h"
#include "../include/hpipm_d_ocp_qp_sol.h"
#include "../include/hpipm_d_dense_qp.h"
#include "../include/hpipm_d_dense_qp_sol.h"
#include "../include/hpipm_d_cond.h"



#define AXPY_LIBSTR daxpy_libstr
#define COND_QP_OCP2DENSE_WORKSPACE d_cond_qp_ocp2dense_workspace
#define DENSE_QP_SOL d_dense_qp_sol
#define GEAD_LIBSTR dgead_libstr
#define GECP_LIBSTR dgecp_libstr
#define GEEX1_LIBSTR dgeex1_libstr
#define GESE_LIBSTR dgese_libstr
#define GEMM_NN_LIBSTR dgemm_nn_libstr
#define GEMV_T_LIBSTR dgemv_t_libstr
#define GEMV_N_LIBSTR dgemv_n_libstr
#define OCP_QP d_ocp_qp
#define OCP_QP_SOL d_ocp_qp_sol
#define POTRF_L_MN_LIBSTR dpotrf_l_mn_libstr
#define REAL double
#define ROWEX_LIBSTR drowex_libstr
#define STRMAT d_strmat
#define STRVEC d_strvec
#define SYMV_L_LIBSTR dsymv_l_libstr
#define SYRK_LN_MN_LIBSTR dsyrk_ln_mn_libstr
#define TRCP_L_LIBSTR dtrcp_l_libstr
#define TRMM_RLNN_LIBSTR dtrmm_rlnn_libstr
#define VECCP_LIBSTR dveccp_libstr

#define COMPUTE_GAMMA d_compute_Gamma
#define COND_BABT d_cond_BAbt
#define COND_RSQRQ_N2NX3 d_cond_RSQrq_N2nx3
#define COND_DCTD d_cond_DCtd
#define EXPAND_SOL d_expand_sol



#include "x_cond_aux.c"
