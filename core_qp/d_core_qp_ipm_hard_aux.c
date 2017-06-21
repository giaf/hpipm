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



#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_d_aux.h>
#include <blasfeo_d_blas.h>

#include "../include/hpipm_d_core_qp_ipm_hard.h"



#define IPM_HARD_CORE_QP_WORKSPACE d_ipm_hard_core_qp_workspace
#define REAL double

#define COMPUTE_QX_QX_HARD_QP d_compute_Qx_qx_hard_qp
#define COMPUTE_LAM_T_HARD_QP d_compute_lam_t_hard_qp
#define COMPUTE_ALPHA_HARD_QP d_compute_alpha_hard_qp
#define UPDATE_VAR_HARD_QP d_update_var_hard_qp
#define COMPUTE_MU_AFF_HARD_QP d_compute_mu_aff_hard_qp
#define COMPUTE_CENTERING_CORRECTION_HARD_QP d_compute_centering_correction_hard_qp
#define COMPUTE_QX_HARD_QP d_compute_qx_hard_qp



#include "x_core_qp_ipm_hard_aux.c"
