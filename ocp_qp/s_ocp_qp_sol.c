/**************************************************************************************************
*                                                                                                 *
* This file is part of HPIPM.                                                                     *
*                                                                                                 *
* HPIPM -- High-Performance Interior Point Method.                                                *
* Copyright (C) 2017-2018 by Gianluca Frison.                                                     *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
*                                                                                                 *
* This program is free software: you can redistribute it and/or modify                            *
* it under the terms of the GNU General Public License as published by                            *
* the Free Software Foundation, either version 3 of the License, or                               *
* (at your option) any later version                                                              *.
*                                                                                                 *
* This program is distributed in the hope that it will be useful,                                 *
* but WITHOUT ANY WARRANTY; without even the implied warranty of                                  *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                   *
* GNU General Public License for more details.                                                    *
*                                                                                                 *
* You should have received a copy of the GNU General Public License                               *
* along with this program.  If not, see <https://www.gnu.org/licenses/>.                          *
*                                                                                                 *
* The authors designate this particular file as subject to the "Classpath" exception              *
* as provided by the authors in the LICENSE file that accompained this code.                      *
*                                                                                                 *
* Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             *
*                                                                                                 *
**************************************************************************************************/



#include <stdlib.h>
#include <stdio.h>

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_s_aux.h>

#include <hpipm_s_ocp_qp_dim.h>
#include <hpipm_s_ocp_qp.h>
#include <hpipm_s_ocp_qp_sol.h>
#include <hpipm_aux_string.h>



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
#define VECSE blasfeo_svecse

#define OCP_QP_SOL_STRSIZE s_ocp_qp_sol_strsize
#define OCP_QP_SOL_MEMSIZE s_ocp_qp_sol_memsize
#define OCP_QP_SOL_CREATE s_ocp_qp_sol_create
#define OCP_QP_SOL_GET_ALL s_ocp_qp_sol_get_all
#define OCP_QP_SOL_SET_ALL s_ocp_qp_sol_set_all
#define OCP_QP_SOL_GET s_ocp_qp_sol_get
#define OCP_QP_SOL_GET_U s_ocp_qp_sol_get_u
#define OCP_QP_SOL_GET_X s_ocp_qp_sol_get_x
#define OCP_QP_SOL_GET_PI s_ocp_qp_sol_get_pi
#define OCP_QP_SOL_GET_LAM_LB s_ocp_qp_sol_get_lam_lb
#define OCP_QP_SOL_GET_LAM_UB s_ocp_qp_sol_get_lam_ub
#define OCP_QP_SOL_GET_LAM_LG s_ocp_qp_sol_get_lam_lg
#define OCP_QP_SOL_GET_LAM_UG s_ocp_qp_sol_get_lam_ug
#define OCP_QP_SOL_SET s_ocp_qp_sol_set
#define OCP_QP_SOL_SET_U s_ocp_qp_sol_set_u
#define OCP_QP_SOL_SET_X s_ocp_qp_sol_set_x
#define OCP_QP_SOL_SET_SL s_ocp_qp_sol_set_sl
#define OCP_QP_SOL_SET_SU s_ocp_qp_sol_set_su


#include "x_ocp_qp_sol.c"

