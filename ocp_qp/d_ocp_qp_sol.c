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
#include <blasfeo_d_aux.h>

#include <hpipm_d_ocp_qp_dim.h>
#include <hpipm_d_ocp_qp.h>
#include <hpipm_d_ocp_qp_sol.h>
#include <hpipm_aux_string.h>



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
#define VECSE blasfeo_dvecse

#define OCP_QP_SOL_STRSIZE d_ocp_qp_sol_strsize
#define OCP_QP_SOL_MEMSIZE d_ocp_qp_sol_memsize
#define OCP_QP_SOL_CREATE d_ocp_qp_sol_create
#define OCP_QP_SOL_GET_ALL d_ocp_qp_sol_get_all
#define OCP_QP_SOL_SET_ALL d_ocp_qp_sol_set_all
#define OCP_QP_SOL_GET d_ocp_qp_sol_get
#define OCP_QP_SOL_GET_U d_ocp_qp_sol_get_u
#define OCP_QP_SOL_GET_X d_ocp_qp_sol_get_x
#define OCP_QP_SOL_GET_PI d_ocp_qp_sol_get_pi
#define OCP_QP_SOL_GET_LAM_LB d_ocp_qp_sol_get_lam_lb
#define OCP_QP_SOL_GET_LAM_UB d_ocp_qp_sol_get_lam_ub
#define OCP_QP_SOL_GET_LAM_LG d_ocp_qp_sol_get_lam_lg
#define OCP_QP_SOL_GET_LAM_UG d_ocp_qp_sol_get_lam_ug
#define OCP_QP_SOL_SET d_ocp_qp_sol_set
#define OCP_QP_SOL_SET_U d_ocp_qp_sol_set_u
#define OCP_QP_SOL_SET_X d_ocp_qp_sol_set_x
#define OCP_QP_SOL_SET_SL d_ocp_qp_sol_set_sl
#define OCP_QP_SOL_SET_SU d_ocp_qp_sol_set_su


#include "x_ocp_qp_sol.c"
