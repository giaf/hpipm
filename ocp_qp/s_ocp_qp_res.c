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
#include <hpipm_s_ocp_qp_res.h>



#define CREATE_STRVEC blasfeo_create_svec
#define CVT_STRVEC2VEC blasfeo_unpack_svec
#define OCP_QP_DIM s_ocp_qp_dim
#define OCP_QP_RES s_ocp_qp_res
#define OCP_QP_RES_WS s_ocp_qp_res_ws
#define REAL float
#define SIZE_STRVEC blasfeo_memsize_svec
#define STRVEC blasfeo_svec



#define OCP_QP_RES_MEMSIZE s_ocp_qp_res_memsize
#define OCP_QP_RES_CREATE s_ocp_qp_res_create
#define OCP_QP_RES_WS_MEMSIZE s_ocp_qp_res_ws_memsize
#define OCP_QP_RES_WS_CREATE s_ocp_qp_res_ws_create
#define OCP_QP_RES_GET_ALL s_ocp_qp_res_get_all



#include "x_ocp_qp_res.c"

