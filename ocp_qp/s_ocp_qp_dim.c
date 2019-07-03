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

#include <hpipm_s_ocp_qp_dim.h>
#include <hpipm_aux_string.h>



#define OCP_QP_DIM s_ocp_qp_dim



#define OCP_QP_DIM_STRSIZE s_ocp_qp_dim_strsize
#define OCP_QP_DIM_MEMSIZE s_ocp_qp_dim_memsize
#define OCP_QP_DIM_CREATE s_ocp_qp_dim_create
#define OCP_QP_DIM_SET_ALL s_ocp_qp_dim_set_all
#define OCP_QP_DIM_SET s_ocp_qp_dim_set
#define OCP_QP_DIM_SET_NX s_ocp_qp_dim_set_nx
#define OCP_QP_DIM_SET_NU s_ocp_qp_dim_set_nu
#define OCP_QP_DIM_SET_NBX s_ocp_qp_dim_set_nbx
#define OCP_QP_DIM_SET_NBU s_ocp_qp_dim_set_nbu
#define OCP_QP_DIM_SET_NG s_ocp_qp_dim_set_ng
#define OCP_QP_DIM_SET_NSBX s_ocp_qp_dim_set_nsbx
#define OCP_QP_DIM_SET_NSBU s_ocp_qp_dim_set_nsbu
#define OCP_QP_DIM_SET_NSG s_ocp_qp_dim_set_nsg



#include "x_ocp_qp_dim.c"

