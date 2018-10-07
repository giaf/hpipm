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



#include <hpipm_d_ocp_qp_dim.h>
#include <hpipm_aux_string.h>



#define OCP_QP_DIM d_ocp_qp_dim

#define SIZEOF_OCP_QP_DIM d_sizeof_ocp_qp_dim
#define MEMSIZE_OCP_QP_DIM d_memsize_ocp_qp_dim
#define CREATE_OCP_QP_DIM d_create_ocp_qp_dim
#define CVT_INT_TO_OCP_QP_DIM d_cvt_int_to_ocp_qp_dim
#define SET_OCP_QP_DIM_GF d_set_ocp_qp_dim_gf
#define SET_OCP_QP_DIM_NX d_set_ocp_qp_dim_nx
#define SET_OCP_QP_DIM_NU d_set_ocp_qp_dim_nu
#define SET_OCP_QP_DIM_NBX d_set_ocp_qp_dim_nbx
#define SET_OCP_QP_DIM_NBU d_set_ocp_qp_dim_nbu
#define SET_OCP_QP_DIM_NG d_set_ocp_qp_dim_ng
#define SET_OCP_QP_DIM_NS d_set_ocp_qp_dim_ns



#include "x_ocp_qp_dim.c"
