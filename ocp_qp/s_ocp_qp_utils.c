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
* Author: Gianluca Frison  {gianluca.frison} (at) imtek.uni-freiburg.de                           *
*       , Andrea Zanelli {andrea.zanelli} (at) imtek.uni-freiburg.de                              *
*                                                                                                 *
**************************************************************************************************/

#include <stdlib.h>
#include <stdio.h>

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_s_aux.h>
#include <blasfeo_s_aux_ext_dep.h>
#include <blasfeo_i_aux_ext_dep.h>

#include <hpipm_s_ocp_qp_dim.h>
#include <hpipm_s_ocp_qp.h>
#include <hpipm_s_ocp_qp_sol.h>



#define BLASFEO_PRINT_MAT blasfeo_print_smat
#define BLASFEO_PRINT_TRAN_VEC blasfeo_print_tran_svec
#define OCP_QP s_ocp_qp
#define OCP_QP_SOL s_ocp_qp_sol
#define OCP_QP_DIM s_ocp_qp_dim



#define PRINT_OCP_QP s_print_ocp_qp
#define PRINT_OCP_QP_SOL s_print_ocp_qp_sol



#include "x_ocp_qp_utils.c"

