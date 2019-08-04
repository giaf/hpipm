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

#ifndef HPIPM_S_OCP_QP_UTILS_H_
#define HPIPM_S_OCP_QP_UTILS_H_



#include <blasfeo_target.h>
#include <blasfeo_common.h>

#include "hpipm_s_ocp_qp_dim.h"
#include "hpipm_s_ocp_qp.h"
#include "hpipm_s_ocp_qp_sol.h"
#include "hpipm_s_ocp_qp_ipm.h"



#ifdef __cplusplus
extern "C" {
#endif



//
void s_ocp_qp_dim_print(struct s_ocp_qp_dim *qp_dim);
//
void s_ocp_qp_dim_codegen(char *file_name, char *mode, char *target, struct s_ocp_qp_dim *qp_dim);
//
void s_ocp_qp_print(struct d_ocp_qp_dim *qp_dim, struct s_ocp_qp *qp);
//
void s_ocp_qp_codegen(char *file_name, char *mode, struct s_ocp_qp_dim *qp_dim, struct s_ocp_qp *qp);
//
void s_ocp_qp_sol_print(struct s_ocp_qp_dim *qp_dim, struct s_ocp_qp_sol *ocp_qp_sol);
//
void s_ocp_qp_ipm_arg_codegen(char *file_name, char *mode, struct s_ocp_qp_dim *qp_dim, struct s_ocp_qp_ipm_arg *arg);



#ifdef __cplusplus
}	// #extern "C"
#endif



#endif // HPIPM_D_OCP_QP_UTILS_H_

