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

#ifndef HPIPM_D_TREE_OCP_QP_SOL_H_
#define HPIPM_D_TREE_OCP_QP_SOL_H_



#include <blasfeo_target.h>
#include <blasfeo_common.h>

#include "hpipm_d_tree_ocp_qp_dim.h"



#ifdef __cplusplus
extern "C" {
#endif

struct d_tree_ocp_qp_sol
	{
	struct d_tree_ocp_qp_dim *dim;
	struct blasfeo_dvec *ux;
	struct blasfeo_dvec *pi;
	struct blasfeo_dvec *lam;
	struct blasfeo_dvec *t;
	void *misc;
	int memsize; // memory size in bytes
	};



//
int d_memsize_tree_ocp_qp_sol(struct d_tree_ocp_qp_dim *dim);
//
void d_create_tree_ocp_qp_sol(struct d_tree_ocp_qp_dim *dim, struct d_tree_ocp_qp_sol *qp_sol, void *memory);
//
void d_cvt_tree_ocp_qp_sol_to_colmaj(struct d_tree_ocp_qp *qp, struct d_tree_ocp_qp_sol *qp_sol, double **u, double **x, double **ls, double **us, double **pi, double **lam_lb, double **lam_ub, double **lam_lg, double **lam_ug, double **lam_ls, double **lam_us);
//
void d_cvt_tree_ocp_qp_sol_to_rowmaj(struct d_tree_ocp_qp *qp, struct d_tree_ocp_qp_sol *qp_sol, double **u, double **x, double **ls, double **us, double **pi, double **lam_lb, double **lam_ub, double **lam_lg, double **lam_ug, double **lam_ls, double **lam_us);

#ifdef __cplusplus
} /* extern "C" */
#endif



#endif // HPIPM_D_TREE_OCP_QP_SOL_H_
