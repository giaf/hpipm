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

#ifndef HPIPM_S_TREE_OCP_QP_SOL_H_
#define HPIPM_S_TREE_OCP_QP_SOL_H_



#include <blasfeo_target.h>
#include <blasfeo_common.h>

#include "hpipm_s_tree_ocp_qp_dim.h"



#ifdef __cplusplus
extern "C" {
#endif

struct s_tree_ocp_qp_sol
	{
	struct s_tree_ocp_qp_dim *dim;
	struct blasfeo_svec *ux;
	struct blasfeo_svec *pi;
	struct blasfeo_svec *lam;
	struct blasfeo_svec *t;
	void *misc;
	int memsize; // memory size in bytes
	};



//
int s_memsize_tree_ocp_qp_sol(struct s_tree_ocp_qp_dim *dim);
//
void s_create_tree_ocp_qp_sol(struct s_tree_ocp_qp_dim *dim, struct s_tree_ocp_qp_sol *qp_sol, void *memory);
//
void s_cvt_tree_ocp_qp_sol_to_colmaj(struct s_tree_ocp_qp *qp, struct s_tree_ocp_qp_sol *qp_sol, float **u, float **x, float **ls, float **us, float **pi, float **lam_lb, float **lam_ub, float **lam_lg, float **lam_ug, float **lam_ls, float **lam_us);
//
void s_cvt_tree_ocp_qp_sol_to_rowmaj(struct s_tree_ocp_qp *qp, struct s_tree_ocp_qp_sol *qp_sol, float **u, float **x, float **ls, float **us, float **pi, float **lam_lb, float **lam_ub, float **lam_lg, float **lam_ug, float **lam_ls, float **lam_us);

#ifdef __cplusplus
} /* extern "C" */
#endif



#endif // HPIPM_S_TREE_OCP_QP_SOL_H_
