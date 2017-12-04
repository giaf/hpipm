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

#ifndef HPIPM_D_TREE_OCP_QP_H_
#define HPIPM_D_TREE_OCP_QP_H_



#include <blasfeo_target.h>
#include <blasfeo_common.h>

#include "hpipm_d_tree_ocp_qp_dim.h"



#ifdef __cplusplus
extern "C" {
#endif



struct d_tree_ocp_qp
	{
	struct d_tree_ocp_qp_dim *dim;
	struct d_strmat *BAbt; // Nn-1
	struct d_strvec *b; // Nn-1
	struct d_strmat *RSQrq; // Nn
	struct d_strvec *rq; // Nn
	struct d_strmat *DCt; // Nn
	struct d_strvec *d; // Nn
	struct d_strvec *Z; // Nn
	struct d_strvec *z; // Nn
	int **idxb; // index of box constraints // Nn
	int **idxs; // index of soft constraints
	int memsize; // memory size in bytes
	};



//
int d_memsize_tree_ocp_qp(struct d_tree_ocp_qp_dim *dim);
//
void d_create_tree_ocp_qp(struct d_tree_ocp_qp_dim *dim, struct d_tree_ocp_qp *qp, void *memory);
//
void d_cvt_colmaj_to_tree_ocp_qp(double **A, double **B, double **b, double **Q, double **S, double **R, double **q, double **r, int **idxb, double **d_lb, double **d_ub, double **C, double **D, double **d_lg, double **d_ug, double **Zl, double **Zu, double **zl, double **zu, int **idxs, struct d_tree_ocp_qp *qp);



#ifdef __cplusplus
} /* extern "C" */
#endif



#endif // HPIPM_D_TREE_OCP_QP_H_
