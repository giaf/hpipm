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



#include <blasfeo_target.h>
#include <blasfeo_common.h>



struct d_tree_ocp_qp
	{
	struct d_strmat *BAbt; // Nn-1
	struct d_strvec *b; // Nn-1
	struct d_strmat *RSQrq; // Nn
	struct d_strvec *rq; // Nn
	struct d_strmat *DCt; // Nn
	struct d_strvec *d_lb; // Nn
	struct d_strvec *d_ub; // Nn
	struct d_strvec *d_lg; // Nn
	struct d_strvec *d_ug; // Nn
	int *nx; // number of states // Nn
	int *nu; // number of inputs // Nn
	int *nb; // number of box constraints // Nn
	int **idxb; // index of box constraints // Nn
	int *ng; // number of general constraints // Nn
	struct tree *ttree;
	int Nn; // number of nodes
	int memsize; // memory size in bytes
	};



//
int d_memsize_tree_ocp_qp(struct tree *ttree, int *nx, int *nu, int *nb, int *ng);
//
void d_create_tree_ocp_qp(struct tree *ttree, int *nx, int *nu, int *nb, int *ng, struct d_tree_ocp_qp *qp, void *memory);
//
void d_cvt_colmaj_to_tree_ocp_qp(double **A, double **B, double **b, double **Q, double **S, double **R, double **q, double **r, int **idxb, double **d_lb, double **d_ub, double **C, double **D, double **d_lg, double **d_ug, struct d_tree_ocp_qp *qp);
