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



struct s_tree_ocp_qp
	{
	struct s_strmat *BAbt; // Nn-1
	struct s_strvec *b; // Nn-1
	struct s_strmat *RSQrq; // Nn
	struct s_strvec *rq; // Nn
	struct s_strmat *DCt; // Nn
	struct s_strvec *d_lb; // Nn
	struct s_strvec *d_ub; // Nn
	struct s_strvec *d_lg; // Nn
	struct s_strvec *d_ug; // Nn
	struct tree *ttree; // tree describing node conndection
	int *nx; // number of states // Nn
	int *nu; // number of inputs // Nn
	int *nb; // number of box constraints // Nn
	int **idxb; // index of box constraints // Nn
	int *ng; // number of general constraints // Nn
	int Nn; // number of nodes
	int memsize; // memory size in bytes
	};



//
int s_memsize_tree_ocp_qp(struct tree *ttree, int *nx, int *nu, int *nb, int *ng);
//
void s_create_tree_ocp_qp(struct tree *ttree, int *nx, int *nu, int *nb, int *ng, struct s_tree_ocp_qp *qp, void *memory);
//
void s_cvt_colmaj_to_tree_ocp_qp(float **A, float **B, float **b, float **Q, float **S, float **R, float **q, float **r, int **idxb, float **d_lb, float **d_ub, float **C, float **D, float **d_lg, float **d_ug, struct s_tree_ocp_qp *qp);

