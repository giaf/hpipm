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
* Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             *
*                                                                                                 *
**************************************************************************************************/



#include <blasfeo_target.h>
#include <blasfeo_common.h>



struct s_dense_qp
	{
	struct s_strmat *H; // hessian
	struct s_strmat *A; // dynamics matrix
	struct s_strmat *Ct; // constraints matrix
	struct s_strvec *g; // gradient
	struct s_strvec *be; // dynamics vector
	struct s_strvec *lb; // lower bound
	struct s_strvec *ub; // upper bound
	struct s_strvec *lg; // lower constraint
	struct s_strvec *ug; // upper constraint
	int *idxb; // index of box constraints
	int nv; // number of variables
	int ne; // number of equality constraints
	int nb; // number of box constraints
	int ng; // number of general constraints
	int mem_size; // memory size in bytes
	};



//
int s_size_dense_qp(int nv, int ne, int nb, int ng);
//
void s_create_dense_qp(int nv, int ne, int nb, int ng, struct s_dense_qp *qp, void *memory);
//
void s_init_dense_qp(struct s_strmat *H, struct s_strmat *A, struct s_strmat *C, struct s_strvec *g, struct s_strvec *be, struct s_strvec *lb, struct s_strvec *ub, struct s_strvec *lg, struct s_strvec *ug, int *idxb, struct s_dense_qp *qp);
//

