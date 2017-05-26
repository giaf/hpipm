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



struct d_dense_qp
	{
	struct d_strmat *H; // hessian
	struct d_strmat *A; // dynamics matrix
	struct d_strmat *Ct; // constraints matrix
	struct d_strvec *g; // gradient
	struct d_strvec *be; // dynamics vector
	struct d_strvec *lb; // lower bound
	struct d_strvec *ub; // upper bound
	struct d_strvec *lg; // lower constraint
	struct d_strvec *ug; // upper constraint
	int *idxb; // index of box constraints
	int nv; // number of variables
	int ne; // number of equality constraints
	int nb; // number of box constraints
	int ng; // number of general constraints
	int mem_size; // memory size in bytes
	};



//
int d_size_dense_qp(int nv, int ne, int nb, int ng);
//
void d_create_dense_qp(int nv, int ne, int nb, int ng, struct d_dense_qp *qp, void *memory);
//
void d_init_dense_qp(struct d_strmat *H, struct d_strmat *A, struct d_strmat *C, struct d_strvec *g, struct d_strvec *be, struct d_strvec *lb, struct d_strvec *ub, struct d_strvec *lg, struct d_strvec *ug, int *idxb, struct d_dense_qp *qp);
//

