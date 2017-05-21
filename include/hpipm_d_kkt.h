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



struct d_qp_dim
	{
	int nv; // number of variables
	int ne; // number of equality constraints
	int nb; // number of box constraints
	int nc; // number of general constraints
	};

struct d_qp_vec
	{
	struct d_strvec *g; // gradient
	struct d_strvec *be; // dynamics vector
	struct d_strvec *lb; // lower bound
	struct d_strvec *ub; // upper bound
	struct d_strvec *lc; // lower constraint
	struct d_strvec *uc; // upper constraint
	int *idxb; // index of box constraints
	};

struct d_qp_mat
	{
	};

struct d_qp
	{
	struct d_qp_dim *dim;
	struct d_qp_vec *vec;
	struct d_qp_mat *mat;
	void (*d_compute_Ct) (void *qp_dim, void *qp_mat, struct d_strvec *v, struct d_strvec *Ctv); // computes Ct * v
	};

