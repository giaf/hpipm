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



struct d_ocp_nlp
	{
	void (**expl_vde)(int t, double *x, double *p, void *vde_args, double *xdot); // explicit vde
	struct d_strmat *RSQ;
	struct d_strvec *ux_ref;
	struct d_strmat *DCt;
	struct d_strvec *d;
	struct d_strvec *Z;
	struct d_strvec *z;
	int *nx; // number of states
	int *nu; // number of inputs
	int *nb; // number of box constraints
	int **idxb; // index of box constraints
	int *ng; // number of general constraints
	int *ns; // number of soft constraints
	int **idxs; // index of soft constraints
	int N; // hotizon lenght
	int memsize; // memory size in bytes
	};



//
int d_memsize_ocp_nlp(int N, int *nx, int *nu, int *nb, int *ng, int *ns);

