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
	int nv; // number of variables
	int ne; // number of equality constraints
	int nb; // number of box constraints
	int *idxb; // index of box constraints
	int ng; // number of general constraints
	struct s_strmat sA;
	struct s_strvec sb;
	struct s_strmat sQ;
	struct s_strvec sq;
	struct s_strmat sCt;
	struct s_strvec slb;
	struct s_strvec sub;
	struct s_strvec slg;
	struct s_strvec sug;
	};



//
int s_size_dense_qp(int nv, int ne, int nb, int ng);
//
void s_create_dense_qp(int nv, int ne, int nb, int ng, struct s_dense_qp *str_out, void *memory);

