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

#ifndef HPIPM_S_DENSE_QP_DIM_H_
#define HPIPM_S_DENSE_QP_DIM_H_



#ifdef __cplusplus
extern "C" {
#endif



struct s_dense_qp_dim
	{
	int nv; // number of variables
	int ne; // number of equality constraints
	int nb; // number of box constraints
	int ng; // number of general constraints
	int ns; // number of softened constraints
	int memsize;
	};



//
int s_memsize_dense_qp_dim();
//
void s_create_dense_qp_dim(struct s_dense_qp_dim *qp_dim, void *memory);
//
void s_cvt_int_to_dense_qp_dim(int nv, int ne, int nb, int ng, int ns, struct s_dense_qp_dim *dim);



#ifdef __cplusplus
}	// #extern "C"
#endif



#endif // HPIPM_S_DENSE_QP_DIM_H_


