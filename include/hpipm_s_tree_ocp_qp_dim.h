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

#ifndef HPIPM_S_TREE_OCP_QP_DIM_H_
#define HPIPM_S_TREE_OCP_QP_DIM_H_



#ifdef __cplusplus
extern "C" {
#endif



struct s_tree_ocp_qp_dim
	{
	struct tree *ttree; // tree describing node conndection
	int *nx; // number of states // Nn
	int *nu; // number of inputs // Nn
	int *nb; // number of box constraints // Nn
	int *nbx; // number of state box constraints // Nn
	int *nbu; // number of input box constraints // Nn
	int *ng; // number of general constraints // Nn
	int *ns; // number of soft constraints // Nn
	int Nn; // number of nodes
	int memsize;
	};



//
int s_memsize_tree_ocp_qp_dim(int Nn);
//
void s_create_tree_ocp_qp_dim(int Nn, struct s_tree_ocp_qp_dim *qp_dim, void *memory);
//
void s_cvt_int_to_tree_ocp_qp_dim(struct tree *ttree, int *nx, int *nu, int *nbx, int *nbu, int *ng, int *ns, struct s_tree_ocp_qp_dim *dim);



#ifdef __cplusplus
}	// #extern "C"
#endif



#endif // HPIPM_S_TREE_OCP_QP_DIM_H_
