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



#ifndef HPIPM_D_PART_COND_H_
#define HPIPM_D_PART_COND_H_



#include <blasfeo_target.h>
#include <blasfeo_common.h>



#ifdef __cplusplus
extern "C" {
#endif



struct d_cond_qp_ocp2ocp_workspace
	{
	struct d_cond_qp_ocp2dense_workspace *cond_workspace;
	int memsize;
	};



//
void d_compute_block_size_cond_qp_ocp2ocp(int N, int N2, int *block_size);
//
void d_compute_qp_dim_ocp2ocp(struct d_ocp_qp_dim *ocp_dim, int *block_size, struct d_ocp_qp_dim *part_dense_dim);
//
int d_memsize_cond_qp_ocp2ocp(struct d_ocp_qp_dim *ocp_dim, int *block_size, struct d_ocp_qp_dim *part_dense_dim);
//
void d_create_cond_qp_ocp2ocp(struct d_ocp_qp_dim *ocp_dim, int *block_size, struct d_ocp_qp_dim *part_dense_dim, struct d_cond_qp_ocp2ocp_workspace *cond_ws, void *mem);
//
void d_cond_qp_ocp2ocp(struct d_ocp_qp *ocp_qp, struct d_ocp_qp *part_dense_qp, struct d_cond_qp_ocp2ocp_workspace *cond_ws);
//
void d_cond_rhs_qp_ocp2ocp(struct d_ocp_qp *ocp_qp, struct d_ocp_qp *part_dense_qp, struct d_cond_qp_ocp2ocp_workspace *cond_ws);
//
void d_expand_sol_ocp2ocp(struct d_ocp_qp *ocp_qp, struct d_ocp_qp *part_dense_qp, struct d_ocp_qp_sol *part_dense_qp_sol, struct d_ocp_qp_sol *ocp_qp_sol, struct d_cond_qp_ocp2ocp_workspace *cond_ws);

//
void d_update_cond_qp_ocp2ocp(int *idxc, struct d_ocp_qp *ocp_qp, struct d_ocp_qp *part_dense_qp, struct d_cond_qp_ocp2ocp_workspace *cond_ws);


#ifdef __cplusplus
} /* extern "C" */
#endif



#endif // HPIPM_D_PART_COND_H_
