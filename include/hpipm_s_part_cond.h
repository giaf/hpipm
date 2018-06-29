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



#ifndef HPIPM_S_PART_COND_H_
#define HPIPM_S_PART_COND_H_



#include <blasfeo_target.h>
#include <blasfeo_common.h>



#ifdef __cplusplus
extern "C" {
#endif



struct s_cond_qp_ocp2ocp_arg
	{
	struct s_cond_qp_ocp2dense_arg *cond_arg;
	int memsize;
	};



struct s_cond_qp_ocp2ocp_workspace
	{
	struct s_cond_qp_ocp2dense_workspace *cond_workspace;
	int memsize;
	};



//
int s_memsize_cond_qp_ocp2ocp_arg(int N2);
//
void s_create_cond_qp_ocp2ocp_arg(int N2, struct s_cond_qp_ocp2ocp_arg *cond_arg, void *mem);
//
void s_set_default_cond_qp_ocp2ocp_arg(int N2, struct s_cond_qp_ocp2ocp_arg *cond_arg);

//
void s_compute_block_size_cond_qp_ocp2ocp(int N, int N2, int *block_size);
//
void s_compute_qp_dim_ocp2ocp(struct s_ocp_qp_dim *ocp_dim, int *block_size, struct s_ocp_qp_dim *part_dense_dim);
//
int s_memsize_cond_qp_ocp2ocp(struct s_ocp_qp_dim *ocp_dim, int *block_size, struct s_ocp_qp_dim *part_dense_dim, struct s_cond_qp_ocp2ocp_arg *cond_arg);
//
void s_create_cond_qp_ocp2ocp(struct s_ocp_qp_dim *ocp_dim, int *block_size, struct s_ocp_qp_dim *part_dense_dim, struct s_cond_qp_ocp2ocp_arg *cond_arg, struct s_cond_qp_ocp2ocp_workspace *cond_ws, void *mem);
//
void s_cond_qp_ocp2ocp(struct s_ocp_qp *ocp_qp, struct s_ocp_qp *part_dense_qp, struct s_cond_qp_ocp2ocp_arg *cond_arg, struct s_cond_qp_ocp2ocp_workspace *cond_ws);
//
void s_cond_rhs_qp_ocp2ocp(struct s_ocp_qp *ocp_qp, struct s_ocp_qp *part_dense_qp, struct s_cond_qp_ocp2ocp_arg *cond_arg, struct s_cond_qp_ocp2ocp_workspace *cond_ws);
//
void s_expand_sol_ocp2ocp(struct s_ocp_qp *ocp_qp, struct s_ocp_qp *part_dense_qp, struct s_ocp_qp_sol *part_dense_qp_sol, struct s_ocp_qp_sol *ocp_qp_sol, struct s_cond_qp_ocp2ocp_arg *cond_arg, struct s_cond_qp_ocp2ocp_workspace *cond_ws);

//
void s_update_cond_qp_ocp2ocp(int *idxc, struct s_ocp_qp *ocp_qp, struct s_ocp_qp *part_dense_qp, struct s_cond_qp_ocp2ocp_arg *cond_arg, struct s_cond_qp_ocp2ocp_workspace *cond_ws);


#ifdef __cplusplus
} /* extern "C" */
#endif



#endif // HPIPM_S_PART_COND_H_
