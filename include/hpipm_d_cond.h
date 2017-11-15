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



#ifndef HPIPM_D_COND_H_
#define HPIPM_D_COND_H_



#include <blasfeo_target.h>
#include <blasfeo_common.h>



#ifdef __cplusplus
extern "C" {
#endif



struct d_cond_qp_ocp2dense_workspace
	{
	struct d_strmat *Gamma;
	struct d_strmat *L;
	struct d_strmat *Lx;
	struct d_strmat *AL;
	struct d_strvec *Gammab;
	struct d_strvec *l;
	struct d_strvec *tmp_ngM;
	struct d_strvec *tmp_nuxM;
	int *idxs_rev;
	int cond_last_stage;
	int memsize;
	};



//
void d_compute_qp_dim_ocp2dense(struct d_ocp_qp_dim *ocp_dim, struct d_dense_qp_dim *dense_dim);
//
int d_memsize_cond_qp_ocp2dense(struct d_ocp_qp_dim *ocp_dim); // XXX + args for algorithm type ???
//
void d_create_cond_qp_ocp2dense(struct d_ocp_qp_dim *ocp_dim, struct d_cond_qp_ocp2dense_workspace *cond_ws, void *mem);
//
void d_cond_qp_ocp2dense(struct d_ocp_qp *ocp_qp, struct d_dense_qp *dense_qp, struct d_cond_qp_ocp2dense_workspace *cond_ws);
//
void d_cond_rhs_qp_ocp2dense(struct d_ocp_qp *ocp_qp, struct d_dense_qp *dense_qp, struct d_cond_qp_ocp2dense_workspace *cond_ws);
//
void d_expand_sol_dense2ocp(struct d_ocp_qp *ocp_qp, struct d_dense_qp_sol *dense_qp_sol, struct d_ocp_qp_sol *ocp_qp_sol, struct d_cond_qp_ocp2dense_workspace *cond_ws);



#ifdef __cplusplus
} /* extern "C" */
#endif



#endif // HPIPM_D_COND_H_
