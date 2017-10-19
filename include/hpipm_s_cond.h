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

#ifdef __cplusplus
extern "C" {
#endif

struct s_cond_qp_ocp2dense_workspace
	{
	struct s_strmat *Gamma;
	struct s_strmat *L;
	struct s_strmat *Lx;
	struct s_strmat *AL;
	struct s_strvec *Gammab;
	struct s_strvec *tmp_ngM;
	struct s_strvec *tmp_nuxM;
	int *idxs_rev;
	int cond_last_stage;
	int memsize;
	};



//
void s_compute_qp_size_ocp2dense(int N, int *nx, int *nu, int *nb, int **idxb, int *ng, int *ns, int *nvc, int *nec, int *nbc, int *ngc, int *nsc);
//
int s_memsize_cond_qp_ocp2dense(struct s_ocp_qp *ocp_qp, struct s_dense_qp *dense_qp); // XXX + args for algorithm type ???
//
void s_create_cond_qp_ocp2dense(struct s_ocp_qp *ocp_qp, struct s_dense_qp *dense_qp, struct s_cond_qp_ocp2dense_workspace *cond_ws, void *mem);
//
void s_cond_qp_ocp2dense(struct s_ocp_qp *ocp_qp, struct s_dense_qp *dense_qp, struct s_cond_qp_ocp2dense_workspace *cond_ws);
//
void s_expand_sol_dense2ocp(struct s_ocp_qp *ocp_qp, struct s_dense_qp_sol *dense_qp_sol, struct s_ocp_qp_sol *ocp_qp_sol, struct s_cond_qp_ocp2dense_workspace *cond_ws);

#ifdef __cplusplus
} /* extern "C" */
#endif
