/**************************************************************************************************
*                                                                                                 *
* This file is part of HPIPM.                                                                     *
*                                                                                                 *
* HPIPM -- High-Performance Interior Point Method.                                                *
* Copyright (C) 2017-2018 by Gianluca Frison.                                                     *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
*                                                                                                 *
* This program is free software: you can redistribute it and/or modify                            *
* it under the terms of the GNU General Public License as published by                            *
* the Free Software Foundation, either version 3 of the License, or                               *
* (at your option) any later version                                                              *.
*                                                                                                 *
* This program is distributed in the hope that it will be useful,                                 *
* but WITHOUT ANY WARRANTY; without even the implied warranty of                                  *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                   *
* GNU General Public License for more details.                                                    *
*                                                                                                 *
* You should have received a copy of the GNU General Public License                               *
* along with this program.  If not, see <https://www.gnu.org/licenses/>.                          *
*                                                                                                 *
* The authors designate this particular file as subject to the "Classpath" exception              *
* as provided by the authors in the LICENSE file that accompained this code.                      *
*                                                                                                 *
* Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             *
*                                                                                                 *
**************************************************************************************************/



#ifndef HPIPM_D_COND_H_
#define HPIPM_D_COND_H_



#include <blasfeo_target.h>
#include <blasfeo_common.h>

#include "hpipm_d_dense_qp.h"
#include "hpipm_d_dense_qp_sol.h"
#include "hpipm_d_ocp_qp.h"
#include "hpipm_d_ocp_qp_dim.h"
#include "hpipm_d_ocp_qp_sol.h"

#ifdef __cplusplus
extern "C" {
#endif



struct d_cond_qp_ocp2dense_arg
	{
	int cond_last_stage; // condense last stage
//	int cond_variant; // TODO
	int comp_dual_sol; // dual solution
	int square_root_alg; // square root algorithm (faster but requires RSQ>0)
	int memsize;
	};



struct d_cond_qp_ocp2dense_workspace
	{
	struct blasfeo_dmat *Gamma;
	struct blasfeo_dmat *L;
	struct blasfeo_dmat *Lx;
	struct blasfeo_dmat *AL;
	struct blasfeo_dvec *Gammab;
	struct blasfeo_dvec *l;
	struct blasfeo_dvec *tmp_ngM;
	struct blasfeo_dvec *tmp_nuxM;
	int *idxs_rev;
	int bs; // block size
	int memsize;
	};



//
int d_memsize_cond_qp_ocp2dense_arg();
//
void d_create_cond_qp_ocp2dense_arg(struct d_cond_qp_ocp2dense_arg *cond_arg, void *mem);
//
void d_set_default_cond_qp_ocp2dense_arg(struct d_cond_qp_ocp2dense_arg *cond_arg);
// set riccati-like algorithm: 0 classical, 1 square-root
void d_set_cond_qp_ocp2dense_arg_ric_alg(int ric_alg, struct d_cond_qp_ocp2dense_arg *cond_arg);

//
void d_compute_qp_dim_ocp2dense(struct d_ocp_qp_dim *ocp_dim, struct d_dense_qp_dim *dense_dim);
//
int d_memsize_cond_qp_ocp2dense(struct d_ocp_qp_dim *ocp_dim, struct d_cond_qp_ocp2dense_arg *cond_arg);
//
void d_create_cond_qp_ocp2dense(struct d_ocp_qp_dim *ocp_dim, struct d_cond_qp_ocp2dense_arg *cond_arg, struct d_cond_qp_ocp2dense_workspace *cond_ws, void *mem);
//
void d_cond_qp_ocp2dense(struct d_ocp_qp *ocp_qp, struct d_dense_qp *dense_qp, struct d_cond_qp_ocp2dense_arg *cond_arg, struct d_cond_qp_ocp2dense_workspace *cond_ws);
//
void d_cond_rhs_qp_ocp2dense(struct d_ocp_qp *ocp_qp, struct d_dense_qp *dense_qp, struct d_cond_qp_ocp2dense_arg *cond_arg, struct d_cond_qp_ocp2dense_workspace *cond_ws);
//
void d_expand_sol_dense2ocp(struct d_ocp_qp *ocp_qp, struct d_dense_qp_sol *dense_qp_sol, struct d_ocp_qp_sol *ocp_qp_sol, struct d_cond_qp_ocp2dense_arg *cond_arg, struct d_cond_qp_ocp2dense_workspace *cond_ws);
// TODO remove
void d_expand_primal_sol_dense2ocp(struct d_ocp_qp *ocp_qp, struct d_dense_qp_sol *dense_qp_sol, struct d_ocp_qp_sol *ocp_qp_sol, struct d_cond_qp_ocp2dense_arg *cond_arg, struct d_cond_qp_ocp2dense_workspace *cond_ws);

//
void d_update_cond_qp_ocp2dense(int *idxc, struct d_ocp_qp *ocp_qp, struct d_dense_qp *dense_qp, struct d_cond_qp_ocp2dense_arg *cond_arg, struct d_cond_qp_ocp2dense_workspace *cond_ws);


#ifdef __cplusplus
} /* extern "C" */
#endif



#endif // HPIPM_D_COND_H_
