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



#ifndef HPIPM_S_COND_H_
#define HPIPM_S_COND_H_



#include "blasfeo_common.h"

#include "hpipm_s_dense_qp.h"
#include "hpipm_s_dense_qp_dim.h"
#include "hpipm_s_dense_qp_sol.h"
#include "hpipm_s_ocp_qp.h"
#include "hpipm_s_ocp_qp_dim.h"
#include "hpipm_s_ocp_qp_sol.h"


#ifdef __cplusplus
extern "C" {
#endif



struct s_cond_qp_ocp2dense_arg
	{
	int cond_last_stage; // condense last stage
//	int cond_variant; // TODO
	int comp_dual_sol; // dual solution
	int memsize;
	};



struct s_cond_qp_ocp2dense_workspace
	{
	struct blasfeo_smat *Gamma;
	struct blasfeo_smat *L;
	struct blasfeo_smat *Lx;
	struct blasfeo_smat *AL;
	struct blasfeo_svec *Gammab;
	struct blasfeo_svec *l;
	struct blasfeo_svec *tmp_ngM;
	struct blasfeo_svec *tmp_nuxM;
	int *idxs_rev;
	int bs; // block size
	int memsize;
	};



//
int s_memsize_cond_qp_ocp2dense_arg();
//
void s_create_cond_qp_ocp2dense_arg(struct s_cond_qp_ocp2dense_arg *cond_arg, void *mem);
//
void s_set_default_cond_qp_ocp2dense_arg(struct s_cond_qp_ocp2dense_arg *cond_arg);

//
void s_compute_qp_dim_ocp2dense(struct s_ocp_qp_dim *ocp_dim, struct s_dense_qp_dim *dense_dim);
//
int s_memsize_cond_qp_ocp2dense(struct s_ocp_qp_dim *ocp_dim, struct s_cond_qp_ocp2dense_arg *cond_arg);
//
void s_create_cond_qp_ocp2dense(struct s_ocp_qp_dim *ocp_dim, struct s_cond_qp_ocp2dense_arg *cond_arg, struct s_cond_qp_ocp2dense_workspace *cond_ws, void *mem);
//
void s_cond_qp_ocp2dense(struct s_ocp_qp *ocp_qp, struct s_dense_qp *dense_qp, struct s_cond_qp_ocp2dense_arg *cond_arg, struct s_cond_qp_ocp2dense_workspace *cond_ws);
//
void s_cond_rhs_qp_ocp2dense(struct s_ocp_qp *ocp_qp, struct s_dense_qp *dense_qp, struct s_cond_qp_ocp2dense_arg *cond_arg, struct s_cond_qp_ocp2dense_workspace *cond_ws);
//
void s_expand_sol_dense2ocp(struct s_ocp_qp *ocp_qp, struct s_dense_qp_sol *dense_qp_sol, struct s_ocp_qp_sol *ocp_qp_sol, struct s_cond_qp_ocp2dense_arg *cond_arg, struct s_cond_qp_ocp2dense_workspace *cond_ws);
// TODO remove
void s_expand_primal_sol_dense2ocp(struct s_ocp_qp *ocp_qp, struct s_dense_qp_sol *dense_qp_sol, struct s_ocp_qp_sol *ocp_qp_sol, struct s_cond_qp_ocp2dense_arg *cond_arg, struct s_cond_qp_ocp2dense_workspace *cond_ws);

//
void s_update_cond_qp_ocp2dense(int *idxc, struct s_ocp_qp *ocp_qp, struct s_dense_qp *dense_qp, struct s_cond_qp_ocp2dense_arg *cond_arg, struct s_cond_qp_ocp2dense_workspace *cond_ws);


#ifdef __cplusplus
} /* extern "C" */
#endif



#endif // HPIPM_S_COND_H_
