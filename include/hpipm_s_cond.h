/**************************************************************************************************
*                                                                                                 *
* This file is part of HPIPM.                                                                     *
*                                                                                                 *
* HPIPM -- High-Performance Interior Point Method.                                                *
* Copyright (C) 2019 by Gianluca Frison.                                                          *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
*                                                                                                 *
* The 2-Clause BSD License                                                                        *
*                                                                                                 *
* Redistribution and use in source and binary forms, with or without                              *
* modification, are permitted provided that the following conditions are met:                     *
*                                                                                                 *
* 1. Redistributions of source code must retain the above copyright notice, this                  *
*    list of conditions and the following disclaimer.                                             *
* 2. Redistributions in binary form must reproduce the above copyright notice,                    *
*    this list of conditions and the following disclaimer in the documentation                    *
*    and/or other materials provided with the distribution.                                       *
*                                                                                                 *
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND                 *
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED                   *
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE                          *
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR                 *
* ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES                  *
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;                    *
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND                     *
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT                      *
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS                   *
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                                    *
*                                                                                                 *
* Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             *
*                                                                                                 *
**************************************************************************************************/



#ifndef HPIPM_S_COND_H_
#define HPIPM_S_COND_H_



#include <blasfeo_target.h>
#include <blasfeo_common.h>



#ifdef __cplusplus
extern "C" {
#endif



struct s_cond_qp_ocp2dense_arg
	{
	int cond_last_stage; // condense last stage
//	int cond_variant; // TODO
	int comp_dual_sol; // dual solution
	int square_root_alg; // square root algorithm (faster but requires RSQ>0)
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
	struct blasfeo_svec *tmp_nbgM;
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
// set riccati-like algorithm: 0 classical, 1 square-root
void s_set_cond_qp_ocp2dense_arg_ric_alg(int ric_alg, struct s_cond_qp_ocp2dense_arg *cond_arg);

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
