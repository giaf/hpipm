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



#ifndef HPIPM_S_TREE_OCP_QP_IPM_H_
#define HPIPM_S_TREE_OCP_QP_IPM_H_



#include <blasfeo_target.h>
#include <blasfeo_common.h>

#include <hpipm_common.h>
#include <hpipm_s_tree_ocp_qp_dim.h>
#include <hpipm_s_tree_ocp_qp.h>
#include <hpipm_s_tree_ocp_qp_res.h>
#include <hpipm_s_tree_ocp_qp_sol.h>



#ifdef __cplusplus
extern "C" {
#endif



struct s_tree_ocp_qp_ipm_arg
	{
	float mu0; // initial value for duality measure
	float alpha_min; // exit cond on step length
	float res_g_max; // exit cond on inf norm of residuals
	float res_b_max; // exit cond on inf norm of residuals
	float res_d_max; // exit cond on inf norm of residuals
	float res_m_max; // exit cond on inf norm of residuals
	float reg_prim; // reg of primal hessian
	float lam_min; // min value in lam vector
	float t_min; // min value in t vector
	int iter_max; // exit cond in iter number
	int stat_max; // iterations saved in stat
	int pred_corr; // use Mehrotra's predictor-corrector IPM algirthm
	int cond_pred_corr; // conditional Mehrotra's predictor-corrector
	int itref_pred_max; // max number of iterative refinement steps for predictor step
	int itref_corr_max; // max number of iterative refinement steps for corrector step
	int warm_start; // 0 no warm start, 1 warm start primal sol, 2 warm start primal and dual sol
	int lq_fact; // 0 syrk+potrf, 1 mix, 2 lq
	int abs_form; // absolute IPM formulation
	int comp_dual_sol; // dual solution (only for abs_form==1)
	int comp_res_exit; // compute residuals on exit (only for abs_form==1 and comp_dual_sol==1)
	int memsize;
	};



struct s_tree_ocp_qp_ipm_workspace
	{
	struct s_core_qp_ipm_workspace *core_workspace;
	struct s_tree_ocp_qp_res_workspace *res_workspace;
	struct s_tree_ocp_qp_sol *sol_step;
	struct s_tree_ocp_qp_sol *sol_itref;
	struct s_tree_ocp_qp *qp_step;
	struct s_tree_ocp_qp *qp_itref;
	struct s_tree_ocp_qp_res *res_itref;
	struct s_tree_ocp_qp_res *res;
	struct blasfeo_svec *Gamma; // hessian update
	struct blasfeo_svec *gamma; // hessian update
	struct blasfeo_svec *tmp_nxM; // work space of size nxM
	struct blasfeo_svec *tmp_nbgM; // work space of size nbgM
	struct blasfeo_svec *tmp_nsM; // work space of size nsM
	struct blasfeo_svec *Pb; // Pb
	struct blasfeo_svec *Zs_inv;
	struct blasfeo_smat *L;
	struct blasfeo_smat *Lh;
	struct blasfeo_smat *AL;
	struct blasfeo_smat *lq0;
	struct blasfeo_svec *tmp_m;
	float *stat; // convergence statistics
	int *use_hess_fact;
	void *lq_work0;
	float qp_res[4]; // infinity norm of residuals
	int iter; // iteration number
	int stat_max; // iterations saved in stat
	int use_Pb;
	int lq_fact; // cache from arg
	int memsize;
	};



//
int s_memsize_tree_ocp_qp_ipm_arg(struct s_tree_ocp_qp_dim *dim);
//
void s_create_tree_ocp_qp_ipm_arg(struct s_tree_ocp_qp_dim *dim, struct s_tree_ocp_qp_ipm_arg *arg, void *mem);
//
void s_set_default_tree_ocp_qp_ipm_arg(enum hpipm_mode mode, struct s_tree_ocp_qp_ipm_arg *arg);
//
void s_set_tree_ocp_qp_ipm_arg_iter_max(int iter_max, struct s_tree_ocp_qp_ipm_arg *arg);
//
void s_set_tree_ocp_qp_ipm_arg_mu0(float mu0, struct s_tree_ocp_qp_ipm_arg *arg);
//
void s_set_tree_ocp_qp_ipm_arg_tol_stat(float tol_stat, struct s_tree_ocp_qp_ipm_arg *arg);
//
void s_set_tree_ocp_qp_ipm_arg_tol_eq(float tol_eq, struct s_tree_ocp_qp_ipm_arg *arg);
//
void s_set_tree_ocp_qp_ipm_arg_tol_ineq(float tol_ineq, struct s_tree_ocp_qp_ipm_arg *arg);
//
void s_set_tree_ocp_qp_ipm_arg_tol_comp(float tol_comp, struct s_tree_ocp_qp_ipm_arg *arg);
//
void s_set_tree_ocp_qp_ipm_arg_reg_prim(float reg, struct s_tree_ocp_qp_ipm_arg *arg);

//
int s_memsize_tree_ocp_qp_ipm(struct s_tree_ocp_qp_dim *dim, struct s_tree_ocp_qp_ipm_arg *arg);
//
void s_create_tree_ocp_qp_ipm(struct s_tree_ocp_qp_dim *dim, struct s_tree_ocp_qp_ipm_arg *arg, struct s_tree_ocp_qp_ipm_workspace *ws, void *mem);
//
int s_get_tree_ocp_qp_ipm_iter(struct s_tree_ocp_qp_ipm_workspace *ws);
//
float s_get_tree_ocp_qp_ipm_res_stat(struct s_tree_ocp_qp_ipm_workspace *ws);
//
float s_get_tree_ocp_qp_ipm_res_eq(struct s_tree_ocp_qp_ipm_workspace *ws);
//
float s_get_tree_ocp_qp_ipm_res_ineq(struct s_tree_ocp_qp_ipm_workspace *ws);
//
float s_get_tree_ocp_qp_ipm_res_comp(struct s_tree_ocp_qp_ipm_workspace *ws);
//
float *s_get_tree_ocp_qp_ipm_stat(struct s_tree_ocp_qp_ipm_workspace *ws);
//
int s_solve_tree_ocp_qp_ipm(struct s_tree_ocp_qp *qp, struct s_tree_ocp_qp_sol *qp_sol, struct s_tree_ocp_qp_ipm_arg *arg, struct s_tree_ocp_qp_ipm_workspace *ws);

#ifdef __cplusplus
} /* extern "C" */
#endif



#endif // HPIPM_S_TREE_OCP_QP_IPM_H_
