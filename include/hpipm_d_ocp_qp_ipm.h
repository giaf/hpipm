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

#ifndef HPIPM_D_OCP_QP_IPM_H_
#define HPIPM_D_OCP_QP_IPM_H_



#include <blasfeo_target.h>
#include <blasfeo_common.h>

#include <hpipm_common.h>
#include <hpipm_d_ocp_qp_dim.h>
#include <hpipm_d_ocp_qp.h>
#include <hpipm_d_ocp_qp_res.h>
#include <hpipm_d_ocp_qp_sol.h>



#ifdef __cplusplus
extern "C" {
#endif



struct d_ocp_qp_ipm_arg
	{
	double mu0; // initial value for complementarity slackness
	double alpha_min; // exit cond on step length
	double res_g_max; // exit cond on inf norm of residuals
	double res_b_max; // exit cond on inf norm of residuals
	double res_d_max; // exit cond on inf norm of residuals
	double res_m_max; // exit cond on inf norm of residuals
	double reg_prim; // reg of primal hessian
	double lam_min; // min value in lam vector
	double t_min; // min value in t vector
	int iter_max; // exit cond in iter number
	int stat_max; // iterations saved in stat
	int pred_corr; // use Mehrotra's predictor-corrector IPM algirthm
	int cond_pred_corr; // conditional Mehrotra's predictor-corrector
	int itref_pred_max; // max number of iterative refinement steps for predictor step
	int itref_corr_max; // max number of iterative refinement steps for corrector step
	int warm_start; // 0 no warm start, 1 warm start primal sol, 2 warm start primal and dual sol
	int square_root_alg; // 0 classical Riccati, 1 square-root Riccati
	int lq_fact; // 0 syrk+potrf, 1 mix, 2 lq (for square_root_alg==1)
	int abs_form; // absolute IPM formulation
	int comp_dual_sol; // dual solution (only for abs_form==1)
	int comp_res_exit; // compute residuals on exit (only for abs_form==1 and comp_dual_sol==1)
	int comp_res_pred; // compute residuals of prediction
	int mode;
	int memsize;
	};



struct d_ocp_qp_ipm_ws
	{
	double qp_res[4]; // infinity norm of residuals
	struct d_core_qp_ipm_workspace *core_workspace;
	struct d_ocp_qp_res_ws *res_workspace;
	struct d_ocp_qp_sol *sol_step;
	struct d_ocp_qp_sol *sol_itref;
	struct d_ocp_qp *qp_step;
	struct d_ocp_qp *qp_itref;
	struct d_ocp_qp_res *res_itref;
	struct d_ocp_qp_res *res;
	struct blasfeo_dvec *Gamma; // hessian update
	struct blasfeo_dvec *gamma; // hessian update
	struct blasfeo_dvec *tmp_nxM; // work space of size nxM
	struct blasfeo_dvec *tmp_nbgM; // work space of size nbM+ngM
	struct blasfeo_dvec *tmp_nsM; // work space of size nsM
	struct blasfeo_dvec *Pb; // Pb
	struct blasfeo_dvec *Zs_inv;
	struct blasfeo_dmat *L;
	struct blasfeo_dmat *Ls; // TODO
	struct blasfeo_dmat *P; // TODO
	struct blasfeo_dmat *Lh;
	struct blasfeo_dmat *AL;
	struct blasfeo_dmat *lq0;
	struct blasfeo_dvec *tmp_m;
	double *stat; // convergence statistics
	int *use_hess_fact;
	void *lq_work0;
	int iter; // iteration number
	int stat_max; // iterations saved in stat
	int stat_m; // number of recorded stat per IPM iter
	int use_Pb;
	int status; // solver status
	int memsize;
	};



//
int d_ocp_qp_ipm_arg_strseize();
//
int d_ocp_qp_ipm_arg_memsize(struct d_ocp_qp_dim *ocp_dim);
//
void d_ocp_qp_ipm_arg_create(struct d_ocp_qp_dim *ocp_dim, struct d_ocp_qp_ipm_arg *arg, void *mem);
//
void d_ocp_qp_ipm_arg_set_default(enum hpipm_mode mode, struct d_ocp_qp_ipm_arg *arg);
//
void d_ocp_qp_ipm_arg_set(char *field, void *value, struct d_ocp_qp_ipm_arg *arg);
// set maximum number of iterations
void d_ocp_qp_ipm_arg_set_iter_max(int *iter_max, struct d_ocp_qp_ipm_arg *arg);
// set minimum step lenght
void d_ocp_qp_ipm_arg_set_alpha_min(double *alpha_min, struct d_ocp_qp_ipm_arg *arg);
// set initial value of barrier parameter
void d_ocp_qp_ipm_arg_set_mu0(double *mu0, struct d_ocp_qp_ipm_arg *arg);
// set exit tolerance on stationarity condition
void d_ocp_qp_ipm_arg_set_tol_stat(double *tol_stat, struct d_ocp_qp_ipm_arg *arg);
// set exit tolerance on equality constr
void d_ocp_qp_ipm_arg_set_tol_eq(double *tol_eq, struct d_ocp_qp_ipm_arg *arg);
// set exit tolerance on inequality constr
void d_ocp_qp_ipm_arg_set_tol_ineq(double *tol_ineq, struct d_ocp_qp_ipm_arg *arg);
// set exit tolerance on complementarity condition
void d_ocp_qp_ipm_arg_set_tol_comp(double *tol_comp, struct d_ocp_qp_ipm_arg *arg);
// set regularization of primal variables
void d_ocp_qp_ipm_arg_set_reg_prim(double *tol_comp, struct d_ocp_qp_ipm_arg *arg);
// set warm start: 0 no warm start, 1 primal var
void d_ocp_qp_ipm_arg_set_warm_start(int *warm_start, struct d_ocp_qp_ipm_arg *arg);
// Mehrotra's predictor-corrector IPM algorithm: 0 no predictor-corrector, 1 use predictor-corrector
void d_ocp_qp_ipm_arg_set_pred_corr(int *pred_corr, struct d_ocp_qp_ipm_arg *arg);
// set riccati algorithm: 0 classic, 1 square-root
void d_ocp_qp_ipm_arg_set_ric_alg(int *alg, struct d_ocp_qp_ipm_arg *arg);
// compute residuals of prediction
void d_ocp_qp_ipm_arg_set_comp_res_pred(int *alg, struct d_ocp_qp_ipm_arg *arg);

//
int d_ocp_qp_ipm_ws_strsize();
//
int d_ocp_qp_ipm_ws_memsize(struct d_ocp_qp_dim *ocp_dim, struct d_ocp_qp_ipm_arg *arg);
//
void d_ocp_qp_ipm_ws_create(struct d_ocp_qp_dim *ocp_dim, struct d_ocp_qp_ipm_arg *arg, struct d_ocp_qp_ipm_ws *ws, void *mem);
//
void d_ocp_qp_ipm_get(char *field, struct d_ocp_qp_ipm_ws *ws, void *value);
//
void d_ocp_qp_ipm_get_status(struct d_ocp_qp_ipm_ws *ws, int *status);
//
void d_ocp_qp_ipm_get_iter(struct d_ocp_qp_ipm_ws *ws, int *iter);
//
void d_ocp_qp_ipm_get_max_res_stat(struct d_ocp_qp_ipm_ws *ws, double *res_stat);
//
void d_ocp_qp_ipm_get_max_res_eq(struct d_ocp_qp_ipm_ws *ws, double *res_eq);
//
void d_ocp_qp_ipm_get_max_res_ineq(struct d_ocp_qp_ipm_ws *ws, double *res_ineq);
//
void d_ocp_qp_ipm_get_max_res_comp(struct d_ocp_qp_ipm_ws *ws, double *res_comp);
//
void d_ocp_qp_ipm_get_stat(struct d_ocp_qp_ipm_ws *ws, double **stat);
//
void d_ocp_qp_ipm_get_stat_m(struct d_ocp_qp_ipm_ws *ws, int *stat_m);
//
void d_ocp_qp_ipm_solve(struct d_ocp_qp *qp, struct d_ocp_qp_sol *qp_sol, struct d_ocp_qp_ipm_arg *arg, struct d_ocp_qp_ipm_ws *ws);
//
void d_ocp_qp_ipm_predict(struct d_ocp_qp *qp, struct d_ocp_qp_sol *qp_sol, struct d_ocp_qp_ipm_arg *arg, struct d_ocp_qp_ipm_ws *ws);
//
void d_ocp_qp_ipm_sens(struct d_ocp_qp *qp, struct d_ocp_qp_sol *qp_sol, struct d_ocp_qp_ipm_arg *arg, struct d_ocp_qp_ipm_ws *ws);




#ifdef __cplusplus
}	// #extern "C"
#endif


#endif // HPIPM_D_OCP_QP_IPM_H_
