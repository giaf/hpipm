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



#ifndef HPIPM_S_OCP_QP_IPM_H_
#define HPIPM_S_OCP_QP_IPM_H_



#include <blasfeo_target.h>
#include <blasfeo_common.h>

#include <hpipm_common.h>
#include <hpipm_s_ocp_qp_dim.h>
#include <hpipm_s_ocp_qp.h>
#include <hpipm_s_ocp_qp_res.h>
#include <hpipm_s_ocp_qp_sol.h>



#ifdef __cplusplus
extern "C" {
#endif



struct s_ocp_qp_ipm_arg
	{
	float mu0; // initial value for complementarity slackness
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
	int square_root_alg; // 0 classical Riccati, 1 square-root Riccati
	int lq_fact; // 0 syrk+potrf, 1 mix, 2 lq (for square_root_alg==1)
	int abs_form; // absolute IPM formulation
	int comp_dual_sol; // dual solution (only for abs_form==1)
	int comp_res_exit; // compute residuals on exit (only for abs_form==1 and comp_dual_sol==1)
	int mode;
	int memsize;
	};



struct s_ocp_qp_ipm_ws
	{
	struct s_core_qp_ipm_workspace *core_workspace;
	struct s_ocp_qp_res_ws *res_workspace;
	struct s_ocp_qp_sol *sol_step;
	struct s_ocp_qp_sol *sol_itref;
	struct s_ocp_qp *qp_step;
	struct s_ocp_qp *qp_itref;
	struct s_ocp_qp_res *res;
	struct s_ocp_qp_res *res_itref;
	struct blasfeo_svec *Gamma; // hessian update
	struct blasfeo_svec *gamma; // hessian update
	struct blasfeo_svec *tmp_nxM; // work space of size nxM
	struct blasfeo_svec *tmp_nbgM; // work space of size nbM+ngM
	struct blasfeo_svec *tmp_nsM; // work space of size nsM
	struct blasfeo_svec *Pb; // Pb
	struct blasfeo_svec *Zs_inv;
	struct blasfeo_smat *L;
	struct blasfeo_smat *Ls; // TODO
	struct blasfeo_smat *P; // TODO
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
	int stat_m; // number of recorded stat per IPM iter
	int use_Pb;
	int status; // solver status
	int memsize;
	};



//
int s_ocp_qp_ipm_arg_strsize();
//
int s_ocp_qp_ipm_arg_memsize(struct s_ocp_qp_dim *ocp_dim);
//
void s_ocp_qp_ipm_arg_create(struct s_ocp_qp_dim *ocp_dim, struct s_ocp_qp_ipm_arg *arg, void *mem);
//
void s_ocp_qp_ipm_arg_set_default(enum hpipm_mode mode, struct s_ocp_qp_ipm_arg *arg);
//
void s_ocp_qp_ipm_arg_set(char *field, void *value, struct s_ocp_qp_ipm_arg *arg);
// set maximum number of iterations
void s_ocp_qp_ipm_arg_set_iter_max(int *iter_max, struct s_ocp_qp_ipm_arg *arg);
// set minimum step lenght
void s_ocp_qp_ipm_arg_set_alpha_min(float *alpha_min, struct s_ocp_qp_ipm_arg *arg);
// set initial value of barrier parameter
void s_ocp_qp_ipm_arg_set_mu0(float *mu0, struct s_ocp_qp_ipm_arg *arg);
// set exit tolerance on stationarity condition
void s_ocp_qp_ipm_arg_set_tol_stat(float *tol_stat, struct s_ocp_qp_ipm_arg *arg);
// set exit tolerance on equality constr
void s_ocp_qp_ipm_arg_set_tol_eq(float *tol_eq, struct s_ocp_qp_ipm_arg *arg);
// set exit tolerance on inequality constr
void s_ocp_qp_ipm_arg_set_tol_ineq(float *tol_ineq, struct s_ocp_qp_ipm_arg *arg);
// set exit tolerance on complementarity condition
void s_ocp_qp_ipm_arg_set_tol_comp(float *tol_comp, struct s_ocp_qp_ipm_arg *arg);
// set regularization of primal variables
void s_ocp_qp_ipm_arg_set_reg_prim(float *tol_comp, struct s_ocp_qp_ipm_arg *arg);
// set warm start: 0 no warm start, 1 primal var
void s_ocp_qp_ipm_arg_set_warm_start(int *warm_start, struct s_ocp_qp_ipm_arg *arg);
// Mehrotra's predictor-corrector IPM algorithm: 0 no predictor-corrector, 1 use predictor-corrector
void s_ocp_qp_ipm_arg_set_pred_corr(int *pred_corr, struct s_ocp_qp_ipm_arg *arg);
// set riccati algorithm: 0 classic, 1 square-root
void s_ocp_qp_ipm_arg_set_ric_alg(int *alg, struct s_ocp_qp_ipm_arg *arg);

//
int s_ocp_qp_ipm_ws_strsize();
//
int s_ocp_qp_ipm_ws_memsize(struct s_ocp_qp_dim *ocp_dim, struct s_ocp_qp_ipm_arg *arg);
//
void s_ocp_qp_ipm_ws_create(struct s_ocp_qp_dim *ocp_dim, struct s_ocp_qp_ipm_arg *arg, struct s_ocp_qp_ipm_ws *ws, void *mem);
//
void s_ocp_qp_ipm_get(char *field, struct s_ocp_qp_ipm_ws *ws, void *value);
//
void s_ocp_qp_ipm_get_status(struct s_ocp_qp_ipm_ws *ws, int *status);
//
void s_ocp_qp_ipm_get_iter(struct s_ocp_qp_ipm_ws *ws, int *iter);
//
void s_ocp_qp_ipm_get_res_stat(struct s_ocp_qp_ipm_ws *ws, float *res_stat);
//
void s_ocp_qp_ipm_get_res_eq(struct s_ocp_qp_ipm_ws *ws, float *res_eq);
//
void s_ocp_qp_ipm_get_res_ineq(struct s_ocp_qp_ipm_ws *ws, float *res_ineq);
//
void s_ocp_qp_ipm_get_res_comp(struct s_ocp_qp_ipm_ws *ws, float *res_comp);
//
void s_ocp_qp_ipm_get_stat(struct s_ocp_qp_ipm_ws *ws, float **stat);
//
void s_ocp_qp_ipm_get_stat_m(struct s_ocp_qp_ipm_ws *ws, int *stat_m);
//
void s_ocp_qp_ipm_solve(struct s_ocp_qp *qp, struct s_ocp_qp_sol *qp_sol, struct s_ocp_qp_ipm_arg *arg, struct s_ocp_qp_ipm_ws *ws);



#ifdef __cplusplus
} /* extern "C" */
#endif



#endif // HPIPM_S_OCP_QP_IPM_H_
