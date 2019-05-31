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
	int memsize;
	};



struct d_ocp_qp_ipm_workspace
	{
	struct d_core_qp_ipm_workspace *core_workspace;
	struct d_ocp_qp_res_workspace *res_workspace;
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
	double qp_res[4]; // infinity norm of residuals
	int iter; // iteration number
	int stat_max; // iterations saved in stat
	int use_Pb;
	int memsize;
	};



//
int d_sizeof_ocp_qp_ipm_arg();
//
int d_memsize_ocp_qp_ipm_arg(struct d_ocp_qp_dim *ocp_dim);
//
void d_create_ocp_qp_ipm_arg(struct d_ocp_qp_dim *ocp_dim, struct d_ocp_qp_ipm_arg *arg, void *mem);
//
void d_set_default_ocp_qp_ipm_arg(enum hpipm_mode mode, struct d_ocp_qp_ipm_arg *arg);
// set maximum number of iterations
void d_set_ocp_qp_ipm_arg_iter_max(int iter_max, struct d_ocp_qp_ipm_arg *arg);
// set initial value of barrier parameter
void d_set_ocp_qp_ipm_arg_mu0(double mu0, struct d_ocp_qp_ipm_arg *arg);
// set exit tolerance on stationarity condition
void d_set_ocp_qp_ipm_arg_tol_stat(double tol_stat, struct d_ocp_qp_ipm_arg *arg);
// set exit tolerance on equality constr
void d_set_ocp_qp_ipm_arg_tol_eq(double tol_eq, struct d_ocp_qp_ipm_arg *arg);
// set exit tolerance on inequality constr
void d_set_ocp_qp_ipm_arg_tol_ineq(double tol_ineq, struct d_ocp_qp_ipm_arg *arg);
// set exit tolerance on complementarity condition
void d_set_ocp_qp_ipm_arg_tol_comp(double tol_comp, struct d_ocp_qp_ipm_arg *arg);
// set regularization of primal variables
void d_set_ocp_qp_ipm_arg_reg_prim(double tol_comp, struct d_ocp_qp_ipm_arg *arg);
// set warm start: 0 no warm start, 1 primal var
void d_set_ocp_qp_ipm_arg_warm_start(int warm_start, struct d_ocp_qp_ipm_arg *arg);
// set riccati algorithm: 0 classic, 1 square-root
void d_set_ocp_qp_ipm_arg_ric_alg(int alg, struct d_ocp_qp_ipm_arg *arg);

//
int d_sizeof_ocp_qp_ipm_workspace();
//
int d_memsize_ocp_qp_ipm(struct d_ocp_qp_dim *ocp_dim, struct d_ocp_qp_ipm_arg *arg);
//
void d_create_ocp_qp_ipm(struct d_ocp_qp_dim *ocp_dim, struct d_ocp_qp_ipm_arg *arg, struct d_ocp_qp_ipm_workspace *ws, void *mem);
//
int d_get_ocp_qp_ipm_iter(struct d_ocp_qp_ipm_workspace *ws);
//
double d_get_ocp_qp_ipm_res_stat(struct d_ocp_qp_ipm_workspace *ws);
//
double d_get_ocp_qp_ipm_res_eq(struct d_ocp_qp_ipm_workspace *ws);
//
double d_get_ocp_qp_ipm_res_ineq(struct d_ocp_qp_ipm_workspace *ws);
//
double d_get_ocp_qp_ipm_res_comp(struct d_ocp_qp_ipm_workspace *ws);
//
double *d_get_ocp_qp_ipm_stat(struct d_ocp_qp_ipm_workspace *ws);
//
int d_solve_ocp_qp_ipm(struct d_ocp_qp *qp, struct d_ocp_qp_sol *qp_sol, struct d_ocp_qp_ipm_arg *arg, struct d_ocp_qp_ipm_workspace *ws);




#ifdef __cplusplus
}	// #extern "C"
#endif


#endif // HPIPM_D_OCP_QP_IPM_H_
