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



#ifndef HPIPM_D_DENSE_QP_IPM_H_
#define HPIPM_D_DENSE_QP_IPM_H_



#include <blasfeo_target.h>
#include <blasfeo_common.h>

#include "hpipm_x_dense_qp_ipm.h"


#ifdef __cplusplus
extern "C" {
#endif



struct d_dense_qp_ipm_arg
	{
	double mu0; // initial value for duality measure
	double alpha_min; // exit cond on step length
	double res_g_max; // exit cond on inf norm of residuals
	double res_b_max; // exit cond on inf norm of residuals
	double res_d_max; // exit cond on inf norm of residuals
	double res_m_max; // exit cond on inf norm of residuals
	double reg_prim; // reg of primal hessian
	double reg_dual; // reg of dual hessian
	double lam_min; // min value in lam vector
	double t_min; // min value in t vector
	int iter_max; // exit cond in iter number
	int stat_max; // iterations saved in stat
	int pred_corr; // Mehrotra's predictor-corrector IPM algirthm
	int cond_pred_corr; // conditional Mehrotra's predictor-corrector
	int scale; // scale hessian
	int itref_pred_max; // max number of iterative refinement steps for predictor step
	int itref_corr_max; // max number of iterative refinement steps for corrector step
	int warm_start; // 0 no warm start, 1 warm start primal sol
	int lq_fact; // 0 syrk+potrf, 1 mix, 2 lq
	int abs_form; // absolute IPM formulation
	int comp_res_exit; // compute residuals on exit (only for abs_form==1)
	int memsize;
	};



struct d_dense_qp_ipm_workspace
	{
	struct d_core_qp_ipm_workspace *core_workspace;
	struct d_dense_qp_res_workspace *res_workspace;
	struct d_dense_qp_sol *sol_step;
	struct d_dense_qp_sol *sol_itref;
	struct d_dense_qp *qp_step;
	struct d_dense_qp *qp_itref;
	struct d_dense_qp_res *res_itref;
	struct d_dense_qp_res *res;
	struct blasfeo_dvec *Gamma; //
	struct blasfeo_dvec *gamma; //
	struct blasfeo_dvec *Zs_inv; //
	struct blasfeo_dmat *Lv; //
	struct blasfeo_dmat *AL; //
	struct blasfeo_dmat *Le; //
	struct blasfeo_dmat *Ctx; //
	struct blasfeo_dvec *lv; //
	struct blasfeo_dvec *sv; // scale for Lv
	struct blasfeo_dvec *se; // scale for Le
	struct blasfeo_dvec *tmp_nbg; // work space of size nb+ng
	struct blasfeo_dvec *tmp_ns; // work space of size ns
	struct blasfeo_dmat *lq0;
	struct blasfeo_dmat *lq1;
	struct blasfeo_dvec *tmp_m;
	double *stat; // convergence statistics
//	int *ipiv_v;
//	int *ipiv_e;
	void *lq_work0;
	void *lq_work1;
	double qp_res[4]; // infinity norm of residuals
	int iter; // iteration number
	int stat_max; // iterations saved in stat
	int scale;
	int use_hess_fact;
	int memsize; // memory size (in bytes) of workspace
	};



//
int d_memsize_dense_qp_ipm_arg(struct d_dense_qp_dim *dim);
//
void d_create_dense_qp_ipm_arg(struct d_dense_qp_dim *dim, struct d_dense_qp_ipm_arg *arg, void *mem);
//
void d_set_default_dense_qp_ipm_arg(enum dense_qp_ipm_mode mode, struct d_dense_qp_ipm_arg *arg);

//
int d_memsize_dense_qp_ipm(struct d_dense_qp_dim *qp_dim, struct d_dense_qp_ipm_arg *arg);
//
void d_create_dense_qp_ipm(struct d_dense_qp_dim *qp_dim, struct d_dense_qp_ipm_arg *arg, struct d_dense_qp_ipm_workspace *ws, void *mem);
//
int d_solve_dense_qp_ipm(struct d_dense_qp *qp, struct d_dense_qp_sol *qp_sol, struct d_dense_qp_ipm_arg *arg, struct d_dense_qp_ipm_workspace *ws);



#ifdef __cplusplus
} /* extern "C" */
#endif



#endif // HPIPM_D_DENSE_QP_IPM_H_
