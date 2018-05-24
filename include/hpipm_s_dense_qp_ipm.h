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



#ifndef HPIPM_S_DENSE_QP_IPM_H_
#define HPIPM_S_DENSE_QP_IPM_H_



#ifdef __cplusplus
extern "C" {
#endif



struct s_dense_qp_ipm_arg
	{
	float mu0; // initial value for duality measure
	float alpha_min; // exit cond on step length
	float res_g_max; // exit cond on inf norm of residuals
	float res_b_max; // exit cond on inf norm of residuals
	float res_d_max; // exit cond on inf norm of residuals
	float res_m_max; // exit cond on inf norm of residuals
	float reg_prim; // reg of primal hessian
	float reg_dual; // reg of dual hessian
	float lam_min; // min value in lam vector
	float t_min; // min value in t vector
	int iter_max; // exit cond in iter number
	int stat_max; // iterations saved in stat
	int pred_corr; // Mehrotra's predictor-corrector IPM algirthm
	int cond_pred_corr; // conditional Mehrotra's predictor-corrector
	int scale; // scale hessian
	int itref_pred_max; // max number of iterative refinement steps for predictor step
	int itref_corr_max; // max number of iterative refinement steps for corrector step
	int warm_start; // 0 no warm start, 1 warm start primal sol
	int lq_fact; // 0 syrk+potrf, 1 mix, 2 lq
	int memsize;
	};



struct s_dense_qp_ipm_workspace
	{
	struct s_core_qp_ipm_workspace *core_workspace;
	struct s_dense_qp_res *res;
	struct s_dense_qp_res_workspace *res_workspace;
	struct s_dense_qp_sol *sol_step;
	struct s_dense_qp_sol *sol_itref;
	struct s_dense_qp *qp_step;
	struct s_dense_qp *qp_itref;
	struct s_dense_qp_res *res_itref;
	struct blasfeo_svec *Gamma; //
	struct blasfeo_svec *gamma; //
	struct blasfeo_svec *Zs_inv; //
	struct blasfeo_smat *Lv; //
	struct blasfeo_smat *AL; //
	struct blasfeo_smat *Le; //
	struct blasfeo_smat *Ctx; //
	struct blasfeo_svec *lv; //
	struct blasfeo_svec *sv; // scale for Lv
	struct blasfeo_svec *se; // scale for Le
	struct blasfeo_svec *tmp_nbg; // work space of size nb+ng
	struct blasfeo_svec *tmp_ns; // work space of size ns
	struct blasfeo_smat *lq0;
	struct blasfeo_smat *lq1;
	float *stat; // convergence statistics
	int *ipiv_v;
	int *ipiv_e;
	void *lq_work0;
	void *lq_work1;
	float qp_res[4]; // infinity norm of residuals
	float mu0; // mu0
	int iter; // iteration number
	int stat_max; // iterations saved in stat
	int warm_start; // 0 no warm start, 1 warm start primal sol
	int scale;
	int use_hess_fact;
	int memsize; // memory size (in bytes) of workspace
	};



//
int s_memsize_dense_qp_ipm_arg(struct s_dense_qp_dim *qp_dim);
//
void s_create_dense_qp_ipm_arg(struct s_dense_qp_dim *qp_dim, struct s_dense_qp_ipm_arg *arg, void *mem);
//
void s_set_default_dense_qp_ipm_arg(struct s_dense_qp_ipm_arg *arg);

//
int s_memsize_dense_qp_ipm(struct s_dense_qp_dim *qp_dim, struct s_dense_qp_ipm_arg *arg);
//
void s_create_dense_qp_ipm(struct s_dense_qp_dim *qp_dim, struct s_dense_qp_ipm_arg *arg, struct s_dense_qp_ipm_workspace *ws, void *mem);
//
int s_solve_dense_qp_ipm(struct s_dense_qp *qp, struct s_dense_qp_sol *qp_sol, struct s_dense_qp_ipm_arg *arg, struct s_dense_qp_ipm_workspace *ws);



#ifdef __cplusplus
} /* extern "C" */
#endif



#endif // HPIPM_S_DENSE_QP_IPM_H_
