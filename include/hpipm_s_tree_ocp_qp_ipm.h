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



#ifndef HPIPM_S_TREE_OCP_QP_IPM_H_
#define HPIPM_S_TREE_OCP_QP_IPM_H_



#include <blasfeo_target.h>
#include <blasfeo_common.h>

#include "hpipm_x_tree_ocp_qp_ipm.h"



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
	double reg_prim; // reg of primal hessian
	float lam_min; // min value in lam vector
	float t_min; // min value in t vector
	int iter_max; // exit cond in iter number
	int stat_max; // iterations saved in stat
	int pred_corr; // use Mehrotra's predictor-corrector IPM algirthm
	int cond_pred_corr; // conditional Mehrotra's predictor-corrector
	int itref_pred_max; // max number of iterative refinement steps for predictor step
	int itref_corr_max; // max number of iterative refinement steps for corrector step
	int warm_start; // 0 no warm start, 1 warm start primal sol
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
	int memsize;
	};



//
int s_memsize_tree_ocp_qp_ipm_arg(struct s_tree_ocp_qp_dim *dim);
//
void s_create_tree_ocp_qp_ipm_arg(struct s_tree_ocp_qp_dim *dim, struct s_tree_ocp_qp_ipm_arg *arg, void *mem);
//
void s_set_default_tree_ocp_qp_ipm_arg(enum tree_ocp_qp_ipm_mode mode, struct s_tree_ocp_qp_ipm_arg *arg);

//
int s_memsize_tree_ocp_qp_ipm(struct s_tree_ocp_qp_dim *dim, struct s_tree_ocp_qp_ipm_arg *arg);
//
void s_create_tree_ocp_qp_ipm(struct s_tree_ocp_qp_dim *dim, struct s_tree_ocp_qp_ipm_arg *arg, struct s_tree_ocp_qp_ipm_workspace *ws, void *mem);
//
int s_solve_tree_ocp_qp_ipm(struct s_tree_ocp_qp *qp, struct s_tree_ocp_qp_sol *qp_sol, struct s_tree_ocp_qp_ipm_arg *arg, struct s_tree_ocp_qp_ipm_workspace *ws);

#ifdef __cplusplus
} /* extern "C" */
#endif



#endif // HPIPM_S_TREE_OCP_QP_IPM_H_
