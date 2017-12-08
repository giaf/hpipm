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
	int iter_max; // exit cond in iter number
	int stat_max; // iterations saved in stat
	int pred_corr; // Mehrotra's predictor-corrector IPM algirthm
	int cond_pred_corr; // conditional Mehrotra's predictor-corrector
	int warm_start; // 0 no warm start, 1 warm start primal sol
	int memsize;
	};



struct d_dense_qp_ipm_workspace
	{
	struct d_core_qp_ipm_workspace *core_workspace;
	struct d_dense_qp_res *res;
	struct d_dense_qp_res_workspace *res_workspace;
	struct d_dense_qp_sol *step;
	struct d_dense_qp *itref_qp;
	struct d_dense_qp_res *itref_res;
	struct d_strvec *Gamma; //
	struct d_strvec *gamma; //
	struct d_strvec *Zs_inv; //
	struct d_strmat *Lv; //
	struct d_strmat *AL; //
	struct d_strmat *Le; //
	struct d_strmat *Ctx; //
	struct d_strvec *lv; //
	struct d_strvec *tmp_nbg; // work space of size nb+ng
	struct d_strvec *tmp_ns; // work space of size ns
	double *stat; // convergence statistics
	int *ipiv;
	double qp_res[4]; // infinity norm of residuals
	double mu0; // mu0
	int iter; // iteration number
	int stat_max; // iterations saved in stat
	int warm_start; // 0 no warm start, 1 warm start primal sol
	int memsize; // memory size (in bytes) of workspace
	};



//
int d_memsize_dense_qp_ipm_arg(struct d_dense_qp_dim *dim);
//
void d_create_dense_qp_ipm_arg(struct d_dense_qp_dim *dim, struct d_dense_qp_ipm_arg *arg, void *mem);
//
void d_set_default_dense_qp_ipm_arg(struct d_dense_qp_ipm_arg *arg);

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
