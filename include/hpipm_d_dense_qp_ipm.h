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



struct d_dense_qp_ipm_workspace
	{
	struct d_core_qp_ipm_workspace *core_workspace;
	struct d_strvec *dv; // step in v
	struct d_strvec *dpi; // step in pi
	struct d_strvec *dlam; // step in lam XXX needed ???
	struct d_strvec *dt; // step in t
	struct d_strvec *res_g; // q-residuals
	struct d_strvec *res_b; // b-residuals
	struct d_strvec *res_d; // d-residuals
	struct d_strvec *res_m; // m-residuals
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
	double mu0; // mu0
	double res_mu; // mu-residual
	int iter; // iteration number
	int stat_max; // iterations saved in stat
	int memsize; // memory size (in bytes) of workspace
	};



struct d_dense_qp_ipm_arg
	{
	double alpha_min; // exit cond on step length
	double mu_max; // exit cond on duality measure
	double mu0; // initial value for duality measure
	int iter_max; // exit cond in iter number
	int stat_max; // iterations saved in stat
	int pred_corr; // use Mehrotra's predictor-corrector IPM algirthm
	};



//
int d_memsize_dense_qp_ipm(struct d_dense_qp *qp, struct d_dense_qp_ipm_arg *arg);
//
void d_create_dense_qp_ipm(struct d_dense_qp *qp, struct d_dense_qp_ipm_arg *arg, struct d_dense_qp_ipm_workspace *ws, void *mem);
//
int d_solve_dense_qp_ipm(struct d_dense_qp *qp, struct d_dense_qp_sol *qp_sol, struct d_dense_qp_ipm_arg *arg, struct d_dense_qp_ipm_workspace *ws);
