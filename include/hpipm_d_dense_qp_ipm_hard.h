/**************************************************************************************************
*                                                                                                 *
* This file is part of HPIPM.                                                                     *
*                                                                                                 *
* HPIPM -- High Performance Interior Point Method.                                                *
* Copyright (C) 2017 by Gianluca Frison.                                                          *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
*                                                                                                 *
* HPMPC is free software; you can redistribute it and/or                                          *
* modify it under the terms of the GNU Lesser General Public                                      *
* License as published by the Free Software Foundation; either                                    *
* version 2.1 of the License, or (at your option) any later version.                              *
*                                                                                                 *
* HPMPC is distributed in the hope that it will be useful,                                        *
* but WITHOUT ANY WARRANTY; without even the implied warranty of                                  *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                                            *
* See the GNU Lesser General Public License for more details.                                     *
*                                                                                                 *
* You should have received a copy of the GNU Lesser General Public                                *
* License along with HPMPC; if not, write to the Free Software                                    *
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA                  *
*                                                                                                 *
* Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             *
*                                                                                                 *
**************************************************************************************************/



#include <blasfeo_target.h>
#include <blasfeo_common.h>



struct d_ipm_hard_dense_qp_workspace
	{
	struct d_ipm_hard_core_qp_workspace *core_workspace;
	struct d_strvec *dv; // step in v
	struct d_strvec *dpi; // step in pi
	struct d_strvec *dlam; // step in lam XXX needed ???
	struct d_strvec *dt; // step in t XXX needed ???
	struct d_strvec *dt_lb; //
	struct d_strvec *dt_ub; // XXX needed ???
	struct d_strvec *dt_lg; //
	struct d_strvec *dt_ug; // XXX needed ???
	struct d_strvec *res_g; // q-residuals
	struct d_strvec *res_b; // b-residuals
	struct d_strvec *res_d; // d-residuals
	struct d_strvec *res_d_lb; // d-residuals
	struct d_strvec *res_d_ub; // d-residuals
	struct d_strvec *res_d_lg; // d-residuals
	struct d_strvec *res_d_ug; // d-residuals
	struct d_strvec *res_m; // m-residuals
	struct d_strvec *Qx; // hessian update
	struct d_strvec *qx; // gradient update
	struct d_strmat *Lv; //
	struct d_strmat *AL; //
	struct d_strmat *Le; //
	struct d_strmat *Ctx; //
	struct d_strvec *lv; //
	struct d_strvec *tmp_nb; // work space of size nb
	struct d_strvec *tmp_ng0; // work space of size nb
	struct d_strvec *tmp_ng1; // work space of size nb
	double *stat; // convergence statistics
	double res_mu; // mu-residual
	int iter; // iteration number
	int memsize; // memory size (in bytes) of workspace
	};



struct d_ipm_hard_dense_qp_arg
	{
	double alpha_min; // exit cond on step length
	double mu_max; // exit cond on duality measure
	double mu0; // initial value for duality measure
	int iter_max; // exit cond in iter number
	};



//
int d_memsize_ipm_hard_dense_qp(struct d_dense_qp *qp, struct d_ipm_hard_dense_qp_arg *arg);
//
void d_create_ipm_hard_dense_qp(struct d_dense_qp *qp, struct d_ipm_hard_dense_qp_arg *arg, struct d_ipm_hard_dense_qp_workspace *ws, void *mem);
//
void d_solve_ipm_hard_dense_qp(struct d_dense_qp *qp, struct d_dense_qp_sol *qp_sol, struct d_ipm_hard_dense_qp_workspace *ws);
//
void d_solve_ipm2_hard_dense_qp(struct d_dense_qp *qp, struct d_dense_qp_sol *qp_sol, struct d_ipm_hard_dense_qp_workspace *ws);
