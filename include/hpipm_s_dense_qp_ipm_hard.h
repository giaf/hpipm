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



struct s_ipm_hard_dense_qp_workspace
	{
	struct s_ipm_hard_core_qp_workspace *core_workspace;
	struct s_strvec *dv; // step in v
	struct s_strvec *dpi; // step in pi
	struct s_strvec *dlam; // step in lam XXX needed ???
	struct s_strvec *dt; // step in t XXX needed ???
	struct s_strvec *dt_lb; //
	struct s_strvec *dt_ub; // XXX needed ???
	struct s_strvec *dt_lg; //
	struct s_strvec *dt_ug; // XXX needed ???
	struct s_strvec *res_g; // q-residuals
	struct s_strvec *res_b; // b-residuals
	struct s_strvec *res_d; // d-residuals
	struct s_strvec *res_d_lb; // d-residuals
	struct s_strvec *res_d_ub; // d-residuals
	struct s_strvec *res_d_lg; // d-residuals
	struct s_strvec *res_d_ug; // d-residuals
	struct s_strvec *res_m; // m-residuals
	struct s_strvec *Qx; // hessian update
	struct s_strvec *qx; // gradient update
	struct s_strmat *Lv; //
	struct s_strmat *AL; //
	struct s_strmat *Le; //
	struct s_strmat *Ctx; //
	struct s_strvec *lv; //
	struct s_strvec *tmp_nb; // work space of size nb
	struct s_strvec *tmp_ng0; // work space of size nb
	struct s_strvec *tmp_ng1; // work space of size nb
	float *stat; // convergence statistics
	float res_mu; // mu-residual
	int iter; // iteration number
	int memsize; // memory size (in bytes) of workspace
	};



struct s_ipm_hard_dense_qp_arg
	{
	float alpha_min; // exit cond on step length
	float mu_max; // exit cond on duality measure
	float mu0; // initial value for duality measure
	int iter_max; // exit cond in iter number
	};



//
int s_memsize_ipm_hard_dense_qp(struct s_dense_qp *qp, struct s_ipm_hard_dense_qp_arg *arg);
//
void s_create_ipm_hard_dense_qp(struct s_dense_qp *qp, struct s_ipm_hard_dense_qp_arg *arg, struct s_ipm_hard_dense_qp_workspace *ws, void *mem);
//
void s_solve_ipm_hard_dense_qp(struct s_dense_qp *qp, struct s_dense_qp_sol *qp_sol, struct s_ipm_hard_dense_qp_workspace *ws);
//
void s_solve_ipm2_hard_dense_qp(struct s_dense_qp *qp, struct s_dense_qp_sol *qp_sol, struct s_ipm_hard_dense_qp_workspace *ws);

