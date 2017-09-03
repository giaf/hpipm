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




struct d_tree_ocp_qp_ipm_workspace
	{
	struct d_core_qp_ipm_workspace *core_workspace;
	struct d_strvec *dux;
	struct d_strvec *dpi;
	struct d_strvec *dt;
	struct d_strvec *res_g; // q-residuals
	struct d_strvec *res_b; // b-residuals
	struct d_strvec *res_d; // d-residuals
	struct d_strvec *res_m; // m-residuals
	struct d_strvec *Gamma; // hessian update
	struct d_strvec *gamma; // hessian update
	struct d_strvec *tmp_nxM; // work space of size nxM
	struct d_strvec *tmp_nbgM; // work space of size nbgM
	struct d_strvec *tmp_nsM; // work space of size nsM
	struct d_strvec *Pb; // Pb
	struct d_strvec *Zs_inv;
	struct d_strmat *L;
	struct d_strmat *AL;
	double *stat; // convergence statistics
	double mu0; // mu0
	double res_mu; // mu-residual
	int iter; // iteration number
	int memsize;
	};



struct d_tree_ocp_qp_ipm_arg
	{
	double alpha_min; // exit cond on step length
	double mu_max; // exit cond on duality measure
	double mu0; // initial value for duality measure
	int iter_max; // exit cond in iter number
	int stat_max; // iterations saved in stat
	int pred_corr; // use Mehrotra's predictor-corrector IPM algirthm
	};



//
int d_memsize_tree_ocp_qp_ipm(struct d_tree_ocp_qp *qp, struct d_tree_ocp_qp_ipm_arg *arg);
//
void d_create_tree_ocp_qp_ipm(struct d_tree_ocp_qp *qp, struct d_tree_ocp_qp_ipm_arg *arg, struct d_tree_ocp_qp_ipm_workspace *ws, void *mem);
//
int d_solve_tree_ocp_qp_ipm(struct d_tree_ocp_qp *qp, struct d_tree_ocp_qp_sol *qp_sol, struct d_tree_ocp_qp_ipm_arg *arg, struct d_tree_ocp_qp_ipm_workspace *ws);
