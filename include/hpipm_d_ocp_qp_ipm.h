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




struct d_ocp_qp_ipm_workspace
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
	struct d_strvec *tmp_nbgM; // work space of size nbM+ngM
	struct d_strvec *tmp_nsM; // work space of size nsM
	struct d_strvec *Pb; // Pb
	struct d_strvec *Zs_inv;
	struct d_strmat *L;
	struct d_strmat *AL;
	double *stat; // convergence statistics
	double mu0; // mu0
	double res_mu; // mu-residual
	int iter; // iteration number
	int stat_max; // iterations saved in stat
	int memsize;
	};



struct d_ocp_qp_ipm_arg
	{
	double mu0; // initial value for duality measure
	double alpha_min; // exit cond on step length
	double res_g_max; // exit cond on inf norm of residuals
	double res_b_max; // exit cond on inf norm of residuals
	double res_d_max; // exit cond on inf norm of residuals
	double res_m_max; // exit cond on inf norm of residuals
	int iter_max; // exit cond in iter number
	int stat_max; // iterations saved in stat
	int pred_corr; // use Mehrotra's predictor-corrector IPM algirthm
	};



//
int d_memsize_ocp_qp_ipm(struct d_ocp_qp *qp, struct d_ocp_qp_ipm_arg *arg);
//
void d_create_ocp_qp_ipm(struct d_ocp_qp *qp, struct d_ocp_qp_ipm_arg *arg, struct d_ocp_qp_ipm_workspace *ws, void *mem);
//
int d_solve_ocp_qp_ipm(struct d_ocp_qp *qp, struct d_ocp_qp_sol *qp_sol, struct d_ocp_qp_ipm_arg *arg, struct d_ocp_qp_ipm_workspace *ws);
//
void d_cvt_ocp_qp_res_to_colmaj(struct d_ocp_qp *qp, struct d_ocp_qp_ipm_workspace *ws, double **res_r, double **res_q, double **res_ls, double **res_us, double **res_b, double **res_d_lb, double **res_d_ub, double **res_d_lg, double **res_d_ug, double **res_d_ls, double **res_d_us, double **res_m_lb, double **res_m_ub, double **res_m_lg, double **res_m_ug, double **res_m_ls, double **res_m_us);
//
void d_cvt_ocp_qp_res_to_rowmaj(struct d_ocp_qp *qp, struct d_ocp_qp_ipm_workspace *ws, double **res_r, double **res_q, double **res_ls, double **res_us, double **res_b, double **res_d_lb, double **res_d_ub, double **res_d_lg, double **res_d_ug, double **res_d_ls, double **res_d_us, double **res_m_lb, double **res_m_ub, double **res_m_lg, double **res_m_ug, double **res_m_ls, double **res_m_us);

