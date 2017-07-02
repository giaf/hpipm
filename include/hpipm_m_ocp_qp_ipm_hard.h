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




struct m_ipm_hard_ocp_qp_workspace
	{
	struct d_ipm_hard_core_qp_workspace *core_workspace;
	struct d_strvec *dux;
	struct d_strvec *dpi;
	struct d_strvec *dt_lb;
	struct d_strvec *dt_lg;
	struct s_strvec *sdux; // XXX
	struct s_strvec *sdpi; // XXX
	struct d_strvec *res_g; // q-residuals
	struct d_strvec *res_b; // b-residuals
	struct d_strvec *res_d; // d-residuals XXX remove ???
	struct d_strvec *res_d_lb; // d-residuals
	struct d_strvec *res_d_ub; // d-residuals
	struct d_strvec *res_d_lg; // d-residuals
	struct d_strvec *res_d_ug; // d-residuals
	struct d_strvec *res_m; // m-residuals
	struct d_strvec *res_m_lb; // m-residuals
	struct d_strvec *res_m_ub; // m-residuals
	struct d_strvec *res_m_lg; // m-residuals
	struct d_strvec *res_m_ug; // m-residuals
	struct s_strvec *sres_g; // q-residuals // XXX
	struct s_strvec *sres_b; // b-residuals // XXX
	struct d_strvec *Qx_lb; // hessian update
	struct d_strvec *Qx_lg; // hessian update
	struct d_strvec *qx_lb; // gradient update
	struct d_strvec *qx_lg; // gradient update
	struct s_strvec *sQx_lb; // hessian update // XXX
	struct s_strvec *sQx_lg; // hessian update // XXX
	struct s_strvec *sqx_lb; // gradient update // XXX
	struct s_strvec *sqx_lg; // gradient update // XXX
	struct d_strvec *tmp_nbM; // work space of size nbM
	struct s_strvec *tmp_nxM; // work space of size nxM // XXX
	struct d_strvec *tmp_ngM; // work space of size ngM
	struct s_strvec *Pb; // Pb // XXX
	struct s_strmat *L; // XXX
	struct s_strmat *AL; // XXX
	double *stat; // convergence statistics
	double res_mu; // mu-residual
	int iter; // iteration number
	int compute_Pb;
	};



struct m_ipm_hard_ocp_qp_arg
	{
	double alpha_min; // exit cond on step length
	double mu_max; // exit cond on duality measure
	double mu0; // initial value for duality measure
	int iter_max; // exit cond in iter number
	};



//
int m_memsize_ipm_hard_ocp_qp(struct d_ocp_qp *d_qp, struct s_ocp_qp *s_qp, struct m_ipm_hard_ocp_qp_arg *arg);
//
void m_create_ipm_hard_ocp_qp(struct d_ocp_qp *d_qp, struct s_ocp_qp *s_qp, struct m_ipm_hard_ocp_qp_arg *arg, struct m_ipm_hard_ocp_qp_workspace *ws, void *mem);
//
void m_solve_ipm_hard_ocp_qp(struct d_ocp_qp *d_qp, struct s_ocp_qp *s_qp, struct d_ocp_qp_sol *qp_sol, struct m_ipm_hard_ocp_qp_workspace *ws);
//
void m_solve_ipm2_hard_ocp_qp(struct d_ocp_qp *d_qp, struct s_ocp_qp *s_qp, struct d_ocp_qp_sol *qp_sol, struct m_ipm_hard_ocp_qp_workspace *ws);

