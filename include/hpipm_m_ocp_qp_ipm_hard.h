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

#ifdef __cplusplus
extern "C" {
#endif


struct m_ipm_hard_ocp_qp_workspace
	{
	struct d_ipm_hard_core_qp_workspace *core_workspace;
	struct blasfeo_dvec *dux;
	struct blasfeo_dvec *dpi;
	struct blasfeo_dvec *dt_lb;
	struct blasfeo_dvec *dt_lg;
	struct blasfeo_svec *sdux; // XXX
	struct blasfeo_svec *sdpi; // XXX
	struct blasfeo_dvec *res_g; // q-residuals
	struct blasfeo_dvec *res_b; // b-residuals
	struct blasfeo_dvec *res_d; // d-residuals XXX remove ???
	struct blasfeo_dvec *res_d_lb; // d-residuals
	struct blasfeo_dvec *res_d_ub; // d-residuals
	struct blasfeo_dvec *res_d_lg; // d-residuals
	struct blasfeo_dvec *res_d_ug; // d-residuals
	struct blasfeo_dvec *res_m; // m-residuals
	struct blasfeo_dvec *res_m_lb; // m-residuals
	struct blasfeo_dvec *res_m_ub; // m-residuals
	struct blasfeo_dvec *res_m_lg; // m-residuals
	struct blasfeo_dvec *res_m_ug; // m-residuals
	struct blasfeo_svec *sres_g; // q-residuals // XXX
	struct blasfeo_svec *sres_b; // b-residuals // XXX
	struct blasfeo_dvec *Qx_lb; // hessian update
	struct blasfeo_dvec *Qx_lg; // hessian update
	struct blasfeo_dvec *qx_lb; // gradient update
	struct blasfeo_dvec *qx_lg; // gradient update
	struct blasfeo_svec *sQx_lb; // hessian update // XXX
	struct blasfeo_svec *sQx_lg; // hessian update // XXX
	struct blasfeo_svec *sqx_lb; // gradient update // XXX
	struct blasfeo_svec *sqx_lg; // gradient update // XXX
	struct blasfeo_dvec *tmp_nbM; // work space of size nbM
	struct blasfeo_svec *tmp_nxM; // work space of size nxM // XXX
	struct blasfeo_dvec *tmp_ngM; // work space of size ngM
	struct blasfeo_svec *Pb; // Pb // XXX
	struct blasfeo_smat *L; // XXX
	struct blasfeo_smat *AL; // XXX
	struct blasfeo_svec *sSx; // scaling
	struct blasfeo_svec *sSi; // scaling inverted
	double *stat; // convergence statistics
	double res_mu; // mu-residual
	int iter; // iteration number
	int compute_Pb;
	int scale;
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

#ifdef __cplusplus
} /* extern "C" */
#endif
