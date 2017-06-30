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



struct d_ipm_hard_core_qp_workspace
	{
	double *d; // constraints
	double *d_lb; // lower box constraints
	double *d_ub; // upper box constraints
	double *d_lg; // lower general constraints
	double *d_ug; // upper general constraints
	double *v; // primal variables
	double *pi; // equality constraints multipliers
	double *lam; // inequality constraints multipliers
	double *lam_lb; // lower bounds multipliers
	double *lam_lg; // lower general constraints multipliers
	double *lam_ub; // upper bounds multipliers
	double *lam_ug; // upper general constraints multipliers
	double *t; // inequality constraints slacks
	double *t_lb; // lower box constraints slacks
	double *t_lg; // lower general constraints slacks
	double *t_ub; // upper box constraints slacks
	double *t_ug; // upper general constraints slacks
	double *t_inv; // inverse of t
	double *t_inv_lb; // inverse of t
	double *t_inv_ub; // inverse of t
	double *t_inv_lg; // inverse of t
	double *t_inv_ug; // inverse of t
	double *dv; // step in v
	double *dpi; // step in pi
	double *dlam; // step in lam
	double *dlam_lb; //
	double *dlam_lg; //
	double *dlam_ub; //
	double *dlam_ug; //
	double *dt; // step in t
	double *dt_lb; // step in t_lb
	double *dt_ub; // step in t_ub
	double *dt_lg; // step in t_lg
	double *dt_ug; // step in t_ug
	double *res_g; // q-residuals
	double *res_b; // b-residuals
	double *res_d; // d-residuals
	double *res_d_lb; // d-residuals
	double *res_d_ub; // d-residuals
	double *res_d_lg; // d-residuals
	double *res_d_ug; // d-residuals
	double *res_m; // m-residuals
	double *res_m_lb; // m-residuals
	double *res_m_ub; // m-residuals
	double *res_m_lg; // m-residuals
	double *res_m_ug; // m-residuals
	double *Qx; // Hessian update
	double *Qx_lb; // Hessian update
	double *Qx_lg; // Hessian update
	double *qx; // gradient update
	double *qx_lb; // gradient update
	double *qx_lg; // gradient update
	double *stat; // convergence statistics
	double alpha; // step length
	double alpha_min; // exit cond on step lenght
	double sigma; // centering XXX
	double mu; // duality measuere
	double mu_aff; // affine duality measuere
	double mu0; // initial duality measuere
	double mu_max; // exit cond on mu
	double nt_inv; // 1.0/nt, where nt is the total number of constraints
	int nv; // number of primal variables
	int ne; // number of equality constraints
	int nb; // number of two-sized bounds
	int ng; // number of two-sized constraints
	int iter_max; // exit cond on iter mumber
	int memsize; // memory size (in bytes) of workspace
	};



//
int d_memsize_ipm_hard_core_qp(int nv, int ne, int nb, int ng, int iter_max);
//

void d_create_ipm_hard_core_qp(struct d_ipm_hard_core_qp_workspace *workspace, void *mem);
//
void d_ipm_hard_core_qp(struct d_ipm_hard_core_qp_workspace *workspace);
