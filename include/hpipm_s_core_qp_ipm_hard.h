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



struct s_ipm_hard_core_qp_workspace
	{
	float *d; // constraints
	float *d_lb; // lower box constraints
	float *d_ub; // upper box constraints
	float *d_lg; // lower general constraints
	float *d_ug; // upper general constraints
	float *v; // primal variables
	float *pi; // equality constraints multipliers
	float *lam; // inequality constraints multipliers
	float *lam_lb; // lower bounds multipliers
	float *lam_lg; // lower general constraints multipliers
	float *lam_ub; // upper bounds multipliers
	float *lam_ug; // upper general constraints multipliers
	float *t; // inequality constraints slacks
	float *t_lb; // lower box constraints slacks
	float *t_lg; // lower general constraints slacks
	float *t_ub; // upper box constraints slacks
	float *t_ug; // upper general constraints slacks
	float *t_inv; // inverse of t
	float *t_inv_lb; // inverse of t
	float *t_inv_ub; // inverse of t
	float *t_inv_lg; // inverse of t
	float *t_inv_ug; // inverse of t
	float *dv; // step in v
	float *dpi; // step in pi
	float *dlam; // step in lam
	float *dlam_lb; //
	float *dlam_lg; //
	float *dlam_ub; //
	float *dlam_ug; //
	float *dt; // step in t
	float *dt_lb; // step in t_lb
	float *dt_ub; // step in t_ub
	float *dt_lg; // step in t_lg
	float *dt_ug; // step in t_ug
	float *res_g; // q-residuals
	float *res_b; // b-residuals
	float *res_d; // d-residuals
	float *res_d_lb; // d-residuals
	float *res_d_ub; // d-residuals
	float *res_d_lg; // d-residuals
	float *res_d_ug; // d-residuals
	float *res_m; // m-residuals
	float *res_m_lb; // m-residuals
	float *res_m_ub; // m-residuals
	float *res_m_lg; // m-residuals
	float *res_m_ug; // m-residuals
	float *Qx; // Hessian update
	float *Qx_lb; // Hessian update
	float *Qx_lg; // Hessian update
	float *qx; // gradient update
	float *qx_lb; // gradient update
	float *qx_lg; // gradient update
	float *stat; // convergence statistics
	float alpha; // step length
	float alpha_min; // exit cond on step lenght
	float sigma; // centering XXX
	float mu; // duality measuere
	float mu_aff; // affine duality measuere
	float mu0; // initial duality measuere
	float mu_max; // exit cond on mu
	float nt_inv; // 1.0/nt, where nt is the total number of constraints
	int nv; // number of primal variables
	int ne; // number of equality constraints
	int nb; // number of two-sized bounds
	int ng; // number of two-sized constraints
	int memsize; // memory size (in bytes) of workspace
	int iter_max; // exit cond on iter mumber
	};



//
int s_memsize_ipm_hard_core_qp(int nv, int ne, int nb, int ng, int iter_max);
//

void s_create_ipm_hard_core_qp(struct s_ipm_hard_core_qp_workspace *workspace, void *mem);
//
void s_ipm_hard_core_qp(struct s_ipm_hard_core_qp_workspace *workspace);

