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



struct s_ipm_hard_core_qp_workspace
	{
	float *v; // primal variables
	float *pi; // equality constraints multipliers
	float *lam; // inequality constraints multipliers
	float *t; // inequality constraints slacks
	float *t_inv; // inverse of t
	float *dv; // step in v
	float *dpi; // step in pi
	float *dlam; // step in lam
	float *dt; // step in t
	float *res_g; // q-residuals
	float *res_b; // b-residuals
	float *res_d; // d-residuals
	float *res_m; // m-residuals
	float *Gamma; // Hessian update
	float *gamma; // gradient update
	float *Qx; // Hessian update
	float *qx; // gradient update
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
	int nc; // number of (two-sided) inequality constraints
	int iter_max; // exit cond on iter mumber
	int memsize; // memory size (in bytes) of workspace
	};



//
int s_memsize_ipm_hard_core_qp(int nv, int ne, int nc, int iter_max);
//
void s_create_ipm_hard_core_qp(struct s_ipm_hard_core_qp_workspace *workspace, void *mem);
//
void s_ipm_hard_core_qp(struct s_ipm_hard_core_qp_workspace *workspace);

