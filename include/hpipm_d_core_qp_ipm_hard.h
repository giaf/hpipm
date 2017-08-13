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



struct d_ipm_hard_core_qp_workspace
	{
	double *v; // primal variables
	double *pi; // equality constraints multipliers
	double *lam; // inequality constraints multipliers
	double *t; // inequality constraints slacks
	double *t_inv; // inverse of t
	double *dv; // step in v
	double *dpi; // step in pi
	double *dlam; // step in lam
	double *dt; // step in t
	double *res_g; // q-residuals
	double *res_b; // b-residuals
	double *res_d; // d-residuals
	double *res_m; // m-residuals
	double *Qx; // Hessian update
	double *qx; // gradient update
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
