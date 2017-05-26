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



struct d_ipm2_hard_revcom_workspace
	{
	double *v; // primal variables
	double *pi; // equality constraints multipliers
	double *lam; // inequality constraints multipliers
	double *lam_lb; // lower bounds multipliers
	double *lam_lg; // lower general constraints multipliers
	double *lam_ub; // upper bounds multipliers
	double *lam_ug; // upper general constraints multipliers
	double *t; // inequality constraints slacks
	double *t_lb; // lower bound slacks
	double *t_lg; // lower general constraints slacks
	double *t_ub; // upper bound slacks
	double *t_ug; // upper general constraints slacks
	double *t_inv; // inverse of inequality constraints slacks
	double *dv; // step in v
	double *dlam; // step in lam
	double *dt; // step in t
	double *Qx; // Hessian update
	double *qx; // gradient update
	double *res_m; // m-residuals
	double *res_d; // d-residuals
	double *Dv; // holds the product D*v
	double alpha; // step length
	int nv; // number of primal variables
	int ne; // number of equality constraints
	int nb; // number of two-sized bounds
	int ng; // number of two-sized constraints
	int status; // status
	int memsize; // memory size (in bytes) of workspace
	};
	
