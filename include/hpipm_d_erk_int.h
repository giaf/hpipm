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



struct d_erk_arg
	{
	double h; // step size
	int steps; // number of steps
	int adj_sens; // compute adjoint sensitivities
	};



struct d_erk_workspace
	{
	void (*ode)(int t, double *x, double *p, void *ode_args, double *xdot); // function pointer to vde
	void *ode_args; // pointer to ode args
	struct d_rk_data *rk_data; // integrator data
	struct d_erk_arg *erk_arg; // erk arg
	double *K; // internal variables
	double *x; // states and forward sensitivities
	double *p; // parameter
	double *xt; // temporary states and forward sensitivities
	int nx; // number of states
	int nf; // number of forward sensitivities
	int np; // number of parameters
	int memsize; // TODO
	};



//
int d_memsize_erk_int(struct d_rk_data *rk_data, struct d_erk_arg *erk_arg, int nx, int nf, int np);
//
void d_create_erk_int(struct d_rk_data *rk_data, struct d_erk_arg *erk_arg, int nx, int nf, int np, struct d_erk_workspace *workspace, void *memory);
//
void d_init_erk_int(double *x0, double *fs0, double *p0, void (*ode)(int t, double *x, double *p, void *ode_args, double *xdot), void *ode_args, struct d_erk_workspace *ws);
//
void d_update_p_erk_int(double *p0, struct d_erk_workspace *ws);
//
void d_erk_int(struct d_erk_workspace *workspace);
