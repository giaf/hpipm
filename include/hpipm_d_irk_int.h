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



struct d_irk_workspace
	{
	void (*d_res_impl_vde)(int t, double *xdot, double *x, double *p, void *ode_args, double *res); // function pointer to residuals of implicit vde
	void (*d_jac_impl_ode)(int t, double *xdot, double *x, double *p, void *ode_args, double *jac); // function pointer to jacobian of implicit ode
	void *ode_args; // pointer to ode args
	struct d_rk_data *rk_data; // integrator data
	struct d_strmat *JG; // jacobian of G
	struct d_strmat *rG; // residuals of G
	struct d_strmat *K; // internal variables
	double *x; // states and forward sensitivities
	double *p; // parameter
	double *xt0; // temporary states and forward sensitivities
	double *xt1; // temporary states and forward sensitivities
	double *Kt; // temporary internal variables
	double *rGt; // temporary residuals of G
	double *Jt0; // temporary Jacobian of ode
	double *Jt1; // temporary Jacobian of ode
	int *ipiv; // index of pivot vector
	int nx; // number of states
	int nf; // number of forward sensitivities
	int np; // number of parameters
	int memsize; // TODO
	};



struct d_irk_args
	{
	double h; // step size
	int steps; // number of steps
	int newton_iter;
	};



//
int d_memsize_irk_int(struct d_rk_data *rk_data, int nx, int nf, int np);
//
void d_create_irk_int(struct d_rk_data *rk_data, int nx, int nf, int np, struct d_irk_workspace *workspace, void *memory);
//
void d_init_irk_int(double *x0, double *fs0, double *p0, void (*d_res_impl_vde)(int t, double *xdot, double *x, double *p, void *ode_args, double *res), void (*d_jac_impl_ode)(int t, double *xdot, double *x, double *p, void *ode_args, double *jac), void *ode_args, struct d_irk_workspace *ws);
//
void d_update_p_irk_int(double *p0, struct d_irk_workspace *ws);
//
void d_irk_int(struct d_irk_args *irk_args, struct d_irk_workspace *workspace);

