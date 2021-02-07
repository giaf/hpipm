/**************************************************************************************************
*                                                                                                 *
* This file is part of HPIPM.                                                                     *
*                                                                                                 *
* HPIPM -- High-Performance Interior Point Method.                                                *
* Copyright (C) 2017-2018 by Gianluca Frison.                                                     *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
*                                                                                                 *
* This program is free software: you can redistribute it and/or modify                            *
* it under the terms of the GNU General Public License as published by                            *
* the Free Software Foundation, either version 3 of the License, or                               *
* (at your option) any later version                                                              *.
*                                                                                                 *
* This program is distributed in the hope that it will be useful,                                 *
* but WITHOUT ANY WARRANTY; without even the implied warranty of                                  *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                   *
* GNU General Public License for more details.                                                    *
*                                                                                                 *
* You should have received a copy of the GNU General Public License                               *
* along with this program.  If not, see <https://www.gnu.org/licenses/>.                          *
*                                                                                                 *
* The authors designate this particular file as subject to the "Classpath" exception              *
* as provided by the authors in the LICENSE file that accompained this code.                      *
*                                                                                                 *
* Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             *
*                                                                                                 *
**************************************************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

struct d_irk_workspace
	{
	void (*d_res_impl_vde)(int t, double *xdot, double *x, double *p, void *ode_args, double *res); // function pointer to residuals of implicit vde
	void (*d_jac_impl_ode)(int t, double *xdot, double *x, double *p, void *ode_args, double *jac); // function pointer to jacobian of implicit ode
	void *ode_args; // pointer to ode args
	struct d_rk_data *rk_data; // integrator data
	struct blasfeo_dmat *JG; // jacobian of G
	struct blasfeo_dmat *rG; // residuals of G
	struct blasfeo_dmat *K; // internal variables
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
	hpipm_size_t memsize; // TODO
	};



struct d_irk_args
	{
	double h; // step size
	int steps; // number of steps
	int newton_iter;
	};



//
hpipm_size_t d_memsize_irk_int(struct d_rk_data *rk_data, int nx, int nf, int np);
//
void d_create_irk_int(struct d_rk_data *rk_data, int nx, int nf, int np, struct d_irk_workspace *workspace, void *memory);
//
void d_init_irk_int(double *x0, double *fs0, double *p0, void (*d_res_impl_vde)(int t, double *xdot, double *x, double *p, void *ode_args, double *res), void (*d_jac_impl_ode)(int t, double *xdot, double *x, double *p, void *ode_args, double *jac), void *ode_args, struct d_irk_workspace *ws);
//
void d_update_p_irk_int(double *p0, struct d_irk_workspace *ws);
//
void d_irk_int(struct d_irk_args *irk_args, struct d_irk_workspace *workspace);

#ifdef __cplusplus
} /* extern "C" */
#endif
