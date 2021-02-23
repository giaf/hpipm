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



#include <blasfeo_target.h>
#include <blasfeo_common.h>

#ifdef __cplusplus
extern "C" {
#endif

struct d_ocp_nlp_model
	{
	void (*expl_ode)(int t, double *x, double *p, void *vde_arg, double *xdot); // explicit ode
	void (*expl_vde_for)(int t, double *x, double *p, void *vde_arg, double *xdot); // explicit forward vde
	void (*expl_vde_adj)(int t, double *adj_in, void *vde_arg, double *adj_out); // explicit adjoint vde
	double *forward_seed; // seeds to initialize forward sensitivities
	void *arg; // arguments to vde
	};



struct d_ocp_nlp
	{
	struct d_ocp_nlp_model *model;
	struct blasfeo_dmat *RSQ; // hessian
	struct blasfeo_dvec *rq; // gradient
	struct blasfeo_dmat *DCt; // constraints matrix
	struct blasfeo_dvec *d; // constraints rhs
	struct blasfeo_dvec *Z; // soft constr slack diag hessian
	struct blasfeo_dvec *z; // soft constr slack grad
	struct blasfeo_dvec *tmp_nuxM; // work space of size max(nu+nx)
	int *nx; // number of states
	int *nu; // number of inputs
	int *nb; // number of box constraints
	int **idxb; // index of box constraints
	int *ng; // number of general constraints
	int *ns; // number of soft constraints
	int **idxs; // index of soft constraints
	int N; // hotizon lenght
	hpipm_size_t memsize; // memory size in bytes
	};



//
hpipm_size_t d_memsize_ocp_nlp(int N, int *nx, int *nu, int *nb, int *ng, int *ns);
//
void d_create_ocp_nlp(int N, int *nx, int *nu, int *nb, int *ng, int *ns, struct d_ocp_nlp *nlp, void *memory);
//
void d_cvt_colmaj_to_ocp_nlp(struct d_ocp_nlp_model *model, double **Q, double **S, double **R, double **x_ref, double **u_ref, int **idxb, double **d_lb, double **d_ub, double **C, double **D, double **d_lg, double **d_ug, double **Zl, double **Zu, double **zl, double **zu, int **idxs, struct d_ocp_nlp *nlp);

#ifdef __cplusplus
} /* extern "C" */
#endif
