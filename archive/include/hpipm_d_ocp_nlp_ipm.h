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


struct d_ocp_nlp_ipm_arg
	{
	struct d_erk_arg *erk_arg;
	double mu0; // initial value for duality measure
	double alpha_min; // exit cond on step length
	double nlp_res_g_max; // exit cond on inf norm of residuals
	double nlp_res_b_max; // exit cond on inf norm of residuals
	double nlp_res_d_max; // exit cond on inf norm of residuals
	double nlp_res_m_max; // exit cond on inf norm of residuals
	int nlp_iter_max; // exit cond in iter number
	int stat_max; // iterations saved in stat
	int pred_corr; // use Mehrotra's predictor-corrector IPM algirthm
	hpipm_size_t memsize;
	};



struct d_ocp_nlp_ipm_workspace
	{
	struct d_ocp_qp *qp;
	struct d_ocp_qp_sol *qp_sol;
	struct d_ocp_qp_ipm_workspace *ipm_workspace;
	struct d_erk_workspace *erk_workspace;
	double nlp_res_g; // exit inf norm of residuals
	double nlp_res_b; // exit inf norm of residuals
	double nlp_res_d; // exit inf norm of residuals
	double nlp_res_m; // exit inf norm of residuals
	int iter; // iteration number
	hpipm_size_t memsize;
	};



//
hpipm_size_t d_memsize_ocp_nlp_ipm_arg(struct d_ocp_nlp *nlp);
//
void d_create_ocp_nlp_ipm_arg(struct d_ocp_nlp *nlp, struct d_ocp_nlp_ipm_arg *arg, void *mem);
//
void d_set_default_ocp_nlp_ipm_arg(struct d_ocp_nlp_ipm_arg *arg);

//
hpipm_size_t d_memsize_ocp_nlp_ipm(struct d_ocp_nlp *nlp, struct d_ocp_nlp_ipm_arg *arg);
//
void d_create_ocp_nlp_ipm(struct d_ocp_nlp *nlp, struct d_ocp_nlp_ipm_arg *arg, struct d_ocp_nlp_ipm_workspace *ws, void *mem);
//
int d_solve_ocp_nlp_ipm(struct d_ocp_nlp *nlp, struct d_ocp_nlp_sol *nlp_sol, struct d_ocp_nlp_ipm_arg *arg, struct d_ocp_nlp_ipm_workspace *ws);

#ifdef __cplusplus
} /* extern "C" */
#endif
