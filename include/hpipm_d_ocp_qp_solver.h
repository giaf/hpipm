/**************************************************************************************************
*                                                                                                 *
* This file is part of HPIPM.                                                                     *
*                                                                                                 *
* HPIPM -- High-Performance Interior Point Method.                                                *
* Copyright (C) 2019 by Gianluca Frison.                                                          *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
*                                                                                                 *
* The 2-Clause BSD License                                                                        *
*                                                                                                 *
* Redistribution and use in source and binary forms, with or without                              *
* modification, are permitted provided that the following conditions are met:                     *
*                                                                                                 *
* 1. Redistributions of source code must retain the above copyright notice, this                  *
*    list of conditions and the following disclaimer.                                             *
* 2. Redistributions in binary form must reproduce the above copyright notice,                    *
*    this list of conditions and the following disclaimer in the documentation                    *
*    and/or other materials provided with the distribution.                                       *
*                                                                                                 *
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND                 *
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED                   *
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE                          *
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR                 *
* ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES                  *
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;                    *
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND                     *
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT                      *
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS                   *
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                                    *
*                                                                                                 *
* Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             *
*                                                                                                 *
**************************************************************************************************/

#ifndef HPIPM_D_OCP_QP_SOLVER_H_
#define HPIPM_D_OCP_QP_SOLVER_H_



#include <blasfeo_target.h>
#include <blasfeo_common.h>

#include <hpipm_common.h>
#include <hpipm_d_ocp_qp_dim.h>
#include <hpipm_d_ocp_qp.h>
#include <hpipm_d_ocp_qp_res.h>
#include <hpipm_d_ocp_qp_sol.h>
#include <hpipm_d_ocp_qp_ipm.h>
#include <hpipm_d_ocp_qp_red.h>
#include <hpipm_d_ocp_qp_seed.h>



#ifdef __cplusplus
extern "C" {
#endif



struct d_ocp_qp_solver_arg
	{
	struct d_ocp_qp_ipm_arg *ipm_arg;
	// dof removal
	struct d_ocp_qp_reduce_eq_dof_arg *red_arg;
	struct d_ocp_qp_dim *red_dim;
	// full condensing
	// partial condensing
	int reduce_eq_dof;
	int valid_red_dim;
	hpipm_size_t memsize;
	};



struct d_ocp_qp_solver_ws
	{
	struct d_ocp_qp_ipm_ws *ipm_ws;
	struct d_ocp_qp_solver_arg *arg; // args to avoid on the fly changes to arbitrary fields
	struct d_ocp_qp_reduce_eq_dof_ws *red_ws;
	struct d_ocp_qp *red_qp;
	struct d_ocp_qp_sol *red_sol;
	struct d_ocp_qp_sol *red_sens;
	struct d_ocp_qp_seed *red_seed;
	hpipm_size_t memsize;
	};



//
hpipm_size_t d_ocp_qp_solver_arg_strsize();
//
hpipm_size_t d_ocp_qp_solver_arg_memsize(struct d_ocp_qp_dim *ocp_dim);
//
void d_ocp_qp_solver_arg_create(struct d_ocp_qp_dim *ocp_dim, struct d_ocp_qp_solver_arg *arg, void *mem);
//
void d_ocp_qp_solver_arg_set_default(enum hpipm_mode mode, struct d_ocp_qp_dim *ocp_dim, struct d_ocp_qp_solver_arg *arg);
//
void d_ocp_qp_solver_arg_set(char *field, void *value, struct d_ocp_qp_solver_arg *arg);
//
void d_ocp_qp_solver_arg_set_reduce_eq_dof(int *value, struct d_ocp_qp_solver_arg *arg);
//
void d_ocp_qp_solver_arg_deepcopy(struct d_ocp_qp_solver_arg *arg_s, struct d_ocp_qp_solver_arg *arg_d);
//
void d_ocp_qp_solver_arg_get_reduce_eq_dof(struct d_ocp_qp_solver_arg *arg, int *value);

//
hpipm_size_t d_ocp_qp_solver_ws_strsize();
//
hpipm_size_t d_ocp_qp_solver_ws_memsize(struct d_ocp_qp_dim *ocp_dim, struct d_ocp_qp_solver_arg *arg);
//
void d_ocp_qp_solver_ws_create(struct d_ocp_qp_dim *ocp_dim, struct d_ocp_qp_solver_arg *arg, struct d_ocp_qp_solver_ws *ws, void *mem);
//
void d_ocp_qp_solver_set(char *field, void *value, struct d_ocp_qp_solver_ws *ws); // XXX set selected args that can safely change
//
void d_ocp_qp_solver_set_iter_max(int *value, struct d_ocp_qp_solver_ws *ws);
//
void d_ocp_qp_solver_set_alpha_min(double *value, struct d_ocp_qp_solver_ws *ws);
//
void d_ocp_qp_solver_set_mu0(double *value, struct d_ocp_qp_solver_ws *ws);
//
void d_ocp_qp_solver_set_tol_stat(double *value, struct d_ocp_qp_solver_ws *ws);
//
void d_ocp_qp_solver_set_tol_eq(double *value, struct d_ocp_qp_solver_ws *ws);
//
void d_ocp_qp_solver_set_tol_ineq(double *value, struct d_ocp_qp_solver_ws *ws);
//
void d_ocp_qp_solver_set_tol_comp(double *value, struct d_ocp_qp_solver_ws *ws);
//
void d_ocp_qp_solver_set_pred_corr(int *value, struct d_ocp_qp_solver_ws *ws);
//
void d_ocp_qp_solver_set_split_step(int *value, struct d_ocp_qp_solver_ws *ws);
//
void d_ocp_qp_solver_set_reg_prim(double *value, struct d_ocp_qp_solver_ws *ws);
//
void d_ocp_qp_solver_get(char *field, struct d_ocp_qp_solver_ws *ws, void *value);
//
void d_ocp_qp_solver_get_status(struct d_ocp_qp_solver_ws *ws, int *value);
//
void d_ocp_qp_solver_get_iter(struct d_ocp_qp_solver_ws *ws, int *value);
//
void d_ocp_qp_solver_get_stat_m(struct d_ocp_qp_solver_ws *ws, int *value);
//
void d_ocp_qp_solver_get_stat(struct d_ocp_qp_solver_ws *ws, double **value);
//
void d_ocp_qp_solver_get_reduce_eq_dof(struct d_ocp_qp_solver_ws *ws, int *value);
//
void d_ocp_qp_solver_solve(struct d_ocp_qp *qp, struct d_ocp_qp_sol *qp_sol, struct d_ocp_qp_solver_ws *ws); // XXX no arg
//
void d_ocp_qp_solver_sens_frw(struct d_ocp_qp *qp, struct d_ocp_qp_seed *qp_seed, struct d_ocp_qp_sol *qp_sens, struct d_ocp_qp_solver_ws *ws);


#ifdef __cplusplus
}	// #extern "C"
#endif


#endif // HPIPM_D_OCP_QP_SOLVER_H_
