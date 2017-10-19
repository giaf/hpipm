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



#include <blasfeo_target.h>
#include <blasfeo_common.h>

#ifdef __cplusplus
extern "C" {
#endif


struct s_ocp_qp_ipm_arg
	{
	float mu0; // initial value for duality measure
	float alpha_min; // exit cond on step length
	float res_g_max; // exit cond on inf norm of residuals
	float res_b_max; // exit cond on inf norm of residuals
	float res_d_max; // exit cond on inf norm of residuals
	float res_m_max; // exit cond on inf norm of residuals
	int iter_max; // exit cond in iter number
	int stat_max; // iterations saved in stat
	int pred_corr; // use Mehrotra's predictor-corrector IPM algirthm
	int warm_start; // 0 no warm start, 1 warm start primal sol
	int memsize;
	};



struct s_ocp_qp_res
	{
	struct s_strvec *res_g; // q-residuals
	struct s_strvec *res_b; // b-residuals
	struct s_strvec *res_d; // d-residuals
	struct s_strvec *res_m; // m-residuals
	struct s_strvec *tmp_nbgM; // work space of size nbM+ngM
	struct s_strvec *tmp_nsM; // work space of size nsM
	float res_mu; // mu-residual
	int memsize;
	};



struct s_ocp_qp_ipm_workspace
	{
	struct s_core_qp_ipm_workspace *core_workspace;
	struct s_ocp_qp_res *res_workspace;
	struct s_strvec *dux;
	struct s_strvec *dpi;
	struct s_strvec *dt;
	struct s_strvec *Gamma; // hessian update
	struct s_strvec *gamma; // hessian update
	struct s_strvec *tmp_nxM; // work space of size nxM
	struct s_strvec *tmp_nbgM; // work space of size nbM+ngM
	struct s_strvec *tmp_nsM; // work space of size nsM
	struct s_strvec *Pb; // Pb
	struct s_strvec *Zs_inv;
	struct s_strmat *L;
	struct s_strmat *AL;
	float *stat; // convergence statistics
	float qp_res[4]; // infinity norm of residuals
	float mu0; // mu0
	int iter; // iteration number
	int stat_max; // iterations saved in stat
	int warm_start; // 0 no warm start, 1 warm start primal sol
	int memsize;
	};



//
int s_memsize_ocp_qp_ipm_arg(struct s_ocp_qp *qp);
//
void s_create_ocp_qp_ipm_arg(struct s_ocp_qp *qp, struct s_ocp_qp_ipm_arg *arg, void *mem);
//
void s_set_default_ocp_qp_ipm_arg(struct s_ocp_qp_ipm_arg *arg);

//
int s_memsize_ocp_qp_res(struct s_ocp_qp *qp);
//
void s_create_ocp_qp_res(struct s_ocp_qp *qp, struct s_ocp_qp_res *workspace, void *mem);

//
int s_memsize_ocp_qp_ipm(struct s_ocp_qp *qp, struct s_ocp_qp_ipm_arg *arg);
//
void s_create_ocp_qp_ipm(struct s_ocp_qp *qp, struct s_ocp_qp_ipm_arg *arg, struct s_ocp_qp_ipm_workspace *ws, void *mem);
//
int s_solve_ocp_qp_ipm(struct s_ocp_qp *qp, struct s_ocp_qp_sol *qp_sol, struct s_ocp_qp_ipm_arg *arg, struct s_ocp_qp_ipm_workspace *ws);
//
//
void s_cvt_ocp_qp_res_to_colmaj(struct s_ocp_qp *qp, struct s_ocp_qp_ipm_workspace *ws, float **res_r, float **res_q, float **res_ls, float **res_us, float **res_b, float **res_d_lb, float **res_d_ub, float **res_d_lg, float **res_d_ug, float **res_d_ls, float **res_d_us, float **res_m_lb, float **res_m_ub, float **res_m_lg, float **res_m_ug, float **res_m_ls, float **res_m_us);
//
void s_cvt_ocp_qp_res_to_rowmaj(struct s_ocp_qp *qp, struct s_ocp_qp_ipm_workspace *ws, float **res_r, float **res_q, float **res_ls, float **res_us, float **res_b, float **res_d_lb, float **res_d_ub, float **res_d_lg, float **res_d_ug, float **res_d_ls, float **res_d_us, float **res_m_lb, float **res_m_ub, float **res_m_lg, float **res_m_ug, float **res_m_ls, float **res_m_us);

#ifdef __cplusplus
} /* extern "C" */
#endif

