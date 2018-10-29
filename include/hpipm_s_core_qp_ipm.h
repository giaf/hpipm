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

#ifndef HPIPM_S_CORE_QP_IPM_
#define HPIPM_S_CORE_QP_IPM_

#ifdef __cplusplus
extern "C" {
#endif

struct s_core_qp_ipm_workspace
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
	float *res_m_bkp; // m-residuals
	float *Gamma; // Hessian update
	float *gamma; // gradient update
	float alpha_prim; // step length
	float alpha_dual; // step length
	float alpha; // step length
	float sigma; // centering XXX
	float mu; // duality measuere
	float mu_aff; // affine duality measuere
	float nc_inv; // 1.0/nt, where nt is the total number of constraints
	float lam_min; // min value in t vector
	float t_min; // min value in lam vector
	int nv; // number of primal variables
	int ne; // number of equality constraints
	int nc; // number of (two-sided) inequality constraints
	int memsize; // memory size (in bytes) of workspace
	};



//
int s_memsize_core_qp_ipm(int nv, int ne, int nc);
//
void s_create_core_qp_ipm(int nv, int ne, int nc, struct s_core_qp_ipm_workspace *workspace, void *mem);
//
void s_core_qp_ipm(struct s_core_qp_ipm_workspace *workspace);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // HPIPM_S_CORE_QP_IPM_
