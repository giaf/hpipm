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

#ifdef __cplusplus
extern "C" {
#endif

struct d_core_qp_ipm_workspace
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
	double *res_m_bkp; // m-residuals
	double *Gamma; // Hessian update
	double *gamma; // gradient update
	double alpha; // step length
	double alpha_prim; // step length
	double alpha_dual; // step length
	double sigma; // centering XXX
	double mu; // duality measuere
	double mu_aff; // affine duality measuere
	double nc_inv; // 1.0/nt, where nt is the total number of constraints
	double lam_min; // min value in lam vector
	double t_min; // min value in t vector
	int nv; // number of primal variables
	int ne; // number of equality constraints
	int nc; // number of (two-sided) inequality constraints
	int memsize; // memory size (in bytes) of workspace
	};



//
int d_memsize_core_qp_ipm(int nv, int ne, int nc);
//
void d_create_core_qp_ipm(int nv, int ne, int nc, struct d_core_qp_ipm_workspace *workspace, void *mem);
//
void d_core_qp_ipm(struct d_core_qp_ipm_workspace *workspace);

#ifdef __cplusplus
} /* extern "C" */
#endif
