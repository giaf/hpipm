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

#ifndef HPIPM_D_OCP_QP_SOL_H_
#define HPIPM_D_OCP_QP_SOL_H_



#include <blasfeo_target.h>
#include <blasfeo_common.h>

#include "hpipm_d_ocp_qp_dim.h"



#ifdef __cplusplus
extern "C" {
#endif



struct d_ocp_qp_sol
	{
	struct d_ocp_qp_dim *dim;
	struct blasfeo_dvec *ux;
	struct blasfeo_dvec *pi;
	struct blasfeo_dvec *lam;
	struct blasfeo_dvec *t;
	void *misc;
	int memsize; // memory size in bytes
	};



//
int d_memsize_ocp_qp_sol(struct d_ocp_qp_dim *dim);
//
void d_create_ocp_qp_sol(struct d_ocp_qp_dim *dim, struct d_ocp_qp_sol *qp_sol, void *memory);
//
void d_cvt_ocp_qp_sol_to_colmaj(struct d_ocp_qp_sol *qp_sol, double **u, double **x, double **ls, double **us, double **pi, double **lam_lb, double **lam_ub, double **lam_lg, double **lam_ug, double **lam_ls, double **lam_us);
//
void d_cvt_colmaj_to_ocp_qp_sol(double **u, double **x, double **ls, double **us, double **pi, double **lam_lb, double **lam_ub, double **lam_lg, double **lam_ug, double **lam_ls, double **lam_us, struct d_ocp_qp_sol *qp_sol);
//
void d_cvt_ocp_qp_sol_to_rowmaj(struct d_ocp_qp_sol *qp_sol, double **u, double **x, double **ls, double **us, double **pi, double **lam_lb, double **lam_ub, double **lam_lg, double **lam_ug, double **lam_ls, double **lam_us);
//
void d_cvt_ocp_qp_sol_to_libstr(struct d_ocp_qp_sol *qp_sol, struct blasfeo_dvec *u, struct blasfeo_dvec *ls, struct blasfeo_dvec *us, struct blasfeo_dvec *x, struct blasfeo_dvec *pi, struct blasfeo_dvec *lam_lb, struct blasfeo_dvec *lam_ub, struct blasfeo_dvec *lam_lg, struct blasfeo_dvec *lam_ug, struct blasfeo_dvec *lam_ls, struct blasfeo_dvec *lam_us);
//
void d_cvt_ocp_qp_sol_to_colmaj_x(struct d_ocp_qp_sol *qp_sol, double *vec, int stage);
//
void d_cvt_ocp_qp_sol_to_colmaj_u(struct d_ocp_qp_sol *qp_sol, double *vec, int stage);
//
void d_cvt_ocp_qp_sol_to_colmaj_pi(struct d_ocp_qp_sol *qp_sol, double *vec, int stage);
//
void d_cvt_ocp_qp_sol_to_colmaj_lam_lb(struct d_ocp_qp_sol *qp_sol, double *vec, int stage);
//
void d_cvt_ocp_qp_sol_to_colmaj_lam_ub(struct d_ocp_qp_sol *qp_sol, double *vec, int stage);
//
void d_cvt_ocp_qp_sol_to_colmaj_lam_lg(struct d_ocp_qp_sol *qp_sol, double *vec, int stage);
//
void d_cvt_ocp_qp_sol_to_colmaj_lam_ug(struct d_ocp_qp_sol *qp_sol, double *vec, int stage);

#ifdef __cplusplus
}	// #extern "C"
#endif



#endif // HPIPM_D_OCP_QP_SOL_H_
