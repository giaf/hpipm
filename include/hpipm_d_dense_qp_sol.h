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



#ifndef HPIPM_D_DENSE_QP_SOL_H_
#define HPIPM_D_DENSE_QP_SOL_H_



#include <blasfeo_target.h>
#include <blasfeo_common.h>

#include "hpipm_d_dense_qp_dim.h"



#ifdef __cplusplus
extern "C" {
#endif



struct d_dense_qp_sol
	{
	struct d_dense_qp_dim *dim;
	struct d_strvec *v;
	struct d_strvec *pi;
	struct d_strvec *lam;
	struct d_strvec *t;
	void *misc;
	int memsize;
	};



//
int d_memsize_dense_qp_sol(struct d_dense_qp_dim *dim);
//
void d_create_dense_qp_sol(struct d_dense_qp_dim *dim, struct d_dense_qp_sol *qp_sol, void *memory);
//
void d_cvt_dense_qp_sol_to_colmaj(struct d_dense_qp_sol *qp_sol, double *v, double *ls, double *us, double *pi, double *lam_lb, double *lam_ub, double *lam_lg, double *lam_ug, double *lam_ls, double *lam_us);
//
void d_cvt_dense_qp_sol_to_rowmaj(struct d_dense_qp_sol *qp_sol, double *v, double *ls, double *us, double *pi, double *lam_lb, double *lam_ub, double *lam_lg, double *lam_ug, double *lam_ls, double *lam_us);
//
void d_cvt_dense_qp_sol_to_libstr(struct d_dense_qp_sol *qp_sol, struct d_strvec *v, struct d_strvec *ls, struct d_strvec *us, struct d_strvec *pi, struct d_strvec *lam_lb, struct d_strvec *lam_ub, struct d_strvec *lam_lg, struct d_strvec *lam_ug, struct d_strvec *lam_ls, struct d_strvec *lam_us);



#ifdef __cplusplus
} /* extern "C" */
#endif



#endif // HPIPM_D_DENSE_QP_SOL_H_
