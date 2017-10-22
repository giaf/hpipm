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



#ifndef HPIPM_S_DENSE_QP_SOL_H_
#define HPIPM_S_DENSE_QP_SOL_H_



#include <blasfeo_target.h>
#include <blasfeo_common.h>



#ifdef __cplusplus
extern "C" {
#endif



struct s_dense_qp_sol
	{
	struct s_strvec *v;
	struct s_strvec *pi;
	struct s_strvec *lam;
	struct s_strvec *t;
	int memsize;
	};



//
int s_memsize_dense_qp_sol(int nv, int ne, int nb, int ng, int ns);
//
void s_create_dense_qp_sol(int nv, int ne, int nb, int ng, int ns, struct s_dense_qp_sol *qp_sol, void *memory);
//
void s_cvt_dense_qp_sol_to_colmaj(struct s_dense_qp *qp, struct s_dense_qp_sol *qp_sol, float *v, float *ls, float *us, float *pi, float *lam_lb, float *lam_ub, float *lam_lg, float *lam_ug, float *lam_ls, float *lam_us);
//
void s_cvt_dense_qp_sol_to_rowmaj(struct s_dense_qp *qp, struct s_dense_qp_sol *qp_sol, float *v, float *ls, float *us, float *pi, float *lam_lb, float *lam_ub, float *lam_lg, float *lam_ug, float *lam_ls, float *lam_us);
//
void s_cvt_dense_qp_sol_to_libstr(struct s_dense_qp *qp, struct s_dense_qp_sol *qp_sol, struct s_strvec *v, struct s_strvec *ls, struct s_strvec *us, struct s_strvec *pi, struct s_strvec *lam_lb, struct s_strvec *lam_ub, struct s_strvec *lam_lg, struct s_strvec *lam_ug, struct s_strvec *lam_ls, struct s_strvec *lam_us);



#ifdef __cplusplus
} /* extern "C" */
#endif



#endif // HPIPM_S_DENSE_QP_SOL_H_
