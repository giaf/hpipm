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



#ifndef HPIPM_D_DENSE_QP_H_
#define HPIPM_D_DENSE_QP_H_



#include <blasfeo_target.h>
#include <blasfeo_common.h>

#include "hpipm_d_dense_qp_dim.h"



#ifdef __cplusplus
extern "C" {
#endif



struct d_dense_qp
	{
	struct d_dense_qp_dim *dim;
	struct d_strmat *Hg; // hessian & gradient
	struct d_strmat *A; // equality constraint matrix
	struct d_strmat *Ct; // inequality constraints matrix
	struct d_strvec *g; // gradient
	struct d_strvec *b; // equality constraint vector
	struct d_strvec *d; // inequality constraints vector
	struct d_strvec *Z; // (diagonal) hessian of slacks
	struct d_strvec *z; // gradient of slacks
	int *idxb; // index of box constraints
	int *idxs; // index of soft constraints
	int memsize; // memory size in bytes
	};



//
int d_memsize_dense_qp(struct d_dense_qp_dim *dim);
//
void d_create_dense_qp(struct d_dense_qp_dim *dim, struct d_dense_qp *qp, void *memory);
//
void d_cvt_colmaj_to_dense_qp(double *H, double *g, double *A, double *b, int *idxb, double *d_lb, double *d_ub, double *C, double *d_lg, double *d_ug, double *Zl, double *Zu, double *zl, double *zu, int *idxs, struct d_dense_qp *qp);
//
void d_cvt_dense_qp_to_colmaj(struct d_dense_qp *qp, double *H, double *g, double *A, double *b, int *idxb, double *d_lb, double *d_ub, double *C, double *d_lg, double *d_ug, double *Zl, double *Zu, double *zl, double *zu, int *idxs);
//
void d_cvt_rowmaj_to_dense_qp(double *H, double *g, double *A, double *b, int *idxb, double *d_lb, double *d_ub, double *C, double *d_lg, double *d_ug, double *Zl, double *Zu, double *zl, double *zu, int *idxs, struct d_dense_qp *qp);
//
void d_cvt_dense_qp_to_rowmaj(struct d_dense_qp *qp, double *H, double *g, double *A, double *b, int *idxb, double *d_lb, double *d_ub, double *C, double *d_lg, double *d_ug, double *Zl, double *Zu, double *zl, double *zu, int *idxs);
//
void d_cvt_libstr_to_dense_qp(struct d_strmat *H, struct d_strmat *A, struct d_strmat *C, struct d_strvec *g, struct d_strvec *b, struct d_strvec *d_lb, struct d_strvec *d_ub, struct d_strvec *d_lg, struct d_strvec *d_ug, int *idxb, struct d_strvec *Zl, struct d_strvec *Zu, struct d_strvec *zl, struct d_strvec *zu, int *idxs, struct d_dense_qp *qp);
//
void d_cvt_dense_qp_to_libstr(struct d_dense_qp *qp, struct d_strmat *H, struct d_strmat *A, struct d_strmat *C, struct d_strvec *g, struct d_strvec *b, struct d_strvec *d_lb, struct d_strvec *d_ub, struct d_strvec *d_lg, struct d_strvec *d_ug, int *idxb, struct d_strvec *Zl, struct d_strvec *Zu, struct d_strvec *zl, struct d_strvec *zu, int *idxs);



#ifdef __cplusplus
} /* extern "C" */
#endif



#endif // HPIPM_D_DENSE_QP_H_
