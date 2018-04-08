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
	struct blasfeo_dmat *Hv; // hessian & gradient
	struct blasfeo_dmat *A; // equality constraint matrix
	struct blasfeo_dmat *Ct; // inequality constraints matrix
	struct blasfeo_dvec *g; // gradient
	struct blasfeo_dvec *b; // equality constraint vector
	struct blasfeo_dvec *d; // inequality constraints vector
	struct blasfeo_dvec *m; // rhs of complementarity condition
	struct blasfeo_dvec *Z; // (diagonal) hessian of slacks
	struct blasfeo_dvec *z; // gradient of slacks
	int *idxb; // index of box constraints
	int *idxs; // index of soft constraints
	int memsize; // memory size in bytes
	};



//
int d_memsize_dense_qp(struct d_dense_qp_dim *dim);
//
void d_create_dense_qp(struct d_dense_qp_dim *dim, struct d_dense_qp *qp, void *memory);
//
void d_cvt_colmaj_to_dense_qp(double *H, double *g, double *A, double *b, int *idxb, double *d_lb, double *d_ub, double *C, double *d_lg, double *d_ug, double *Zl, double *Zu, double *zl, double *zu, int *idxs, double *d_ls, double *d_us, struct d_dense_qp *qp);
//
void d_cvt_dense_qp_to_colmaj(struct d_dense_qp *qp, double *H, double *g, double *A, double *b, int *idxb, double *d_lb, double *d_ub, double *C, double *d_lg, double *d_ug, double *Zl, double *Zu, double *zl, double *zu, int *idxs, double *d_ls, double *d_us);
//
void d_cvt_rowmaj_to_dense_qp(double *H, double *g, double *A, double *b, int *idxb, double *d_lb, double *d_ub, double *C, double *d_lg, double *d_ug, double *Zl, double *Zu, double *zl, double *zu, int *idxs, double *d_ls, double *d_us, struct d_dense_qp *qp);
//
void d_cvt_dense_qp_to_rowmaj(struct d_dense_qp *qp, double *H, double *g, double *A, double *b, int *idxb, double *d_lb, double *d_ub, double *C, double *d_lg, double *d_ug, double *Zl, double *Zu, double *zl, double *zu, int *idxs, double *d_ls, double *d_us);



#ifdef __cplusplus
} /* extern "C" */
#endif



#endif // HPIPM_D_DENSE_QP_H_
