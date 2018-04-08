
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



#ifndef HPIPM_S_DENSE_QP_H_
#define HPIPM_S_DENSE_QP_H_



#include <blasfeo_target.h>
#include <blasfeo_common.h>

#include "hpipm_s_dense_qp_dim.h"



#ifdef __cplusplus
extern "C" {
#endif



struct s_dense_qp
	{
	struct s_dense_qp_dim *dim;
	struct blasfeo_smat *Hv; // hessian & gradient
	struct blasfeo_smat *A; // dynamics matrix
	struct blasfeo_smat *Ct; // constraints matrix
	struct blasfeo_svec *g; // gradient
	struct blasfeo_svec *b; // dynamics vector
	struct blasfeo_svec *d; // constraints vector
	struct blasfeo_svec *m; // rhs of complementarity condition
	struct blasfeo_svec *Z; // (diagonal) hessian of slacks
	struct blasfeo_svec *z; // gradient of slacks
	int *idxb; // index of box constraints
	int *idxs; // index of soft constraints
	int memsize; // memory size in bytes
	};



//
int s_memsize_dense_qp(struct s_dense_qp_dim *dim);
//
void s_create_dense_qp(struct s_dense_qp_dim *dim, struct s_dense_qp *qp, void *memory);
//
void s_cvt_colmaj_to_dense_qp(float *H, float *g, float *A, float *b, int *idxb, float *d_lb, float *d_ub, float *C, float *d_lg, float *d_ug, float *Zl, float *Zu, float *zl, float *zu, int *idxs, float *d_ls, float *d_us, struct s_dense_qp *qp);
//
void s_cvt_dense_qp_to_colmaj(struct s_dense_qp *qp, float *H, float *g, float *A, float *b, int *idxb, float *d_lb, float *d_ub, float *C, float *d_lg, float *d_ug, float *Zl, float *Zu, float *zl, float *zu, int *idxs, float *d_ls, float *d_us);
//
void s_cvt_rowmaj_to_dense_qp(float *H, float *g, float *A, float *b, int *idxb, float *d_lb, float *d_ub, float *C, float *d_lg, float *d_ug, float *Zl, float *Zu, float *zl, float *zu, int *idxs, float *d_ls, float *d_us, struct s_dense_qp *qp);
//
void s_cvt_dense_qp_to_rowmaj(struct s_dense_qp *qp, float *H, float *g, float *A, float *b, int *idxb, float *d_lb, float *d_ub, float *C, float *d_lg, float *d_ug, float *Zl, float *Zu, float *zl, float *zu, int *idxs, float *d_ls, float *d_us);



#ifdef __cplusplus
} /* extern "C" */
#endif



#endif // HPIPM_S_DENSE_QP_H_
