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
	struct blasfeo_svec *gz; // gradient & gradient of slacks
	struct blasfeo_svec *b; // dynamics vector
	struct blasfeo_svec *d; // constraints vector
	struct blasfeo_svec *m; // rhs of complementarity condition
	struct blasfeo_svec *Z; // (diagonal) hessian of slacks
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


// setters (COLMAJ)

void s_dense_qp_set_H(float *H, struct s_dense_qp *qp);
//
void s_dense_qp_set_g(float *g, struct s_dense_qp *qp);
//
void s_dense_qp_set_A(float *A, struct s_dense_qp *qp);
//
void s_dense_qp_set_b(float *b, struct s_dense_qp *qp);
//
void s_dense_qp_set_idxb(int *idxb, struct s_dense_qp *qp);
//
void s_dense_qp_set_lb(float *lb, struct s_dense_qp *qp);
//
void s_dense_qp_set_ub(float *ub, struct s_dense_qp *qp);
//
void s_dense_qp_set_C(float *C, struct s_dense_qp *qp);
//
void s_dense_qp_set_lg(float *lg, struct s_dense_qp *qp);
//
void s_dense_qp_set_ug(float *ug, struct s_dense_qp *qp);
//
void s_dense_qp_set_idxs(int *idxs, struct s_dense_qp *qp);
//
void s_dense_qp_set_Zl(float *Zl, struct s_dense_qp *qp);
//
void s_dense_qp_set_Zu(float *Zu, struct s_dense_qp *qp);
//
void s_dense_qp_set_zl(float *zl, struct s_dense_qp *qp);
//
void s_dense_qp_set_zu(float *zu, struct s_dense_qp *qp);
//
void s_dense_qp_set_ls(float *ls, struct s_dense_qp *qp);
//
void s_dense_qp_set_us(float *us, struct s_dense_qp *qp);

// getters (COLMAJ)

void s_dense_qp_get_H(struct s_dense_qp *qp, float *H);
//
void s_dense_qp_get_g(struct s_dense_qp *qp, float *g);
//
void s_dense_qp_get_A(struct s_dense_qp *qp, float *A);
//
void s_dense_qp_get_b(struct s_dense_qp *qp, float *b);
//
void s_dense_qp_get_idxb(struct s_dense_qp *qp, int *idxb);
//
void s_dense_qp_get_lb(struct s_dense_qp *qp, float *lb);
//
void s_dense_qp_get_ub(struct s_dense_qp *qp, float *ub);
//
void s_dense_qp_get_C(struct s_dense_qp *qp, float *C);
//
void s_dense_qp_get_lg(struct s_dense_qp *qp, float *lg);
//
void s_dense_qp_get_ug(struct s_dense_qp *qp, float *ug);
//
void s_dense_qp_get_idxs(struct s_dense_qp *qp, int *idxs);
//
void s_dense_qp_get_Zl(struct s_dense_qp *qp, float *Zl);
//
void s_dense_qp_get_Zu(struct s_dense_qp *qp, float *Zu);
//
void s_dense_qp_get_zl(struct s_dense_qp *qp, float *zl);
//
void s_dense_qp_get_zu(struct s_dense_qp *qp, float *zu);
//
void s_dense_qp_get_ls(struct s_dense_qp *qp, float *ls);
//
void s_dense_qp_get_us(struct s_dense_qp *qp, float *us);


#ifdef __cplusplus
} /* extern "C" */
#endif



#endif // HPIPM_S_DENSE_QP_H_
