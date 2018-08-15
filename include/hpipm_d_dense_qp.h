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
	struct blasfeo_dmat *Hv; // hessian of cost & vector work space
	struct blasfeo_dmat *A; // equality constraint matrix
	struct blasfeo_dmat *Ct; // inequality constraints matrix
	struct blasfeo_dvec *gz; // gradient of cost & gradient of slacks
	struct blasfeo_dvec *b; // equality constraint vector
	struct blasfeo_dvec *d; // inequality constraints vector
	struct blasfeo_dvec *m; // rhs of complementarity condition
	struct blasfeo_dvec *Z; // (diagonal) hessian of slacks
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

// setters (COLMAJ)

void d_dense_qp_set_H(double *H, struct d_dense_qp *qp);
//
void d_dense_qp_set_g(double *g, struct d_dense_qp *qp);
//
void d_dense_qp_set_A(double *A, struct d_dense_qp *qp);
//
void d_dense_qp_set_b(double *b, struct d_dense_qp *qp);
//
void d_dense_qp_set_idxb(int *idxb, struct d_dense_qp *qp);
//
void d_dense_qp_set_lb(double *lb, struct d_dense_qp *qp);
//
void d_dense_qp_set_ub(double *ub, struct d_dense_qp *qp);
//
void d_dense_qp_set_C(double *C, struct d_dense_qp *qp);
//
void d_dense_qp_set_lg(double *lg, struct d_dense_qp *qp);
//
void d_dense_qp_set_ug(double *ug, struct d_dense_qp *qp);
//
void d_dense_qp_set_idxs(int *idxs, struct d_dense_qp *qp);
//
void d_dense_qp_set_Zl(double *Zl, struct d_dense_qp *qp);
//
void d_dense_qp_set_Zu(double *Zu, struct d_dense_qp *qp);
//
void d_dense_qp_set_zl(double *zl, struct d_dense_qp *qp);
//
void d_dense_qp_set_zu(double *zu, struct d_dense_qp *qp);
//
void d_dense_qp_set_ls(double *ls, struct d_dense_qp *qp);
//
void d_dense_qp_set_us(double *us, struct d_dense_qp *qp);

// getters (COLMAJ)

void d_dense_qp_get_H(struct d_dense_qp *qp, double *H);
//
void d_dense_qp_get_g(struct d_dense_qp *qp, double *g);
//
void d_dense_qp_get_A(struct d_dense_qp *qp, double *A);
//
void d_dense_qp_get_b(struct d_dense_qp *qp, double *b);
//
void d_dense_qp_get_idxb(struct d_dense_qp *qp, int *idxb);
//
void d_dense_qp_get_lb(struct d_dense_qp *qp, double *lb);
//
void d_dense_qp_get_ub(struct d_dense_qp *qp, double *ub);
//
void d_dense_qp_get_C(struct d_dense_qp *qp, double *C);
//
void d_dense_qp_get_lg(struct d_dense_qp *qp, double *lg);
//
void d_dense_qp_get_ug(struct d_dense_qp *qp, double *ug);
//
void d_dense_qp_get_idxs(struct d_dense_qp *qp, int *idxs);
//
void d_dense_qp_get_Zl(struct d_dense_qp *qp, double *Zl);
//
void d_dense_qp_get_Zu(struct d_dense_qp *qp, double *Zu);
//
void d_dense_qp_get_zl(struct d_dense_qp *qp, double *zl);
//
void d_dense_qp_get_zu(struct d_dense_qp *qp, double *zu);
//
void d_dense_qp_get_ls(struct d_dense_qp *qp, double *ls);
//
void d_dense_qp_get_us(struct d_dense_qp *qp, double *us);


#ifdef __cplusplus
} /* extern "C" */
#endif



#endif // HPIPM_D_DENSE_QP_H_
