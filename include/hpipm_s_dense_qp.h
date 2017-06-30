
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



#include <blasfeo_target.h>
#include <blasfeo_common.h>



struct s_dense_qp
	{
	struct s_strmat *Hg; // hessian & gradient
	struct s_strmat *A; // dynamics matrix
	struct s_strmat *Ct; // constraints matrix
	struct s_strvec *g; // gradient
	struct s_strvec *b; // dynamics vector
	struct s_strvec *d; // constraint
	struct s_strvec *d_lb; // lower box constraint
	struct s_strvec *d_ub; // upper box constraint
	struct s_strvec *d_lg; // lower general constraint
	struct s_strvec *d_ug; // upper general constraint
	int *idxb; // index of box constraints
	int nv; // number of variables
	int ne; // number of equality constraints
	int nb; // number of box constraints
	int ng; // number of general constraints
	int memsize; // memory size in bytes
	};



//
int s_memsize_dense_qp(int nv, int ne, int nb, int ng);
//
void s_create_dense_qp(int nv, int ne, int nb, int ng, struct s_dense_qp *qp, void *memory);
//
void s_cvt_colmaj_to_dense_qp(float *H, float *g, float *A, float *b, int *idxb, float *d_lb, float *d_ub, float *C, float *d_lg, float *d_ug, struct s_dense_qp *qp);
//
void s_cvt_dense_qp_to_colmaj(struct s_dense_qp *qp, float *H, float *g, float *A, float *b, int *idxb, float *d_lb, float *d_ub, float *C, float *d_lg, float *d_ug);
//
void s_cvt_rowmaj_to_dense_qp(float *H, float *g, float *A, float *b, int *idxb, float *d_lb, float *d_ub, float *C, float *d_lg, float *d_ug, struct s_dense_qp *qp);
//
void s_cvt_dense_qp_to_rowmaj(struct s_dense_qp *qp, float *H, float *g, float *A, float *b, int *idxb, float *d_lb, float *d_ub, float *C, float *d_lg, float *d_ug);
//
void s_cvt_libstr_to_dense_qp(struct s_strmat *H, struct s_strmat *A, struct s_strmat *C, struct s_strvec *g, struct s_strvec *b, struct s_strvec *d_lb, struct s_strvec *d_ub, struct s_strvec *d_lg, struct s_strvec *d_ug, int *idxb, struct s_dense_qp *qp);
//
void s_cvt_dense_qp_to_libstr(struct s_dense_qp *qp, struct s_strmat *H, struct s_strmat *A, struct s_strmat *C, struct s_strvec *g, struct s_strvec *b, struct s_strvec *d_lb, struct s_strvec *d_ub, struct s_strvec *d_lg, struct s_strvec *d_ug, int *idxb);


