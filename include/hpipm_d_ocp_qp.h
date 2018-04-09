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

#ifndef HPIPM_D_OCP_QP_H_
#define HPIPM_D_OCP_QP_H_



#include <blasfeo_target.h>
#include <blasfeo_common.h>

#include "hpipm_d_ocp_qp_dim.h"



#ifdef __cplusplus
extern "C" {
#endif



struct d_ocp_qp
	{
	struct d_ocp_qp_dim *dim;
	struct blasfeo_dmat *BAbt;
	struct blasfeo_dvec *b;
	struct blasfeo_dmat *RSQrq;
	struct blasfeo_dvec *rqz;
	struct blasfeo_dmat *DCt;
	struct blasfeo_dvec *d;
	struct blasfeo_dvec *Z;
	int **idxb; // index of box constraints
	int **idxs; // index of soft constraints
	int memsize; // memory size in bytes
	};



//
int d_memsize_ocp_qp(struct d_ocp_qp_dim *dim);
//
void d_create_ocp_qp(struct d_ocp_qp_dim *dim, struct d_ocp_qp *qp, void *memory);
//
void d_change_bounds_dimensions_ocp_qp(int *nbu, int *nbx, struct d_ocp_qp *qp);
//
void d_cvt_colmaj_to_ocp_qp(double **A, double **B, double **b, double **Q, double **S, double **R, double **q, double **r, int **idxb, double **lb, double **ub, double **C, double **D, double **lg, double **ug, double **Zl, double **Zu, double **zl, double **zu, int **idxs, double **ls, double **us, struct d_ocp_qp *qp);
//
void d_cvt_rowmaj_to_ocp_qp(double **A, double **B, double **b, double **Q, double **S, double **R, double **q, double **r, int **idxb, double **lb, double **ub, double **C, double **D, double **lg, double **ug, double **Zl, double **Zu, double **zl, double **zu, int **idxs, double **ls, double **us, struct d_ocp_qp *qp);
//
void d_cvt_colmaj_to_ocp_qp_Q(int stage, double *mat, struct d_ocp_qp *qp);
//
void d_cvt_ocp_qp_to_colmaj_Q(int stage, struct d_ocp_qp *qp, double *mat);
//
void d_cvt_colmaj_to_ocp_qp_S(int stage, double *mat, struct d_ocp_qp *qp);
//
void d_cvt_ocp_qp_to_colmaj_S(int stage, struct d_ocp_qp *qp, double *mat);
//
void d_cvt_colmaj_to_ocp_qp_R(int stage, double *mat, struct d_ocp_qp *qp);
//
void d_cvt_ocp_qp_to_colmaj_R(int stage, struct d_ocp_qp *qp, double *mat);
//
void d_cvt_colmaj_to_ocp_qp_q(int stage, double *vec, struct d_ocp_qp *qp);
//
void d_cvt_ocp_qp_to_colmaj_q(int stage, struct d_ocp_qp *qp, double *vec);
//
void d_cvt_colmaj_to_ocp_qp_r(int stage, double *vec, struct d_ocp_qp *qp);
//
void d_cvt_ocp_qp_to_colmaj_r(int stage, struct d_ocp_qp *qp, double *vec);
//
void d_cvt_colmaj_to_ocp_qp_A(int stage, double *mat, struct d_ocp_qp *qp);
//
void d_cvt_ocp_qp_to_colmaj_A(int stage, struct d_ocp_qp *qp, double *mat);
//
void d_cvt_colmaj_to_ocp_qp_B(int stage, double *mat, struct d_ocp_qp *qp);
//
void d_cvt_ocp_qp_to_colmaj_B(int stage, struct d_ocp_qp *qp, double *mat);
//
void d_cvt_colmaj_to_ocp_qp_b(int stage, double *vec, struct d_ocp_qp *qp);
//
void d_cvt_ocp_qp_to_colmaj_b(int stage, struct d_ocp_qp *qp, double *vec);
//
void d_cvt_colmaj_to_ocp_qp_lbx(int stage, double *vec, struct d_ocp_qp *qp);
//
void d_cvt_ocp_qp_to_colmaj_lbx(int stage, struct d_ocp_qp *qp, double *vec);
//
void d_cvt_colmaj_to_ocp_qp_ubx(int stage, double *vec, struct d_ocp_qp *qp);
//
void d_cvt_ocp_qp_to_colmaj_ubx(int stage, struct d_ocp_qp *qp, double *vec);
//
void d_cvt_colmaj_to_ocp_qp_lbu(int stage, double *vec, struct d_ocp_qp *qp);
//
void d_cvt_ocp_qp_to_colmaj_lbu(int stage, struct d_ocp_qp *qp, double *vec);
//
void d_cvt_colmaj_to_ocp_qp_ubu(int stage, double *vec, struct d_ocp_qp *qp);
//
void d_cvt_ocp_qp_to_colmaj_ubu(int stage, struct d_ocp_qp *qp, double *vec);
//
void d_cvt_colmaj_to_ocp_qp_C(int stage, double *mat, struct d_ocp_qp *qp);
//
void d_cvt_ocp_qp_to_colmaj_C(int stage, struct d_ocp_qp *qp, double *mat);
//
void d_cvt_colmaj_to_ocp_qp_D(int stage, double *mat, struct d_ocp_qp *qp);
//
void d_cvt_ocp_qp_to_colmaj_D(int stage, struct d_ocp_qp *qp, double *mat);
//
void d_cvt_colmaj_to_ocp_qp_lg(int stage, double *vec, struct d_ocp_qp *qp);
//
void d_cvt_ocp_qp_to_colmaj_lg(int stage, struct d_ocp_qp *qp, double *vec);
//
void d_cvt_colmaj_to_ocp_qp_ug(int stage, double *vec, struct d_ocp_qp *qp);
//
void d_cvt_ocp_qp_to_colmaj_ug(int stage, struct d_ocp_qp *qp, double *vec);


#ifdef __cplusplus
}	// #extern "C"
#endif



#endif // HPIPM_D_OCP_QP_H_
