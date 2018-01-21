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

#ifndef HPIPM_S_OCP_QP_H_
#define HPIPM_S_OCP_QP_H_



#include <blasfeo_target.h>
#include <blasfeo_common.h>

#include "hpipm_s_ocp_qp_dim.h"



#ifdef __cplusplus
extern "C" {
#endif



struct s_ocp_qp
	{
	struct s_ocp_qp_dim *dim;
	struct blasfeo_smat *BAbt;
	struct blasfeo_svec *b;
	struct blasfeo_smat *RSQrq;
	struct blasfeo_svec *rq;
	struct blasfeo_smat *DCt;
	struct blasfeo_svec *d;
	struct blasfeo_svec *Z;
	struct blasfeo_svec *z;
	int **idxb; // index of box constraints
	int **idxs; // index of soft constraints
	int memsize; // memory size in bytes
	};



//
int s_memsize_ocp_qp(struct s_ocp_qp_dim *dim);
//
void s_create_ocp_qp(struct s_ocp_qp_dim *dim, struct s_ocp_qp *qp, void *memory);
//
void s_cvt_colmaj_to_ocp_qp(float **A, float **B, float **b, float **Q, float **S, float **R, float **q, float **r, int **idxb, float **lb, float **ub, float **C, float **D, float **lg, float **ug, float **Zl, float **Zu, float **zl, float **zu, int **idxs, struct s_ocp_qp *qp);
//
void s_cvt_rowmaj_to_ocp_qp(float **A, float **B, float **b, float **Q, float **S, float **R, float **q, float **r, int **idxb, float **lb, float **ub, float **C, float **D, float **lg, float **ug, float **Zl, float **Zu, float **zl, float **zu, int **idxs, struct s_ocp_qp *qp);
//
void s_update_Q(int stage, float *mat, struct s_ocp_qp *qp);
//
void s_update_S(int stage, float *mat, struct s_ocp_qp *qp);
//
void s_update_R(int stage, float *mat, struct s_ocp_qp *qp);
//
void s_update_q(int stage, float *vec, struct s_ocp_qp *qp);
//
void s_update_r(int stage, float *vec, struct s_ocp_qp *qp);
//
void s_update_A(int stage, float *mat, struct s_ocp_qp *qp);
//
void s_update_B(int stage, float *mat, struct s_ocp_qp *qp);
//
void s_update_b(int stage, float *vec, struct s_ocp_qp *qp);
//
void s_update_lbx(int stage, float *vec, struct s_ocp_qp *qp);
//
void s_update_ubx(int stage, float *vec, struct s_ocp_qp *qp);
//
void s_update_lbu(int stage, float *vec, struct s_ocp_qp *qp);
//
void s_update_ubu(int stage, float *vec, struct s_ocp_qp *qp);
//
void s_update_C(int stage, float *mat, struct s_ocp_qp *qp);
//
void s_update_D(int stage, float *mat, struct s_ocp_qp *qp);
//
void s_update_lg(int stage, float *vec, struct s_ocp_qp *qp);
//
void s_update_ug(int stage, float *vec, struct s_ocp_qp *qp);

#ifdef __cplusplus
} /* extern "C" */
#endif



#endif // HPIPM_S_OCP_QP_H_
