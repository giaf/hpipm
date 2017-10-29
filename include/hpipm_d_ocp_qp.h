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

#ifdef __cplusplus
extern "C" {
#endif

struct d_ocp_qp
	{
	struct d_strmat *BAbt;
	struct d_strvec *b;
	struct d_strmat *RSQrq;
	struct d_strvec *rq;
	struct d_strmat *DCt;
	struct d_strvec *d;
	struct d_strvec *Z;
	struct d_strvec *z;
	int *nx; // number of states
	int *nu; // number of inputs
	int *nb; // number of box constraints
	int **idxb; // index of box constraints
	int *ng; // number of general constraints
	int *ns; // number of soft constraints
	int **idxs; // index of soft constraints
	int N; // hotizon lenght
	int memsize; // memory size in bytes
	};



//
int d_memsize_ocp_qp(int N, int *nx, int *nu, int *nb, int *ng, int *ns);
//
void d_create_ocp_qp(int N, int *nx, int *nu, int *nb, int *ng, int *ns, struct d_ocp_qp *qp, void *memory);
//
void d_cvt_colmaj_to_ocp_qp(double **A, double **B, double **b, double **Q, double **S, double **R, double **q, double **r, int **idxb, double **lb, double **ub, double **C, double **D, double **lg, double **ug, double **Zl, double **Zu, double **zl, double **zu, int **idxs, struct d_ocp_qp *qp);
//
void d_cvt_rowmaj_to_ocp_qp(double **A, double **B, double **b, double **Q, double **S, double **R, double **q, double **r, int **idxb, double **lb, double **ub, double **C, double **D, double **lg, double **ug, double **Zl, double **Zu, double **zl, double **zu, int **idxs, struct d_ocp_qp *qp);

#ifdef __cplusplus
}	// #extern "C"
#endif



#endif // HPIPM_D_OCP_QP_H_
