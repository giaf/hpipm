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
	struct blasfeo_dmat *BAbt; // dynamics matrix & vector work space
	struct blasfeo_dmat *RSQrq; // hessian of cost & vector work space
	struct blasfeo_dmat *DCt; // inequality constraints matrix
	struct blasfeo_dvec *b; // dynamics vector
	struct blasfeo_dvec *rqz; // gradient of cost & gradient of slacks
	struct blasfeo_dvec *d; // inequality constraints vector
	struct blasfeo_dvec *m; // rhs of complementarity condition
	struct blasfeo_dvec *Z; // (diagonal) hessian of slacks
	int **idxb; // index of box constraints
	int **idxs; // index of soft constraints
	int memsize; // memory size in bytes
	};



//
int d_ocp_qp_strsize();
//
int d_ocp_qp_memsize(struct d_ocp_qp_dim *dim);
//
void d_ocp_qp_create(struct d_ocp_qp_dim *dim, struct d_ocp_qp *qp, void *memory);
//
void d_ocp_qp_set_all(double **A, double **B, double **b, double **Q, double **S, double **R, double **q, double **r, int **idxbx, double **lbx, double **ubx, int **idxbu, double **lbu, double **ubu, double **C, double **D, double **lg, double **ug, double **Zl, double **Zu, double **zl, double **zu, int **idxs, double **ls, double **us, struct d_ocp_qp *qp);
//
void d_ocp_qp_set(char *field_name, int stage, void *value, struct d_ocp_qp *qp);
//
void d_ocp_qp_set_A(int stage, double *mat, struct d_ocp_qp *qp);
//
void d_ocp_qp_get_A(int stage, struct d_ocp_qp *qp, double *mat);
//
void d_ocp_qp_set_B(int stage, double *mat, struct d_ocp_qp *qp);
//
void d_ocp_qp_get_B(int stage, struct d_ocp_qp *qp, double *mat);
//
void d_ocp_qp_set_b(int stage, double *vec, struct d_ocp_qp *qp);
//
void d_ocp_qp_get_b(int stage, struct d_ocp_qp *qp, double *vec);
//
void d_ocp_qp_set_Q(int stage, double *mat, struct d_ocp_qp *qp);
//
void d_ocp_qp_get_Q(int stage, struct d_ocp_qp *qp, double *mat);
//
void d_ocp_qp_set_S(int stage, double *mat, struct d_ocp_qp *qp);
//
void d_ocp_qp_get_S(int stage, struct d_ocp_qp *qp, double *mat);
//
void d_ocp_qp_set_R(int stage, double *mat, struct d_ocp_qp *qp);
//
void d_ocp_qp_get_R(int stage, struct d_ocp_qp *qp, double *mat);
//
void d_ocp_qp_set_q(int stage, double *vec, struct d_ocp_qp *qp);
//
void d_ocp_qp_get_q(int stage, struct d_ocp_qp *qp, double *vec);
//
void d_ocp_qp_set_r(int stage, double *vec, struct d_ocp_qp *qp);
//
void d_ocp_qp_get_r(int stage, struct d_ocp_qp *qp, double *vec);
//
void d_ocp_qp_set_lbx(int stage, double *vec, struct d_ocp_qp *qp);
//
void d_ocp_qp_get_lbx(int stage, struct d_ocp_qp *qp, double *vec);
//
void d_ocp_qp_set_ubx(int stage, double *vec, struct d_ocp_qp *qp);
//
void d_ocp_qp_get_ubx(int stage, struct d_ocp_qp *qp, double *vec);
//
void d_ocp_qp_set_lbu(int stage, double *vec, struct d_ocp_qp *qp);
//
void d_ocp_qp_get_lbu(int stage, struct d_ocp_qp *qp, double *vec);
//
void d_ocp_qp_set_ubu(int stage, double *vec, struct d_ocp_qp *qp);
//
void d_ocp_qp_get_ubu(int stage, struct d_ocp_qp *qp, double *vec);
//
void d_ocp_qp_set_idxb(int stage, int *vec, struct d_ocp_qp *qp);
//
void d_ocp_qp_get_idxb(int stage, struct d_ocp_qp *qp, int *vec);
//
void d_ocp_qp_set_idxbx(int stage, int *vec, struct d_ocp_qp *qp);
//
//void d_ocp_qp_get_idxbx(int stage, struct d_ocp_qp *qp, int *vec);
//
void d_ocp_qp_set_Jbx(int stage, double *vec, struct d_ocp_qp *qp);
//
//void d_ocp_qp_get_Jbx(int stage, struct d_ocp_qp *qp, double *vec);
//
void d_ocp_qp_set_idxbu(int stage, int *vec, struct d_ocp_qp *qp);
//
//void d_ocp_qp_get_idxbu(int stage, struct d_ocp_qp *qp, int *vec);
//
void d_ocp_qp_set_Jbu(int stage, double *vec, struct d_ocp_qp *qp);
//
//void d_ocp_qp_get_Jbu(int stage, struct d_ocp_qp *qp, double *vec);
//
void d_ocp_qp_set_lb(int stage, double *vec, struct d_ocp_qp *qp);
//
void d_ocp_qp_get_lb(int stage, struct d_ocp_qp *qp, double *vec);
//
void d_ocp_qp_set_ub(int stage, double *vec, struct d_ocp_qp *qp);
//
void d_ocp_qp_get_ub(int stage, struct d_ocp_qp *qp, double *vec);
//
void d_ocp_qp_set_C(int stage, double *mat, struct d_ocp_qp *qp);
//
void d_ocp_qp_get_C(int stage, struct d_ocp_qp *qp, double *mat);
//
void d_ocp_qp_set_D(int stage, double *mat, struct d_ocp_qp *qp);
//
void d_ocp_qp_get_D(int stage, struct d_ocp_qp *qp, double *mat);
//
void d_ocp_qp_set_lg(int stage, double *vec, struct d_ocp_qp *qp);
//
void d_ocp_qp_get_lg(int stage, struct d_ocp_qp *qp, double *vec);
//
void d_ocp_qp_set_ug(int stage, double *vec, struct d_ocp_qp *qp);
//
void d_ocp_qp_get_ug(int stage, struct d_ocp_qp *qp, double *vec);
//
void d_ocp_qp_set_Zl(int stage, double *vec, struct d_ocp_qp *qp);
//
void d_ocp_qp_get_Zl(int stage, struct d_ocp_qp *qp, double *vec);
//
void d_ocp_qp_set_Zu(int stage, double *vec, struct d_ocp_qp *qp);
//
void d_ocp_qp_get_Zu(int stage, struct d_ocp_qp *qp, double *vec);
//
void d_ocp_qp_set_zl(int stage, double *vec, struct d_ocp_qp *qp);
//
void d_ocp_qp_get_zl(int stage, struct d_ocp_qp *qp, double *vec);
//
void d_ocp_qp_set_zu(int stage, double *vec, struct d_ocp_qp *qp);
//
void d_ocp_qp_get_zu(int stage, struct d_ocp_qp *qp, double *vec);
//
void d_ocp_qp_set_idxs(int stage, int *vec, struct d_ocp_qp *qp);
//
void d_ocp_qp_get_idxs(int stage, struct d_ocp_qp *qp, int *vec);
//
void d_ocp_qp_set_lls(int stage, double *vec, struct d_ocp_qp *qp);
//
void d_ocp_qp_get_lls(int stage, struct d_ocp_qp *qp, double *vec);
//
void d_ocp_qp_set_lus(int stage, double *vec, struct d_ocp_qp *qp);
//
void d_ocp_qp_get_lus(int stage, struct d_ocp_qp *qp, double *vec);
// TODO remove ???
void d_change_bounds_dimensions_ocp_qp(int *nbu, int *nbx, struct d_ocp_qp *qp);



#ifdef __cplusplus
}	// #extern "C"
#endif



#endif // HPIPM_D_OCP_QP_H_
