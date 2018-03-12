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



#ifndef HPIPM_D_COND_AUX_H_
#define HPIPM_D_COND_AUX_H_



#include <blasfeo_target.h>
#include <blasfeo_common.h>



#ifdef __cplusplus
extern "C" {
#endif



//
void d_cond_BAbt(struct d_ocp_qp *ocp_qp, struct blasfeo_dmat *BAbt2, struct blasfeo_dvec *b2, struct d_cond_qp_ocp2dense_workspace *cond_ws);
//
void d_cond_b(struct d_ocp_qp *ocp_qp, struct blasfeo_dvec *b2, struct d_cond_qp_ocp2dense_workspace *cond_ws);
//
void d_cond_RSQrq_N2nx3(struct d_ocp_qp *ocp_qp, struct blasfeo_dmat *RSQrq2, struct blasfeo_dvec *rq2, struct d_cond_qp_ocp2dense_workspace *cond_ws);
//
void d_cond_rq_N2nx3(struct d_ocp_qp *ocp_qp, struct blasfeo_dvec *rq2, struct d_cond_qp_ocp2dense_workspace *cond_ws);
//
void d_cond_DCtd(struct d_ocp_qp *ocp_qp, int *idxb2, struct blasfeo_dmat *DCt2, struct blasfeo_dvec *d2, int *idxs2, struct blasfeo_dvec *Z2, struct blasfeo_dvec *z2, struct d_cond_qp_ocp2dense_workspace *cond_ws);
//
void d_cond_d(struct d_ocp_qp *ocp_qp, struct blasfeo_dvec *d2, struct blasfeo_dvec *z2, struct d_cond_qp_ocp2dense_workspace *cond_ws);
//
void d_expand_sol(struct d_ocp_qp *ocp_qp, struct d_dense_qp_sol *dense_qp_sol, struct d_ocp_qp_sol *ocp_qp_sol, struct d_cond_qp_ocp2dense_workspace *cond_ws);
//
void d_expand_primal_sol(struct d_ocp_qp *ocp_qp, struct d_dense_qp_sol *dense_qp_sol, struct d_ocp_qp_sol *ocp_qp_sol, struct d_cond_qp_ocp2dense_workspace *cond_ws);



#ifdef __cplusplus
} /* extern "C" */
#endif



#endif // HPIPM_D_COND_AUX_H_
