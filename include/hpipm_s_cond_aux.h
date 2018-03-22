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



#ifndef HPIPM_S_COND_AUX_H_
#define HPIPM_S_COND_AUX_H_



#include <blasfeo_target.h>
#include <blasfeo_common.h>



#ifdef __cplusplus
extern "C" {
#endif



//
void s_cond_BAbt(struct s_ocp_qp *ocp_qp, struct blasfeo_smat *BAbt2, struct blasfeo_svec *b2, struct s_cond_qp_ocp2dense_arg *cond_arg, struct s_cond_qp_ocp2dense_workspace *cond_ws);
//
void s_cond_b(struct s_ocp_qp *ocp_qp, struct blasfeo_svec *b2, struct s_cond_qp_ocp2dense_arg *cond_arg, struct s_cond_qp_ocp2dense_workspace *cond_ws);
//
void s_cond_RSQrq_N2nx3(struct s_ocp_qp *ocp_qp, struct blasfeo_smat *RSQrq2, struct blasfeo_svec *rq2, struct s_cond_qp_ocp2dense_arg *cond_arg, struct s_cond_qp_ocp2dense_workspace *cond_ws);
//
void s_cond_rq_N2nx3(struct s_ocp_qp *ocp_qp, struct blasfeo_svec *rq2, struct s_cond_qp_ocp2dense_arg *cond_arg, struct s_cond_qp_ocp2dense_workspace *cond_ws);
//
void s_cond_DCtd(struct s_ocp_qp *ocp_qp, int *idxb2, struct blasfeo_smat *DCt2, struct blasfeo_svec *d2, int *idxs2, struct blasfeo_svec *Z2, struct blasfeo_svec *z2, struct s_cond_qp_ocp2dense_arg *cond_arg, struct s_cond_qp_ocp2dense_workspace *cond_ws);
//
void s_cond_d(struct s_ocp_qp *ocp_qp, struct blasfeo_svec *d2, struct blasfeo_svec *z2, struct s_cond_qp_ocp2dense_arg *cond_arg, struct s_cond_qp_ocp2dense_workspace *cond_ws);
//
void s_expand_sol(struct s_ocp_qp *ocp_qp, struct s_dense_qp_sol *dense_qp_sol, struct s_ocp_qp_sol *ocp_qp_sol, struct s_cond_qp_ocp2dense_arg *cond_arg, struct s_cond_qp_ocp2dense_workspace *cond_ws);
//
void s_expand_primal_sol(struct s_ocp_qp *ocp_qp, struct s_dense_qp_sol *dense_qp_sol, struct s_ocp_qp_sol *ocp_qp_sol, struct s_cond_qp_ocp2dense_arg *cond_arg, struct s_cond_qp_ocp2dense_workspace *cond_ws);

//
void s_update_cond_BAbt(int *idxc, struct s_ocp_qp *ocp_qp, struct blasfeo_smat *BAbt2, struct blasfeo_svec *b2, struct s_cond_qp_ocp2dense_arg *cond_arg, struct s_cond_qp_ocp2dense_workspace *cond_ws);
//
void s_update_cond_RSQrq_N2nx3(int *idxc, struct s_ocp_qp *ocp_qp, struct blasfeo_smat *RSQrq2, struct blasfeo_svec *rq2, struct s_cond_qp_ocp2dense_arg *cond_arg, struct s_cond_qp_ocp2dense_workspace *cond_ws);
//
void s_update_cond_DCtd(int *idxc, struct s_ocp_qp *ocp_qp, int *idxb2, struct blasfeo_smat *DCt2, struct blasfeo_svec *d2, int *idxs2, struct blasfeo_svec *Z2, struct blasfeo_svec *z2, struct s_cond_qp_ocp2dense_arg *cond_arg, struct s_cond_qp_ocp2dense_workspace *cond_ws);


#ifdef __cplusplus
} /* extern "C" */
#endif



#endif // HPIPM_S_COND_AUX_H_
