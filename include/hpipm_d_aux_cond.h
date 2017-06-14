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



//
void d_compute_Gamma(struct d_ocp_qp *ocp_qp, struct d_cond_qp_ocp2dense_workspace *cond_ws);
//
void d_cond_RSQrq_N2nx3(struct d_ocp_qp *ocp_qp, struct d_strmat *RSQrq2, struct d_strvec *rq2, struct d_cond_qp_ocp2dense_workspace *cond_ws);
//
void d_cond_DCtd(struct d_ocp_qp *ocp_qp, int *idxb2, struct d_strvec *d_lb2, struct d_strvec *d_ub2, struct d_strmat *DCt2, struct d_strvec *d_lg2, struct d_strvec *d_ug2, struct d_cond_qp_ocp2dense_workspace *cond_ws);
//
void d_expand_sol(struct d_ocp_qp *ocp_qp, struct d_dense_qp_sol *dense_qp_sol, struct d_strvec *ux, struct d_strvec *pi, struct d_strvec *lam_lb, struct d_strvec *lam_ub, struct d_strvec *lam_lg, struct d_strvec *lam_ug, struct d_strvec *t_lb, struct d_strvec *t_ub, struct d_strvec *t_lg, struct d_strvec *t_ug, struct d_strvec *tmp_nuxM, struct d_strvec *tmp_ngM);
