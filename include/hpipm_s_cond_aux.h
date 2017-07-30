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



#include <blasfeo_target.h>
#include <blasfeo_common.h>



//
void s_compute_Gamma(struct s_ocp_qp *ocp_qp, struct s_cond_qp_ocp2dense_workspace *cond_ws);
//
void s_cond_BAbt(struct s_ocp_qp *ocp_qp, struct s_strmat *BAbt2, struct s_strvec *b2, struct s_cond_qp_ocp2dense_workspace *cond_ws);
//
void s_cond_RSQrq_N2nx3(struct s_ocp_qp *ocp_qp, struct s_strmat *RSQrq2, struct s_strvec *rq2, struct s_cond_qp_ocp2dense_workspace *cond_ws);
//
void s_cond_DCtd(struct s_ocp_qp *ocp_qp, int *idxb2, struct s_strvec *d_lb2, struct s_strvec *d_ub2, struct s_strmat *DCt2, struct s_strvec *d_lg2, struct s_strvec *d_ug2, struct s_cond_qp_ocp2dense_workspace *cond_ws);
//
void s_expand_sol(struct s_ocp_qp *ocp_qp, struct s_dense_qp_sol *dense_qp_sol, struct s_strvec *ux, struct s_strvec *pi, struct s_strvec *lam_lb, struct s_strvec *lam_ub, struct s_strvec *lam_lg, struct s_strvec *lam_ug, struct s_strvec *t_lb, struct s_strvec *t_ub, struct s_strvec *t_lg, struct s_strvec *t_ug, struct s_strvec *tmp_nuxM, struct s_strvec *tmp_ngM);

