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




struct d_cond_qp_ocp2ocp_workspace
	{
	struct d_strmat *Gamma;
	struct d_strmat *L;
	struct d_strmat *Lx;
	struct d_strmat *AL;
	struct d_strvec *Gammab;
	struct d_strvec *tmp_ngM;
	struct d_strvec *tmp_nuxM;
	int memsize;
	};



//
void d_compute_qp_size_ocp2ocp(int N, int *nx, int *nu, int *nb, int **idxb, int *ng, int N2, int *nx2, int *nu2, int *nb2, int *ng2);
//
int d_memsize_cond_qp_ocp2ocp(struct d_ocp_qp *ocp_qp, struct d_ocp_qp *part_dense_qp);
//
void d_create_cond_qp_ocp2ocp(struct d_ocp_qp *ocp_qp, struct d_ocp_qp *part_dense_qp, struct d_cond_qp_ocp2ocp_workspace *cond_ws, void *mem);
//
void d_cond_qp_ocp2ocp(struct d_ocp_qp *ocp_qp, struct d_ocp_qp *part_dense_qp, struct d_cond_qp_ocp2ocp_workspace *cond_ws);
