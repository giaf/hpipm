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



struct d_ocp_qp_sol
	{
	struct d_strvec *ux;
	struct d_strvec *pi;
	struct d_strvec *lam_lb;
	struct d_strvec *lam_ub;
	struct d_strvec *lam_lg;
	struct d_strvec *lam_ug;
	struct d_strvec *t_lb;
	struct d_strvec *t_ub;
	struct d_strvec *t_lg;
	struct d_strvec *t_ug;
	int memsize; // memory size in bytes
	};



//
int d_memsize_ocp_qp_sol(int N, int *nx, int *nu, int *nb, int *ng);
//
void d_create_ocp_qp_sol(int N, int *nx, int *nu, int *nb, int *ng, struct d_ocp_qp_sol *qp_sol, void *memory);
//
void d_cvt_ocp_qp_sol_to_colmaj(struct d_ocp_qp *qp, struct d_ocp_qp_sol *qp_sol, double **u, double **x, double **pi, double **lam_lb, double **lam_ub, double **lam_lg, double **lam_ug);
//
void d_cvt_ocp_qp_sol_to_rowmaj(struct d_ocp_qp *qp, struct d_ocp_qp_sol *qp_sol, double **u, double **x, double **pi, double **lam_lb, double **lam_ub, double **lam_lg, double **lam_ug);
//
void d_cvt_ocp_qp_sol_to_libstr(struct d_ocp_qp *qp, struct d_ocp_qp_sol *qp_sol, struct d_strvec *u, struct d_strvec *x, struct d_strvec *pi, struct d_strvec *lam_lb, struct d_strvec *lam_ub, struct d_strvec *lam_lg, struct d_strvec *lam_ug);
