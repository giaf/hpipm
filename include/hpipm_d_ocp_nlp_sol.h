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



struct d_ocp_nlp_sol
	{
	struct d_strvec *ux;
	struct d_strvec *pi;
	struct d_strvec *lam;
	struct d_strvec *t;
	struct d_strvec *eta0;
	int memsize; // memory size in bytes
	};



//
int d_memsize_ocp_nlp_sol(int N, int *nx, int *nu, int *nb, int *ng, int *ns, int ne0);
//
void d_create_ocp_nlp_sol(int N, int *nx, int *nu, int *nb, int *ng, int *ns, int ne0, struct d_ocp_nlp_sol *qp_sol, void *memory);
//
void d_cvt_ocp_nlp_sol_to_colmaj(struct d_ocp_nlp *qp, struct d_ocp_nlp_sol *qp_sol, double **u, double **x, double **ls, double **us, double **pi, double **lam_lb, double **lam_ub, double **lam_lg, double **lam_ug, double **lam_ls, double **lam_us, double *eta0);
//
void d_cvt_ocp_nlp_sol_to_rowmaj(struct d_ocp_nlp *qp, struct d_ocp_nlp_sol *qp_sol, double **u, double **x, double **ls, double **us, double **pi, double **lam_lb, double **lam_ub, double **lam_lg, double **lam_ug, double **lam_ls, double **lam_us, double *eta0);
//
void d_cvt_ocp_nlp_sol_to_libstr(struct d_ocp_nlp *qp, struct d_ocp_nlp_sol *qp_sol, struct d_strvec *u, struct d_strvec *ls, struct d_strvec *us, struct d_strvec *x, struct d_strvec *pi, struct d_strvec *lam_lb, struct d_strvec *lam_ub, struct d_strvec *lam_lg, struct d_strvec *lam_ug, struct d_strvec *lam_ls, struct d_strvec *lam_us, struct d_strvec *eta0);
