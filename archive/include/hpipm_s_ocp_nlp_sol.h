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

#ifdef __cplusplus
extern "C" {
#endif

struct s_ocp_nlp_sol
	{
	struct blasfeo_svec *ux;
	struct blasfeo_svec *pi;
	struct blasfeo_svec *lam;
	struct blasfeo_svec *t;
	struct blasfeo_svec *eta0;
	int memsize; // memory size in bytes
	};



//
int s_memsize_ocp_nlp_sol(int N, int *nx, int *nu, int *nb, int *ng, int *ns, int ne0);
//
void s_create_ocp_nlp_sol(int N, int *nx, int *nu, int *nb, int *ng, int *ns, int ne0, struct s_ocp_nlp_sol *qp_sol, void *memory);
//
void s_cvt_ocp_nlp_sol_to_colmaj(struct s_ocp_nlp *qp, struct s_ocp_nlp_sol *qp_sol, float **u, float **x, float **ls, float **us, float **pi, float **lam_lb, float **lam_ub, float **lam_lg, float **lam_ug, float **lam_ls, float **lam_us, float *eta0);
//
void s_cvt_ocp_nlp_sol_to_rowmaj(struct s_ocp_nlp *qp, struct s_ocp_nlp_sol *qp_sol, float **u, float **x, float **ls, float **us, float **pi, float **lam_lb, float **lam_ub, float **lam_lg, float **lam_ug, float **lam_ls, float **lam_us, float *eta0);
//
void s_cvt_ocp_nlp_sol_to_libstr(struct s_ocp_nlp *qp, struct s_ocp_nlp_sol *qp_sol, struct blasfeo_svec *u, struct blasfeo_svec *ls, struct blasfeo_svec *us, struct blasfeo_svec *x, struct blasfeo_svec *pi, struct blasfeo_svec *lam_lb, struct blasfeo_svec *lam_ub, struct blasfeo_svec *lam_lg, struct blasfeo_svec *lam_ug, struct blasfeo_svec *lam_ls, struct blasfeo_svec *lam_us, struct s_stevec *eta0);

#ifdef __cplusplus
} /* extern "C" */
#endif
