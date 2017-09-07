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
#include <blasfeo_d_aux.h>
#include <blasfeo_d_blas.h>

#include "../include/hpipm_d_erk_int.h"
#include "../include/hpipm_d_ocp_qp.h"
#include "../include/hpipm_d_ocp_nlp.h"
#include "../include/hpipm_d_ocp_nlp_sol.h"



void d_cvt_erk_int_to_ocp_qp(int n, struct d_erk_workspace *erk_ws, double *xn, struct d_ocp_qp *qp, struct d_ocp_nlp_sol *nlp_sol)
	{

	int ii;

//	int *nx = qp->nx+n;
//	int *nu = qp->nu+n;
	struct d_strmat *BAbt = qp->BAbt+n;
	struct d_strvec *b = qp->b+n;

	struct d_strvec *ux = nlp_sol->ux+n;

	int nx = erk_ws->nx;
	int nf = erk_ws->nf;

	double *x = erk_ws->x;
	double *xt = b->pa;

	double *tmp;

//	d_cvt_tran_mat2strmat(nx[1], nu[0]+nx[0], x, nx[1], BAbt, 0, 0);
	d_cvt_tran_mat2strmat(nx, nf, x, nx, BAbt, 0, 0);

	// XXX not compute this again in residuals !!!
//	tmp = x+nx*nf;
//	for(ii=0; ii<nx; ii++)
//		xt[ii] = tmp[ii] - xn[ii];
//printf("\n%d\n", nf);
	d_cvt_vec2strvec(nx, x+nf*nx, b, 0);
	dgemv_t_libstr(nf, nx, -1.0, BAbt, 0, 0, ux, 0, 1.0, b, 0, b, 0);
//d_print_strvec(nx, b, 0);
	
	drowin_libstr(nx, 1.0, b, 0, BAbt, nf, 0);

	return;

	}
