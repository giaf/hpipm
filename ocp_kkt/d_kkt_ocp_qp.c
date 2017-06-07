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



#include <math.h>

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_d_aux.h>
#include <blasfeo_d_blas.h>

#include "../include/hpipm_d_ocp_qp.h"
#include "../include/hpipm_d_ipm_hard_ocp_qp.h"
#include "../include/hpipm_d_ipm_hard_core_qp.h"
#include "../include/hpipm_d_aux_ipm_hard.h"



void d_init_var_hard_ocp_qp(struct d_ocp_qp *qp, struct d_ipm_hard_ocp_qp_workspace *ws)
	{

//	struct d_ipm_hard_core_qp_workspace *rws = ws->core_workspace;
	
	// loop index
	int ii, jj;

	//
	int N = qp->N;
	int *nx = qp->nx;
	int *nu = qp->nu;
	int *nb = qp->nb;
	int *ng = qp->ng;

	//
	double *ux, *pi, *d_lb, *d_ub, *d_lg, *d_ug, *lam_lb, *lam_ub, *lam_lg, *lam_ug, *t_lb, *t_ub, *t_lg, *t_ug;
	int *idxb;

	double thr0 = 0.1;
	double mu0 = ws->mu0;

	// warm start TODO

	// cold start

	// ux
	for(ii=0; ii<=N; ii++)
		{
		ux = ws->ux[ii].pa;
		for(jj=0; jj<nu[ii]+nx[ii]; jj++)
			{
			ux[jj] = 0.0;
			}
		}
	
	// pi
	for(ii=0; ii<N; ii++)
		{
		pi = ws->pi[ii].pa;
		for(jj=0; jj<nx[ii+1]; jj++)
			{
			pi[jj] = 0.0;
			}
		}
	
	// box constraints
	for(ii=0; ii<=N; ii++)
		{
		ux = ws->ux[ii].pa;
		d_lb = qp->lb[ii].pa;
		d_ub = qp->ub[ii].pa;
		lam_lb = ws->lam_lb[ii].pa;
		lam_ub = ws->lam_ub[ii].pa;
		t_lb = ws->t_lb[ii].pa;
		t_ub = ws->t_ub[ii].pa;
		idxb = qp->idxb[ii];
		for(jj=0; jj<nb[ii]; jj++)
			{
			t_lb[jj] = - d_lb[jj] + ux[idxb[jj]];
			t_ub[jj] =   d_ub[jj] - ux[idxb[jj]];
			if(t_lb[jj]<thr0)
				{
				if(t_ub[jj]<thr0)
					{
					ux[idxb[jj]] = 0.5*(d_lb[jj]-d_ub[jj]);
					t_lb[jj] = thr0;
					t_ub[jj] = thr0;
					}
				else
					{
					t_lb[jj] = thr0;
					ux[idxb[jj]] = d_lb[jj] + thr0;
					}
				}
			else if(t_ub[jj]<thr0)
				{
				t_ub[jj] = thr0;
				ux[idxb[jj]] = d_ub[jj] - thr0;
				}
			lam_lb[jj] = mu0/t_lb[jj];
			lam_ub[jj] = mu0/t_ub[jj];
			}
		}
	
	// general constraints
	for(ii=0; ii<=N; ii++)
		{
		t_lg = ws->t_lg[ii].pa;
		t_ug = ws->t_ug[ii].pa;
		lam_lg = ws->lam_lg[ii].pa;
		lam_ug = ws->lam_ug[ii].pa;
		d_lg = qp->lg[ii].pa;
		d_ug = qp->ug[ii].pa;
		ux = ws->ux[ii].pa;
		dgemv_t_libstr(nu[ii]+nx[ii], ng[ii], 1.0, qp->DCt+ii, 0, 0, ws->ux+ii, 0, 0.0, ws->t_lg+ii, 0, ws->t_lg+ii, 0);
		for(jj=0; jj<ng[ii]; jj++)
			{
			t_ug[jj] = - t_lg[jj];
			t_lg[jj] -= d_lg[jj];
			t_ug[jj] += d_ug[jj];
			t_lg[jj] = fmax(thr0, t_lg[jj]);
			t_ug[jj] = fmax(thr0, t_ug[jj]);
			lam_lg[jj] = mu0/t_lg[jj];
			lam_ug[jj] = mu0/t_ug[jj];
			}
		}


	return;

	}


