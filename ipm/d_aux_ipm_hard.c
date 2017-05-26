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
#include <blasfeo_d_aux.h>
#include <blasfeo_d_blas.h>

#include "../include/hpipm_d_ipm2_hard_revcom.h"



void d_update_hessian_gradient_res_hard(struct d_ipm2_hard_revcom_workspace *workspace)
	{

	// extract workspace members
	int nb = workspace->nb;
	int ng = workspace->ng;
	double *lam = workspace->lam;
	double *t = workspace->t;
	double *res_m = workspace->res_m;
	double *res_d = workspace->res_d;
	double *t_inv = workspace->t_inv;
	double *Qx = workspace->Qx;
	double *qx = workspace->qx;

	// local variables
	int nt = nb+ng;
	int ii;

	for(ii=0; ii<nt; ii++)
		{

		t_inv[ii]    = 1.0/t[ii];
		t_inv[ii+nt] = 1.0/t[ii+nt];
		// TODO mask out unconstrained components for one-sided
		Qx[ii] = t_inv[ii]*lam[ii] \
		       + t_inv[ii+nt]*lam[ii+nt];
		qx[ii] = t_inv[ii]*(res_m[ii]-lam[ii]*res_d[ii]) \
		       - t_inv[ii+nt]*(res_m[ii+nt]+lam[ii+nt]*res_d[ii+nt]);

		}
	
	return;

	}



void d_compute_alpha_res_hard(struct d_ipm2_hard_revcom_workspace *workspace)
	{
	
	// extract workspace members
	int nb = workspace->nb;
	int ng = workspace->ng;
	int *idxb = workspace->idxb;
	double *lam = workspace->lam;
	double *t = workspace->t;
	double *dv = workspace->dv;
	double *dlam = workspace->dlam;
	double *dt = workspace->dt;
	double *dt_lb = workspace->dt_lb;
	double *dt_lg = workspace->dt_lg;
	double *res_d = workspace->res_d;
	double *res_m = workspace->res_m;
	double *t_inv = workspace->t_inv;
	double *Dv = workspace->Dv;
	double alpha = workspace->alpha;
	
	// local variables
	int nt = nb+ng;
	int ii;

	// box constraints
	for(ii=0; ii<nb; ii++)
		dt_lb[ii] = dv[idxb[ii]];

	// general constraints TODO call back for that
//	dgemv_t_libstr(nx0+nu0, ng0, 1.0, &hsDCt[jj], 0, 0, &hsdux[jj], 0, 0.0, &hsdt[jj], nb0, &hsdt[jj], nb0);
	for(ii=0; ii<ng; ii++)
		dt_lg[ii] = Dv[ii];

	for(ii=0; ii<nt; ii++)
		{

		dt[ii+nt] = - dt[ii];

		dt[ii+0]  -= res_d[ii+0];
		dt[ii+nt] += res_d[ii+nt];

		dlam[ii+0]  = - t_inv[ii+0]  * ( lam[ii+0]*dt[ii+0]   + res_m[ii+0] );
		dlam[ii+nt] = - t_inv[ii+nt] * ( lam[ii+nt]*dt[ii+nt] + res_m[ii+nt] );

		if( -alpha*dlam[ii+0]>lam[ii+0] )
			{
			alpha = - lam[ii+0] / dlam[ii+0];
			}
		if( -alpha*dlam[ii+nt]>lam[ii+nt] )
			{
			alpha = - lam[ii+nt] / dlam[ii+nt];
			}
		if( -alpha*dt[ii+0]>t[ii+0] )
			{
			alpha = - t[ii+0] / dt[ii+0];
			}
		if( -alpha*dt[ii+nt]>t[ii+nt] )
			{
			alpha = - t[ii+nt] / dt[ii+nt];
			}

		}

	// store alpha
	workspace->alpha = alpha;

	return;
	
	}



void d_update_var_res_hard(struct d_ipm2_hard_revcom_workspace *workspace)
	{
	
	// extract workspace members
	int nv = workspace->nv;
	int ne = workspace->ne;
	int nb = workspace->nb;
	int ng = workspace->ng;
	double *v = workspace->v;
	double *pi = workspace->pi;
	double *lam = workspace->lam;
	double *t = workspace->t;
	double *dv = workspace->dv;
	double *dpi = workspace->dpi;
	double *dlam = workspace->dlam;
	double *dt = workspace->dt;
	double alpha = workspace->alpha;

	// local variables
	int nt = nb+ng;
	int ii;

	// update v
	for(ii=0; ii<nv; ii++)
		{
		v[ii] += alpha * dv[ii];
		}

	// update pi
	for(ii=0; ii<ne; ii++)
		{
		pi[ii] += alpha * dpi[ii];
		}

	// update lam
	for(ii=0; ii<2*nt; ii++)
		{
		lam[ii] += alpha * dlam[ii];
		}

	// update t
	for(ii=0; ii<2*nt; ii++)
		{
		t[ii] += alpha * dt[ii];
		}
	
	return;

	}

