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

#include "../include/hpipm_d_ipm_hard_revcom_qp.h"



void d_compute_Qx_qx_hard_qp(struct d_ipm_hard_revcom_qp_workspace *rws)
	{

	int nv = rws->nv;
	int ne = rws->ne;
	int nb = rws->nb;
	int ng = rws->ng;

	double *lam_lb = rws->lam_lb;
	double *lam_ub = rws->lam_ub;
	double *t_lb = rws->t_lb;
	double *t_ub = rws->t_ub;
	double *res_m_lb = rws->res_m_lb;
	double *res_m_ub = rws->res_m_ub;
	double *res_d_lb = rws->res_d_lb;
	double *res_d_ub = rws->res_d_ub;
	double *t_inv_lb = rws->t_inv_lb;
	double *t_inv_ub = rws->t_inv_ub;
	double *Qx = rws->Qx;
	double *qx = rws->qx;

	// local variables
	int nt = nb+ng;
	int ii;

	for(ii=0; ii<nt; ii++)
		{

		t_inv_lb[ii] = 1.0/t_lb[ii];
		t_inv_ub[ii] = 1.0/t_ub[ii];
		// TODO mask out unconstrained components for one-sided
		Qx[ii] = t_inv_lb[ii]*lam_lb[ii] \
		       + t_inv_ub[ii]*lam_ub[ii];
		qx[ii] = t_inv_lb[ii]*(res_m_lb[ii]-lam_lb[ii]*res_d_lb[ii]) \
		       - t_inv_ub[ii]*(res_m_ub[ii]+lam_ub[ii]*res_d_ub[ii]);

		}
		return;

	}



void d_compute_lam_t_hard_qp(struct d_ipm_hard_revcom_qp_workspace *rws)
	{

	int nv = rws->nv;
	int ne = rws->ne;
	int nb = rws->nb;
	int ng = rws->ng;

	double *lam_lb = rws->lam_lb;
	double *lam_ub = rws->lam_ub;
	double *dlam_lb = rws->dlam_lb;
	double *dlam_ub = rws->dlam_ub;
	double *dt_lb = rws->dt_lb;
	double *dt_ub = rws->dt_ub;
	double *res_d_lb = rws->res_d_lb;
	double *res_d_ub = rws->res_d_ub;
	double *res_m_lb = rws->res_m_lb;
	double *res_m_ub = rws->res_m_ub;
	double *t_inv_lb = rws->t_inv_lb;
	double *t_inv_ub = rws->t_inv_ub;

	// local variables
	int ii;
	int nt = nb+ng;

	for(ii=0; ii<nt; ii++)
		{

		dt_ub[ii] = - dt_lb[ii];

		dt_lb[ii] -= res_d_lb[ii];
		dt_ub[ii] += res_d_ub[ii];

		// TODO compute lamda alone ???
		dlam_lb[ii] = - t_inv_lb[ii] * (lam_lb[ii]*dt_lb[ii] + res_m_lb[ii]);
		dlam_ub[ii] = - t_inv_ub[ii] * (lam_ub[ii]*dt_ub[ii] + res_m_ub[ii]);

		}
	
	return;

	}



void d_compute_alpha_hard_qp(struct d_ipm_hard_revcom_qp_workspace *rws)
	{
	
	// extract workspace members
	int nb = rws->nb;
	int ng = rws->ng;

	double *lam_lb = rws->lam_lb;
	double *lam_ub = rws->lam_ub;
	double *t_lb = rws->t_lb;
	double *t_ub = rws->t_ub;
	double *dlam_lb = rws->dlam_lb;
	double *dlam_ub = rws->dlam_ub;
	double *dt_lb = rws->dt_lb;
	double *dt_ub = rws->dt_ub;

	double alpha = 1.0;
	
	// local variables
	int nt = nb+ng;
	int ii;

	for(ii=0; ii<nt; ii++)
		{

		if( -alpha*dlam_lb[ii+0]>lam_lb[ii+0] )
			{
			alpha = - lam_lb[ii+0] / dlam_lb[ii+0];
			}
		if( -alpha*dlam_ub[ii]>lam_ub[ii] )
			{
			alpha = - lam_ub[ii] / dlam_ub[ii];
			}
		if( -alpha*dt_lb[ii+0]>t_lb[ii+0] )
			{
			alpha = - t_lb[ii+0] / dt_lb[ii+0];
			}
		if( -alpha*dt_ub[ii]>t_ub[ii] )
			{
			alpha = - t_ub[ii] / dt_ub[ii];
			}

		}

	// store alpha
	rws->alpha = alpha;

	return;

	}
	


void d_update_var_hard_qp(struct d_ipm_hard_revcom_qp_workspace *rws)
	{
	
	// extract workspace members
	int nv = rws->nv;
	int ne = rws->ne;
	int nb = rws->nb;
	int ng = rws->ng;

	double *v = rws->v;
	double *pi = rws->pi;
	double *lam = rws->lam;
	double *t = rws->t;
	double *dv = rws->dv;
	double *dpi = rws->dpi;
	double *dlam = rws->dlam;
	double *dt = rws->dt;
	double alpha = 0.995*rws->alpha;

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


