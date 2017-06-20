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



void COMPUTE_QX_QX_HARD_QP(struct IPM_HARD_CORE_QP_WORKSPACE *rws)
	{

	int nb = rws->nb;
	int ng = rws->ng;

	REAL *lam_lb = rws->lam_lb;
	REAL *lam_ub = rws->lam_ub;
	REAL *t_lb = rws->t_lb;
	REAL *t_ub = rws->t_ub;
	REAL *res_m_lb = rws->res_m_lb;
	REAL *res_m_ub = rws->res_m_ub;
	REAL *res_d_lb = rws->res_d_lb;
	REAL *res_d_ub = rws->res_d_ub;
	REAL *t_inv_lb = rws->t_inv_lb;
	REAL *t_inv_ub = rws->t_inv_ub;
	REAL *Qx = rws->Qx;
	REAL *qx = rws->qx;

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



void COMPUTE_LAM_T_HARD_QP(struct IPM_HARD_CORE_QP_WORKSPACE *rws)
	{

	int nb = rws->nb;
	int ng = rws->ng;

	REAL *lam_lb = rws->lam_lb;
	REAL *lam_ub = rws->lam_ub;
	REAL *dlam_lb = rws->dlam_lb;
	REAL *dlam_ub = rws->dlam_ub;
	REAL *dt_lb = rws->dt_lb;
	REAL *dt_ub = rws->dt_ub;
	REAL *res_d_lb = rws->res_d_lb;
	REAL *res_d_ub = rws->res_d_ub;
	REAL *res_m_lb = rws->res_m_lb;
	REAL *res_m_ub = rws->res_m_ub;
	REAL *t_inv_lb = rws->t_inv_lb;
	REAL *t_inv_ub = rws->t_inv_ub;

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



void COMPUTE_ALPHA_HARD_QP(struct IPM_HARD_CORE_QP_WORKSPACE *rws)
	{
	
	// extract workspace members
	int nb = rws->nb;
	int ng = rws->ng;

	REAL *lam_lb = rws->lam_lb;
	REAL *t_lb = rws->t_lb;
	REAL *dlam_lb = rws->dlam_lb;
	REAL *dt_lb = rws->dt_lb;

	REAL alpha = - 1.0;
	
	// local variables
	int nt = nb+ng;
	int ii;

	for(ii=0; ii<2*nt; ii++)
		{

		if( alpha*dlam_lb[ii+0]>lam_lb[ii+0] )
			{
			alpha = lam_lb[ii+0] / dlam_lb[ii+0];
			}
		if( alpha*dt_lb[ii+0]>t_lb[ii+0] )
			{
			alpha = t_lb[ii+0] / dt_lb[ii+0];
			}

		}

	// store alpha
	rws->alpha = - alpha;

	return;

	}
	


void UPDATE_VAR_HARD_QP(struct IPM_HARD_CORE_QP_WORKSPACE *rws)
	{
	
	// extract workspace members
	int nv = rws->nv;
	int ne = rws->ne;
	int nb = rws->nb;
	int ng = rws->ng;

	REAL *v = rws->v;
	REAL *pi = rws->pi;
	REAL *lam = rws->lam;
	REAL *t = rws->t;
	REAL *dv = rws->dv;
	REAL *dpi = rws->dpi;
	REAL *dlam = rws->dlam;
	REAL *dt = rws->dt;
	REAL alpha = rws->alpha;
#if 0
	if(alpha<1.0)
		alpha *= 0.995;
#else
	alpha = alpha * ((1.0-alpha)*0.99 + alpha*0.9999999);
#endif

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



void COMPUTE_MU_AFF_HARD_QP(struct IPM_HARD_CORE_QP_WORKSPACE *rws)
	{

	int ii;

	// extract workspace members
	int nb = rws->nb;
	int ng = rws->ng;
	int nt = nb+ng;

	REAL *ptr_lam = rws->lam_lb;
	REAL *ptr_t = rws->t_lb;
	REAL *ptr_dlam = rws->dlam_lb;
	REAL *ptr_dt = rws->dt_lb;
	REAL alpha = rws->alpha;
	// this affects the minimum value of signa !!!
//		alpha *= 0.99;

	REAL mu = 0;

	for(ii=0; ii<2*nt; ii++)
		{
		mu += (ptr_lam[ii+0] + alpha*ptr_dlam[ii+0]) * (ptr_t[ii+0] + alpha*ptr_dt[ii+0]);
		}
	
	rws->mu_aff = mu*rws->nt_inv;

	return;

	}



void COMPUTE_CENTERING_CORRECTION_HARD_QP(struct IPM_HARD_CORE_QP_WORKSPACE *rws)
	{

	int ii;

	// extract workspace members
	int nb = rws->nb;
	int ng = rws->ng;
	int nt = nb+ng;

	REAL *ptr_dlam = rws->dlam_lb;
	REAL *ptr_dt = rws->dt_lb;
	REAL *ptr_res_m = rws->res_m_lb;

	REAL sigma_mu = rws->sigma*rws->mu;

	for(ii=0; ii<2*nt; ii++)
		{
		ptr_res_m[ii+0] += ptr_dt[ii+0] * ptr_dlam[ii+0] - sigma_mu;
		}

	return;

	}



void COMPUTE_QX_HARD_QP(struct IPM_HARD_CORE_QP_WORKSPACE *rws)
	{

	int nb = rws->nb;
	int ng = rws->ng;

	REAL *lam_lb = rws->lam_lb;
	REAL *lam_ub = rws->lam_ub;
	REAL *res_m_lb = rws->res_m_lb;
	REAL *res_m_ub = rws->res_m_ub;
	REAL *res_d_lb = rws->res_d_lb;
	REAL *res_d_ub = rws->res_d_ub;
	REAL *t_inv_lb = rws->t_inv_lb;
	REAL *t_inv_ub = rws->t_inv_ub;
	REAL *qx = rws->qx;

	// local variables
	int nt = nb+ng;
	int ii;

	for(ii=0; ii<nt; ii++)
		{

		// TODO mask out unconstrained components for one-sided
		qx[ii] = t_inv_lb[ii]*(res_m_lb[ii]-lam_lb[ii]*res_d_lb[ii]) \
		       - t_inv_ub[ii]*(res_m_ub[ii]+lam_ub[ii]*res_d_ub[ii]);

		}

	return;

	}


