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



void COMPUTE_QX_QX_QP(struct IPM_CORE_QP_WORKSPACE *cws)
	{

	int nc = cws->nc;

	REAL *lam = cws->lam;
	REAL *t = cws->t;
	REAL *res_m = cws->res_m;
	REAL *res_d = cws->res_d;
	REAL *t_inv = cws->t_inv;
	REAL *Qx = cws->Qx;
	REAL *qx = cws->qx;
	REAL *Gamma = cws->Gamma;
	REAL *gamma = cws->gamma;

	// local variables
	int ii;

	for(ii=0; ii<2*nc; ii++)
		{
		t_inv[ii] = 1.0/t[ii];
		Gamma[ii] = t_inv[ii]*lam[ii];
		gamma[ii] = t_inv[ii]*(res_m[ii]-lam[ii]*res_d[ii]);
		}

	// TODO move outside core
	for(ii=0; ii<nc; ii++)
		{
		// TODO mask out unconstrained components for one-sided (multiply by zero?)
		Qx[ii] = Gamma[ii] + Gamma[nc+ii];
		qx[ii] = gamma[ii] - gamma[nc+ii];
		}
	
	return;

	}



void COMPUTE_LAM_T_QP(struct IPM_CORE_QP_WORKSPACE *cws)
	{

	int nc = cws->nc;

	REAL *lam = cws->lam;
	REAL *dlam = cws->dlam;
	REAL *dt = cws->dt;
	REAL *res_d = cws->res_d;
	REAL *res_m = cws->res_m;
	REAL *t_inv = cws->t_inv;

	// local variables
	int ii;

	// TODO move outside core
	// XXX only ocp_qp works now !!!!!!!!!!!!!!!!!!!!
//	for(ii=0; ii<nc; ii++)
//		{
//		dt[nc+ii] = - dt[ii];
//		}
	
	for(ii=0; ii<2*nc; ii++)
		{
		dt[ii] -= res_d[ii]; // XXX change sign for upper?
		// TODO compute lamda alone ???
		dlam[ii] = - t_inv[ii] * (lam[ii]*dt[ii] + res_m[ii]);
		}
	
	return;

	}



void COMPUTE_ALPHA_QP(struct IPM_CORE_QP_WORKSPACE *cws)
	{
	
	// extract workspace members
	int nc = cws->nc;

	REAL *lam_lb = cws->lam;
	REAL *t_lb = cws->t;
	REAL *dlam_lb = cws->dlam;
	REAL *dt_lb = cws->dt;

	REAL alpha = - 1.0;
	
	// local variables
	int ii;

	for(ii=0; ii<2*nc; ii++)
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
	cws->alpha = - alpha;

	return;

	}
	


void UPDATE_VAR_QP(struct IPM_CORE_QP_WORKSPACE *cws)
	{
	
	// extract workspace members
	int nv = cws->nv;
	int ne = cws->ne;
	int nc = cws->nc;

	REAL *v = cws->v;
	REAL *pi = cws->pi;
	REAL *lam = cws->lam;
	REAL *t = cws->t;
	REAL *dv = cws->dv;
	REAL *dpi = cws->dpi;
	REAL *dlam = cws->dlam;
	REAL *dt = cws->dt;
	REAL alpha = cws->alpha;
#if 0
	if(alpha<1.0)
		alpha *= 0.995;
#else
	alpha = alpha * ((1.0-alpha)*0.99 + alpha*0.9999999);
#endif

	// local variables
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
	for(ii=0; ii<2*nc; ii++)
		{
		lam[ii] += alpha * dlam[ii];
		}

	// update t
	for(ii=0; ii<2*nc; ii++)
		{
		t[ii] += alpha * dt[ii];
		}
	
	return;

	}



void COMPUTE_MU_AFF_QP(struct IPM_CORE_QP_WORKSPACE *cws)
	{

	int ii;

	// extract workspace members
	int nc = cws->nc;

	REAL *ptr_lam = cws->lam;
	REAL *ptr_t = cws->t;
	REAL *ptr_dlam = cws->dlam;
	REAL *ptr_dt = cws->dt;
	REAL alpha = cws->alpha;
	// this affects the minimum value of signa !!!
//		alpha *= 0.99;

	REAL mu = 0;

	for(ii=0; ii<2*nc; ii++)
		{
		mu += (ptr_lam[ii+0] + alpha*ptr_dlam[ii+0]) * (ptr_t[ii+0] + alpha*ptr_dt[ii+0]);
		}
	
	cws->mu_aff = mu*cws->nt_inv;

	return;

	}



void COMPUTE_CENTERING_CORRECTION_QP(struct IPM_CORE_QP_WORKSPACE *cws)
	{

	int ii;

	// extract workspace members
	int nc = cws->nc;

	REAL *ptr_dlam = cws->dlam;
	REAL *ptr_dt = cws->dt;
	REAL *ptr_res_m = cws->res_m;

	REAL sigma_mu = cws->sigma*cws->mu;

	for(ii=0; ii<2*nc; ii++)
		{
		ptr_res_m[ii+0] += ptr_dt[ii+0] * ptr_dlam[ii+0] - sigma_mu;
		}

	return;

	}



void COMPUTE_QX_QP(struct IPM_CORE_QP_WORKSPACE *cws)
	{

	int nc = cws->nc;

	REAL *lam = cws->lam;
	REAL *res_m = cws->res_m;
	REAL *res_d = cws->res_d;
	REAL *t_inv = cws->t_inv;
	REAL *qx = cws->qx;
	REAL *gamma = cws->gamma;

	// local variables
	int ii;

	for(ii=0; ii<2*nc; ii++)
		{
		gamma[ii] = t_inv[ii]*(res_m[ii]-lam[ii]*res_d[ii]);
		}

	// TODO move outside core
	for(ii=0; ii<nc; ii++)
		{
		// TODO mask out unconstrained components for one-sided
		qx[ii] = gamma[ii] - gamma[nc+ii];
		}

	return;

	}


