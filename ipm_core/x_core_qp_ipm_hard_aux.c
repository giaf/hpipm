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



// TODO rename in compute Gamma_gamma
void COMPUTE_QX_QX_HARD_QP(struct IPM_HARD_CORE_QP_WORKSPACE *rws)
	{

	int nc = rws->nc;

	REAL *lam = rws->lam;
	REAL *t = rws->t;
	REAL *res_m = rws->res_m;
	REAL *res_d = rws->res_d;
	REAL *t_inv = rws->t_inv;
	REAL *Gamma = rws->Gamma;
	REAL *gamma = rws->gamma;
	REAL *Qx = rws->Qx;
	REAL *qx = rws->qx;

	// local variables
	int ii;

	for(ii=0; ii<2*nc; ii++)
		{
		// TODO mask out unconstrained components for one-sided (multiply by zero/one?)
		t_inv[ii] = 1.0/t[ii];
		Gamma[ii] = t_inv[ii]*lam[ii];
		gamma[ii] = t_inv[ii]*(res_m[ii]-lam[ii]*res_d[ii]);
		}

	// TODO remove !!!
	for(ii=0; ii<nc; ii++)
		{
		Qx[ii] = Gamma[0+ii] + Gamma[nc+ii];
		qx[ii] = gamma[0+ii] - gamma[nc+ii];
		}

	return;

	}



void COMPUTE_LAM_T_HARD_QP(struct IPM_HARD_CORE_QP_WORKSPACE *rws)
	{

	int nc = rws->nc;

	REAL *lam = rws->lam;
	REAL *dlam = rws->dlam;
	REAL *dt = rws->dt;
	REAL *res_d = rws->res_d;
	REAL *res_m = rws->res_m;
	REAL *t_inv = rws->t_inv;

	// local variables
	int ii;

	// TODO move outside core !!!
	for(ii=0; ii<nc; ii++)
		{
		dt[nc+ii] = - dt[0 +ii];
		}

	for(ii=0; ii<nc; ii++)
		{

		dt[0 +ii] -= res_d[0 +ii];
		dt[nc+ii] += res_d[nc+ii];

		// TODO compute lamda alone ???
		dlam[0 +ii] = - t_inv[0 +ii] * (lam[0 +ii]*dt[0 +ii] + res_m[0 +ii]);
		dlam[nc+ii] = - t_inv[nc+ii] * (lam[nc+ii]*dt[nc+ii] + res_m[nc+ii]);

		}
	
	return;

	}



void COMPUTE_ALPHA_HARD_QP(struct IPM_HARD_CORE_QP_WORKSPACE *rws)
	{
	
	// extract workspace members
	int nc = rws->nc;

	REAL *lam = rws->lam;
	REAL *t = rws->t;
	REAL *dlam = rws->dlam;
	REAL *dt = rws->dt;

	REAL alpha = - 1.0;
	
	// local variables
	int ii;

	for(ii=0; ii<2*nc; ii++)
		{

		if( alpha*dlam[ii+0]>lam[ii+0] )
			{
			alpha = lam[ii+0] / dlam[ii+0];
			}
		if( alpha*dt[ii+0]>t[ii+0] )
			{
			alpha = t[ii+0] / dt[ii+0];
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
	int nc = rws->nc;

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



void COMPUTE_MU_AFF_HARD_QP(struct IPM_HARD_CORE_QP_WORKSPACE *rws)
	{

	int ii;

	// extract workspace members
	int nc = rws->nc;

	REAL *ptr_lam = rws->lam;
	REAL *ptr_t = rws->t;
	REAL *ptr_dlam = rws->dlam;
	REAL *ptr_dt = rws->dt;
	REAL alpha = rws->alpha;
	// this affects the minimum value of signa !!!
//		alpha *= 0.99;

	REAL mu = 0;

	for(ii=0; ii<2*nc; ii++)
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
	int nc = rws->nc;

	REAL *ptr_dlam = rws->dlam;
	REAL *ptr_dt = rws->dt;
	REAL *ptr_res_m = rws->res_m;

	REAL sigma_mu = rws->sigma*rws->mu;

	for(ii=0; ii<2*nc; ii++)
		{
		ptr_res_m[ii+0] += ptr_dt[ii+0] * ptr_dlam[ii+0] - sigma_mu;
		}

	return;

	}



void COMPUTE_QX_HARD_QP(struct IPM_HARD_CORE_QP_WORKSPACE *rws)
	{

	int nc = rws->nc;

	REAL *lam = rws->lam;
	REAL *res_m = rws->res_m;
	REAL *res_d = rws->res_d;
	REAL *t_inv = rws->t_inv;
	REAL *gamma = rws->gamma;
	REAL *qx = rws->qx;

	// local variables
	int ii;

	for(ii=0; ii<2*nc; ii++)
		{
		// TODO mask out unconstrained components for one-sided (multiply by zero/one?)
		gamma[ii] = t_inv[ii]*(res_m[ii]-lam[ii]*res_d[ii]);
		}

	// TODO remove !!!
	for(ii=0; ii<nc; ii++)
		{
		qx[ii] = gamma[0+ii] - gamma[nc+ii];
		}

	return;

	}


