/**************************************************************************************************
*                                                                                                 *
* This file is part of HPIPM.                                                                     *
*                                                                                                 *
* HPIPM -- High-Performance Interior Point Method.                                                *
* Copyright (C) 2017-2018 by Gianluca Frison.                                                     *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
*                                                                                                 *
* This program is free software: you can redistribute it and/or modify                            *
* it under the terms of the GNU General Public License as published by                            *
* the Free Software Foundation, either version 3 of the License, or                               *
* (at your option) any later version                                                              *.
*                                                                                                 *
* This program is distributed in the hope that it will be useful,                                 *
* but WITHOUT ANY WARRANTY; without even the implied warranty of                                  *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                   *
* GNU General Public License for more details.                                                    *
*                                                                                                 *
* You should have received a copy of the GNU General Public License                               *
* along with this program.  If not, see <https://www.gnu.org/licenses/>.                          *
*                                                                                                 *
* The authors designate this particular file as subject to the "Classpath" exception              *
* as provided by the authors in the LICENSE file that accompained this code.                      *
*                                                                                                 *
* Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             *
*                                                                                                 *
**************************************************************************************************/

#include "../include/hpipm_s_core_qp_ipm.h"



void s_compute_Gamma_gamma_qp(float *res_d, float *res_m, struct s_core_qp_ipm_workspace *cws)
	{

	int nc = cws->nc;

	float *lam = cws->lam;
	float *t = cws->t;
//	float *res_d = cws->res_d; // TODO rename d ???
//	float *res_m = cws->res_m; // TODO rename m ???
	float *t_inv = cws->t_inv;
	float *Gamma = cws->Gamma;
	float *gamma = cws->gamma;

	// local variables
	int ii;

	for(ii=0; ii<nc; ii++)
		{
		t_inv[ii] = 1.0/t[ii];
		Gamma[ii] = t_inv[ii]*lam[ii];
		gamma[ii] = t_inv[ii]*(res_m[ii]-lam[ii]*res_d[ii]);
		}

	return;

	}



void s_compute_gamma_qp(float *res_d, float *res_m, struct s_core_qp_ipm_workspace *cws)
	{

	int nc = cws->nc;

	float *lam = cws->lam;
//	float *res_m = cws->res_m;
//	float *res_d = cws->res_d;
	float *t_inv = cws->t_inv;
	float *gamma = cws->gamma;

	// local variables
	int ii;

	for(ii=0; ii<nc; ii++)
		{
		gamma[ii] = t_inv[ii]*(res_m[ii]-lam[ii]*res_d[ii]);
		}

	return;

	}

void s_compute_lam_t_qp(float *res_d, float *res_m, float *dlam, float *dt, struct s_core_qp_ipm_workspace *cws)
	{

	int nc = cws->nc;

	float *lam = cws->lam;
//	float *dlam = cws->dlam;
//	float *dt = cws->dt;
//	float *res_d = cws->res_d; // TODO rename d ???
//	float *res_m = cws->res_m; // TODO rename m ???
	float *t_inv = cws->t_inv;

	// local variables
	int ii;

	for(ii=0; ii<nc; ii++)
		{
		dt[ii] -= res_d[ii];
		// TODO compute lamda alone ???
		dlam[ii] = - t_inv[ii] * (lam[ii]*dt[ii] + res_m[ii]);
		}
	
	return;

	}



void s_compute_alpha_qp(struct s_core_qp_ipm_workspace *cws)
	{
	
	// extract workspace members
	int nc = cws->nc;

	float *lam_lb = cws->lam;
	float *t_lb = cws->t;
	float *dlam_lb = cws->dlam;
	float *dt_lb = cws->dt;

	float alpha_prim = - 1.0;
	float alpha_dual = - 1.0;
	float alpha = - 1.0;

#if 1

	// local variables
	int ii;

	for(ii=0; ii<nc; ii++)
		{

		if( alpha_dual*dlam_lb[ii+0]>lam_lb[ii+0] )
			{
			alpha_dual = lam_lb[ii+0] / dlam_lb[ii+0];
			}
		if( alpha_prim*dt_lb[ii+0]>t_lb[ii+0] )
			{
			alpha_prim = t_lb[ii+0] / dt_lb[ii+0];
			}

		}

#else // fraction to the boundary

	float mu = cws->mu;
	float tau = 1.0;
//	float tau = 0.995;

	tau = tau>(1-mu) ? tau : 1-mu;
	
	// local variables
	int ii;

	for(ii=0; ii<nc; ii++)
		{

		if( alpha_dual*dlam_lb[ii+0]>tau*lam_lb[ii+0] )
			{
			alpha_dual = tau*lam_lb[ii+0] / dlam_lb[ii+0];
			}
		if( alpha_prim*dt_lb[ii+0]>tau*t_lb[ii+0] )
			{
			alpha_prim = tau*t_lb[ii+0] / dt_lb[ii+0];
			}

		}

#endif
	alpha = alpha_prim>alpha_dual ? alpha_prim : alpha_dual;

	// store alpha
	cws->alpha_prim = - alpha_prim;
	cws->alpha_dual = - alpha_dual;
	cws->alpha = - alpha;

	return;

	}
	


void s_update_var_qp(struct s_core_qp_ipm_workspace *cws)
	{
	
	// extract workspace members
	int nv = cws->nv;
	int ne = cws->ne;
	int nc = cws->nc;

	float *v = cws->v;
	float *pi = cws->pi;
	float *lam = cws->lam;
	float *t = cws->t;
	float *dv = cws->dv;
	float *dpi = cws->dpi;
	float *dlam = cws->dlam;
	float *dt = cws->dt;
	float alpha = cws->alpha;
	float alpha_prim = cws->alpha_prim;
	float alpha_dual = cws->alpha_dual;
#if 0
	if(alpha<1.0)
		alpha *= 0.995;
#else
//	alpha_prim = alpha_prim * ((1.0-alpha)*0.99 + alpha*0.9999999);
//	alpha_dual = alpha_dual * ((1.0-alpha)*0.99 + alpha*0.9999999);
	alpha_prim = alpha_prim * ((1.0-alpha_prim)*0.99 + alpha_prim*0.9999999);
	alpha_dual = alpha_dual * ((1.0-alpha_dual)*0.99 + alpha_dual*0.9999999);
	alpha = alpha * ((1.0-alpha)*0.99 + alpha*0.9999999);
#endif

	// local variables
	int ii;

#if 1

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
	for(ii=0; ii<nc; ii++)
		{
		lam[ii] += alpha * dlam[ii];
		}

	// update t
	for(ii=0; ii<nc; ii++)
		{
		t[ii] += alpha * dt[ii];
		}

#else // split step

	// update v
	for(ii=0; ii<nv; ii++)
		{
		v[ii] += alpha_prim * dv[ii];
		}

	// update pi
	for(ii=0; ii<ne; ii++)
		{
		pi[ii] += alpha_dual * dpi[ii];
		}

	// update lam
	for(ii=0; ii<nc; ii++)
		{
		lam[ii] += alpha_dual * dlam[ii];
		}

	// update t
	for(ii=0; ii<nc; ii++)
		{
		t[ii] += alpha_prim * dt[ii];
		}

#endif

	return;

	}



void s_compute_mu_aff_qp(struct s_core_qp_ipm_workspace *cws)
	{

	int ii;

	// extract workspace members
	int nc = cws->nc;

	float *ptr_lam = cws->lam;
	float *ptr_t = cws->t;
	float *ptr_dlam = cws->dlam;
	float *ptr_dt = cws->dt;
	float alpha = cws->alpha;
	// this affects the minimum value of signa !!!
//		alpha *= 0.99;

	float mu = 0;

	for(ii=0; ii<nc; ii++)
		{
		mu += (ptr_lam[ii+0] + alpha*ptr_dlam[ii+0]) * (ptr_t[ii+0] + alpha*ptr_dt[ii+0]);
		}
	
	cws->mu_aff = mu*cws->nc_inv;

	return;

	}



void s_backup_res_m(struct s_core_qp_ipm_workspace *cws)
	{

	int ii;

	// extract workspace members
	int nc = cws->nc;

	float *ptr_res_m = cws->res_m;
	float *ptr_res_m_bkp = cws->res_m_bkp;

	for(ii=0; ii<nc; ii++)
		{
		ptr_res_m_bkp[ii+0] = ptr_res_m[ii+0];
		}

	return;

	}



void s_compute_centering_correction_qp(struct s_core_qp_ipm_workspace *cws)
	{

	int ii;

	// extract workspace members
	int nc = cws->nc;

	float *ptr_dlam = cws->dlam;
	float *ptr_dt = cws->dt;
	float *ptr_res_m = cws->res_m;
	float *ptr_res_m_bkp = cws->res_m_bkp;

	float sigma_mu = cws->sigma*cws->mu;

	for(ii=0; ii<nc; ii++)
		{
		ptr_res_m[ii+0] = ptr_res_m_bkp[ii+0] + ptr_dt[ii+0] * ptr_dlam[ii+0] - sigma_mu;
		}

	return;

	}



void s_compute_centering_qp(struct s_core_qp_ipm_workspace *cws)
	{

	int ii;

	// extract workspace members
	int nc = cws->nc;

	float *ptr_dlam = cws->dlam;
	float *ptr_dt = cws->dt;
	float *ptr_res_m = cws->res_m;
	float *ptr_res_m_bkp = cws->res_m_bkp;

	float sigma_mu = cws->sigma*cws->mu;

	for(ii=0; ii<nc; ii++)
		{
		ptr_res_m[ii+0] = ptr_res_m_bkp[ii+0] - sigma_mu;
		}

	return;

	}




