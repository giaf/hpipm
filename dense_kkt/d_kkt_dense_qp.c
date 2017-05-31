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

#include "../include/hpipm_d_dense_qp.h"
#include "../include/hpipm_d_ipm2_hard_dense_qp.h"
#include "../include/hpipm_d_ipm2_hard_revcom_qp.h"



void d_compute_res_dense_qp(struct d_dense_qp *qp, struct d_ipm2_hard_dense_qp_workspace *workspace)
	{

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;

	double mu;

	// res g
	dsymv_l_libstr(nv, nv, 1.0, qp->H, 0, 0, workspace->v, 0, 1.0, qp->g, 0, workspace->res_g, 0);

	if(nb>0)
		{
		// res_g
		daxpy_libstr(nb, -1.0, workspace->lam_lb, 0, workspace->lam_ub, 0, workspace->tmp_nb, 0);
		dvecad_sp_libstr(nb, 1.0, workspace->tmp_nb, 0, qp->idxb, workspace->res_g, 0);
		// res_d
		dvecex_sp_libstr(nb, -1.0, qp->idxb, workspace->v, 0, workspace->res_d_lb, 0);
		dveccp_libstr(nb, workspace->res_d_lb, 0, workspace->res_d_ub, 0);
		daxpy_libstr(nb, 1.0, qp->d_lb, 0, workspace->res_d_lb, 0, workspace->res_d_lb, 0);
		daxpy_libstr(nb, 1.0, qp->d_ub, 0, workspace->res_d_ub, 0, workspace->res_d_ub, 0);
		daxpy_libstr(nb, 1.0, workspace->t_lb, 0, workspace->res_d_lb, 0, workspace->res_d_lb, 0);
		daxpy_libstr(nb, -1.0, workspace->t_ub, 0, workspace->res_d_ub, 0, workspace->res_d_ub, 0);
		}

	if(ng>0)
		{
		// TODO
		}
	
	// res b, res g
	dgemv_nt_libstr(ne, nv, -1.0, -1.0, qp->A, 0, 0, workspace->v, 0, workspace->pi, 0, 1.0, 1.0, qp->b, 0, workspace->res_g, 0, workspace->res_b, 0, workspace->res_g, 0);

	// res_mu
	mu = dvecmuldot_libstr(2*nb+2*ng, workspace->lam, 0, workspace->t, 0, workspace->res_m, 0);

	workspace->res_mu = mu*workspace->nt_inv;


	return;

	}


void d_compute_Qx_qx_step_dense_qp(struct d_ipm2_hard_revcom_qp_workspace *rws)
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



void d_compute_lam_t_step_dense_qp(struct d_ipm2_hard_revcom_qp_workspace *rws)
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



// range-space (Schur complement) method
void d_fact_solve_kkt_step_dense_qp(struct d_dense_qp *qp, struct d_ipm2_hard_dense_qp_workspace *ws)
	{

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;

	struct d_ipm2_hard_revcom_qp_workspace *rws = ws->revcom_workspace;

	if(nb>0 | ng>0)
		{
		d_compute_Qx_qx_step_dense_qp(rws);
		}

	dtrcp_l_libstr(nv, qp->H, 0, 0, ws->Lv, 0, 0);

	dveccp_libstr(nv, ws->res_g, 0, ws->lv, 0);

	if(nb>0)
		{
		ddiaad_sp_libstr(nb, 1.0, ws->Qx, 0, qp->idxb, ws->Lv, 0, 0);
		dvecad_sp_libstr(nb, 1.0, ws->qx, 0, qp->idxb, ws->lv, 0);
		}

	dveccp_libstr(nv, ws->lv, 0, ws->dv, 0);

	dpotrf_l_libstr(nv, ws->Lv, 0, 0, ws->Lv, 0, 0);

	dgecp_libstr(ne, nv, qp->A, 0, 0, ws->AL, 0, 0);
	dtrsm_rltn_libstr(ne, nv, 1.0, ws->Lv, 0, 0, qp->A, 0, 0, ws->AL, 0, 0);

	dgese_libstr(ne, ne, 0.0, ws->Le, 0, 0);
	dsyrk_dpotrf_ln_libstr(ne, ne, nv, ws->AL, 0, 0, ws->AL, 0, 0, ws->Le, 0, 0, ws->Le, 0, 0);

	dtrsv_lnn_libstr(nv, ws->Lv, 0, 0, ws->lv, 0, ws->lv, 0);

	dgemv_n_libstr(ne, nv, 1.0, ws->AL, 0, 0, ws->lv, 0, 1.0, ws->res_b, 0, ws->dpi, 0);

	dtrsv_lnn_libstr(ne, ws->Le, 0, 0, ws->dpi, 0, ws->dpi, 0);
	dtrsv_ltn_libstr(ne, ws->Le, 0, 0, ws->dpi, 0, ws->dpi, 0);

	dgemv_t_libstr(ne, nv, 1.0, qp->A, 0, 0, ws->dpi, 0, -1.0, ws->dv, 0, ws->dv, 0);
	dtrsv_lnn_libstr(nv, ws->Lv, 0, 0, ws->dv, 0, ws->dv, 0);
	dtrsv_ltn_libstr(nv, ws->Lv, 0, 0, ws->dv, 0, ws->dv, 0);

	if(nb>0)
		{
		dvecex_sp_libstr(nb, 1.0, qp->idxb, ws->dv, 0, ws->dt_lb, 0);
		}

	if(ng>0)
		{
		dgemv_t_libstr(nv, ng, 1.0, qp->Ct, 0, 0, ws->dv, 0, 0.0, ws->dt_lg, 0, ws->dt_lg, 0);
		}

	if(nb>0 | ng>0)
		{
		d_compute_lam_t_step_dense_qp(rws);
		}

	return;

	}



void d_compute_alpha_dense_qp(struct d_ipm2_hard_revcom_qp_workspace *rws)
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
	


void d_update_var_dense_qp(struct d_ipm2_hard_revcom_qp_workspace *rws)
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




