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


void d_compute_Qx_qx_step_dense_qp(struct d_dense_qp *qp, struct d_ipm2_hard_dense_qp_workspace *ws)
	{

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;

	struct d_ipm2_hard_revcom_qp_workspace *rws = ws->revcom_workspace;

	double *lam_lb = rws->lam_lb;
	double *lam_ub = rws->lam_ub;
	double *t_lb = rws->t_lb;
	double *t_ub = rws->t_ub;
	double *res_m_lb = rws->res_m;
	double *res_m_ub = rws->res_m;
	double *res_d_lb = rws->res_d;
	double *res_d_ub = rws->res_d;
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



void d_compute_lam_t_step_dense_qp(struct d_dense_qp *qp, struct d_ipm2_hard_dense_qp_workspace *ws)
	{

	return;

	}



// range-space (Schur complement) method
void d_fact_solve_kkt_step_dense_qp(struct d_dense_qp *qp, struct d_ipm2_hard_dense_qp_workspace *ws)
	{

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;

	if(nb>0 | ng>0)
		{
		d_compute_Qx_qx_step_dense_qp(qp, ws);
		}

	dtrcp_l_libstr(nv, qp->H, 0, 0, ws->Lv, 0, 0);

	dveccp_libstr(nv, qp->g, 0, ws->lv, 0);

	if(nb>0)
		{
		ddiaad_sp_libstr(nb, 1.0, ws->Qx, 0, qp->idxb, ws->Lv, 0, 0);
		dvecad_sp_libstr(nb, 1.0, ws->qx, 0, qp->idxb, ws->lv, 0);
		}

	dpotrf_l_libstr(nv, ws->Lv, 0, 0, ws->Lv, 0, 0);

	dgecp_libstr(ne, nv, qp->A, 0, 0, ws->AL, 0, 0);
	dtrsm_rltn_libstr(ne, nv, 1.0, ws->Lv, 0, 0, qp->A, 0, 0, ws->AL, 0, 0);

	dgese_libstr(ne, ne, 0.0, ws->Le, 0, 0);
	dsyrk_dpotrf_ln_libstr(ne, ne, nv, ws->AL, 0, 0, ws->AL, 0, 0, ws->Le, 0, 0, ws->Le, 0, 0);

	dtrsv_lnn_libstr(nv, ws->Lv, 0, 0, ws->lv, 0, ws->lv, 0);

	dgemv_n_libstr(ne, nv, 1.0, ws->AL, 0, 0, ws->lv, 0, 1.0, qp->b, 0, ws->le, 0);

	dtrsv_lnn_libstr(ne, ws->Le, 0, 0, ws->le, 0, ws->le, 0);
	dtrsv_ltn_libstr(ne, ws->Le, 0, 0, ws->le, 0, ws->pi, 0);

	dgemv_t_libstr(ne, nv, 1.0, qp->A, 0, 0, ws->pi, 0, -1.0, qp->g, 0, ws->lv, 0);
	dtrsv_lnn_libstr(nv, ws->Lv, 0, 0, ws->lv, 0, ws->lv, 0);
	dtrsv_ltn_libstr(nv, ws->Lv, 0, 0, ws->lv, 0, ws->v, 0);

	if(nb>0 | ng>0)
		{
		d_compute_lam_t_step_dense_qp(qp, ws);
		}

	return;

	}
