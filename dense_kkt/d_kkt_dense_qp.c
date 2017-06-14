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

#include "../include/hpipm_d_dense_qp.h"
#include "../include/hpipm_d_dense_qp_sol.h"
#include "../include/hpipm_d_ipm_hard_dense_qp.h"
#include "../include/hpipm_d_ipm_hard_core_qp.h"
#include "../include/hpipm_d_aux_ipm_hard.h"



void d_init_var_hard_dense_qp(struct d_dense_qp *qp, struct d_dense_qp_sol *qp_sol, struct d_ipm_hard_dense_qp_workspace *ws)
	{

	struct d_ipm_hard_core_qp_workspace *rws = ws->core_workspace;

	// extract rws members
	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;

	double *d_lb = qp->d_lb->pa;
	double *d_ub = qp->d_ub->pa;
	double *d_lg = qp->d_lg->pa;
	double *d_ug = qp->d_ug->pa;
	int *idxb = qp->idxb;

	double *v = qp_sol->v->pa;
	double *pi = qp_sol->pi->pa;
	double *lam_lb = qp_sol->lam_lb->pa;
	double *lam_ub = qp_sol->lam_ub->pa;
	double *lam_lg = qp_sol->lam_lg->pa;
	double *lam_ug = qp_sol->lam_ug->pa;
	double *t_lb = qp_sol->t_lb->pa;
	double *t_ub = qp_sol->t_ub->pa;
	double *t_lg = qp_sol->t_lg->pa;
	double *t_ug = qp_sol->t_ug->pa;
	double mu0 = rws->mu0;

	// local variables
	int ii;
	int idxb0;
	double thr0 = 0.1;

	// warm start TODO


	// cold start

	// primal variables
	for(ii=0; ii<nv; ii++)
		{
		v[ii] = 0.0;
		}
	
	// equality constraints
	for(ii=0; ii<ne; ii++)
		{
		pi[ii] = 0.0;
		}
	
	// box constraints
	for(ii=0; ii<nb; ii++)
		{
		idxb0 = idxb[ii];
		t_lb[ii] = - d_lb[ii] + v[idxb0];
		t_ub[ii] =   d_ub[ii] - v[idxb0];
		if(t_lb[ii]<thr0)
			{
			if(t_ub[ii]<thr0)
				{
				v[idxb0] = 0.5*(d_lb[ii] - d_ub[ii]);
				t_lb[ii] = thr0;
				t_ub[ii] = thr0;
				}
			else
				{
				t_lb[ii] = thr0;
				v[idxb0] = d_lb[ii] + thr0;
				}
			}
		else if(t_ub[ii]<thr0)
			{
			t_ub[ii] = thr0;
			v[idxb0] = d_ub[ii] - thr0;
			}
		lam_lb[ii] = mu0/t_lb[ii];
		lam_ub[ii] = mu0/t_ub[ii];
		}
	
	// inequality constraints
	dgemv_t_libstr(nv, ng, 1.0, qp->Ct, 0, 0, qp_sol->v, 0, 0.0, qp_sol->t_lg, 0, qp_sol->t_lg, 0);
	for(ii=0; ii<ng; ii++)
		{
		t_ug[ii] = t_lg[ii];
		t_lg[ii] -= d_lg[ii];
		t_ug[ii] += d_ug[ii];
		t_lg[ii] = fmax( thr0, t_lg[ii] );
		t_ug[ii] = fmax( thr0, t_ug[ii] );
//		t_lg[ii] = thr0>t_lg[ii] ? thr0 : t_lg[ii];
//		t_ug[ii] = thr0>t_ug[ii] ? thr0 : t_ug[ii];
		lam_lg[ii] = mu0/t_lg[ii];
		lam_ug[ii] = mu0/t_ug[ii];
		}

	return;

	}



void d_compute_res_hard_dense_qp(struct d_dense_qp *qp, struct d_dense_qp_sol *qp_sol, struct d_ipm_hard_dense_qp_workspace *ws)
	{

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;

	// TODO extract qp arguments !!!!!

	struct d_strvec *v = qp_sol->v;
	struct d_strvec *pi = qp_sol->pi;
	struct d_strvec *lam_lb = qp_sol->lam_lb;
	struct d_strvec *lam_ub = qp_sol->lam_ub;
	struct d_strvec *lam_lg = qp_sol->lam_lg;
	struct d_strvec *lam_ug = qp_sol->lam_ug;
	struct d_strvec *t_lb = qp_sol->t_lb;
	struct d_strvec *t_ub = qp_sol->t_ub;
	struct d_strvec *t_lg = qp_sol->t_lg;
	struct d_strvec *t_ug = qp_sol->t_ug;

	double mu;

	// res g
	dsymv_l_libstr(nv, nv, 1.0, qp->Hg, 0, 0, v, 0, 1.0, qp->g, 0, ws->res_g, 0);

	if(nb>0)
		{
		// res_g
		daxpy_libstr(nb, -1.0, lam_lb, 0, lam_ub, 0, ws->tmp_nb, 0);
		dvecad_sp_libstr(nb, 1.0, ws->tmp_nb, 0, qp->idxb, ws->res_g, 0);
		// res_d
		dvecex_sp_libstr(nb, -1.0, qp->idxb, v, 0, ws->res_d_lb, 0);
		dveccp_libstr(nb, ws->res_d_lb, 0, ws->res_d_ub, 0);
		daxpy_libstr(nb, 1.0, qp->d_lb, 0, ws->res_d_lb, 0, ws->res_d_lb, 0);
		daxpy_libstr(nb, 1.0, qp->d_ub, 0, ws->res_d_ub, 0, ws->res_d_ub, 0);
		daxpy_libstr(nb, 1.0, t_lb, 0, ws->res_d_lb, 0, ws->res_d_lb, 0);
		daxpy_libstr(nb, -1.0, t_ub, 0, ws->res_d_ub, 0, ws->res_d_ub, 0);
		}

	if(ng>0)
		{
		daxpy_libstr(ng, -1.0, lam_lg, 0, lam_ug, 0, ws->tmp_ng0, 0);
		daxpy_libstr(ng, 1.0, t_lg, 0, qp->d_lg, 0, ws->res_d_lg, 0);
		daxpy_libstr(ng, -1.0, t_ug, 0, qp->d_ug, 0, ws->res_d_ug, 0);
		dgemv_nt_libstr(nv, ng, 1.0, 1.0, qp->Ct, 0, 0, ws->tmp_ng0, 0, v, 0, 1.0, 0.0, ws->res_g, 0, ws->tmp_ng1, 0, ws->res_g, 0, ws->tmp_ng1, 0);
		daxpy_libstr(ng, -1.0, ws->tmp_ng1, 0, ws->res_d_lg, 0, ws->res_d_lg, 0);
		daxpy_libstr(ng, -1.0, ws->tmp_ng1, 0, ws->res_d_ug, 0, ws->res_d_ug, 0);
		}
	
	// res b, res g
	dgemv_nt_libstr(ne, nv, -1.0, -1.0, qp->A, 0, 0, v, 0, pi, 0, 1.0, 1.0, qp->b, 0, ws->res_g, 0, ws->res_b, 0, ws->res_g, 0);

	// res_mu
	mu = dvecmuldot_libstr(2*nb+2*ng, lam_lb, 0, t_lb, 0, ws->res_m, 0);

	ws->res_mu = mu*ws->nt_inv;


	return;

	}


// range-space (Schur complement) method
void d_fact_solve_kkt_step_hard_dense_qp(struct d_dense_qp *qp, struct d_ipm_hard_dense_qp_workspace *ws)
	{

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;

	struct d_ipm_hard_core_qp_workspace *rws = ws->core_workspace;

	if(nb>0 | ng>0)
		{
		d_compute_Qx_qx_hard_qp(rws);
		}

	dtrcp_l_libstr(nv, qp->Hg, 0, 0, ws->Lv, 0, 0);

	dveccp_libstr(nv, ws->res_g, 0, ws->lv, 0);

	if(nb>0)
		{
		ddiaad_sp_libstr(nb, 1.0, ws->Qx, 0, qp->idxb, ws->Lv, 0, 0);
		dvecad_sp_libstr(nb, 1.0, ws->qx, 0, qp->idxb, ws->lv, 0);
		}

	if(ng>0)
		{
		dgemv_n_libstr(nv, ng, 1.0, qp->Ct, 0, 0, ws->qx, nb, 1.0, ws->lv, 0, ws->lv, 0);
		dgemm_r_diag_libstr(nv, ng, 1.0, qp->Ct, 0, 0, ws->Qx, nb, 0.0, ws->Ctx, 0, 0, ws->Ctx, 0, 0);
		dsyrk_dpotrf_ln_libstr(nv, nv, ng, ws->Ctx, 0, 0, qp->Ct, 0, 0, ws->Lv, 0, 0, ws->Lv, 0, 0);
		}
	else
		{
		dpotrf_l_libstr(nv, ws->Lv, 0, 0, ws->Lv, 0, 0);
		}

	dveccp_libstr(nv, ws->lv, 0, ws->dv, 0);

	if(ne>0)
		{
		dgecp_libstr(ne, nv, qp->A, 0, 0, ws->AL, 0, 0);
		dtrsm_rltn_libstr(ne, nv, 1.0, ws->Lv, 0, 0, qp->A, 0, 0, ws->AL, 0, 0);

		dgese_libstr(ne, ne, 0.0, ws->Le, 0, 0);
		dsyrk_dpotrf_ln_libstr(ne, ne, nv, ws->AL, 0, 0, ws->AL, 0, 0, ws->Le, 0, 0, ws->Le, 0, 0);

		dtrsv_lnn_libstr(nv, ws->Lv, 0, 0, ws->lv, 0, ws->lv, 0);

		dgemv_n_libstr(ne, nv, 1.0, ws->AL, 0, 0, ws->lv, 0, 1.0, ws->res_b, 0, ws->dpi, 0);

		dtrsv_lnn_libstr(ne, ws->Le, 0, 0, ws->dpi, 0, ws->dpi, 0);
		dtrsv_ltn_libstr(ne, ws->Le, 0, 0, ws->dpi, 0, ws->dpi, 0);

		dgemv_t_libstr(ne, nv, 1.0, qp->A, 0, 0, ws->dpi, 0, -1.0, ws->dv, 0, ws->dv, 0);
		}
	else
		{
		dvecsc_libstr(nv, -1.0, ws->dv, 0);
		}

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
		d_compute_lam_t_hard_qp(rws);
		}

	return;

	}

