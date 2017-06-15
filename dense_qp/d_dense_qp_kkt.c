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
#include "../include/hpipm_d_dense_qp_ipm_hard.h"
#include "../include/hpipm_d_core_qp_ipm_hard.h"
#include "../include/hpipm_d_core_qp_ipm_hard_aux.h"



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

	struct d_ipm_hard_core_qp_workspace *cws = ws->core_workspace;

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

	ws->res_mu = mu*cws->nt_inv;


	return;

	}


// range-space (Schur complement) method
void d_fact_solve_kkt_step_hard_dense_qp(struct d_dense_qp *qp, struct d_ipm_hard_dense_qp_workspace *ws)
	{

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;
	struct d_strmat *Hg = qp->Hg;
	struct d_strmat *A = qp->A;
	struct d_strmat *Ct = qp->Ct;
	int *idxb = qp->idxb;

	struct d_strmat *Lv = ws->Lv;
	struct d_strmat *Le = ws->Le;
	struct d_strmat *Ctx = ws->Ctx;
	struct d_strmat *AL = ws->AL;
	struct d_strvec *lv = ws->lv;
	struct d_strvec *dv = ws->dv;
	struct d_strvec *dpi = ws->dpi;
	struct d_strvec *dt_lb = ws->dt_lb;
	struct d_strvec *dt_lg = ws->dt_lg;
	struct d_strvec *res_g = ws->res_g;
	struct d_strvec *res_b = ws->res_b;
	struct d_strvec *Qx = ws->Qx;
	struct d_strvec *qx = ws->qx;

	struct d_ipm_hard_core_qp_workspace *rws = ws->core_workspace;

	if(nb>0 | ng>0)
		{
		d_compute_Qx_qx_hard_qp(rws);
		}

	if(ne>0)
		{
		dtrcp_l_libstr(nv, Hg, 0, 0, Lv, 0, 0);

		dveccp_libstr(nv, res_g, 0, lv, 0);

		if(nb>0)
			{
			ddiaad_sp_libstr(nb, 1.0, Qx, 0, idxb, Lv, 0, 0);
			dvecad_sp_libstr(nb, 1.0, qx, 0, idxb, lv, 0);
			}

		if(ng>0)
			{
			dgemv_n_libstr(nv, ng, 1.0, Ct, 0, 0, qx, nb, 1.0, lv, 0, lv, 0);
			dgemm_r_diag_libstr(nv, ng, 1.0, Ct, 0, 0, Qx, nb, 0.0, Ctx, 0, 0, Ctx, 0, 0);
			dsyrk_dpotrf_ln_libstr(nv, nv, ng, Ctx, 0, 0, Ct, 0, 0, Lv, 0, 0, Lv, 0, 0);
			}
		else
			{
			dpotrf_l_libstr(nv, Lv, 0, 0, Lv, 0, 0);
			}

		dveccp_libstr(nv, lv, 0, dv, 0);

		dgecp_libstr(ne, nv, A, 0, 0, AL, 0, 0);
		dtrsm_rltn_libstr(ne, nv, 1.0, Lv, 0, 0, A, 0, 0, AL, 0, 0);

		dgese_libstr(ne, ne, 0.0, Le, 0, 0);
		dsyrk_dpotrf_ln_libstr(ne, ne, nv, AL, 0, 0, AL, 0, 0, Le, 0, 0, Le, 0, 0);

		dtrsv_lnn_libstr(nv, Lv, 0, 0, lv, 0, lv, 0);

		dgemv_n_libstr(ne, nv, 1.0, AL, 0, 0, lv, 0, 1.0, res_b, 0, dpi, 0);

		dtrsv_lnn_libstr(ne, Le, 0, 0, dpi, 0, dpi, 0);
		dtrsv_ltn_libstr(ne, Le, 0, 0, dpi, 0, dpi, 0);

		dgemv_t_libstr(ne, nv, 1.0, A, 0, 0, dpi, 0, -1.0, dv, 0, dv, 0);

		dtrsv_lnn_libstr(nv, Lv, 0, 0, dv, 0, dv, 0);
		dtrsv_ltn_libstr(nv, Lv, 0, 0, dv, 0, dv, 0);
		}
	else
		{
#if 0
		dtrcp_l_libstr(nv, Hg, 0, 0, Lv, 0, 0);
		dveccp_libstr(nv, res_g, 0, lv, 0);

		if(nb>0)
			{
			ddiaad_sp_libstr(nb, 1.0, Qx, 0, idxb, Lv, 0, 0);
			dvecad_sp_libstr(nb, 1.0, qx, 0, idxb, lv, 0);
			}

		if(ng>0)
			{
			dgemm_r_diag_libstr(nv, ng, 1.0, Ct, 0, 0, Qx, nb, 0.0, Ctx, 0, 0, Ctx, 0, 0);
			dgemv_n_libstr(nv, ng, 1.0, Ct, 0, 0, qx, nb, 1.0, lv, 0, lv, 0);
			dsyrk_dpotrf_ln_libstr(nv, nv, ng, Ctx, 0, 0, Ct, 0, 0, Lv, 0, 0, Lv, 0, 0); // TODO _mn_ routine in BLASFEO !!!
			}
		else
			{
			dpotrf_l_libstr(nv, Lv, 0, 0, Lv, 0, 0);
			}

		dveccp_libstr(nv, lv, 0, dv, 0);
		dvecsc_libstr(nv, -1.0, dv, 0);

		dtrsv_lnn_libstr(nv, Lv, 0, 0, dv, 0, dv, 0);
		dtrsv_ltn_libstr(nv, Lv, 0, 0, dv, 0, dv, 0);
#else
		dtrcp_l_libstr(nv, Hg, 0, 0, Lv, 0, 0);
		drowin_libstr(nv, 1.0, res_g, 0, Lv, nv, 0);

		if(nb>0)
			{
			ddiaad_sp_libstr(nb, 1.0, Qx, 0, idxb, Lv, 0, 0);
			drowad_sp_libstr(nb, 1.0, qx, 0, idxb, Lv, nv, 0);
			}

		if(ng>0)
			{
			dgemm_r_diag_libstr(nv, ng, 1.0, Ct, 0, 0, Qx, nb, 0.0, Ctx, 0, 0, Ctx, 0, 0);
			drowin_libstr(ng, 1.0, qx, nb, Ctx, nv, 0);
			dsyrk_dpotrf_ln_libstr(nv+1, nv, ng, Ctx, 0, 0, Ct, 0, 0, Lv, 0, 0, Lv, 0, 0); // TODO _mn_ routine in BLASFEO !!!
			}
		else
			{
			dpotrf_l_mn_libstr(nv+1, nv, Lv, 0, 0, Lv, 0, 0);
			}

		drowex_libstr(nv, -1.0, Lv, nv, 0, dv, 0);
		dtrsv_ltn_libstr(nv, Lv, 0, 0, dv, 0, dv, 0);
#endif
		}

	if(nb>0)
		{
		dvecex_sp_libstr(nb, 1.0, idxb, dv, 0, dt_lb, 0);
		}

	if(ng>0)
		{
		dgemv_t_libstr(nv, ng, 1.0, Ct, 0, 0, dv, 0, 0.0, dt_lg, 0, dt_lg, 0);
		}

	if(nb>0 | ng>0)
		{
		d_compute_lam_t_hard_qp(rws);
		}

	return;

	}



// range-space (Schur complement) method
void d_solve_kkt_step_hard_dense_qp(struct d_dense_qp *qp, struct d_ipm_hard_dense_qp_workspace *ws)
	{

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;
	struct d_strmat *A = qp->A;
	struct d_strmat *Ct = qp->Ct;
	int *idxb = qp->idxb;

	struct d_strmat *Lv = ws->Lv;
	struct d_strmat *Le = ws->Le;
	struct d_strmat *Ctx = ws->Ctx;
	struct d_strmat *AL = ws->AL;
	struct d_strvec *lv = ws->lv;
	struct d_strvec *dv = ws->dv;
	struct d_strvec *dpi = ws->dpi;
	struct d_strvec *dt_lb = ws->dt_lb;
	struct d_strvec *dt_lg = ws->dt_lg;
	struct d_strvec *res_g = ws->res_g;
	struct d_strvec *res_b = ws->res_b;
	struct d_strvec *qx = ws->qx;

	struct d_ipm_hard_core_qp_workspace *rws = ws->core_workspace;

	if(nb>0 | ng>0)
		{
		d_compute_qx_hard_qp(rws);
		}

	if(ne>0)
		{
		dveccp_libstr(nv, res_g, 0, lv, 0);

		if(nb>0)
			{
			dvecad_sp_libstr(nb, 1.0, qx, 0, idxb, lv, 0);
			}

		if(ng>0)
			{
			dgemv_n_libstr(nv, ng, 1.0, Ct, 0, 0, qx, nb, 1.0, lv, 0, lv, 0);
			}

		dveccp_libstr(nv, lv, 0, dv, 0);

		dtrsv_lnn_libstr(nv, Lv, 0, 0, lv, 0, lv, 0);

		dgemv_n_libstr(ne, nv, 1.0, AL, 0, 0, lv, 0, 1.0, res_b, 0, dpi, 0);

		dtrsv_lnn_libstr(ne, Le, 0, 0, dpi, 0, dpi, 0);
		dtrsv_ltn_libstr(ne, Le, 0, 0, dpi, 0, dpi, 0);

		dgemv_t_libstr(ne, nv, 1.0, A, 0, 0, dpi, 0, -1.0, dv, 0, dv, 0);

		dtrsv_lnn_libstr(nv, Lv, 0, 0, dv, 0, dv, 0);
		dtrsv_ltn_libstr(nv, Lv, 0, 0, dv, 0, dv, 0);
		}
	else
		{
		dveccp_libstr(nv, res_g, 0, lv, 0);

		if(nb>0)
			{
			dvecad_sp_libstr(nb, 1.0, qx, 0, idxb, lv, 0);
			}

		if(ng>0)
			{
			dgemv_n_libstr(nv, ng, 1.0, Ct, 0, 0, qx, nb, 1.0, lv, 0, lv, 0);
			}

		dveccp_libstr(nv, lv, 0, dv, 0);
		dvecsc_libstr(nv, -1.0, dv, 0);

		dtrsv_lnn_libstr(nv, Lv, 0, 0, dv, 0, dv, 0);
		dtrsv_ltn_libstr(nv, Lv, 0, 0, dv, 0, dv, 0);
		}

	if(nb>0)
		{
		dvecex_sp_libstr(nb, 1.0, idxb, dv, 0, dt_lb, 0);
		}

	if(ng>0)
		{
		dgemv_t_libstr(nv, ng, 1.0, Ct, 0, 0, dv, 0, 0.0, dt_lg, 0, dt_lg, 0);
		}

	if(nb>0 | ng>0)
		{
		d_compute_lam_t_hard_qp(rws);
		}

	return;

	}

