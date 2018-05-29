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



void INIT_VAR_OCP_QP(struct OCP_QP *qp, struct OCP_QP_SOL *qp_sol, struct OCP_QP_IPM_WORKSPACE *ws)
	{

//	struct CORE_QP_IPM_WORKSPACE *cws = ws->core_workspace;
	
	// loop index
	int ii, jj;

	//
	int N = qp->dim->N;
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;
	int *nb = qp->dim->nb;
	int *ng = qp->dim->ng;
	int *ns = qp->dim->ns;

	REAL mu0 = ws->mu0;

	//
	REAL *ux, *pi, *d_lb, *d_ub, *d_lg, *d_ug, *lam_lb, *lam_ub, *lam_lg, *lam_ug, *t_lb, *t_ub, *t_lg, *t_ug;
	int *idxb;

	REAL thr0 = 0.1;

	// ux
	if(ws->warm_start==0)
		{
		// cold start
		for(ii=0; ii<=N; ii++)
			{
			ux = qp_sol->ux[ii].pa;
			for(jj=0; jj<nu[ii]+nx[ii]+2*ns[ii]; jj++)
				{
				ux[jj] = 0.0;
				}
			}
		}
	else
		{
		// warm start (keep u and x in solution)
		for(ii=0; ii<=N; ii++)
			{
			ux = qp_sol->ux[ii].pa;
			for(jj=nu[ii]+nx[ii]; jj<nu[ii]+nx[ii]+2*ns[ii]; jj++)
				{
				ux[jj] = 0.0;
				}
			}
		}
	
	// pi
	for(ii=0; ii<N; ii++)
		{
		pi = qp_sol->pi[ii].pa;
		for(jj=0; jj<nx[ii+1]; jj++)
			{
			pi[jj] = 0.0;
			}
		}

	// box constraints
	for(ii=0; ii<=N; ii++)
		{
		ux = qp_sol->ux[ii].pa;
		d_lb = qp->d[ii].pa+0;
		d_ub = qp->d[ii].pa+nb[ii]+ng[ii];
		lam_lb = qp_sol->lam[ii].pa+0;
		lam_ub = qp_sol->lam[ii].pa+nb[ii]+ng[ii];
		t_lb = qp_sol->t[ii].pa+0;
		t_ub = qp_sol->t[ii].pa+nb[ii]+ng[ii];
		idxb = qp->idxb[ii];
		for(jj=0; jj<nb[ii]; jj++)
			{
#if 1
			t_lb[jj] = - d_lb[jj] + ux[idxb[jj]];
			t_ub[jj] = - d_ub[jj] - ux[idxb[jj]];
			if(t_lb[jj]<thr0)
				{
				if(t_ub[jj]<thr0)
					{
					ux[idxb[jj]] = 0.5*(d_lb[jj] + d_ub[jj]);
					t_lb[jj] = thr0;
					t_ub[jj] = thr0;
					}
				else
					{
					t_lb[jj] = thr0;
					ux[idxb[jj]] = d_lb[jj] + thr0;
					}
				}
			else if(t_ub[jj]<thr0)
				{
				t_ub[jj] = thr0;
				ux[idxb[jj]] = - d_ub[jj] - thr0;
				}
#else
			t_lb[jj] = 1.0;
			t_ub[jj] = 1.0;
#endif
			lam_lb[jj] = mu0/t_lb[jj];
			lam_ub[jj] = mu0/t_ub[jj];
			}
		}
	
	// general constraints
	for(ii=0; ii<=N; ii++)
		{
		t_lg = qp_sol->t[ii].pa+nb[ii];
		t_ug = qp_sol->t[ii].pa+2*nb[ii]+ng[ii];
		lam_lg = qp_sol->lam[ii].pa+nb[ii];
		lam_ug = qp_sol->lam[ii].pa+2*nb[ii]+ng[ii];
		d_lg = qp->d[ii].pa+nb[ii];
		d_ug = qp->d[ii].pa+2*nb[ii]+ng[ii];
		ux = qp_sol->ux[ii].pa;
		GEMV_T(nu[ii]+nx[ii], ng[ii], 1.0, qp->DCt+ii, 0, 0, qp_sol->ux+ii, 0, 0.0, qp_sol->t+ii, nb[ii], qp_sol->t+ii, nb[ii]);
		for(jj=0; jj<ng[ii]; jj++)
			{
#if 1
			t_ug[jj] = - t_lg[jj];
			t_lg[jj] -= d_lg[jj];
			t_ug[jj] -= d_ug[jj];
//			t_lg[jj] = fmax(thr0, t_lg[jj]);
//			t_ug[jj] = fmax(thr0, t_ug[jj]);
			t_lg[jj] = thr0>t_lg[jj] ? thr0 : t_lg[jj];
			t_ug[jj] = thr0>t_ug[jj] ? thr0 : t_ug[jj];
#else
			t_lg[jj] = 1.0;
			t_ug[jj] = 1.0;
#endif
			lam_lg[jj] = mu0/t_lg[jj];
			lam_ug[jj] = mu0/t_ug[jj];
			}
		}

	// soft constraints
	for(ii=0; ii<=N; ii++)
		{
		lam_lb = qp_sol->lam[ii].pa+2*nb[ii]+2*ng[ii];
		lam_ub = qp_sol->lam[ii].pa+2*nb[ii]+2*ng[ii]+ns[ii];
		t_lb = qp_sol->t[ii].pa+2*nb[ii]+2*ng[ii];
		t_ub = qp_sol->t[ii].pa+2*nb[ii]+2*ng[ii]+ns[ii];
		for(jj=0; jj<ns[ii]; jj++)
			{
			t_lb[jj] = 1.0; // thr0;
			t_ub[jj] = 1.0; // thr0;
			lam_lb[jj] = mu0/t_lb[jj];
			lam_ub[jj] = mu0/t_ub[jj];
			}
		}

	return;

	}



void COMPUTE_RES_OCP_QP(struct OCP_QP *qp, struct OCP_QP_SOL *qp_sol, struct OCP_QP_RES *res, struct OCP_QP_RES_WORKSPACE *ws)
	{

	// loop index
	int ii;

	//
	int N = qp->dim->N;
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;
	int *nb = qp->dim->nb;
	int *ng = qp->dim->ng;
	int *ns = qp->dim->ns;

	int nct = 0;
	for(ii=0; ii<=N; ii++)
		nct += 2*nb[ii]+2*ng[ii]+2*ns[ii];
	
	REAL nct_inv = 1.0/nct;

	struct STRMAT *BAbt = qp->BAbt;
	struct STRMAT *RSQrq = qp->RSQrq;
	struct STRMAT *DCt = qp->DCt;
	struct STRVEC *b = qp->b;
	struct STRVEC *rqz = qp->rqz;
	struct STRVEC *d = qp->d;
	struct STRVEC *m = qp->m;
	int **idxb = qp->idxb;
	struct STRVEC *Z = qp->Z;
	int **idxs = qp->idxs;

	struct STRVEC *ux = qp_sol->ux;
	struct STRVEC *pi = qp_sol->pi;
	struct STRVEC *lam = qp_sol->lam;
	struct STRVEC *t = qp_sol->t;

	struct STRVEC *res_g = res->res_g;
	struct STRVEC *res_b = res->res_b;
	struct STRVEC *res_d = res->res_d;
	struct STRVEC *res_m = res->res_m;

	struct STRVEC *tmp_nbgM = ws->tmp_nbgM;
	struct STRVEC *tmp_nsM = ws->tmp_nsM;

	int nx0, nx1, nu0, nu1, nb0, ng0, ns0;

	//
	REAL mu = 0.0;

	// loop over stages
	for(ii=0; ii<=N; ii++)
		{

		nx0 = nx[ii];
		nu0 = nu[ii];
		nb0 = nb[ii];
		ng0 = ng[ii];
		ns0 = ns[ii];

		SYMV_L(nu0+nx0, nu0+nx0, 1.0, RSQrq+ii, 0, 0, ux+ii, 0, 1.0, rqz+ii, 0, res_g+ii, 0);

		if(ii>0)
			AXPY(nx0, -1.0, pi+(ii-1), 0, res_g+ii, nu0, res_g+ii, nu0);

		if(nb0+ng0>0)
			{
			AXPY(nb0+ng0, -1.0, lam+ii, 0, lam+ii, nb[ii]+ng[ii], tmp_nbgM+0, 0);
//			AXPY(nb0+ng0,  1.0, d+ii, 0, t+ii, 0, res_d+ii, 0);
//			AXPY(nb0+ng0,  1.0, d+ii, nb0+ng0, t+ii, nb0+ng0, res_d+ii, nb0+ng0);
			AXPY(2*nb0+2*ng0,  1.0, d+ii, 0, t+ii, 0, res_d+ii, 0);
			// box
			if(nb0>0)
				{
				VECAD_SP(nb0, 1.0, tmp_nbgM+0, 0, idxb[ii], res_g+ii, 0);
				VECEX_SP(nb0, 1.0, idxb[ii], ux+ii, 0, tmp_nbgM+1, 0);
				}
			// general
			if(ng0>0)
				{
				GEMV_NT(nu0+nx0, ng0, 1.0, 1.0, DCt+ii, 0, 0, tmp_nbgM+0, nb[ii], ux+ii, 0, 1.0, 0.0, res_g+ii, 0, tmp_nbgM+1, nb0, res_g+ii, 0, tmp_nbgM+1, nb0);
				}

			AXPY(nb0+ng0, -1.0, tmp_nbgM+1, 0, res_d+ii, 0, res_d+ii, 0);
			AXPY(nb0+ng0,  1.0, tmp_nbgM+1, 0, res_d+ii, nb0+ng0, res_d+ii, nb0+ng0);
			}
		if(ns0>0)
			{
			// res_g
			GEMV_DIAG(2*ns0, 1.0, Z+ii, 0, ux+ii, nu0+nx0, 1.0, rqz+ii, nu0+nx0, res_g+ii, nu0+nx0);
			AXPY(2*ns0, -1.0, lam+ii, 2*nb0+2*ng0, res_g+ii, nu0+nx0, res_g+ii, nu0+nx0);
			VECEX_SP(ns0, 1.0, idxs[ii], lam+ii, 0, tmp_nsM, 0);
			AXPY(ns0, -1.0, tmp_nsM, 0, res_g+ii, nu0+nx0, res_g+ii, nu0+nx0);
			VECEX_SP(ns0, 1.0, idxs[ii], lam+ii, nb0+ng0, tmp_nsM, 0);
			AXPY(ns0, -1.0, tmp_nsM, 0, res_g+ii, nu0+nx0+ns0, res_g+ii, nu0+nx0+ns0);
			// res_d
			VECAD_SP(ns0, -1.0, ux+ii, nu0+nx0, idxs[ii], res_d+ii, 0);
			VECAD_SP(ns0, -1.0, ux+ii, nu0+nx0+ns0, idxs[ii], res_d+ii, nb0+ng0);
			AXPY(2*ns0, -1.0, ux+ii, nu0+nx0, t+ii, 2*nb0+2*ng0, res_d+ii, 2*nb0+2*ng0);
			AXPY(2*ns0, 1.0, d+ii, 2*nb0+2*ng0, res_d+ii, 2*nb0+2*ng0, res_d+ii, 2*nb0+2*ng0);
			}

		if(ii<N)
			{

			nu1 = nu[ii+1];
			nx1 = nx[ii+1];

			AXPY(nx1, -1.0, ux+(ii+1), nu1, b+ii, 0, res_b+ii, 0);

			GEMV_NT(nu0+nx0, nx1, 1.0, 1.0, BAbt+ii, 0, 0, pi+ii, 0, ux+ii, 0, 1.0, 1.0, res_g+ii, 0, res_b+ii, 0, res_g+ii, 0, res_b+ii, 0);

			}

		mu += VECMULDOT(2*nb0+2*ng0+2*ns0, lam+ii, 0, t+ii, 0, res_m+ii, 0);
		AXPY(2*nb0+2*ng0+2*ns0, -1.0, m+ii, 0, res_m+ii, 0, res_m+ii, 0);

		}

	res->res_mu = mu*nct_inv;

	return;

	}



void COMPUTE_LIN_RES_OCP_QP(struct OCP_QP *qp, struct OCP_QP_SOL *qp_sol, struct OCP_QP_SOL *qp_step, struct OCP_QP_RES *res, struct OCP_QP_RES_WORKSPACE *ws)
	{

	// loop index
	int ii;

	//
	int N = qp->dim->N;
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;
	int *nb = qp->dim->nb;
	int *ng = qp->dim->ng;
	int *ns = qp->dim->ns;

	int nct = 0;
	for(ii=0; ii<=N; ii++)
		nct += 2*nb[ii]+2*ng[ii]+2*ns[ii];
	
	REAL nct_inv = 1.0/nct;

	struct STRMAT *BAbt = qp->BAbt;
	struct STRMAT *RSQrq = qp->RSQrq;
	struct STRMAT *DCt = qp->DCt;
	struct STRVEC *b = qp->b;
	struct STRVEC *rqz = qp->rqz;
	struct STRVEC *d = qp->d;
	struct STRVEC *m = qp->m;
	int **idxb = qp->idxb;
	struct STRVEC *Z = qp->Z;
	int **idxs = qp->idxs;

	struct STRVEC *ux = qp_step->ux;
	struct STRVEC *pi = qp_step->pi;
	struct STRVEC *lam = qp_step->lam;
	struct STRVEC *t = qp_step->t;

	struct STRVEC *Lam = qp_sol->lam;
	struct STRVEC *T = qp_sol->t;

	struct STRVEC *res_g = res->res_g;
	struct STRVEC *res_b = res->res_b;
	struct STRVEC *res_d = res->res_d;
	struct STRVEC *res_m = res->res_m;

	struct STRVEC *tmp_nbgM = ws->tmp_nbgM;
	struct STRVEC *tmp_nsM = ws->tmp_nsM;

	int nx0, nx1, nu0, nu1, nb0, ng0, ns0;

	//
	REAL mu = 0.0;

	// loop over stages
	for(ii=0; ii<=N; ii++)
		{

		nx0 = nx[ii];
		nu0 = nu[ii];
		nb0 = nb[ii];
		ng0 = ng[ii];
		ns0 = ns[ii];

		SYMV_L(nu0+nx0, nu0+nx0, 1.0, RSQrq+ii, 0, 0, ux+ii, 0, 1.0, rqz+ii, 0, res_g+ii, 0);

		if(ii>0)
			AXPY(nx0, -1.0, pi+(ii-1), 0, res_g+ii, nu0, res_g+ii, nu0);

		if(nb0+ng0>0)
			{
			AXPY(nb0+ng0, -1.0, lam+ii, 0, lam+ii, nb[ii]+ng[ii], tmp_nbgM+0, 0);
//			AXPY(nb0+ng0,  1.0, d+ii, 0, t+ii, 0, res_d+ii, 0);
//			AXPY(nb0+ng0,  1.0, d+ii, nb0+ng0, t+ii, nb0+ng0, res_d+ii, nb0+ng0);
			AXPY(2*nb0+2*ng0,  1.0, d+ii, 0, t+ii, 0, res_d+ii, 0);
			// box
			if(nb0>0)
				{
				VECAD_SP(nb0, 1.0, tmp_nbgM+0, 0, idxb[ii], res_g+ii, 0);
				VECEX_SP(nb0, 1.0, idxb[ii], ux+ii, 0, tmp_nbgM+1, 0);
				}
			// general
			if(ng0>0)
				{
				GEMV_NT(nu0+nx0, ng0, 1.0, 1.0, DCt+ii, 0, 0, tmp_nbgM+0, nb[ii], ux+ii, 0, 1.0, 0.0, res_g+ii, 0, tmp_nbgM+1, nb0, res_g+ii, 0, tmp_nbgM+1, nb0);
				}

			AXPY(nb0+ng0, -1.0, tmp_nbgM+1, 0, res_d+ii, 0, res_d+ii, 0);
			AXPY(nb0+ng0,  1.0, tmp_nbgM+1, 0, res_d+ii, nb0+ng0, res_d+ii, nb0+ng0);
			}
		if(ns0>0)
			{
			// res_g
			GEMV_DIAG(2*ns0, 1.0, Z+ii, 0, ux+ii, nu0+nx0, 1.0, rqz+ii, nu0+nx0, res_g+ii, nu0+nx0);
			AXPY(2*ns0, -1.0, lam+ii, 2*nb0+2*ng0, res_g+ii, nu0+nx0, res_g+ii, nu0+nx0);
			VECEX_SP(ns0, 1.0, idxs[ii], lam+ii, 0, tmp_nsM, 0);
			AXPY(ns0, -1.0, tmp_nsM, 0, res_g+ii, nu0+nx0, res_g+ii, nu0+nx0);
			VECEX_SP(ns0, 1.0, idxs[ii], lam+ii, nb0+ng0, tmp_nsM, 0);
			AXPY(ns0, -1.0, tmp_nsM, 0, res_g+ii, nu0+nx0+ns0, res_g+ii, nu0+nx0+ns0);
			// res_d
			VECAD_SP(ns0, -1.0, ux+ii, nu0+nx0, idxs[ii], res_d+ii, 0);
			VECAD_SP(ns0, -1.0, ux+ii, nu0+nx0+ns0, idxs[ii], res_d+ii, nb0+ng0);
			AXPY(2*ns0, -1.0, ux+ii, nu0+nx0, t+ii, 2*nb0+2*ng0, res_d+ii, 2*nb0+2*ng0);
			AXPY(2*ns0, 1.0, d+ii, 2*nb0+2*ng0, res_d+ii, 2*nb0+2*ng0, res_d+ii, 2*nb0+2*ng0);
			}

		if(ii<N)
			{

			nu1 = nu[ii+1];
			nx1 = nx[ii+1];

			AXPY(nx1, -1.0, ux+(ii+1), nu1, b+ii, 0, res_b+ii, 0);

			GEMV_NT(nu0+nx0, nx1, 1.0, 1.0, BAbt+ii, 0, 0, pi+ii, 0, ux+ii, 0, 1.0, 1.0, res_g+ii, 0, res_b+ii, 0, res_g+ii, 0, res_b+ii, 0);

			}

//		mu += VECMULDOT(2*nb0+2*ng0+2*ns0, lam+ii, 0, t+ii, 0, res_m+ii, 0);
		VECCP(2*nb0+2*ng0+2*ns0, m+ii, 0, res_m+ii, 0);
		VECMULACC(2*nb0+2*ng0+2*ns0, Lam+ii, 0, t+ii, 0, res_m+ii, 0);
		VECMULACC(2*nb0+2*ng0+2*ns0, lam+ii, 0, T+ii, 0, res_m+ii, 0);

		}

//	res->res_mu = mu*nct_inv;

	return;

	}



// backward Riccati recursion
void FACT_SOLVE_KKT_UNCONSTR_OCP_QP(struct OCP_QP *qp, struct OCP_QP_SOL *qp_sol, struct OCP_QP_IPM_ARG *arg, struct OCP_QP_IPM_WORKSPACE *ws)
	{

	int N = qp->dim->N;
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;
	int *nb = qp->dim->nb;
	int *ng = qp->dim->ng;

	struct STRMAT *BAbt = qp->BAbt;
	struct STRMAT *RSQrq = qp->RSQrq;
	struct STRVEC *b = qp->b;
	struct STRVEC *rqz = qp->rqz;

	struct STRVEC *ux = qp_sol->ux;
	struct STRVEC *pi = qp_sol->pi;

	struct STRMAT *L = ws->L;
	struct STRMAT *AL = ws->AL;
	struct STRVEC *tmp_nxM = ws->tmp_nxM;

	//
	int ii;

	// factorization and backward substitution

	// last stage
	POTRF_L_MN(nu[N]+nx[N]+1, nu[N]+nx[N], RSQrq+N, 0, 0, L+N, 0, 0);

	// middle stages
	for(ii=0; ii<N; ii++)
		{
		ROWIN(nx[N-ii], 1.0, b+N-ii-1, 0, BAbt+N-ii-1, nu[N-ii-1]+nx[N-ii-1], 0);
		TRMM_RLNN(nu[N-ii-1]+nx[N-ii-1]+1, nx[N-ii], 1.0, L+(N-ii), nu[N-ii], nu[N-ii], BAbt+(N-ii-1), 0, 0, AL, 0, 0);
		GEAD(1, nx[N-ii], 1.0, L+(N-ii), nu[N-ii]+nx[N-ii], nu[N-ii], AL, nu[N-ii-1]+nx[N-ii-1], 0);

		ROWIN(nu[N-ii-1]+nx[N-ii-1], 1.0, rqz+N-ii-1, 0, RSQrq+N-ii-1, nu[N-ii-1]+nx[N-ii-1], 0);
		SYRK_POTRF_LN(nu[N-ii-1]+nx[N-ii-1]+1, nu[N-ii-1]+nx[N-ii-1], nx[N-ii], AL, 0, 0, AL, 0, 0, RSQrq+(N-ii-1), 0, 0, L+(N-ii-1), 0, 0);
		}

	// forward substitution

	// first stage
	ii = 0;
	ROWEX(nu[ii]+nx[ii], -1.0, L+(ii), nu[ii]+nx[ii], 0, ux+ii, 0);
	TRSV_LTN(nu[ii]+nx[ii], L+ii, 0, 0, ux+ii, 0, ux+ii, 0);
	GEMV_T(nu[ii]+nx[ii], nx[ii+1], 1.0, BAbt+ii, 0, 0, ux+ii, 0, 1.0, b+ii, 0, ux+(ii+1), nu[ii+1]);
	ROWEX(nx[ii+1], 1.0, L+(ii+1), nu[ii+1]+nx[ii+1], nu[ii+1], tmp_nxM, 0);
	TRMV_LTN(nx[ii+1], nx[ii+1], L+(ii+1), nu[ii+1], nu[ii+1], ux+(ii+1), nu[ii+1], pi+ii, 0);
	AXPY(nx[ii+1], 1.0, tmp_nxM, 0, pi+ii, 0, pi+ii, 0);
	TRMV_LNN(nx[ii+1], nx[ii+1], L+(ii+1), nu[ii+1], nu[ii+1], pi+ii, 0, pi+ii, 0);

	// middle stages
	for(ii=1; ii<N; ii++)
		{
		ROWEX(nu[ii], -1.0, L+(ii), nu[ii]+nx[ii], 0, ux+ii, 0);
		TRSV_LTN_MN(nu[ii]+nx[ii], nu[ii], L+ii, 0, 0, ux+ii, 0, ux+ii, 0);
		GEMV_T(nu[ii]+nx[ii], nx[ii+1], 1.0, BAbt+ii, 0, 0, ux+ii, 0, 1.0, b+ii, 0, ux+(ii+1), nu[ii+1]);
		ROWEX(nx[ii+1], 1.0, L+(ii+1), nu[ii+1]+nx[ii+1], nu[ii+1], tmp_nxM, 0);
		TRMV_LTN(nx[ii+1], nx[ii+1], L+(ii+1), nu[ii+1], nu[ii+1], ux+(ii+1), nu[ii+1], pi+ii, 0);
		AXPY(nx[ii+1], 1.0, tmp_nxM, 0, pi+ii, 0, pi+ii, 0);
		TRMV_LNN(nx[ii+1], nx[ii+1], L+(ii+1), nu[ii+1], nu[ii+1], pi+ii, 0, pi+ii, 0);
		}
	
	ii = N;
	ROWEX(nu[ii], -1.0, L+(ii), nu[ii]+nx[ii], 0, ux+ii, 0);
	TRSV_LTN_MN(nu[ii]+nx[ii], nu[ii], L+ii, 0, 0, ux+ii, 0, ux+ii, 0);

	return;

	}



static void COND_SLACKS_FACT_SOLVE(int ss, struct OCP_QP *qp, struct OCP_QP_IPM_WORKSPACE *ws)
	{

	int ii, idx;

	int nx0 = qp->dim->nx[ss];
	int nu0 = qp->dim->nu[ss];
	int nb0 = qp->dim->nb[ss];
	int ng0 = qp->dim->ng[ss];
	int ns0 = qp->dim->ns[ss];

	struct STRVEC *Z = qp->Z+ss;
	int *idxs0 = qp->idxs[ss];

	struct STRVEC *dux = ws->sol_step->ux+ss;
	struct STRVEC *res_g = ws->res->res_g+ss;
	struct STRVEC *Gamma = ws->Gamma+ss;
	struct STRVEC *gamma = ws->gamma+ss;
	struct STRVEC *Zs_inv = ws->Zs_inv+ss;
	struct STRVEC *tmp_nbgM = ws->tmp_nbgM;

	REAL *ptr_Gamma = Gamma->pa;
	REAL *ptr_gamma = gamma->pa;
	REAL *ptr_Z = Z->pa;
	REAL *ptr_Zs_inv = Zs_inv->pa;
	REAL *ptr_dux = dux->pa;
	REAL *ptr_res_g = res_g->pa;
	REAL *ptr_tmp0 = (tmp_nbgM+0)->pa;
	REAL *ptr_tmp1 = (tmp_nbgM+1)->pa;
	REAL *ptr_tmp2 = (tmp_nbgM+2)->pa;
	REAL *ptr_tmp3 = (tmp_nbgM+3)->pa;

	REAL tmp0, tmp1;

	VECCP(nb0+ng0, Gamma, 0, tmp_nbgM+0, 0);
	VECCP(nb0+ng0, Gamma, nb0+ng0, tmp_nbgM+1, 0);
	VECCP(nb0+ng0, gamma, 0, tmp_nbgM+2, 0);
	VECCP(nb0+ng0, gamma, nb0+ng0, tmp_nbgM+3, 0);

	for(ii=0; ii<ns0; ii++)
		{
		idx = idxs0[ii];
		ptr_Zs_inv[0+ii]   = ptr_Z[0+ii]   + ptr_Gamma[0+idx]       + ptr_Gamma[2*nb0+2*ng0+ii];
		ptr_Zs_inv[ns0+ii] = ptr_Z[ns0+ii] + ptr_Gamma[nb0+ng0+idx] + ptr_Gamma[2*nb0+2*ng0+ns0+ii];
		ptr_dux[nu0+nx0+ii]      = ptr_res_g[nu0+nx0+ii]     + ptr_gamma[0+idx]   + ptr_gamma[2*nb0+2*ng0+ii];
		ptr_dux[nu0+nx0+ns0+ii]  = ptr_res_g[nu0+nx0+ns0+ii] + ptr_gamma[nb0+ng0+idx] + ptr_gamma[2*nb0+2*ng0+ns0+ii];
		ptr_Zs_inv[0+ii]   = 1.0/ptr_Zs_inv[0+ii];
		ptr_Zs_inv[ns0+ii] = 1.0/ptr_Zs_inv[ns0+ii];
		tmp0 = ptr_dux[nu0+nx0+ii]*ptr_Zs_inv[0+ii];
		tmp1 = ptr_dux[nu0+nx0+ns0+ii]*ptr_Zs_inv[ns0+ii];
		ptr_tmp0[idx] = ptr_tmp0[idx] - ptr_tmp0[idx]*ptr_Zs_inv[0+ii]*ptr_tmp0[idx];
		ptr_tmp1[idx] = ptr_tmp1[idx] - ptr_tmp1[idx]*ptr_Zs_inv[ns0+ii]*ptr_tmp1[idx];
		ptr_tmp2[idx] = ptr_tmp2[idx] - ptr_Gamma[0+idx]*tmp0;
		ptr_tmp3[idx] = ptr_tmp3[idx] - ptr_Gamma[nb0+ng0+idx]*tmp1;
		}
	
	AXPY(nb0+ng0,  1.0, tmp_nbgM+1, 0, tmp_nbgM+0, 0, tmp_nbgM+0, 0);
	AXPY(nb0+ng0, -1.0, tmp_nbgM+3, 0, tmp_nbgM+2, 0, tmp_nbgM+1, 0);

	return;

	}



static void COND_SLACKS_SOLVE(int ss, struct OCP_QP *qp, struct OCP_QP_IPM_WORKSPACE *ws)
	{

	int ii, idx;

	int nx0 = qp->dim->nx[ss];
	int nu0 = qp->dim->nu[ss];
	int nb0 = qp->dim->nb[ss];
	int ng0 = qp->dim->ng[ss];
	int ns0 = qp->dim->ns[ss];

	int *idxs0 = qp->idxs[ss];

	struct STRVEC *dux = ws->sol_step->ux+ss;
	struct STRVEC *res_g = ws->res->res_g+ss;
	struct STRVEC *Gamma = ws->Gamma+ss;
	struct STRVEC *gamma = ws->gamma+ss;
	struct STRVEC *Zs_inv = ws->Zs_inv+ss;
	struct STRVEC *tmp_nbgM = ws->tmp_nbgM;

	REAL *ptr_Gamma = Gamma->pa;
	REAL *ptr_gamma = gamma->pa;
	REAL *ptr_Zs_inv = Zs_inv->pa;
	REAL *ptr_dux = dux->pa;
	REAL *ptr_res_g = res_g->pa;
	REAL *ptr_tmp2 = (tmp_nbgM+2)->pa;
	REAL *ptr_tmp3 = (tmp_nbgM+3)->pa;

	REAL tmp0, tmp1;

	VECCP(nb0+ng0, gamma, 0, tmp_nbgM+2, 0);
	VECCP(nb0+ng0, gamma, nb0+ng0, tmp_nbgM+3, 0);

	for(ii=0; ii<ns0; ii++)
		{
		idx = idxs0[ii];
		ptr_dux[nu0+nx0+ii]      = ptr_res_g[nu0+nx0+ii]     + ptr_gamma[0+idx]       + ptr_gamma[2*nb0+2*ng0+ii];
		ptr_dux[nu0+nx0+ns0+ii]  = ptr_res_g[nu0+nx0+ns0+ii] + ptr_gamma[nb0+ng0+idx] + ptr_gamma[2*nb0+2*ng0+ns0+ii];
		tmp0 = ptr_dux[nu0+nx0+ii]*ptr_Zs_inv[0+ii];
		tmp1 = ptr_dux[nu0+nx0+ns0+ii]*ptr_Zs_inv[ns0+ii];
		ptr_tmp2[idx] = ptr_tmp2[idx] - ptr_Gamma[0+idx]*tmp0;
		ptr_tmp3[idx] = ptr_tmp3[idx] - ptr_Gamma[nb0+ng0+idx]*tmp1;
		}
	
	AXPY(nb0+ng0, -1.0, tmp_nbgM+3, 0, tmp_nbgM+2, 0, tmp_nbgM+1, 0);

	return;

	}



static void EXPAND_SLACKS(int ss, struct OCP_QP *qp, struct OCP_QP_SOL *qp_sol, struct OCP_QP_IPM_WORKSPACE *ws)
	{

	int ii, idx;

	int nx0 = qp->dim->nx[ss];
	int nu0 = qp->dim->nu[ss];
	int nb0 = qp->dim->nb[ss];
	int ng0 = qp->dim->ng[ss];
	int ns0 = qp->dim->ns[ss];

	int *idxs0 = qp->idxs[ss];

	struct STRVEC *dux = qp_sol->ux+ss;
	struct STRVEC *dt = qp_sol->t+ss;

	struct STRVEC *Gamma = ws->Gamma+ss;
	struct STRVEC *Zs_inv = ws->Zs_inv+ss;

	REAL *ptr_Gamma = Gamma->pa;
	REAL *ptr_dux = dux->pa;
	REAL *ptr_dt = dt->pa;
	REAL *ptr_Zs_inv = Zs_inv->pa;

	for(ii=0; ii<ns0; ii++)
		{
		idx = idxs0[ii];
		ptr_dux[nu0+nx0+ii]     = - ptr_Zs_inv[0+ii]   * (ptr_dux[nu0+nx0+ii]     + ptr_dt[idx]*ptr_Gamma[idx]);
		ptr_dux[nu0+nx0+ns0+ii] = - ptr_Zs_inv[ns0+ii] * (ptr_dux[nu0+nx0+ns0+ii] + ptr_dt[nb0+ng0+idx]*ptr_Gamma[nb0+ng0+idx]);
		ptr_dt[2*nb0+2*ng0+ii]     = ptr_dux[nu0+nx0+ii];
		ptr_dt[2*nb0+2*ng0+ns0+ii] = ptr_dux[nu0+nx0+ns0+ii];
		ptr_dt[0+idx]       = ptr_dt[0+idx]   + ptr_dux[nu0+nx0+ii];
		ptr_dt[nb0+ng0+idx] = ptr_dt[nb0+ng0+idx] + ptr_dux[nu0+nx0+ns0+ii];

		}

	return;

	}



// backward Riccati recursion
void FACT_SOLVE_KKT_STEP_OCP_QP(struct OCP_QP *qp, struct OCP_QP_SOL *qp_sol, struct OCP_QP_IPM_ARG *arg, struct OCP_QP_IPM_WORKSPACE *ws)
	{

	int N = qp->dim->N;
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;
	int *nb = qp->dim->nb;
	int *ng = qp->dim->ng;
	int *ns = qp->dim->ns;

	struct STRMAT *BAbt = qp->BAbt;
	struct STRMAT *RSQrq = qp->RSQrq;
	struct STRMAT *DCt = qp->DCt;
	struct STRVEC *Z = qp->Z;
	struct STRVEC *res_g = qp->rqz;
	struct STRVEC *res_b = qp->b;
	int **idxb = qp->idxb;
	int **idxs = qp->idxs;

	struct STRVEC *dux = qp_sol->ux;
	struct STRVEC *dpi = qp_sol->pi;
	struct STRVEC *dt = qp_sol->t;

	struct STRMAT *L = ws->L;
	struct STRMAT *AL = ws->AL;
//	struct STRVEC *res_b = ws->res->res_b;
//	struct STRVEC *res_g = ws->res->res_g;
	struct STRVEC *Gamma = ws->Gamma;
	struct STRVEC *gamma = ws->gamma;
	struct STRVEC *Pb = ws->Pb;
	struct STRVEC *Zs_inv = ws->Zs_inv;
	struct STRVEC *tmp_nxM = ws->tmp_nxM;
	struct STRVEC *tmp_nbgM = ws->tmp_nbgM;

	REAL *ptr0, *ptr1, *ptr2, *ptr3;

	//
	int ii, nn, ss, idx;

	struct CORE_QP_IPM_WORKSPACE *cws = ws->core_workspace;

	COMPUTE_GAMMA_GAMMA_QP(ws->res->res_d[0].pa, ws->res->res_m[0].pa, cws);

	// factorization and backward substitution

	// last stage
	ss = N;
#if defined(DOUBLE_PRECISION)
	TRCP_L(nu[ss]+nx[ss], RSQrq+ss, 0, 0, L+ss, 0, 0); // TODO blasfeo_dtrcp_l with m and n, for m>=n
#else
	GECP(nu[ss]+nx[ss], nu[ss]+nx[ss], RSQrq+ss, 0, 0, L+ss, 0, 0); // TODO blasfeo_dtrcp_l with m and n, for m>=n
#endif
	DIARE(nu[ss]+nx[ss], arg->reg_prim, L+ss, 0, 0);
	ROWIN(nu[ss]+nx[ss], 1.0, res_g+ss, 0, L+ss, nu[ss]+nx[ss], 0);

	if(ns[ss]>0)
		{
		COND_SLACKS_FACT_SOLVE(ss, qp, ws);
		}
	else if(nb[ss]+ng[ss]>0)
		{
		AXPY(nb[ss]+ng[ss],  1.0, Gamma+ss, nb[ss]+ng[ss], Gamma+ss, 0, tmp_nbgM+0, 0);
		AXPY(nb[ss]+ng[ss], -1.0, gamma+ss, nb[ss]+ng[ss], gamma+ss, 0, tmp_nbgM+1, 0);
		}
	if(nb[ss]>0)
		{
		DIAAD_SP(nb[ss], 1.0, tmp_nbgM+0, 0, idxb[ss], L+ss, 0, 0);
		ROWAD_SP(nb[ss], 1.0, tmp_nbgM+1, 0, idxb[ss], L+ss, nu[ss]+nx[ss], 0);
		}
	if(ng[ss]>0)
		{
		GEMM_R_DIAG(nu[ss]+nx[ss], ng[ss], 1.0, DCt+ss, 0, 0, tmp_nbgM+0, nb[ss], 0.0, AL+0, 0, 0, AL+0, 0, 0);
		ROWIN(ng[ss], 1.0, tmp_nbgM+1, nb[ss], AL+0, nu[ss]+nx[ss], 0);
		SYRK_POTRF_LN(nu[ss]+nx[ss]+1, nu[ss]+nx[ss], ng[ss], AL+0, 0, 0, DCt+ss, 0, 0, L+ss, 0, 0, L+ss, 0, 0);
		}
	else
		{
		POTRF_L_MN(nu[ss]+nx[ss]+1, nu[ss]+nx[ss], L+ss, 0, 0, L+ss, 0, 0);
		}
	
	// middle stages
	for(nn=0; nn<N; nn++)
		{
		ss = N-nn-1;
		ROWIN(nx[ss+1], 1.0, res_b+ss, 0, BAbt+ss, nu[ss]+nx[ss], 0);
		TRMM_RLNN(nu[ss]+nx[ss]+1, nx[ss+1], 1.0, L+ss+1, nu[ss+1], nu[ss+1], BAbt+ss, 0, 0, AL, 0, 0);
		ROWEX(nx[ss+1], 1.0, AL, nu[ss]+nx[ss], 0, tmp_nxM, 0);
		TRMV_LNN(nx[ss+1], nx[ss+1], L+ss+1, nu[ss+1], nu[ss+1], tmp_nxM, 0, Pb+ss, 0);
		GEAD(1, nx[ss+1], 1.0, L+ss+1, nu[ss+1]+nx[ss+1], nu[ss+1], AL, nu[ss]+nx[ss], 0);

#if defined(DOUBLE_PRECISION)
		TRCP_L(nu[ss]+nx[ss], RSQrq+ss, 0, 0, L+ss, 0, 0);
#else
		GECP(nu[ss]+nx[ss], nu[ss]+nx[ss], RSQrq+ss, 0, 0, L+ss, 0, 0);
#endif
		DIARE(nu[ss]+nx[ss], arg->reg_prim, L+ss, 0, 0);
		ROWIN(nu[ss]+nx[ss], 1.0, res_g+ss, 0, L+ss, nu[ss]+nx[ss], 0);

		if(ns[ss]>0)
			{
			COND_SLACKS_FACT_SOLVE(ss, qp, ws);
			}
		else if(nb[ss]+ng[ss]>0)
			{
			AXPY(nb[ss]+ng[ss],  1.0, Gamma+ss, nb[ss]+ng[ss], Gamma+ss, 0, tmp_nbgM+0, 0);
			AXPY(nb[ss]+ng[ss], -1.0, gamma+ss, nb[ss]+ng[ss], gamma+ss, 0, tmp_nbgM+1, 0);
			}
		if(nb[ss]>0)
			{
			DIAAD_SP(nb[ss], 1.0, tmp_nbgM+0, 0, idxb[ss], L+ss, 0, 0);
			ROWAD_SP(nb[ss], 1.0, tmp_nbgM+1, 0, idxb[ss], L+ss, nu[ss]+nx[ss], 0);
			}
		if(ng[ss]>0)
			{
			GEMM_R_DIAG(nu[ss]+nx[ss], ng[ss], 1.0, DCt+ss, 0, 0, tmp_nbgM+0, nb[ss], 0.0, AL+0, 0, nx[ss+1], AL+0, 0, nx[ss+1]);
			ROWIN(ng[ss], 1.0, tmp_nbgM+1, nb[ss], AL+0, nu[ss]+nx[ss], nx[ss+1]);
			GECP(nu[ss]+nx[ss], nx[ss+1], AL+0, 0, 0, AL+1, 0, 0);
			GECP(nu[ss]+nx[ss], ng[ss], DCt+ss, 0, 0, AL+1, 0, nx[ss+1]);
			SYRK_POTRF_LN(nu[ss]+nx[ss]+1, nu[ss]+nx[ss], nx[ss+1]+ng[ss], AL+0, 0, 0, AL+1, 0, 0, L+ss, 0, 0, L+ss, 0, 0);
			}
		else
			{
			SYRK_POTRF_LN(nu[ss]+nx[ss]+1, nu[ss]+nx[ss], nx[ss+1], AL, 0, 0, AL, 0, 0, L+ss, 0, 0, L+ss, 0, 0);
			}

		}

	// forward substitution

	// first stage
	ss = 0;
	ROWEX(nu[ss]+nx[ss], -1.0, L+ss, nu[ss]+nx[ss], 0, dux+ss, 0);
	TRSV_LTN(nu[ss]+nx[ss], L+ss, 0, 0, dux+ss, 0, dux+ss, 0);
	GEMV_T(nu[ss]+nx[ss], nx[ss+1], 1.0, BAbt+ss, 0, 0, dux+ss, 0, 1.0, res_b+ss, 0, dux+ss+1, nu[ss+1]);
	ROWEX(nx[ss+1], 1.0, L+ss+1, nu[ss+1]+nx[ss+1], nu[ss+1], tmp_nxM, 0);
	TRMV_LTN(nx[ss+1], nx[ss+1], L+ss+1, nu[ss+1], nu[ss+1], dux+ss+1, nu[ss+1], dpi+ss, 0);
	AXPY(nx[ss+1], 1.0, tmp_nxM, 0, dpi+ss, 0, dpi+ss, 0);
	TRMV_LNN(nx[ss+1], nx[ss+1], L+ss+1, nu[ss+1], nu[ss+1], dpi+ss, 0, dpi+ss, 0);

	// middle stages
	for(ss=1; ss<N; ss++)
		{
		ROWEX(nu[ss], -1.0, L+ss, nu[ss]+nx[ss], 0, dux+ss, 0);
		TRSV_LTN_MN(nu[ss]+nx[ss], nu[ss], L+ss, 0, 0, dux+ss, 0, dux+ss, 0);
		GEMV_T(nu[ss]+nx[ss], nx[ss+1], 1.0, BAbt+ss, 0, 0, dux+ss, 0, 1.0, res_b+ss, 0, dux+(ss+1), nu[ss+1]);
		ROWEX(nx[ss+1], 1.0, L+ss+1, nu[ss+1]+nx[ss+1], nu[ss+1], tmp_nxM, 0);
		TRMV_LTN(nx[ss+1], nx[ss+1], L+ss+1, nu[ss+1], nu[ss+1], dux+ss+1, nu[ss+1], dpi+ss, 0);
		AXPY(nx[ss+1], 1.0, tmp_nxM, 0, dpi+ss, 0, dpi+ss, 0);
		TRMV_LNN(nx[ss+1], nx[ss+1], L+ss+1, nu[ss+1], nu[ss+1], dpi+ss, 0, dpi+ss, 0);
		}

	ss = N;
	ROWEX(nu[ss], -1.0, L+ss, nu[ss]+nx[ss], 0, dux+ss, 0);
	TRSV_LTN_MN(nu[ss]+nx[ss], nu[ss], L+ss, 0, 0, dux+ss, 0, dux+ss, 0);


	for(ss=0; ss<=N; ss++)
		VECEX_SP(nb[ss], 1.0, idxb[ss], dux+ss, 0, dt+ss, 0);
	for(ss=0; ss<=N; ss++)
		GEMV_T(nu[ss]+nx[ss], ng[ss], 1.0, DCt+ss, 0, 0, dux+ss, 0, 0.0, dt+ss, nb[ss], dt+ss, nb[ss]);
	
	for(ss=0; ss<=N; ss++)
		{
		VECCP(nb[ss]+ng[ss], dt+ss, 0, dt+ss, nb[ss]+ng[ss]);
		VECSC(nb[ss]+ng[ss], -1.0, dt+ss, nb[ss]+ng[ss]);
		}

	for(ss=0; ss<=N; ss++)
		{
		if(ns[ss]>0)
			EXPAND_SLACKS(ss, qp, qp_sol, ws);
		}

	COMPUTE_LAM_T_QP(ws->res->res_d[0].pa, ws->res->res_m[0].pa, cws->dlam, cws->dt, cws);

	return;

	}



void FACT_SOLVE_LQ_KKT_STEP_OCP_QP(struct OCP_QP *qp, struct OCP_QP_SOL *qp_sol, struct OCP_QP_IPM_ARG *arg, struct OCP_QP_IPM_WORKSPACE *ws)
	{

	int N = qp->dim->N;
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;
	int *nb = qp->dim->nb;
	int *ng = qp->dim->ng;
	int *ns = qp->dim->ns;

	struct STRMAT *BAbt = qp->BAbt;
	struct STRMAT *RSQrq = qp->RSQrq;
	struct STRMAT *DCt = qp->DCt;
	struct STRVEC *Z = qp->Z;
	struct STRVEC *res_g = qp->rqz;
	struct STRVEC *res_b = qp->b;
	int **idxb = qp->idxb;
	int **idxs = qp->idxs;

	struct STRVEC *dux = qp_sol->ux;
	struct STRVEC *dpi = qp_sol->pi;
	struct STRVEC *dt = qp_sol->t;

	struct STRMAT *L = ws->L;
	struct STRMAT *Lh = ws->Lh;
	struct STRMAT *AL = ws->AL;
//	struct STRVEC *res_b = ws->res->res_b;
//	struct STRVEC *res_g = ws->res->res_g;
	struct STRVEC *Gamma = ws->Gamma;
	struct STRVEC *gamma = ws->gamma;
	struct STRVEC *Pb = ws->Pb;
	struct STRVEC *Zs_inv = ws->Zs_inv;
	struct STRVEC *tmp_nxM = ws->tmp_nxM;
	struct STRVEC *tmp_nbgM = ws->tmp_nbgM;
	struct STRMAT *lq0 = ws->lq0;
	void *lq_work0 = ws->lq_work0;

	REAL *ptr0, *ptr1, *ptr2, *ptr3;

	REAL tmp;

	//
	int ii, nn, ss, idx;

	struct CORE_QP_IPM_WORKSPACE *cws = ws->core_workspace;

	COMPUTE_GAMMA_GAMMA_QP(ws->res->res_d[0].pa, ws->res->res_m[0].pa, cws);

	// factorization and backward substitution

	// last stage
	ss = N;

	VECCP(nu[ss]+nx[ss], res_g+ss, 0, dux+ss, 0);

//	GESE(nu[ss]+nx[ss], 2*nu[ss]+2*nx[ss]+ng[ss], 0.0, lq0, 0, 0);
	GESE(nu[ss]+nx[ss], nu[ss]+nx[ss]+ng[ss], 0.0, lq0, 0, nu[ss]+nx[ss]);

	if(ns[ss]>0)
		{
		COND_SLACKS_FACT_SOLVE(ss, qp, ws);
		}
	else if(nb[ss]+ng[ss]>0)
		{
		AXPY(nb[ss]+ng[ss],  1.0, Gamma+ss, nb[ss]+ng[ss], Gamma+ss, 0, tmp_nbgM+0, 0);
		AXPY(nb[ss]+ng[ss], -1.0, gamma+ss, nb[ss]+ng[ss], gamma+ss, 0, tmp_nbgM+1, 0);
		}
	if(nb[ss]>0)
		{
		for(ii=0; ii<nb[ss]; ii++)
			{
			tmp = BLASFEO_DVECEL(tmp_nbgM+0, ii);
			tmp = tmp>=0.0 ? tmp : 0.0;
			tmp = sqrt( tmp );
			BLASFEO_DMATEL(lq0, idxb[ss][ii], nu[ss]+nx[ss]+idxb[ss][ii]) = tmp>0.0 ? tmp : 0.0;
			}
		VECAD_SP(nb[ss], 1.0, tmp_nbgM+1, 0, idxb[ss], dux+ss, 0);
		}
	if(ng[ss]>0)
		{
		for(ii=0; ii<ng[ss]; ii++)
			{
			tmp = BLASFEO_DVECEL(tmp_nbgM+0, nb[ss]+ii);
			tmp = tmp>=0.0 ? tmp : 0.0;
			tmp = sqrt( tmp );
			BLASFEO_DVECEL(tmp_nbgM+0, nb[ss]+ii) = tmp;
			}
		GEMM_R_DIAG(nu[ss]+nx[ss], ng[ss], 1.0, DCt+ss, 0, 0, tmp_nbgM+0, nb[ss], 0.0, lq0, 0, 2*nu[ss]+2*nx[ss], lq0, 0, 2*nu[ss]+2*nx[ss]);
		GEMV_N(nu[ss]+nx[ss], ng[ss], 1.0, DCt+ss, 0, 0, tmp_nbgM+1, nb[ss], 1.0, dux+ss, 0, dux+ss, 0);
		}

	if(ws->use_hess_fact[ss]==0)
		{
		POTRF_L(nu[ss]+nx[ss], RSQrq+ss, 0, 0, Lh+ss, 0, 0);
		ws->use_hess_fact[ss]==1;
		}

	DIARE(nu[ss]+nx[ss], arg->reg_prim, lq0, 0, nu[ss]+nx[ss]);
#if defined(LA_HIGH_PERFORMANCE) | defined(LA_REFERENCE)
	TRCP_L(nu[ss]+nx[ss], Lh+ss, 0, 0, L+ss, 0, 0);
	GELQF_PD_LLA(nu[ss]+nx[ss], ng[ss], L+ss, 0, 0, lq0, 0, nu[ss]+nx[ss], lq0, 0, 2*nu[ss]+2*nx[ss], lq_work0); // TODO reduce lq1 size !!!
#else
	TRCP_L(nu[ss]+nx[ss], Lh+ss, 0, 0, lq0, 0, 0);
	GELQF(nu[ss]+nx[ss], 2*nu[ss]+2*nx[ss]+ng[ss], lq0, 0, 0, lq0, 0, 0, lq_work0);
	TRCP_L(nu[ss]+nx[ss], lq0, 0, 0, L+ss, 0, 0);
	for(ii=0; ii<nu[ss]+nx[ss]; ii++)
		if(BLASFEO_DMATEL(L+ss, ii, ii) < 0)
			COLSC(nu[ss]+nx[ss]-ii, -1.0, L+ss, ii, ii);
#endif

	TRSV_LNN_MN(nu[ss]+nx[ss], nu[ss], L+ss, 0, 0, dux+ss, 0, dux+ss, 0);


	// middle stages
	for(nn=0; nn<N-1; nn++)
		{
		ss = N-nn-1;

//		GESE(nu[ss]+nx[ss], 2*nu[ss]+2*nx[ss]+nx[ss+1]+ng[ss], 0.0, lq0, 0, 0);
		GESE(nu[ss]+nx[ss], nu[ss]+nx[ss]+ng[ss], 0.0, lq0, 0, nu[ss]+nx[ss]);

		TRMM_RLNN(nu[ss]+nx[ss], nx[ss+1], 1.0, L+ss+1, nu[ss+1], nu[ss+1], BAbt+ss, 0, 0, lq0, 0, 2*nu[ss]+2*nx[ss]+ng[ss]);
		TRMV_LTN(nx[ss+1], nx[ss+1], L+ss+1, nu[ss+1], nu[ss+1], res_b+ss, 0, Pb+ss, 0);
		TRMV_LNN(nx[ss+1], nx[ss+1], L+ss+1, nu[ss+1], nu[ss+1], Pb+ss, 0, Pb+ss, 0);

		VECCP(nu[ss]+nx[ss], res_g+ss, 0, dux+ss, 0);
		AXPY(nx[ss+1], 1.0, dux+ss+1, nu[ss+1], Pb+ss, 0, tmp_nxM, 0);
		GEMV_N(nu[ss]+nx[ss], nx[ss+1], 1.0, BAbt+ss, 0, 0, tmp_nxM, 0, 1.0, dux+ss, 0, dux+ss, 0);

		if(ns[ss]>0)
			{
			COND_SLACKS_FACT_SOLVE(ss, qp, ws);
			}
		else if(nb[ss]+ng[ss]>0)
			{
			AXPY(nb[ss]+ng[ss],  1.0, Gamma+ss, nb[ss]+ng[ss], Gamma+ss, 0, tmp_nbgM+0, 0);
			AXPY(nb[ss]+ng[ss], -1.0, gamma+ss, nb[ss]+ng[ss], gamma+ss, 0, tmp_nbgM+1, 0);
			}
		if(nb[ss]>0)
			{
			for(ii=0; ii<nb[ss]; ii++)
				{
				tmp = BLASFEO_DVECEL(tmp_nbgM+0, ii);
				tmp = tmp>=0.0 ? tmp : 0.0;
				tmp = sqrt( tmp );
				BLASFEO_DMATEL(lq0, idxb[ss][ii], nu[ss]+nx[ss]+idxb[ss][ii]) = tmp>0.0 ? tmp : 0.0;
				}
			VECAD_SP(nb[ss], 1.0, tmp_nbgM+1, 0, idxb[ss], dux+ss, 0);
			}
		if(ng[ss]>0)
			{
			for(ii=0; ii<ng[ss]; ii++)
				{
				tmp = BLASFEO_DVECEL(tmp_nbgM+0, nb[ss]+ii);
				tmp = tmp>=0.0 ? tmp : 0.0;
				tmp = sqrt( tmp );
				BLASFEO_DVECEL(tmp_nbgM+0, nb[ss]+ii) = tmp;
				}
			GEMM_R_DIAG(nu[ss]+nx[ss], ng[ss], 1.0, DCt+ss, 0, 0, tmp_nbgM+0, nb[ss], 0.0, lq0, 0, 2*nu[ss]+2*nx[ss], lq0, 0, 2*nu[ss]+2*nx[ss]);
			GEMV_N(nu[ss]+nx[ss], ng[ss], 1.0, DCt+ss, 0, 0, tmp_nbgM+1, nb[ss], 1.0, dux+ss, 0, dux+ss, 0);
			}

		if(ws->use_hess_fact[ss]==0)
			{
			POTRF_L(nu[ss]+nx[ss], RSQrq+ss, 0, 0, Lh+ss, 0, 0);
			ws->use_hess_fact[ss]==1;
			}

		DIARE(nu[ss]+nx[ss], arg->reg_prim, lq0, 0, nu[ss]+nx[ss]);
#if defined(LA_HIGH_PERFORMANCE) | defined(LA_REFERENCE)
		TRCP_L(nu[ss]+nx[ss], Lh+ss, 0, 0, L+ss, 0, 0);
		GELQF_PD_LLA(nu[ss]+nx[ss], nx[ss+1]+ng[ss], L+ss, 0, 0, lq0, 0, nu[ss]+nx[ss], lq0, 0, 2*nu[ss]+2*nx[ss], lq_work0); // TODO reduce lq1 size !!!
#else
		TRCP_L(nu[ss]+nx[ss], Lh+ss, 0, 0, lq0, 0, 0);
		GELQF(nu[ss]+nx[ss], 2*nu[ss]+2*nx[ss]+nx[ss+1]+ng[ss], lq0, 0, 0, lq0, 0, 0, lq_work0);
		TRCP_L(nu[ss]+nx[ss], lq0, 0, 0, L+ss, 0, 0);
		for(ii=0; ii<nu[ss]+nx[ss]; ii++)
			if(BLASFEO_DMATEL(L+ss, ii, ii) < 0)
				COLSC(nu[ss]+nx[ss]-ii, -1.0, L+ss, ii, ii);
#endif

		TRSV_LNN_MN(nu[ss]+nx[ss], nu[ss], L+ss, 0, 0, dux+ss, 0, dux+ss, 0);

		}

	// first stage
	nn = N-1;
	ss = N-nn-1;

//	GESE(nu[ss]+nx[ss], 2*nu[ss]+2*nx[ss]+nx[ss+1]+ng[ss], 0.0, lq0, 0, 0);
	GESE(nu[ss]+nx[ss], nu[ss]+nx[ss]+ng[ss], 0.0, lq0, 0, nu[ss]+nx[ss]);

	TRMM_RLNN(nu[ss]+nx[ss], nx[ss+1], 1.0, L+ss+1, nu[ss+1], nu[ss+1], BAbt+ss, 0, 0, lq0, 0, 2*nu[ss]+2*nx[ss]+ng[ss]);
TRMM_RLNN(nu[ss]+nx[ss], nx[ss+1], 1.0, L+ss+1, nu[ss+1], nu[ss+1], BAbt+ss, 0, 0, AL, 0, 0);
	TRMV_LTN(nx[ss+1], nx[ss+1], L+ss+1, nu[ss+1], nu[ss+1], res_b+ss, 0, Pb+ss, 0);
	TRMV_LNN(nx[ss+1], nx[ss+1], L+ss+1, nu[ss+1], nu[ss+1], Pb+ss, 0, Pb+ss, 0);

	VECCP(nu[ss]+nx[ss], res_g+ss, 0, dux+ss, 0);
	AXPY(nx[ss+1], 1.0, dux+ss+1, nu[ss+1], Pb+ss, 0, tmp_nxM, 0);
	GEMV_N(nu[ss]+nx[ss], nx[ss+1], 1.0, BAbt+ss, 0, 0, tmp_nxM, 0, 1.0, dux+ss, 0, dux+ss, 0);

	if(ns[ss]>0)
		{
		COND_SLACKS_FACT_SOLVE(ss, qp, ws);
		}
	else if(nb[ss]+ng[ss]>0)
		{
		AXPY(nb[ss]+ng[ss],  1.0, Gamma+ss, nb[ss]+ng[ss], Gamma+ss, 0, tmp_nbgM+0, 0);
		AXPY(nb[ss]+ng[ss], -1.0, gamma+ss, nb[ss]+ng[ss], gamma+ss, 0, tmp_nbgM+1, 0);
		}
	if(nb[ss]>0)
		{
		for(ii=0; ii<nb[ss]; ii++)
			{
			tmp = BLASFEO_DVECEL(tmp_nbgM+0, ii);
			tmp = tmp>=0.0 ? tmp : 0.0;
			tmp = sqrt( tmp );
			BLASFEO_DMATEL(lq0, idxb[ss][ii], nu[ss]+nx[ss]+idxb[ss][ii]) = tmp>0.0 ? tmp : 0.0;
			}
		VECAD_SP(nb[ss], 1.0, tmp_nbgM+1, 0, idxb[ss], dux+ss, 0);
		}
	if(ng[ss]>0)
		{
		for(ii=0; ii<ng[ss]; ii++)
			{
			tmp = BLASFEO_DVECEL(tmp_nbgM+0, nb[ss]+ii);
			tmp = tmp>=0.0 ? tmp : 0.0;
			tmp = sqrt( tmp );
			BLASFEO_DVECEL(tmp_nbgM+0, nb[ss]+ii) = tmp;
			}
		GEMM_R_DIAG(nu[ss]+nx[ss], ng[ss], 1.0, DCt+ss, 0, 0, tmp_nbgM+0, nb[ss], 0.0, lq0, 0, 2*nu[ss]+2*nx[ss], lq0, 0, 2*nu[ss]+2*nx[ss]);
		GEMV_N(nu[ss]+nx[ss], ng[ss], 1.0, DCt+ss, 0, 0, tmp_nbgM+1, nb[ss], 1.0, dux+ss, 0, dux+ss, 0);
		}

	if(ws->use_hess_fact[ss]==0)
		{
		POTRF_L(nu[ss]+nx[ss], RSQrq+ss, 0, 0, Lh+ss, 0, 0);
		ws->use_hess_fact[ss]==1;
		}

	DIARE(nu[ss]+nx[ss], arg->reg_prim, lq0, 0, nu[ss]+nx[ss]);
#if defined(LA_HIGH_PERFORMANCE) | defined(LA_REFERENCE)
	TRCP_L(nu[ss]+nx[ss], Lh+ss, 0, 0, L+ss, 0, 0);
	GELQF_PD_LLA(nu[ss]+nx[ss], nx[ss+1]+ng[ss], L+ss, 0, 0, lq0, 0, nu[ss]+nx[ss], lq0, 0, 2*nu[ss]+2*nx[ss], lq_work0); // TODO reduce lq1 size !!!
#else
	TRCP_L(nu[ss]+nx[ss], Lh+ss, 0, 0, lq0, 0, 0);
	GELQF(nu[ss]+nx[ss], 2*nu[ss]+2*nx[ss]+nx[ss+1]+ng[ss], lq0, 0, 0, lq0, 0, 0, lq_work0);
	TRCP_L(nu[ss]+nx[ss], lq0, 0, 0, L+ss, 0, 0);
	for(ii=0; ii<nu[ss]+nx[ss]; ii++)
		if(BLASFEO_DMATEL(L+ss, ii, ii) < 0)
			COLSC(nu[ss]+nx[ss]-ii, -1.0, L+ss, ii, ii);
#endif

	TRSV_LNN(nu[ss]+nx[ss], L+ss, 0, 0, dux+ss, 0, dux+ss, 0);

	// forward substitution

	// first stage
	ss = 0;
	VECCP(nx[ss+1], dux+ss+1, nu[ss+1], dpi+ss, 0);
	VECSC(nu[ss]+nx[ss], -1.0, dux+ss, 0);
	TRSV_LTN(nu[ss]+nx[ss], L+ss, 0, 0, dux+ss, 0, dux+ss, 0);
	GEMV_T(nu[ss]+nx[ss], nx[ss+1], 1.0, BAbt+ss, 0, 0, dux+ss, 0, 1.0, res_b+ss, 0, dux+ss+1, nu[ss+1]);
	VECCP(nx[ss+1], dux+ss+1, nu[ss+1], tmp_nxM, 0);
	TRMV_LTN(nx[ss+1], nx[ss+1], L+ss+1, nu[ss+1], nu[ss+1], tmp_nxM, 0, tmp_nxM, 0);
	TRMV_LNN(nx[ss+1], nx[ss+1], L+ss+1, nu[ss+1], nu[ss+1], tmp_nxM, 0, tmp_nxM, 0);
	AXPY(nx[ss+1], 1.0, tmp_nxM, 0, dpi+ss, 0, dpi+ss, 0);

	// middle stages
	for(ss=1; ss<N; ss++)
		{
		VECCP(nx[ss+1], dux+ss+1, nu[ss+1], dpi+ss, 0);
		VECSC(nu[ss], -1.0, dux+ss, 0);
		TRSV_LTN_MN(nu[ss]+nx[ss], nu[ss], L+ss, 0, 0, dux+ss, 0, dux+ss, 0);
		GEMV_T(nu[ss]+nx[ss], nx[ss+1], 1.0, BAbt+ss, 0, 0, dux+ss, 0, 1.0, res_b+ss, 0, dux+ss+1, nu[ss+1]);
		VECCP(nx[ss+1], dux+ss+1, nu[ss+1], tmp_nxM, 0);
		TRMV_LTN(nx[ss+1], nx[ss+1], L+ss+1, nu[ss+1], nu[ss+1], tmp_nxM, 0, tmp_nxM, 0);
		TRMV_LNN(nx[ss+1], nx[ss+1], L+ss+1, nu[ss+1], nu[ss+1], tmp_nxM, 0, tmp_nxM, 0);
		AXPY(nx[ss+1], 1.0, tmp_nxM, 0, dpi+ss, 0, dpi+ss, 0);
		}

	ss = N;
	VECSC(nu[ss], -1.0, dux+ss, 0);
	TRSV_LTN_MN(nu[ss]+nx[ss], nu[ss], L+ss, 0, 0, dux+ss, 0, dux+ss, 0);


	for(ss=0; ss<=N; ss++)
		VECEX_SP(nb[ss], 1.0, idxb[ss], dux+ss, 0, dt+ss, 0);
	for(ss=0; ss<=N; ss++)
		GEMV_T(nu[ss]+nx[ss], ng[ss], 1.0, DCt+ss, 0, 0, dux+ss, 0, 0.0, dt+ss, nb[ss], dt+ss, nb[ss]);
	
	for(ss=0; ss<=N; ss++)
		{
		VECCP(nb[ss]+ng[ss], dt+ss, 0, dt+ss, nb[ss]+ng[ss]);
		VECSC(nb[ss]+ng[ss], -1.0, dt+ss, nb[ss]+ng[ss]);
		}

	for(ss=0; ss<=N; ss++)
		{
		if(ns[ss]>0)
			EXPAND_SLACKS(ss, qp, qp_sol, ws);
		}

	COMPUTE_LAM_T_QP(ws->res->res_d[0].pa, ws->res->res_m[0].pa, cws->dlam, cws->dt, cws);

	return;

	}



// backward Riccati recursion
void SOLVE_KKT_STEP_OCP_QP(struct OCP_QP *qp, struct OCP_QP_SOL *qp_sol, struct OCP_QP_IPM_ARG *arg, struct OCP_QP_IPM_WORKSPACE *ws)
	{

	int N = qp->dim->N;
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;
	int *nb = qp->dim->nb;
	int *ng = qp->dim->ng;
	int *ns = qp->dim->ns;

	struct STRMAT *BAbt = qp->BAbt;
//	struct STRMAT *RSQrq = qp->RSQrq;
	struct STRMAT *DCt = qp->DCt;
	struct STRVEC *res_g = qp->rqz;
	struct STRVEC *res_b = qp->b;
	int **idxb = qp->idxb;
//	int **idxs = qp->idxs;

	struct STRVEC *dux = qp_sol->ux;
	struct STRVEC *dpi = qp_sol->pi;
	struct STRVEC *dt = qp_sol->t;

	struct STRMAT *L = ws->L;
//	struct STRVEC *res_b = ws->res->res_b;
//	struct STRVEC *res_g = ws->res->res_g;
	struct STRVEC *gamma = ws->gamma;
	struct STRVEC *Pb = ws->Pb;
	struct STRVEC *tmp_nxM = ws->tmp_nxM;
	struct STRVEC *tmp_nbgM = ws->tmp_nbgM;

	//
	int ss, nn, ii;

	struct CORE_QP_IPM_WORKSPACE *cws = ws->core_workspace;

	COMPUTE_GAMMA_QP(ws->res->res_d[0].pa, ws->res->res_m[0].pa, cws);

	// backward substitution

	// last stage
	ss = N;
	VECCP(nu[ss]+nx[ss], res_g+ss, 0, dux+ss, 0);
	if(ns[ss]>0)
		{
		COND_SLACKS_SOLVE(ss, qp, ws);
		}
	else if(nb[ss]+ng[ss]>0)
		{
		AXPY(nb[ss]+ng[ss], -1.0, gamma+ss, nb[ss]+ng[ss], gamma+ss, 0, tmp_nbgM+1, 0);
		}
	if(nb[ss]>0)
		{
		VECAD_SP(nb[ss], 1.0, tmp_nbgM+1, 0, idxb[ss], dux+ss, 0);
		}
	if(ng[ss]>0)
		{
		GEMV_N(nu[ss]+nx[ss], ng[ss], 1.0, DCt+ss, 0, 0, tmp_nbgM+1, nb[ss], 1.0, dux+ss, 0, dux+ss, 0);
		}
	TRSV_LNN_MN(nu[ss]+nx[ss], nu[ss], L+ss, 0, 0, dux+ss, 0, dux+ss, 0);

	// middle stages
	for(nn=0; nn<N-1; nn++)
		{
		ss = N-nn-1;
		VECCP(nu[ss]+nx[ss], res_g+ss, 0, dux+ss, 0);
		if(ns[ss]>0)
			{
			COND_SLACKS_SOLVE(ss, qp, ws);
			}
		else if(nb[ss]+ng[ss]>0)
			{
			AXPY(nb[ss]+ng[ss], -1.0, gamma+ss, nb[ss]+ng[ss], gamma+ss, 0, tmp_nbgM+1, 0);
			}
		if(nb[ss]>0)
			{
			VECAD_SP(nb[ss], 1.0, tmp_nbgM+1, 0, idxb[ss], dux+ss, 0);
			}
		if(ng[ss]>0)
			{
			GEMV_N(nu[ss]+nx[ss], ng[ss], 1.0, DCt+ss, 0, 0, tmp_nbgM+1, nb[ss], 1.0, dux+ss, 0, dux+ss, 0);
			}
		AXPY(nx[ss+1], 1.0, dux+ss+1, nu[ss+1], Pb+ss, 0, tmp_nxM, 0);
		GEMV_N(nu[ss]+nx[ss], nx[ss+1], 1.0, BAbt+ss, 0, 0, tmp_nxM, 0, 1.0, dux+ss, 0, dux+ss, 0);
		TRSV_LNN_MN(nu[ss]+nx[ss], nu[ss], L+ss, 0, 0, dux+ss, 0, dux+ss, 0);
		}

	// first stage
	nn = N-1;
	ss = N-nn-1;
	VECCP(nu[ss]+nx[ss], res_g+ss, 0, dux+ss, 0);
	if(ns[ss]>0)
		{
		COND_SLACKS_SOLVE(ss, qp, ws);
		}
	else if(nb[ss]+ng[ss]>0)
		{
		AXPY(nb[ss]+ng[ss], -1.0, gamma+ss, nb[ss]+ng[ss], gamma+ss, 0, tmp_nbgM+1, 0);
		}
	if(nb[ss]>0)
		{
		VECAD_SP(nb[ss], 1.0, tmp_nbgM+1, 0, idxb[ss], dux+ss, 0);
		}
	if(ng[ss]>0)
		{
		GEMV_N(nu[ss]+nx[ss], ng[ss], 1.0, DCt+ss, 0, 0, tmp_nbgM+1, nb[ss], 1.0, dux+ss, 0, dux+ss, 0);
		}
	AXPY(nx[ss+1], 1.0, dux+ss+1, nu[ss+1], Pb+ss, 0, tmp_nxM, 0);
	GEMV_N(nu[ss]+nx[ss], nx[ss+1], 1.0, BAbt+ss, 0, 0, tmp_nxM, 0, 1.0, dux+ss, 0, dux+ss, 0);
	TRSV_LNN(nu[ss]+nx[ss], L+ss, 0, 0, dux+ss, 0, dux+ss, 0);

	// forward substitution

	// first stage
	ss = 0;
	VECCP(nx[ss+1], dux+ss+1, nu[ss+1], dpi+ss, 0);
	VECSC(nu[ss]+nx[ss], -1.0, dux+ss, 0);
	TRSV_LTN(nu[ss]+nx[ss], L+ss, 0, 0, dux+ss, 0, dux+ss, 0);
	GEMV_T(nu[ss]+nx[ss], nx[ss+1], 1.0, BAbt+ss, 0, 0, dux+ss, 0, 1.0, res_b+ss, 0, dux+ss+1, nu[ss+1]);
	VECCP(nx[ss+1], dux+ss+1, nu[ss+1], tmp_nxM, 0);
	TRMV_LTN(nx[ss+1], nx[ss+1], L+ss+1, nu[ss+1], nu[ss+1], tmp_nxM, 0, tmp_nxM, 0);
	TRMV_LNN(nx[ss+1], nx[ss+1], L+ss+1, nu[ss+1], nu[ss+1], tmp_nxM, 0, tmp_nxM, 0);
	AXPY(nx[ss+1], 1.0, tmp_nxM, 0, dpi+ss, 0, dpi+ss, 0);

	// middle stages
	for(ss=1; ss<N; ss++)
		{
		VECCP(nx[ss+1], dux+ss+1, nu[ss+1], dpi+ss, 0);
		VECSC(nu[ss], -1.0, dux+ss, 0);
		TRSV_LTN_MN(nu[ss]+nx[ss], nu[ss], L+ss, 0, 0, dux+ss, 0, dux+ss, 0);
		GEMV_T(nu[ss]+nx[ss], nx[ss+1], 1.0, BAbt+ss, 0, 0, dux+ss, 0, 1.0, res_b+ss, 0, dux+ss+1, nu[ss+1]);
		VECCP(nx[ss+1], dux+ss+1, nu[ss+1], tmp_nxM, 0);
		TRMV_LTN(nx[ss+1], nx[ss+1], L+ss+1, nu[ss+1], nu[ss+1], tmp_nxM, 0, tmp_nxM, 0);
		TRMV_LNN(nx[ss+1], nx[ss+1], L+ss+1, nu[ss+1], nu[ss+1], tmp_nxM, 0, tmp_nxM, 0);
		AXPY(nx[ss+1], 1.0, tmp_nxM, 0, dpi+ss, 0, dpi+ss, 0);
		}

	ss = N;
	VECSC(nu[ss], -1.0, dux+ss, 0);
	TRSV_LTN_MN(nu[ss]+nx[ss], nu[ss], L+ss, 0, 0, dux+ss, 0, dux+ss, 0);



	for(ss=0; ss<=N; ss++)
		VECEX_SP(nb[ss], 1.0, idxb[ss], dux+ss, 0, dt+ss, 0);
	for(ss=0; ss<=N; ss++)
		GEMV_T(nu[ss]+nx[ss], ng[ss], 1.0, DCt+ss, 0, 0, dux+ss, 0, 0.0, dt+ss, nb[ss], dt+ss, nb[ss]);

	for(ss=0; ss<=N; ss++)
		{
		VECCP(nb[ss]+ng[ss], dt+ss, 0, dt+ss, nb[ss]+ng[ss]);
		VECSC(nb[ss]+ng[ss], -1.0, dt+ss, nb[ss]+ng[ss]);
		}

	for(ss=0; ss<=N; ss++)
		{
		if(ns[ss]>0)
			EXPAND_SLACKS(ss, qp, qp_sol, ws);
		}

	COMPUTE_LAM_T_QP(ws->res->res_d[0].pa, ws->res->res_m[0].pa, cws->dlam, cws->dt, cws);

	return;

	}


