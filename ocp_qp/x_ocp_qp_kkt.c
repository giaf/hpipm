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
	int N = qp->N;
	int *nx = qp->nx;
	int *nu = qp->nu;
	int *nb = qp->nb;
	int *ng = qp->ng;
	int *ns = qp->ns;

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
			t_ub[jj] =   d_ub[jj] - ux[idxb[jj]];
			if(t_lb[jj]<thr0)
				{
				if(t_ub[jj]<thr0)
					{
					ux[idxb[jj]] = 0.5*(d_lb[jj]-d_ub[jj]);
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
				ux[idxb[jj]] = d_ub[jj] - thr0;
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
		GEMV_T_LIBSTR(nu[ii]+nx[ii], ng[ii], 1.0, qp->DCt+ii, 0, 0, qp_sol->ux+ii, 0, 0.0, qp_sol->t+ii, nb[ii], qp_sol->t+ii, nb[ii]);
		for(jj=0; jj<ng[ii]; jj++)
			{
#if 1
			t_ug[jj] = - t_lg[jj];
			t_lg[jj] -= d_lg[jj];
			t_ug[jj] += d_ug[jj];
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



void COMPUTE_RES_OCP_QP(struct OCP_QP *qp, struct OCP_QP_SOL *qp_sol, struct OCP_QP_IPM_WORKSPACE *ws)
	{

	struct CORE_QP_IPM_WORKSPACE *cws = ws->core_workspace;
	
	// loop index
	int ii;

	//
	int N = qp->N;
	int *nx = qp->nx;
	int *nu = qp->nu;
	int *nb = qp->nb;
	int *ng = qp->ng;
	int *ns = qp->ns;

	int nct = ws->core_workspace->nc;

	struct STRMAT *BAbt = qp->BAbt;
	struct STRMAT *RSQrq = qp->RSQrq;
	struct STRMAT *DCt = qp->DCt;
	struct STRVEC *b = qp->b;
	struct STRVEC *rq = qp->rq;
	struct STRVEC *d = qp->d;
	int **idxb = qp->idxb;
	struct STRVEC *Z = qp->Z;
	struct STRVEC *z = qp->z;
	int **idxs = qp->idxs;

	struct STRVEC *ux = qp_sol->ux;
	struct STRVEC *pi = qp_sol->pi;
	struct STRVEC *lam = qp_sol->lam;
	struct STRVEC *t = qp_sol->t;

	struct STRVEC *res_g = ws->res_g;
	struct STRVEC *res_b = ws->res_b;
	struct STRVEC *res_d = ws->res_d;
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

		SYMV_L_LIBSTR(nu0+nx0, nu0+nx0, 1.0, RSQrq+ii, 0, 0, ux+ii, 0, 1.0, rq+ii, 0, res_g+ii, 0);

		if(ii>0)
			AXPY_LIBSTR(nx0, -1.0, pi+(ii-1), 0, res_g+ii, nu0, res_g+ii, nu0);

		if(nb0+ng0>0)
			{
			AXPY_LIBSTR(nb0+ng0, -1.0, lam+ii, 0, lam+ii, nb[ii]+ng[ii], tmp_nbgM+0, 0);
			AXPY_LIBSTR(nb0+ng0,  1.0, d+ii, 0, t+ii, 0, res_d+ii, 0);
			AXPY_LIBSTR(nb0+ng0, -1.0, d+ii, nb0+ng0, t+ii, nb0+ng0, res_d+ii, nb0+ng0);
			// box
			if(nb0>0)
				{
				VECAD_SP_LIBSTR(nb0, 1.0, tmp_nbgM+0, 0, idxb[ii], res_g+ii, 0);
				VECEX_SP_LIBSTR(nb0, 1.0, idxb[ii], ux+ii, 0, tmp_nbgM+1, 0);
				}
			// general
			if(ng0>0)
				{
				GEMV_NT_LIBSTR(nu0+nx0, ng0, 1.0, 1.0, DCt+ii, 0, 0, tmp_nbgM+0, nb[ii], ux+ii, 0, 1.0, 0.0, res_g+ii, 0, tmp_nbgM+1, nb0, res_g+ii, 0, tmp_nbgM+1, nb0);
				}

			AXPY_LIBSTR(nb0+ng0, -1.0, tmp_nbgM+1, 0, res_d+ii, 0, res_d+ii, 0);
			AXPY_LIBSTR(nb0+ng0,  1.0, tmp_nbgM+1, 0, res_d+ii, nb0+ng0, res_d+ii, nb0+ng0);
			}
		if(ns0>0)
			{
			// res_g
			GEMV_DIAG_LIBSTR(2*ns0, 1.0, Z+ii, 0, ux+ii, nu0+nx0, 1.0, z+ii, 0, res_g+ii, nu0+nx0);
			AXPY_LIBSTR(2*ns0, -1.0, lam+ii, 2*nb0+2*ng0, res_g+ii, nu0+nx0, res_g+ii, nu0+nx0);
			VECEX_SP_LIBSTR(ns0, 1.0, idxs[ii], lam+ii, 0, tmp_nsM, 0);
			AXPY_LIBSTR(ns0, -1.0, tmp_nsM, 0, res_g+ii, nu0+nx0, res_g+ii, nu0+nx0);
			VECEX_SP_LIBSTR(ns0, 1.0, idxs[ii], lam+ii, nb0+ng0, tmp_nsM, 0);
			AXPY_LIBSTR(ns0, -1.0, tmp_nsM, 0, res_g+ii, nu0+nx0+ns0, res_g+ii, nu0+nx0+ns0);
			// res_d
			VECAD_SP_LIBSTR(ns0, -1.0, ux+ii, nu0+nx0, idxs[ii], res_d+ii, 0);
			VECAD_SP_LIBSTR(ns0, -1.0, ux+ii, nu0+nx0+ns0, idxs[ii], res_d+ii, nb0+ng0);
			AXPY_LIBSTR(2*ns0, -1.0, ux+ii, nu0+nx0, t+ii, 2*nb0+2*ng0, res_d+ii, 2*nb0+2*ng0);
			}

		if(ii<N)
			{

			nu1 = nu[ii+1];
			nx1 = nx[ii+1];

			AXPY_LIBSTR(nx1, -1.0, ux+(ii+1), nu1, b+ii, 0, res_b+ii, 0);

			GEMV_NT_LIBSTR(nu0+nx0, nx1, 1.0, 1.0, BAbt+ii, 0, 0, pi+ii, 0, ux+ii, 0, 1.0, 1.0, res_g+ii, 0, res_b+ii, 0, res_g+ii, 0, res_b+ii, 0);

			}

		}

	mu += VECMULDOT_LIBSTR(nct, lam, 0, t, 0, ws->res_m, 0);

	ws->res_mu = mu*cws->nc_inv;

	return;

	}



// backward Riccati recursion
void FACT_SOLVE_KKT_UNCONSTR_OCP_QP(struct OCP_QP *qp, struct OCP_QP_SOL *qp_sol, struct OCP_QP_IPM_WORKSPACE *ws)
	{

	int N = qp->N;
	int *nx = qp->nx;
	int *nu = qp->nu;
	int *nb = qp->nb;
	int *ng = qp->ng;

	struct STRMAT *BAbt = qp->BAbt;
	struct STRMAT *RSQrq = qp->RSQrq;
	struct STRVEC *b = qp->b;

	struct STRVEC *ux = qp_sol->ux;
	struct STRVEC *pi = qp_sol->pi;

	struct STRMAT *L = ws->L;
	struct STRMAT *AL = ws->AL;
	struct STRVEC *tmp_nxM = ws->tmp_nxM;

	//
	int ii;

	// factorization and backward substitution

	// last stage
	POTRF_L_MN_LIBSTR(nu[N]+nx[N]+1, nu[N]+nx[N], RSQrq+N, 0, 0, L+N, 0, 0);

	// middle stages
	for(ii=0; ii<N; ii++)
		{
		TRMM_RLNN_LIBSTR(nu[N-ii-1]+nx[N-ii-1]+1, nx[N-ii], 1.0, L+(N-ii), nu[N-ii], nu[N-ii], BAbt+(N-ii-1), 0, 0, AL, 0, 0);
		GEAD_LIBSTR(1, nx[N-ii], 1.0, L+(N-ii), nu[N-ii]+nx[N-ii], nu[N-ii], AL, nu[N-ii-1]+nx[N-ii-1], 0);

		SYRK_POTRF_LN_LIBSTR(nu[N-ii-1]+nx[N-ii-1]+1, nu[N-ii-1]+nx[N-ii-1], nx[N-ii], AL, 0, 0, AL, 0, 0, RSQrq+(N-ii-1), 0, 0, L+(N-ii-1), 0, 0);
		}

	// forward substitution

	// first stage
	ii = 0;
	ROWEX_LIBSTR(nu[ii]+nx[ii], -1.0, L+(ii), nu[ii]+nx[ii], 0, ux+ii, 0);
	TRSV_LTN_LIBSTR(nu[ii]+nx[ii], L+ii, 0, 0, ux+ii, 0, ux+ii, 0);
	GEMV_T_LIBSTR(nu[ii]+nx[ii], nx[ii+1], 1.0, BAbt+ii, 0, 0, ux+ii, 0, 1.0, b+ii, 0, ux+(ii+1), nu[ii+1]);
	ROWEX_LIBSTR(nx[ii+1], 1.0, L+(ii+1), nu[ii+1]+nx[ii+1], nu[ii+1], tmp_nxM, 0);
	TRMV_LTN_LIBSTR(nx[ii+1], nx[ii+1], L+(ii+1), nu[ii+1], nu[ii+1], ux+(ii+1), nu[ii+1], pi+ii, 0);
	AXPY_LIBSTR(nx[ii+1], 1.0, tmp_nxM, 0, pi+ii, 0, pi+ii, 0);
	TRMV_LNN_LIBSTR(nx[ii+1], nx[ii+1], L+(ii+1), nu[ii+1], nu[ii+1], pi+ii, 0, pi+ii, 0);

	// middle stages
	for(ii=1; ii<N; ii++)
		{
		ROWEX_LIBSTR(nu[ii], -1.0, L+(ii), nu[ii]+nx[ii], 0, ux+ii, 0);
		TRSV_LTN_MN_LIBSTR(nu[ii]+nx[ii], nu[ii], L+ii, 0, 0, ux+ii, 0, ux+ii, 0);
		GEMV_T_LIBSTR(nu[ii]+nx[ii], nx[ii+1], 1.0, BAbt+ii, 0, 0, ux+ii, 0, 1.0, b+ii, 0, ux+(ii+1), nu[ii+1]);
		ROWEX_LIBSTR(nx[ii+1], 1.0, L+(ii+1), nu[ii+1]+nx[ii+1], nu[ii+1], tmp_nxM, 0);
		TRMV_LTN_LIBSTR(nx[ii+1], nx[ii+1], L+(ii+1), nu[ii+1], nu[ii+1], ux+(ii+1), nu[ii+1], pi+ii, 0);
		AXPY_LIBSTR(nx[ii+1], 1.0, tmp_nxM, 0, pi+ii, 0, pi+ii, 0);
		TRMV_LNN_LIBSTR(nx[ii+1], nx[ii+1], L+(ii+1), nu[ii+1], nu[ii+1], pi+ii, 0, pi+ii, 0);
		}
	
	ii = N;
	ROWEX_LIBSTR(nu[ii], -1.0, L+(ii), nu[ii]+nx[ii], 0, ux+ii, 0);
	TRSV_LTN_MN_LIBSTR(nu[ii]+nx[ii], nu[ii], L+ii, 0, 0, ux+ii, 0, ux+ii, 0);

	return;

	}



static void COND_SLACKS_FACT_SOLVE(int ss, struct OCP_QP *qp, struct OCP_QP_IPM_WORKSPACE *ws)
	{

	int ii, idx;

	int nx0 = qp->nx[ss];
	int nu0 = qp->nu[ss];
	int nb0 = qp->nb[ss];
	int ng0 = qp->ng[ss];
	int ns0 = qp->ns[ss];

	struct STRVEC *Z = qp->Z+ss;
	int *idxs0 = qp->idxs[ss];

	struct STRVEC *dux = ws->dux+ss;
	struct STRVEC *res_g = ws->res_g+ss;
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

	VECCP_LIBSTR(nb0+ng0, Gamma, 0, tmp_nbgM+0, 0);
	VECCP_LIBSTR(nb0+ng0, Gamma, nb0+ng0, tmp_nbgM+1, 0);
	VECCP_LIBSTR(nb0+ng0, gamma, 0, tmp_nbgM+2, 0);
	VECCP_LIBSTR(nb0+ng0, gamma, nb0+ng0, tmp_nbgM+3, 0);

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
	
	AXPY_LIBSTR(nb0+ng0,  1.0, tmp_nbgM+1, 0, tmp_nbgM+0, 0, tmp_nbgM+0, 0);
	AXPY_LIBSTR(nb0+ng0, -1.0, tmp_nbgM+3, 0, tmp_nbgM+2, 0, tmp_nbgM+1, 0);

	return;

	}



static void COND_SLACKS_SOLVE(int ss, struct OCP_QP *qp, struct OCP_QP_IPM_WORKSPACE *ws)
	{

	int ii, idx;

	int nx0 = qp->nx[ss];
	int nu0 = qp->nu[ss];
	int nb0 = qp->nb[ss];
	int ng0 = qp->ng[ss];
	int ns0 = qp->ns[ss];

	int *idxs0 = qp->idxs[ss];

	struct STRVEC *dux = ws->dux+ss;
	struct STRVEC *res_g = ws->res_g+ss;
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

	VECCP_LIBSTR(nb0+ng0, gamma, 0, tmp_nbgM+2, 0);
	VECCP_LIBSTR(nb0+ng0, gamma, nb0+ng0, tmp_nbgM+3, 0);

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
	
	AXPY_LIBSTR(nb0+ng0, -1.0, tmp_nbgM+3, 0, tmp_nbgM+2, 0, tmp_nbgM+1, 0);

	return;

	}



static void EXPAND_SLACKS(int ss, struct OCP_QP *qp, struct OCP_QP_IPM_WORKSPACE *ws)
	{

	int ii, idx;

	int nx0 = qp->nx[ss];
	int nu0 = qp->nu[ss];
	int nb0 = qp->nb[ss];
	int ng0 = qp->ng[ss];
	int ns0 = qp->ns[ss];

	int *idxs0 = qp->idxs[ss];

	struct STRVEC *dux = ws->dux+ss;
	struct STRVEC *dt = ws->dt+ss;
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
void FACT_SOLVE_KKT_STEP_OCP_QP(struct OCP_QP *qp, struct OCP_QP_IPM_WORKSPACE *ws)
	{

	int N = qp->N;
	int *nx = qp->nx;
	int *nu = qp->nu;
	int *nb = qp->nb;
	int *ng = qp->ng;
	int *ns = qp->ns;

	struct STRMAT *BAbt = qp->BAbt;
	struct STRMAT *RSQrq = qp->RSQrq;
	struct STRMAT *DCt = qp->DCt;
	struct STRVEC *Z = qp->Z;
	struct STRVEC *z = qp->z;
	int **idxb = qp->idxb;
	int **idxs = qp->idxs;

	struct STRMAT *L = ws->L;
	struct STRMAT *AL = ws->AL;
	struct STRVEC *res_b = ws->res_b;
	struct STRVEC *res_g = ws->res_g;
	struct STRVEC *dux = ws->dux;
	struct STRVEC *dpi = ws->dpi;
	struct STRVEC *dt = ws->dt;
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

	COMPUTE_GAMMA_GAMMA_QP(cws);

	// factorization and backward substitution

	// last stage
	ss = N;
#if defined(DOUBLE_PRECISION)
	TRCP_L_LIBSTR(nu[ss]+nx[ss], RSQrq+ss, 0, 0, L+ss, 0, 0); // TODO dtrcp_l_libstr with m and n, for m>=n
#else
	GECP_LIBSTR(nu[ss]+nx[ss], nu[ss]+nx[ss], RSQrq+ss, 0, 0, L+ss, 0, 0); // TODO dtrcp_l_libstr with m and n, for m>=n
#endif
	ROWIN_LIBSTR(nu[ss]+nx[ss], 1.0, res_g+ss, 0, L+ss, nu[ss]+nx[ss], 0);

	if(ns[ss]>0)
		{
		COND_SLACKS_FACT_SOLVE(ss, qp, ws);
		}
	else if(nb[ss]+ng[ss]>0)
		{
		AXPY_LIBSTR(nb[ss]+ng[ss],  1.0, Gamma+ss, nb[ss]+ng[ss], Gamma+ss, 0, tmp_nbgM+0, 0);
		AXPY_LIBSTR(nb[ss]+ng[ss], -1.0, gamma+ss, nb[ss]+ng[ss], gamma+ss, 0, tmp_nbgM+1, 0);
		}
	if(nb[ss]>0)
		{
		DIAAD_SP_LIBSTR(nb[ss], 1.0, tmp_nbgM+0, 0, idxb[ss], L+ss, 0, 0);
		ROWAD_SP_LIBSTR(nb[ss], 1.0, tmp_nbgM+1, 0, idxb[ss], L+ss, nu[ss]+nx[ss], 0);
		}
	if(ng[ss]>0)
		{
		GEMM_R_DIAG_LIBSTR(nu[ss]+nx[ss], ng[ss], 1.0, DCt+ss, 0, 0, tmp_nbgM+0, nb[ss], 0.0, AL+0, 0, 0, AL+0, 0, 0);
		ROWIN_LIBSTR(ng[ss], 1.0, tmp_nbgM+1, nb[ss], AL+0, nu[ss]+nx[ss], 0);
		SYRK_POTRF_LN_LIBSTR(nu[ss]+nx[ss]+1, nu[ss]+nx[ss], ng[ss], AL+0, 0, 0, DCt+ss, 0, 0, L+ss, 0, 0, L+ss, 0, 0);
		}
	else
		{
		POTRF_L_MN_LIBSTR(nu[ss]+nx[ss]+1, nu[ss]+nx[ss], L+ss, 0, 0, L+ss, 0, 0);
		}
	
	// middle stages
	for(nn=0; nn<N; nn++)
		{
		ss = N-nn-1;
		GECP_LIBSTR(nu[ss]+nx[ss], nx[ss+1], BAbt+ss, 0, 0, AL, 0, 0);
		ROWIN_LIBSTR(nx[ss+1], 1.0, res_b+ss, 0, AL, nu[ss]+nx[ss], 0);
		TRMM_RLNN_LIBSTR(nu[ss]+nx[ss]+1, nx[ss+1], 1.0, L+ss+1, nu[ss+1], nu[ss+1], AL, 0, 0, AL, 0, 0);
		ROWEX_LIBSTR(nx[ss+1], 1.0, AL, nu[ss]+nx[ss], 0, tmp_nxM, 0);
		TRMV_LNN_LIBSTR(nx[ss+1], nx[ss+1], L+ss+1, nu[ss+1], nu[ss+1], tmp_nxM, 0, Pb+ss, 0);
		GEAD_LIBSTR(1, nx[ss+1], 1.0, L+ss+1, nu[ss+1]+nx[ss+1], nu[ss+1], AL, nu[ss]+nx[ss], 0);

#if defined(DOUBLE_PRECISION)
		TRCP_L_LIBSTR(nu[ss]+nx[ss], RSQrq+ss, 0, 0, L+ss, 0, 0);
#else
		GECP_LIBSTR(nu[ss]+nx[ss], nu[ss]+nx[ss], RSQrq+ss, 0, 0, L+ss, 0, 0);
#endif
		ROWIN_LIBSTR(nu[ss]+nx[ss], 1.0, res_g+ss, 0, L+ss, nu[ss]+nx[ss], 0);

		if(ns[ss]>0)
			{
			COND_SLACKS_FACT_SOLVE(ss, qp, ws);
			}
		else if(nb[ss]+ng[ss]>0)
			{
			AXPY_LIBSTR(nb[ss]+ng[ss],  1.0, Gamma+ss, nb[ss]+ng[ss], Gamma+ss, 0, tmp_nbgM+0, 0);
			AXPY_LIBSTR(nb[ss]+ng[ss], -1.0, gamma+ss, nb[ss]+ng[ss], gamma+ss, 0, tmp_nbgM+1, 0);
			}
		if(nb[ss]>0)
			{
			DIAAD_SP_LIBSTR(nb[ss], 1.0, tmp_nbgM+0, 0, idxb[ss], L+ss, 0, 0);
			ROWAD_SP_LIBSTR(nb[ss], 1.0, tmp_nbgM+1, 0, idxb[ss], L+ss, nu[ss]+nx[ss], 0);
			}
		if(ng[ss]>0)
			{
			GEMM_R_DIAG_LIBSTR(nu[ss]+nx[ss], ng[ss], 1.0, DCt+ss, 0, 0, tmp_nbgM+0, nb[ss], 0.0, AL+0, 0, nx[ss+1], AL+0, 0, nx[ss+1]);
			ROWIN_LIBSTR(ng[ss], 1.0, tmp_nbgM+1, nb[ss], AL+0, nu[ss]+nx[ss], nx[ss+1]);
			GECP_LIBSTR(nu[ss]+nx[ss], nx[ss+1], AL+0, 0, 0, AL+1, 0, 0);
			GECP_LIBSTR(nu[ss]+nx[ss], ng[ss], DCt+ss, 0, 0, AL+1, 0, nx[ss+1]);
			SYRK_POTRF_LN_LIBSTR(nu[ss]+nx[ss]+1, nu[ss]+nx[ss], nx[ss+1]+ng[ss], AL+0, 0, 0, AL+1, 0, 0, L+ss, 0, 0, L+ss, 0, 0);
			}
		else
			{
			SYRK_POTRF_LN_LIBSTR(nu[ss]+nx[ss]+1, nu[ss]+nx[ss], nx[ss+1], AL, 0, 0, AL, 0, 0, L+ss, 0, 0, L+ss, 0, 0);
			}

		}

	// forward substitution

	// first stage
	ss = 0;
	ROWEX_LIBSTR(nu[ss]+nx[ss], -1.0, L+ss, nu[ss]+nx[ss], 0, dux+ss, 0);
	TRSV_LTN_LIBSTR(nu[ss]+nx[ss], L+ss, 0, 0, dux+ss, 0, dux+ss, 0);
	GEMV_T_LIBSTR(nu[ss]+nx[ss], nx[ss+1], 1.0, BAbt+ss, 0, 0, dux+ss, 0, 1.0, res_b+ss, 0, dux+ss+1, nu[ss+1]);
	ROWEX_LIBSTR(nx[ss+1], 1.0, L+ss+1, nu[ss+1]+nx[ss+1], nu[ss+1], tmp_nxM, 0);
	TRMV_LTN_LIBSTR(nx[ss+1], nx[ss+1], L+ss+1, nu[ss+1], nu[ss+1], dux+ss+1, nu[ss+1], dpi+ss, 0);
	AXPY_LIBSTR(nx[ss+1], 1.0, tmp_nxM, 0, dpi+ss, 0, dpi+ss, 0);
	TRMV_LNN_LIBSTR(nx[ss+1], nx[ss+1], L+ss+1, nu[ss+1], nu[ss+1], dpi+ss, 0, dpi+ss, 0);

	// middle stages
	for(ss=1; ss<N; ss++)
		{
		ROWEX_LIBSTR(nu[ss], -1.0, L+ss, nu[ss]+nx[ss], 0, dux+ss, 0);
		TRSV_LTN_MN_LIBSTR(nu[ss]+nx[ss], nu[ss], L+ss, 0, 0, dux+ss, 0, dux+ss, 0);
		GEMV_T_LIBSTR(nu[ss]+nx[ss], nx[ss+1], 1.0, BAbt+ss, 0, 0, dux+ss, 0, 1.0, res_b+ss, 0, dux+(ss+1), nu[ss+1]);
		ROWEX_LIBSTR(nx[ss+1], 1.0, L+ss+1, nu[ss+1]+nx[ss+1], nu[ss+1], tmp_nxM, 0);
		TRMV_LTN_LIBSTR(nx[ss+1], nx[ss+1], L+ss+1, nu[ss+1], nu[ss+1], dux+ss+1, nu[ss+1], dpi+ss, 0);
		AXPY_LIBSTR(nx[ss+1], 1.0, tmp_nxM, 0, dpi+ss, 0, dpi+ss, 0);
		TRMV_LNN_LIBSTR(nx[ss+1], nx[ss+1], L+ss+1, nu[ss+1], nu[ss+1], dpi+ss, 0, dpi+ss, 0);
		}

	ss = N;
	ROWEX_LIBSTR(nu[ss], -1.0, L+ss, nu[ss]+nx[ss], 0, dux+ss, 0);
	TRSV_LTN_MN_LIBSTR(nu[ss]+nx[ss], nu[ss], L+ss, 0, 0, dux+ss, 0, dux+ss, 0);


	for(ss=0; ss<=N; ss++)
		VECEX_SP_LIBSTR(nb[ss], 1.0, idxb[ss], dux+ss, 0, dt+ss, 0);
	for(ss=0; ss<=N; ss++)
		GEMV_T_LIBSTR(nu[ss]+nx[ss], ng[ss], 1.0, DCt+ss, 0, 0, dux+ss, 0, 0.0, dt+ss, nb[ss], dt+ss, nb[ss]);
	
	for(ss=0; ss<=N; ss++)
		{
		VECCP_LIBSTR(nb[ss]+ng[ss], dt+ss, 0, dt+ss, nb[ss]+ng[ss]);
		VECSC_LIBSTR(nb[ss]+ng[ss], -1.0, dt+ss, nb[ss]+ng[ss]);
		}

	for(ss=0; ss<=N; ss++)
		{
		if(ns[ss]>0)
			EXPAND_SLACKS(ss, qp, ws);
		}

	COMPUTE_LAM_T_QP(cws);

	return;

	}



// backward Riccati recursion
void SOLVE_KKT_STEP_OCP_QP(struct OCP_QP *qp, struct OCP_QP_IPM_WORKSPACE *ws)
	{

	int N = qp->N;
	int *nx = qp->nx;
	int *nu = qp->nu;
	int *nb = qp->nb;
	int *ng = qp->ng;
	int *ns = qp->ns;

	struct STRMAT *BAbt = qp->BAbt;
	struct STRMAT *RSQrq = qp->RSQrq;
	struct STRMAT *DCt = qp->DCt;
	int **idxb = qp->idxb;
	int **idxs = qp->idxs;

	struct STRMAT *L = ws->L;
	struct STRMAT *AL = ws->AL;
	struct STRVEC *res_b = ws->res_b;
	struct STRVEC *res_g = ws->res_g;
	struct STRVEC *dux = ws->dux;
	struct STRVEC *dpi = ws->dpi;
	struct STRVEC *dt = ws->dt;
	struct STRVEC *gamma = ws->gamma;
	struct STRVEC *Pb = ws->Pb;
	struct STRVEC *tmp_nxM = ws->tmp_nxM;
	struct STRVEC *tmp_nbgM = ws->tmp_nbgM;

	//
	int ss, nn, ii;

	struct CORE_QP_IPM_WORKSPACE *cws = ws->core_workspace;

	COMPUTE_GAMMA_QP(cws);

	// backward substitution

	// last stage
	ss = N;
	VECCP_LIBSTR(nu[ss]+nx[ss], res_g+ss, 0, dux+ss, 0);
	if(ns[ss]>0)
		{
		COND_SLACKS_SOLVE(ss, qp, ws);
		}
	else if(nb[ss]+ng[ss]>0)
		{
		AXPY_LIBSTR(nb[ss]+ng[ss], -1.0, gamma+ss, nb[ss]+ng[ss], gamma+ss, 0, tmp_nbgM+1, 0);
		}
	if(nb[ss]>0)
		{
		VECAD_SP_LIBSTR(nb[ss], 1.0, tmp_nbgM+1, 0, idxb[ss], dux+ss, 0);
		}
	if(ng[ss]>0)
		{
		GEMV_N_LIBSTR(nu[ss]+nx[ss], ng[ss], 1.0, DCt+ss, 0, 0, tmp_nbgM+1, nb[ss], 1.0, dux+ss, 0, dux+ss, 0);
		}
	TRSV_LNN_MN_LIBSTR(nu[ss]+nx[ss], nu[ss], L+ss, 0, 0, dux+ss, 0, dux+ss, 0);

	// middle stages
	for(nn=0; nn<N-1; nn++)
		{
		ss = N-nn-1;
		VECCP_LIBSTR(nu[ss]+nx[ss], res_g+ss, 0, dux+ss, 0);
		if(ns[ss]>0)
			{
			COND_SLACKS_SOLVE(ss, qp, ws);
			}
		else if(nb[ss]+ng[ss]>0)
			{
			AXPY_LIBSTR(nb[ss]+ng[ss], -1.0, gamma+ss, nb[ss]+ng[ss], gamma+ss, 0, tmp_nbgM+1, 0);
			}
		if(nb[ss]>0)
			{
			VECAD_SP_LIBSTR(nb[ss], 1.0, tmp_nbgM+1, 0, idxb[ss], dux+ss, 0);
			}
		if(ng[ss]>0)
			{
			GEMV_N_LIBSTR(nu[ss]+nx[ss], ng[ss], 1.0, DCt+ss, 0, 0, tmp_nbgM+1, nb[ss], 1.0, dux+ss, 0, dux+ss, 0);
			}
		AXPY_LIBSTR(nx[ss+1], 1.0, dux+ss+1, nu[ss+1], Pb+ss, 0, tmp_nxM, 0);
		GEMV_N_LIBSTR(nu[ss]+nx[ss], nx[ss+1], 1.0, BAbt+ss, 0, 0, tmp_nxM, 0, 1.0, dux+ss, 0, dux+ss, 0);
		TRSV_LNN_MN_LIBSTR(nu[ss]+nx[ss], nu[ss], L+ss, 0, 0, dux+ss, 0, dux+ss, 0);
		}

	// first stage
	nn = N-1;
	ss = N-nn-1;
	VECCP_LIBSTR(nu[ss]+nx[ss], res_g+ss, 0, dux+ss, 0);
	if(ns[ss]>0)
		{
		COND_SLACKS_SOLVE(ss, qp, ws);
		}
	else if(nb[ss]+ng[ss]>0)
		{
		AXPY_LIBSTR(nb[ss]+ng[ss], -1.0, gamma+ss, nb[ss]+ng[ss], gamma+ss, 0, tmp_nbgM+1, 0);
		}
	if(nb[ss]>0)
		{
		VECAD_SP_LIBSTR(nb[ss], 1.0, tmp_nbgM+1, 0, idxb[ss], dux+ss, 0);
		}
	if(ng[ss]>0)
		{
		GEMV_N_LIBSTR(nu[ss]+nx[ss], ng[ss], 1.0, DCt+ss, 0, 0, tmp_nbgM+1, nb[ss], 1.0, dux+ss, 0, dux+ss, 0);
		}
	AXPY_LIBSTR(nx[ss+1], 1.0, dux+ss+1, nu[ss+1], Pb+ss, 0, tmp_nxM, 0);
	GEMV_N_LIBSTR(nu[ss]+nx[ss], nx[ss+1], 1.0, BAbt+ss, 0, 0, tmp_nxM, 0, 1.0, dux+ss, 0, dux+ss, 0);
	TRSV_LNN_LIBSTR(nu[ss]+nx[ss], L+ss, 0, 0, dux+ss, 0, dux+ss, 0);

	// forward substitution

	// first stage
	ss = 0;
	VECCP_LIBSTR(nx[ss+1], dux+ss+1, nu[ss+1], dpi+ss, 0);
	VECSC_LIBSTR(nu[ss]+nx[ss], -1.0, dux+ss, 0);
	TRSV_LTN_LIBSTR(nu[ss]+nx[ss], L+ss, 0, 0, dux+ss, 0, dux+ss, 0);
	GEMV_T_LIBSTR(nu[ss]+nx[ss], nx[ss+1], 1.0, BAbt+ss, 0, 0, dux+ss, 0, 1.0, res_b+ss, 0, dux+ss+1, nu[ss+1]);
	VECCP_LIBSTR(nx[ss+1], dux+ss+1, nu[ss+1], tmp_nxM, 0);
	TRMV_LTN_LIBSTR(nx[ss+1], nx[ss+1], L+ss+1, nu[ss+1], nu[ss+1], tmp_nxM, 0, tmp_nxM, 0);
	TRMV_LNN_LIBSTR(nx[ss+1], nx[ss+1], L+ss+1, nu[ss+1], nu[ss+1], tmp_nxM, 0, tmp_nxM, 0);
	AXPY_LIBSTR(nx[ss+1], 1.0, tmp_nxM, 0, dpi+ss, 0, dpi+ss, 0);

	// middle stages
	for(ss=1; ss<N; ss++)
		{
		VECCP_LIBSTR(nx[ss+1], dux+ss+1, nu[ss+1], dpi+ss, 0);
		VECSC_LIBSTR(nu[ss], -1.0, dux+ss, 0);
		TRSV_LTN_MN_LIBSTR(nu[ss]+nx[ss], nu[ss], L+ss, 0, 0, dux+ss, 0, dux+ss, 0);
		GEMV_T_LIBSTR(nu[ss]+nx[ss], nx[ss+1], 1.0, BAbt+ss, 0, 0, dux+ss, 0, 1.0, res_b+ss, 0, dux+ss+1, nu[ss+1]);
		VECCP_LIBSTR(nx[ss+1], dux+ss+1, nu[ss+1], tmp_nxM, 0);
		TRMV_LTN_LIBSTR(nx[ss+1], nx[ss+1], L+ss+1, nu[ss+1], nu[ss+1], tmp_nxM, 0, tmp_nxM, 0);
		TRMV_LNN_LIBSTR(nx[ss+1], nx[ss+1], L+ss+1, nu[ss+1], nu[ss+1], tmp_nxM, 0, tmp_nxM, 0);
		AXPY_LIBSTR(nx[ss+1], 1.0, tmp_nxM, 0, dpi+ss, 0, dpi+ss, 0);
		}

	ss = N;
	VECSC_LIBSTR(nu[ss], -1.0, dux+ss, 0);
	TRSV_LTN_MN_LIBSTR(nu[ss]+nx[ss], nu[ss], L+ss, 0, 0, dux+ss, 0, dux+ss, 0);



	for(ss=0; ss<=N; ss++)
		VECEX_SP_LIBSTR(nb[ss], 1.0, idxb[ss], dux+ss, 0, dt+ss, 0);
	for(ss=0; ss<=N; ss++)
		GEMV_T_LIBSTR(nu[ss]+nx[ss], ng[ss], 1.0, DCt+ss, 0, 0, dux+ss, 0, 0.0, dt+ss, nb[ss], dt+ss, nb[ss]);

	for(ss=0; ss<=N; ss++)
		{
		VECCP_LIBSTR(nb[ss]+ng[ss], dt+ss, 0, dt+ss, nb[ss]+ng[ss]);
		VECSC_LIBSTR(nb[ss]+ng[ss], -1.0, dt+ss, nb[ss]+ng[ss]);
		}

	for(ss=0; ss<=N; ss++)
		{
		if(ns[ss]>0)
			EXPAND_SLACKS(ss, qp, ws);
		}

	COMPUTE_LAM_T_QP(cws);

	return;

	}


