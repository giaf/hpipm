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



void INIT_VAR_DENSE_QP(struct DENSE_QP *qp, struct DENSE_QP_SOL *qp_sol, struct DENSE_QP_IPM_WORKSPACE *ws)
	{

	struct CORE_QP_IPM_WORKSPACE *cws = ws->core_workspace;

	// extract cws members
	int nv = qp->dim->nv;
	int ne = qp->dim->ne;
	int nb = qp->dim->nb;
	int ng = qp->dim->ng;
	int ns = qp->dim->ns;

	REAL *d = qp->d->pa;
	int *idxb = qp->idxb;

	REAL *v = qp_sol->v->pa;
	REAL *pi = qp_sol->pi->pa;
	REAL *lam = qp_sol->lam->pa;
	REAL *t = qp_sol->t->pa;

	REAL mu0 = ws->mu0;

	// local variables
	int ii;
	int idxb0;
	REAL thr0 = 0.5;

	// primal variables
	if(ws->warm_start==0)
		{
		// cold start
		for(ii=0; ii<nv+2*ns; ii++)
			{
			v[ii] = 0.0;
			}
		}
		
	// equality constraints
	for(ii=0; ii<ne; ii++)
		{
		pi[ii] = 0.0;
		}
	
	// box constraints
	for(ii=0; ii<nb; ii++)
		{
#if 1
		idxb0 = idxb[ii];
		t[0+ii]     = - d[0+ii]     + v[idxb0];
		t[nb+ng+ii] = - d[nb+ng+ii] - v[idxb0];
		if(t[0+ii]<thr0)
			{
			if(t[nb+ng+ii]<thr0)
				{
				v[idxb0] = 0.5*(d[0+ii] + d[nb+ng+ii]);
				t[0+ii]     = thr0;
				t[nb+ng+ii] = thr0;
				}
			else
				{
				t[0+ii] = thr0;
				v[idxb0] = d[0+ii] + thr0;
				}
			}
		else if(t[nb+ng+ii]<thr0)
			{
			t[nb+ng+ii] = thr0;
			v[idxb0] = - d[nb+ng+ii] - thr0;
			}
#else
		t[0+ii]     = 1.0;
		t[nb+ng+ii] = 1.0;
#endif
		lam[0+ii]     = mu0/t[0+ii];
		lam[nb+ng+ii] = mu0/t[nb+ng+ii];
		}
	
	// general constraints
	GEMV_T(nv, ng, 1.0, qp->Ct, 0, 0, qp_sol->v, 0, 0.0, qp_sol->t, nb, qp_sol->t, nb);
	for(ii=0; ii<ng; ii++)
		{
#if 1
		t[2*nb+ng+ii] = t[nb+ii];
		t[nb+ii]      -= d[nb+ii];
		t[2*nb+ng+ii] -= d[2*nb+ng+ii];
//		t[nb+ii]      = fmax( thr0, t[nb+ii] );
//		t[2*nb+ng+ii] = fmax( thr0, t[2*nb+ng+ii] );
		t[nb+ii]      = thr0>t[nb+ii]      ? thr0 : t[nb+ii];
		t[2*nb+ng+ii] = thr0>t[2*nb+ng+ii] ? thr0 : t[2*nb+ng+ii];
#else
		t[nb+ii]      = 1.0;
		t[2*nb+ng+ii] = 1.0;
#endif
		lam[nb+ii]      = mu0/t[nb+ii];
		lam[2*nb+ng+ii] = mu0/t[2*nb+ng+ii];
		}

	// soft constraints
	for(ii=0; ii<ns; ii++)
		{
		t[2*nb+2*ng+ii]    = 1.0; // thr0;
		t[2*nb+2*ng+ns+ii] = 1.0; // thr0;
		lam[2*nb+2*ng+ii]    = mu0/t[2*nb+2*ng+ii];
		lam[2*nb+2*ng+ns+ii] = mu0/t[2*nb+2*ng+ns+ii];
		}

	return;

	}



void COMPUTE_RES_DENSE_QP(struct DENSE_QP *qp, struct DENSE_QP_SOL *qp_sol, struct DENSE_QP_RES *res, struct DENSE_QP_RES_WORKSPACE *ws)
	{

	int nv = qp->dim->nv;
	int ne = qp->dim->ne;
	int nb = qp->dim->nb;
	int ng = qp->dim->ng;
	int ns = qp->dim->ns;

	int nct = 2*nb+2*ng+2*ns;

	REAL nct_inv = 1.0/nct;

	struct STRMAT *Hg = qp->Hv;
	struct STRMAT *A = qp->A;
	struct STRMAT *Ct = qp->Ct;
	struct STRVEC *gz = qp->gz;
	struct STRVEC *b = qp->b;
	struct STRVEC *d = qp->d;
	struct STRVEC *m = qp->m;
	int *idxb = qp->idxb;
	struct STRVEC *Z = qp->Z;
	int *idxs = qp->idxs;

	struct STRVEC *v = qp_sol->v;
	struct STRVEC *pi = qp_sol->pi;
	struct STRVEC *lam = qp_sol->lam;
	struct STRVEC *t = qp_sol->t;

	struct STRVEC *res_g = res->res_g;
	struct STRVEC *res_b = res->res_b;
	struct STRVEC *res_d = res->res_d;
	struct STRVEC *res_m = res->res_m;

	struct STRVEC *tmp_nbg = ws->tmp_nbg;
	struct STRVEC *tmp_ns = ws->tmp_ns;

	REAL mu;

	// res g
	SYMV_L(nv, nv, 1.0, Hg, 0, 0, v, 0, 1.0, gz, 0, res_g, 0);

	if(nb+ng>0)
		{
		AXPY(nb+ng, -1.0, lam, 0, lam, nb+ng, tmp_nbg+0, 0);
//		AXPY(nb+ng,  1.0, d, 0, t, 0, res_d, 0);
//		AXPY(nb+ng,  1.0, d, nb+ng, t, nb+ng, res_d, nb+ng);
		AXPY(2*nb+2*ng,  1.0, d, 0, t, 0, res_d, 0);
		// box
		if(nb>0)
			{
			VECAD_SP(nb, 1.0, tmp_nbg+0, 0, idxb, res_g, 0);
			VECEX_SP(nb, 1.0, idxb, v, 0, tmp_nbg+1, 0);
			}
		// general
		if(ng>0)
			{
			GEMV_NT(nv, ng, 1.0, 1.0, Ct, 0, 0, tmp_nbg+0, nb, v, 0, 1.0, 0.0, res_g, 0, tmp_nbg+1, nb, res_g, 0, tmp_nbg+1, nb);
			}
		AXPY(nb+ng, -1.0, tmp_nbg+1, 0, res_d, 0, res_d, 0);
		AXPY(nb+ng,  1.0, tmp_nbg+1, 0, res_d, nb+ng, res_d, nb+ng);
		}
	if(ns>0)
		{
		// res_g
		GEMV_DIAG(2*ns, 1.0, Z, 0, v, nv, 1.0, gz, nv, res_g, nv);
		AXPY(2*ns, -1.0, lam, 2*nb+2*ng, res_g, nv, res_g, nv);
		VECEX_SP(ns, 1.0, idxs, lam, 0, tmp_ns, 0);
		AXPY(ns, -1.0, tmp_ns, 0, res_g, nv, res_g, nv);
		VECEX_SP(ns, 1.0, idxs, lam, nb+ng, tmp_ns, 0);
		AXPY(ns, -1.0, tmp_ns, 0, res_g, nv+ns, res_g, nv+ns);
		// res_d
		VECAD_SP(ns, -1.0, v, nv, idxs, res_d, 0);
		VECAD_SP(ns, -1.0, v, nv+ns, idxs, res_d, nb+ng);
		AXPY(2*ns, -1.0, v, nv, t, 2*nb+2*ng, res_d, 2*nb+2*ng);
		AXPY(2*ns, 1.0, d, 2*nb+2*ng, res_d, 2*nb+2*ng, res_d, 2*nb+2*ng);
		}
	
	// res b, res g
	GEMV_NT(ne, nv, -1.0, -1.0, A, 0, 0, v, 0, pi, 0, 1.0, 1.0, b, 0, res_g, 0, res_b, 0, res_g, 0);

	// res_m res_mu
	mu = VECMULDOT(nct, lam, 0, t, 0, res_m, 0);
	AXPY(nct, -1.0, m, 0, res_m, 0, res_m, 0);
	res->res_mu = mu*nct_inv;


	return;

	}


void COMPUTE_LIN_RES_DENSE_QP(struct DENSE_QP *qp, struct DENSE_QP_SOL *qp_sol, struct DENSE_QP_SOL *qp_step, struct DENSE_QP_RES *res, struct DENSE_QP_RES_WORKSPACE *ws)
	{

	int nv = qp->dim->nv;
	int ne = qp->dim->ne;
	int nb = qp->dim->nb;
	int ng = qp->dim->ng;
	int ns = qp->dim->ns;

	int nct = 2*nb+2*ng+2*ns;

	REAL nct_inv = 1.0/nct;

	struct STRMAT *Hg = qp->Hv;
	struct STRMAT *A = qp->A;
	struct STRMAT *Ct = qp->Ct;
	struct STRVEC *gz = qp->gz;
	struct STRVEC *b = qp->b;
	struct STRVEC *d = qp->d;
	struct STRVEC *m = qp->m;
	int *idxb = qp->idxb;
	struct STRVEC *Z = qp->Z;
	int *idxs = qp->idxs;

	struct STRVEC *v = qp_step->v;
	struct STRVEC *pi = qp_step->pi;
	struct STRVEC *lam = qp_step->lam;
	struct STRVEC *t = qp_step->t;

	struct STRVEC *Lam = qp_sol->lam;
	struct STRVEC *T = qp_sol->t;

	struct STRVEC *res_g = res->res_g;
	struct STRVEC *res_b = res->res_b;
	struct STRVEC *res_d = res->res_d;
	struct STRVEC *res_m = res->res_m;

	struct STRVEC *tmp_nbg = ws->tmp_nbg;
	struct STRVEC *tmp_ns = ws->tmp_ns;

	REAL mu;

	// res g
	SYMV_L(nv, nv, 1.0, Hg, 0, 0, v, 0, 1.0, gz, 0, res_g, 0);

	if(nb+ng>0)
		{
		AXPY(nb+ng, -1.0, lam, 0, lam, nb+ng, tmp_nbg+0, 0);
//		AXPY(nb+ng,  1.0, d, 0, t, 0, res_d, 0);
//		AXPY(nb+ng,  1.0, d, nb+ng, t, nb+ng, res_d, nb+ng);
		AXPY(2*nb+2*ng,  1.0, d, 0, t, 0, res_d, 0);
		// box
		if(nb>0)
			{
			VECAD_SP(nb, 1.0, tmp_nbg+0, 0, idxb, res_g, 0);
			VECEX_SP(nb, 1.0, idxb, v, 0, tmp_nbg+1, 0);
			}
		// general
		if(ng>0)
			{
			GEMV_NT(nv, ng, 1.0, 1.0, Ct, 0, 0, tmp_nbg+0, nb, v, 0, 1.0, 0.0, res_g, 0, tmp_nbg+1, nb, res_g, 0, tmp_nbg+1, nb);
			}
		AXPY(nb+ng, -1.0, tmp_nbg+1, 0, res_d, 0, res_d, 0);
		AXPY(nb+ng,  1.0, tmp_nbg+1, 0, res_d, nb+ng, res_d, nb+ng);
		}
	if(ns>0)
		{
		// res_g
		GEMV_DIAG(2*ns, 1.0, Z, 0, v, nv, 1.0, gz, nv, res_g, nv);
		AXPY(2*ns, -1.0, lam, 2*nb+2*ng, res_g, nv, res_g, nv);
		VECEX_SP(ns, 1.0, idxs, lam, 0, tmp_ns, 0);
		AXPY(ns, -1.0, tmp_ns, 0, res_g, nv, res_g, nv);
		VECEX_SP(ns, 1.0, idxs, lam, nb+ng, tmp_ns, 0);
		AXPY(ns, -1.0, tmp_ns, 0, res_g, nv+ns, res_g, nv+ns);
		// res_d
		VECAD_SP(ns, -1.0, v, nv, idxs, res_d, 0);
		VECAD_SP(ns, -1.0, v, nv+ns, idxs, res_d, nb+ng);
		AXPY(2*ns, -1.0, v, nv, t, 2*nb+2*ng, res_d, 2*nb+2*ng);
		AXPY(2*ns, 1.0, d, 2*nb+2*ng, res_d, 2*nb+2*ng, res_d, 2*nb+2*ng);
		}
	
	// res b, res g
	GEMV_NT(ne, nv, -1.0, -1.0, A, 0, 0, v, 0, pi, 0, 1.0, 1.0, b, 0, res_g, 0, res_b, 0, res_g, 0);

	// res_m res_mu
//	VECCPSC(nct, 1.0, m, 0, res_m, 0);
	VECCP(nct, m, 0, res_m, 0);
	VECMULACC(nct, Lam, 0, t, 0, res_m, 0);
	VECMULACC(nct, lam, 0, T, 0, res_m, 0);
//	mu = VECMULDOT(nct, lam, 0, t, 0, res_m, 0);
//	AXPY(nct, -1.0, m, 0, res_m, 0, res_m, 0);
//	res->res_mu = mu*nct_inv;


	return;

	}


// range-space (Schur complement) method
void FACT_SOLVE_KKT_UNCONSTR_DENSE_QP(struct DENSE_QP *qp, struct DENSE_QP_SOL *qp_sol, struct DENSE_QP_IPM_WORKSPACE *ws)
	{

	int nv = qp->dim->nv;
	int ne = qp->dim->ne;
	int nb = qp->dim->nb;
	int ng = qp->dim->ng;

	struct STRMAT *Hg = qp->Hv;
	struct STRMAT *A = qp->A;
	struct STRVEC *gz = qp->gz;
	struct STRVEC *b = qp->b;

	struct STRVEC *v = qp_sol->v;
	struct STRVEC *pi = qp_sol->pi;

	struct STRMAT *Lv = ws->Lv;
	struct STRMAT *Le = ws->Le;
	struct STRMAT *Ctx = ws->Ctx;
	struct STRMAT *AL = ws->AL;
	struct STRVEC *lv = ws->lv;

	if(ne>0)
		{
		POTRF_L(nv, Hg, 0, 0, Lv, 0, 0);

//		GECP(ne, nv, A, 0, 0, AL, 0, 0);
		TRSM_RLTN(ne, nv, 1.0, Lv, 0, 0, A, 0, 0, AL, 0, 0);

		GESE(ne, ne, 0.0, Le, 0, 0);
		SYRK_POTRF_LN(ne, ne, nv, AL, 0, 0, AL, 0, 0, Le, 0, 0, Le, 0, 0);

		TRSV_LNN(nv, Lv, 0, 0, gz, 0, lv, 0);

		GEMV_N(ne, nv, 1.0, AL, 0, 0, lv, 0, 1.0, b, 0, pi, 0);

		TRSV_LNN(ne, Le, 0, 0, pi, 0, pi, 0);
		TRSV_LTN(ne, Le, 0, 0, pi, 0, pi, 0);

		GEMV_T(ne, nv, 1.0, A, 0, 0, pi, 0, -1.0, gz, 0, v, 0);

		TRSV_LNN(nv, Lv, 0, 0, v, 0, v, 0);
		TRSV_LTN(nv, Lv, 0, 0, v, 0, v, 0);
		}
	else
		{
#if 0
		POTRF_L(nv, Hg, 0, 0, Lv, 0, 0);

		VECCP(nv, gz, 0, v, 0);
		VECSC(nv, -1.0, v, 0);

		TRSV_LNN(nv, Lv, 0, 0, v, 0, v, 0);
		TRSV_LTN(nv, Lv, 0, 0, v, 0, v, 0);
#else
		ROWIN(nv, 1.0, gz, 0, Hg, nv, 0);
		POTRF_L_MN(nv+1, nv, Hg, 0, 0, Lv, 0, 0);

		ROWEX(nv, -1.0, Lv, nv, 0, v, 0);
		TRSV_LTN(nv, Lv, 0, 0, v, 0, v, 0);
#endif
		}

	return;

	}



static void COND_SLACKS_FACT_SOLVE(struct DENSE_QP *qp, struct DENSE_QP_IPM_WORKSPACE *ws)
	{

	int ii, idx;

	int nv = qp->dim->nv;
	int nb = qp->dim->nb;
	int ng = qp->dim->ng;
	int ns = qp->dim->ns;

	struct STRVEC *Z = qp->Z;
	int *idxs = qp->idxs;

	struct STRVEC *dv = ws->sol_step->v;
	struct STRVEC *res_g = ws->res->res_g;
	struct STRVEC *Gamma = ws->Gamma;
	struct STRVEC *gamma = ws->gamma;
	struct STRVEC *Zs_inv = ws->Zs_inv;
	struct STRVEC *tmp_nbg = ws->tmp_nbg;

	REAL *ptr_Gamma = Gamma->pa;
	REAL *ptr_gamma = gamma->pa;
	REAL *ptr_Z = Z->pa;
	REAL *ptr_Zs_inv = Zs_inv->pa;
	REAL *ptr_dv = dv->pa;
	REAL *ptr_res_g = res_g->pa;
	REAL *ptr_tmp0 = (tmp_nbg+0)->pa;
	REAL *ptr_tmp1 = (tmp_nbg+1)->pa;
	REAL *ptr_tmp2 = (tmp_nbg+2)->pa;
	REAL *ptr_tmp3 = (tmp_nbg+3)->pa;

	REAL tmp0, tmp1;

	VECCP(nb+ng, Gamma, 0, tmp_nbg+0, 0);
	VECCP(nb+ng, Gamma, nb+ng, tmp_nbg+1, 0);
	VECCP(nb+ng, gamma, 0, tmp_nbg+2, 0);
	VECCP(nb+ng, gamma, nb+ng, tmp_nbg+3, 0);

	for(ii=0; ii<ns; ii++)
		{
		idx = idxs[ii];
		ptr_Zs_inv[0+ii]  = ptr_Z[0+ii]  + ptr_Gamma[0+idx]     + ptr_Gamma[2*nb+2*ng+ii];
		ptr_Zs_inv[ns+ii] = ptr_Z[ns+ii] + ptr_Gamma[nb+ng+idx] + ptr_Gamma[2*nb+2*ng+ns+ii];
		ptr_dv[nv+ii]     = ptr_res_g[nv+ii]    + ptr_gamma[0+idx]     + ptr_gamma[2*nb+2*ng+ii];
		ptr_dv[nv+ns+ii]  = ptr_res_g[nv+ns+ii] + ptr_gamma[nb+ng+idx] + ptr_gamma[2*nb+2*ng+ns+ii];
		ptr_Zs_inv[0+ii]  = 1.0/ptr_Zs_inv[0+ii];
		ptr_Zs_inv[ns+ii] = 1.0/ptr_Zs_inv[ns+ii];
		tmp0 = ptr_dv[nv+ii]*ptr_Zs_inv[0+ii];
		tmp1 = ptr_dv[nv+ns+ii]*ptr_Zs_inv[ns+ii];
		ptr_tmp0[idx] = ptr_tmp0[idx] - ptr_tmp0[idx]*ptr_Zs_inv[0+ii]*ptr_tmp0[idx];
		ptr_tmp1[idx] = ptr_tmp1[idx] - ptr_tmp1[idx]*ptr_Zs_inv[ns+ii]*ptr_tmp1[idx];
		ptr_tmp2[idx] = ptr_tmp2[idx] - ptr_Gamma[0+idx]*tmp0;
		ptr_tmp3[idx] = ptr_tmp3[idx] - ptr_Gamma[nb+ng+idx]*tmp1;
		}
	
	AXPY(nb+ng,  1.0, tmp_nbg+1, 0, tmp_nbg+0, 0, tmp_nbg+0, 0);
	AXPY(nb+ng, -1.0, tmp_nbg+3, 0, tmp_nbg+2, 0, tmp_nbg+1, 0);

	return;

	}



static void COND_SLACKS_SOLVE(struct DENSE_QP *qp, struct DENSE_QP_IPM_WORKSPACE *ws)
	{

	int ii, idx;

	int nv = qp->dim->nv;
	int nb = qp->dim->nb;
	int ng = qp->dim->ng;
	int ns = qp->dim->ns;

	int *idxs = qp->idxs;

	struct STRVEC *dv = ws->sol_step->v;
	struct STRVEC *res_g = ws->res->res_g;
	struct STRVEC *Gamma = ws->Gamma;
	struct STRVEC *gamma = ws->gamma;
	struct STRVEC *Zs_inv = ws->Zs_inv;
	struct STRVEC *tmp_nbg = ws->tmp_nbg;

	REAL *ptr_Gamma = Gamma->pa;
	REAL *ptr_gamma = gamma->pa;
	REAL *ptr_Zs_inv = Zs_inv->pa;
	REAL *ptr_dv = dv->pa;
	REAL *ptr_res_g = res_g->pa;
	REAL *ptr_tmp2 = (tmp_nbg+2)->pa;
	REAL *ptr_tmp3 = (tmp_nbg+3)->pa;

	REAL tmp0, tmp1;

	VECCP(nb+ng, gamma, 0, tmp_nbg+2, 0);
	VECCP(nb+ng, gamma, nb+ng, tmp_nbg+3, 0);

	for(ii=0; ii<ns; ii++)
		{
		idx = idxs[ii];
		ptr_dv[nv+ii]     = ptr_res_g[nv+ii]    + ptr_gamma[0+idx]     + ptr_gamma[2*nb+2*ng+ii];
		ptr_dv[nv+ns+ii]  = ptr_res_g[nv+ns+ii] + ptr_gamma[nb+ng+idx] + ptr_gamma[2*nb+2*ng+ns+ii];
		tmp0 = ptr_dv[nv+ii]*ptr_Zs_inv[0+ii];
		tmp1 = ptr_dv[nv+ns+ii]*ptr_Zs_inv[ns+ii];
		ptr_tmp2[idx] = ptr_tmp2[idx] - ptr_Gamma[0+idx]*tmp0;
		ptr_tmp3[idx] = ptr_tmp3[idx] - ptr_Gamma[nb+ng+idx]*tmp1;
		}
	
	AXPY(nb+ng, -1.0, tmp_nbg+3, 0, tmp_nbg+2, 0, tmp_nbg+1, 0);

	return;

	}



static void EXPAND_SLACKS(struct DENSE_QP *qp, struct DENSE_QP_SOL *qp_sol, struct DENSE_QP_IPM_WORKSPACE *ws)
	{

	int ii, idx;

	int nv = qp->dim->nv;
	int nb = qp->dim->nb;
	int ng = qp->dim->ng;
	int ns = qp->dim->ns;

	int *idxs = qp->idxs;

	struct STRVEC *dv = qp_sol->v;
	struct STRVEC *dt = qp_sol->t;

	struct STRVEC *Gamma = ws->Gamma;
	struct STRVEC *Zs_inv = ws->Zs_inv;

	REAL *ptr_Gamma = Gamma->pa;
	REAL *ptr_dv = dv->pa;
	REAL *ptr_dt = dt->pa;
	REAL *ptr_Zs_inv = Zs_inv->pa;

	for(ii=0; ii<ns; ii++)
		{
		idx = idxs[ii];
		ptr_dv[nv+ii]    = - ptr_Zs_inv[0+ii]  * (ptr_dv[nv+ii]    + ptr_dt[idx]*ptr_Gamma[idx]);
		ptr_dv[nv+ns+ii] = - ptr_Zs_inv[ns+ii] * (ptr_dv[nv+ns+ii] + ptr_dt[nb+ng+idx]*ptr_Gamma[nb+ng+idx]);
		ptr_dt[2*nb+2*ng+ii]    = ptr_dv[nv+ii];
		ptr_dt[2*nb+2*ng+ns+ii] = ptr_dv[nv+ns+ii];
		ptr_dt[0+idx]     = ptr_dt[0+idx]     + ptr_dv[nv+ii];
		ptr_dt[nb+ng+idx] = ptr_dt[nb+ng+idx] + ptr_dv[nv+ns+ii];

		}

	return;

	}



// range-space (Schur complement) method
void FACT_SOLVE_KKT_STEP_DENSE_QP(struct DENSE_QP *qp, struct DENSE_QP_SOL *qp_sol, struct DENSE_QP_IPM_ARG *arg, struct DENSE_QP_IPM_WORKSPACE *ws)
	{

	int ii;

	int nv = qp->dim->nv;
	int ne = qp->dim->ne;
	int nb = qp->dim->nb;
	int ng = qp->dim->ng;
	int ns = qp->dim->ns;

	struct STRMAT *Hg = qp->Hv;
	struct STRMAT *A = qp->A;
	struct STRMAT *Ct = qp->Ct;
	int *idxb = qp->idxb;

	struct STRVEC *res_g = qp->gz;
	struct STRVEC *res_b = qp->b;

	struct STRVEC *dv = qp_sol->v;
	struct STRVEC *dpi = qp_sol->pi;
	struct STRVEC *dt = qp_sol->t;

	struct STRMAT *Lv = ws->Lv;
	struct STRMAT *Le = ws->Le;
	struct STRMAT *Ctx = ws->Ctx;
	struct STRMAT *AL = ws->AL;
	struct STRVEC *lv = ws->lv;
	struct STRVEC *sv = ws->sv;
	struct STRVEC *se = ws->se;
	struct STRVEC *Gamma = ws->Gamma;
	struct STRVEC *gamma = ws->gamma;
	struct STRVEC *tmp_nbg = ws->tmp_nbg;

	REAL tmp;

	struct CORE_QP_IPM_WORKSPACE *cws = ws->core_workspace;

	if(nb+ng>0)
		{
		COMPUTE_GAMMA_GAMMA_QP(qp->d->pa, qp->m->pa, cws);
		}

	if(ne>0)
		{

		if(arg->scale)
			{

//			TRCP_L(nv, Hg, 0, 0, Lv, 0, 0);
			GECP(nv, nv, Hg, 0, 0, Lv, 0, 0);

			VECCP(nv, res_g, 0, lv, 0);

			if(ns>0)
				{
				COND_SLACKS_FACT_SOLVE(qp, ws);
				}
			else if(nb+ng>0)
				{
				AXPY(nb+ng,  1.0, Gamma, nb+ng, Gamma, 0, tmp_nbg+0, 0);
				AXPY(nb+ng, -1.0, gamma, nb+ng, gamma, 0, tmp_nbg+1, 0);
				}
			if(nb>0)
				{
				DIAAD_SP(nb, 1.0, tmp_nbg+0, 0, idxb, Lv, 0, 0);
				VECAD_SP(nb, 1.0, tmp_nbg+1, 0, idxb, lv, 0);
				}
			if(ng>0)
				{
				GEMV_N(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+1, nb, 1.0, lv, 0, lv, 0);
				GEMM_R_DIAG(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+0, nb, 0.0, Ctx, 0, 0, Ctx, 0, 0);
				SYRK_LN(nv, ng, 1.0, Ctx, 0, 0, Ct, 0, 0, 1.0, Lv, 0, 0, Lv, 0, 0);
				}

			DIAEX(nv, 1.0, Lv, 0, 0, sv, 0);
			for(ii=0; ii<nv; ii++)
				{
				tmp = sqrt(sv->pa[ii]);
//				tmp = sqrt(tmp);
//				tmp = sqrt(sv->pa[ii]+tmp);
//				tmp = 1.0;
				sv->pa[ii] = tmp==0 ? 1.0 : 1.0/tmp;
				}

			GEMM_L_DIAG(nv, nv, 1.0, sv, 0, Lv, 0, 0, 0.0, Lv, 0, 0, Lv, 0, 0);
			GEMM_R_DIAG(nv, nv, 1.0, Lv, 0, 0, sv, 0, 0.0, Lv, 0, 0, Lv, 0, 0);
			DIARE(nv, arg->reg_prim, Lv, 0, 0);
			POTRF_L(nv, Lv, 0, 0, Lv, 0, 0);

			GEMV_DIAG(nv, 1.0, sv, 0, lv, 0, 0.0, lv, 0, lv, 0);
			VECCP(nv, lv, 0, dv, 0);

			GECP(ne, nv, A, 0, 0, AL, 0, 0);
			GEMM_R_DIAG(ne, nv, 1.0, AL, 0, 0, sv, 0, 0.0, AL, 0, 0, AL, 0, 0);
			TRSM_RLTN(ne, nv, 1.0, Lv, 0, 0, AL, 0, 0, AL, 0, 0);

			TRSV_LNN(nv, Lv, 0, 0, lv, 0, lv, 0);

			GESE(ne, ne, 0.0, Le, 0, 0);
			SYRK_LN(ne, nv, 1.0, AL, 0, 0, AL, 0, 0, 1.0, Le, 0, 0, Le, 0, 0);

			DIAEX(ne, 1.0, Le, 0, 0, se, 0);
			for(ii=0; ii<ne; ii++)
				{
				tmp = sqrt(se->pa[ii]);
//				tmp = sqrt(tmp);
//				tmp = sqrt(se->pa[ii]+tmp);
//				tmp = 1.0;
				se->pa[ii] = tmp==0 ? 1.0 : 1.0/tmp;
				}

			GEMM_L_DIAG(ne, ne, 1.0, se, 0, Le, 0, 0, 0.0, Le, 0, 0, Le, 0, 0);
			GEMM_R_DIAG(ne, ne, 1.0, Le, 0, 0, se, 0, 0.0, Le, 0, 0, Le, 0, 0);
			DIARE(ne, arg->reg_prim, Le, 0, 0);
			POTRF_L(ne, Le, 0, 0, Le, 0, 0);

			GEMV_N(ne, nv, 1.0, AL, 0, 0, lv, 0, 1.0, res_b, 0, dpi, 0);

			GEMV_DIAG(ne, 1.0, se, 0, dpi, 0, 0.0, dpi, 0, dpi, 0);
			TRSV_LNN(ne, Le, 0, 0, dpi, 0, dpi, 0);
			TRSV_LTN(ne, Le, 0, 0, dpi, 0, dpi, 0);
			GEMV_DIAG(ne, 1.0, se, 0, dpi, 0, 0.0, dpi, 0, dpi, 0);

			GEMV_T(ne, nv, 1.0, A, 0, 0, dpi, 0, 0.0, lv, 0, lv, 0);
			GEMV_DIAG(nv, 1.0, sv, 0, lv, 0, -1.0, dv, 0, dv, 0);

			TRSV_LNN(nv, Lv, 0, 0, dv, 0, dv, 0);
			TRSV_LTN(nv, Lv, 0, 0, dv, 0, dv, 0);
			GEMV_DIAG(nv, 1.0, sv, 0, dv, 0, 0.0, dv, 0, dv, 0);

			}
		else // no scale
			{

//			TRCP_L(nv, Hg, 0, 0, Lv, 0, 0);
			GECP(nv, nv, Hg, 0, 0, Lv, 0, 0);

			VECCP(nv, res_g, 0, lv, 0);

			if(ns>0)
				{
				COND_SLACKS_FACT_SOLVE(qp, ws);
				}
			else if(nb+ng>0)
				{
				AXPY(nb+ng,  1.0, Gamma, nb+ng, Gamma, 0, tmp_nbg+0, 0);
				AXPY(nb+ng, -1.0, gamma, nb+ng, gamma, 0, tmp_nbg+1, 0);
				}
			if(nb>0)
				{
				DIAAD_SP(nb, 1.0, tmp_nbg+0, 0, idxb, Lv, 0, 0);
				VECAD_SP(nb, 1.0, tmp_nbg+1, 0, idxb, lv, 0);
				}
			if(ng>0)
				{
				GEMV_N(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+1, nb, 1.0, lv, 0, lv, 0);
				GEMM_R_DIAG(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+0, nb, 0.0, Ctx, 0, 0, Ctx, 0, 0);
				SYRK_POTRF_LN(nv, nv, ng, Ctx, 0, 0, Ct, 0, 0, Lv, 0, 0, Lv, 0, 0);
				}
			else
				{
				POTRF_L(nv, Lv, 0, 0, Lv, 0, 0);
				}

			VECCP(nv, lv, 0, dv, 0);

			TRSM_RLTN(ne, nv, 1.0, Lv, 0, 0, A, 0, 0, AL, 0, 0);

			TRSV_LNN(nv, Lv, 0, 0, lv, 0, lv, 0);

			GEMV_N(ne, nv, 1.0, AL, 0, 0, lv, 0, 1.0, res_b, 0, dpi, 0);

			GESE(ne, ne, 0.0, Le, 0, 0);
			SYRK_POTRF_LN(ne, ne, nv, AL, 0, 0, AL, 0, 0, Le, 0, 0, Le, 0, 0);

			TRSV_LNN(ne, Le, 0, 0, dpi, 0, dpi, 0);
			TRSV_LTN(ne, Le, 0, 0, dpi, 0, dpi, 0);

			GEMV_T(ne, nv, 1.0, A, 0, 0, dpi, 0, -1.0, dv, 0, dv, 0);

			TRSV_LNN(nv, Lv, 0, 0, dv, 0, dv, 0);
			TRSV_LTN(nv, Lv, 0, 0, dv, 0, dv, 0);


			} // scale

		}
	else // ne==0
		{

		if(arg->scale)
			{

			TRCP_L(nv, Hg, 0, 0, Lv, 0, 0);
			VECCP(nv, res_g, 0, lv, 0);

			if(ns>0)
				{
				COND_SLACKS_FACT_SOLVE(qp, ws);
				}
			else if(nb+ng>0)
				{
				AXPY(nb+ng,  1.0, Gamma, nb+ng, Gamma, 0, tmp_nbg+0, 0);
				AXPY(nb+ng, -1.0, gamma, nb+ng, gamma, 0, tmp_nbg+1, 0);
				}
			if(nb>0)
				{
				DIAAD_SP(nb, 1.0, tmp_nbg+0, 0, idxb, Lv, 0, 0);
				VECAD_SP(nb, 1.0, tmp_nbg+1, 0, idxb, lv, 0);
				}
			if(ng>0)
				{
				GEMM_R_DIAG(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+0, nb, 0.0, Ctx, 0, 0, Ctx, 0, 0);
				GEMV_N(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+1, nb, 1.0, lv, 0, lv, 0);
				SYRK_LN(nv, ng, 1.0, Ctx, 0, 0, Ct, 0, 0, 1.0, Lv, 0, 0, Lv, 0, 0);
				}

			DIAEX(nv, 1.0, Lv, 0, 0, sv, 0);
			for(ii=0; ii<nv; ii++)
				{
				tmp = sqrt(sv->pa[ii]);
//				tmp = sqrt(tmp);
//				tmp = sqrt(sv->pa[ii]+tmp);
//				tmp = 1.0;
				sv->pa[ii] = tmp==0 ? 1.0 : 1.0/tmp;
				}

			GEMM_L_DIAG(nv, nv, 1.0, sv, 0, Lv, 0, 0, 0.0, Lv, 0, 0, Lv, 0, 0);
			GEMM_R_DIAG(nv, nv, 1.0, Lv, 0, 0, sv, 0, 0.0, Lv, 0, 0, Lv, 0, 0);
			DIARE(nv, arg->reg_prim, Lv, 0, 0);
			POTRF_L_MN(nv, nv, Lv, 0, 0, Lv, 0, 0);

			VECCP(nv, lv, 0, dv, 0);
			VECSC(nv, -1.0, dv, 0);

			GEMV_DIAG(nv, 1.0, sv, 0, dv, 0, 0.0, dv, 0, dv, 0);
			TRSV_LNN(nv, Lv, 0, 0, dv, 0, dv, 0);
			TRSV_LTN(nv, Lv, 0, 0, dv, 0, dv, 0);
			GEMV_DIAG(nv, 1.0, sv, 0, dv, 0, 0.0, dv, 0, dv, 0);
			}
		else // no scale
			{
	//		TRCP_L(nv, Hg, 0, 0, Lv, 0, 0);
			GECP(nv, nv, Hg, 0, 0, Lv, 0, 0);
			ROWIN(nv, 1.0, res_g, 0, Lv, nv, 0);

			if(ns>0)
				{
				COND_SLACKS_FACT_SOLVE(qp, ws);
				}
			else if(nb+ng>0)
				{
				AXPY(nb+ng,  1.0, Gamma, nb+ng, Gamma, 0, tmp_nbg+0, 0);
				AXPY(nb+ng, -1.0, gamma, nb+ng, gamma, 0, tmp_nbg+1, 0);
				}
			if(nb>0)
				{
				DIAAD_SP(nb, 1.0, tmp_nbg+0, 0, idxb, Lv, 0, 0);
				ROWAD_SP(nb, 1.0, tmp_nbg+1, 0, idxb, Lv, nv, 0);
				}
			if(ng>0)
				{
				GEMM_R_DIAG(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+0, nb, 0.0, Ctx, 0, 0, Ctx, 0, 0);
				ROWIN(ng, 1.0, tmp_nbg+1, nb, Ctx, nv, 0);
				SYRK_POTRF_LN(nv+1, nv, ng, Ctx, 0, 0, Ct, 0, 0, Lv, 0, 0, Lv, 0, 0); // TODO _mn_ routine in BLASFEO !!!
				}
			else
				{
				POTRF_L_MN(nv+1, nv, Lv, 0, 0, Lv, 0, 0);
				}

			ROWEX(nv, -1.0, Lv, nv, 0, dv, 0);
			TRSV_LTN(nv, Lv, 0, 0, dv, 0, dv, 0);

			} // scale

		} // ne>0

	if(nb+ng>0)
		{
		if(nb>0)
			VECEX_SP(nb, 1.0, idxb, dv, 0, dt, 0);

		if(ng>0)
			GEMV_T(nv, ng, 1.0, Ct, 0, 0, dv, 0, 0.0, dt, nb, dt, nb);

		VECCP(nb+ng, dt, 0, dt, nb+ng);
		VECSC(nb+ng, -1.0, dt, nb+ng);

		if(ns>0)
			EXPAND_SLACKS(qp, qp_sol, ws);

		COMPUTE_LAM_T_QP(qp->d->pa, qp->m->pa, qp_sol->lam->pa, qp_sol->t->pa, cws);
		}

	return;

	}



void FACT_SOLVE_LQ_KKT_STEP_DENSE_QP(struct DENSE_QP *qp, struct DENSE_QP_SOL *qp_sol, struct DENSE_QP_IPM_ARG *arg, struct DENSE_QP_IPM_WORKSPACE *ws)
	{

	int ii;

	int nv = qp->dim->nv;
	int ne = qp->dim->ne;
	int nb = qp->dim->nb;
	int ng = qp->dim->ng;
	int ns = qp->dim->ns;

	struct STRMAT *Hg = qp->Hv;
	struct STRMAT *A = qp->A;
	struct STRMAT *Ct = qp->Ct;
	int *idxb = qp->idxb;

	struct STRVEC *res_g = qp->gz;
	struct STRVEC *res_b = qp->b;

	struct STRVEC *dv = qp_sol->v;
	struct STRVEC *dpi = qp_sol->pi;
	struct STRVEC *dt = qp_sol->t;

	struct STRMAT *Lv = ws->Lv;
	struct STRMAT *Le = ws->Le;
	struct STRMAT *Ctx = ws->Ctx;
	struct STRMAT *AL = ws->AL;
	struct STRVEC *lv = ws->lv;
	struct STRVEC *sv = ws->sv;
	struct STRVEC *se = ws->se;
	struct STRVEC *Gamma = ws->Gamma;
	struct STRVEC *gamma = ws->gamma;
	struct STRVEC *tmp_nbg = ws->tmp_nbg;
	void *lq_work0 = ws->lq_work0;
	void *lq_work1 = ws->lq_work1;
	struct STRMAT *lq0 = ws->lq0;
	struct STRMAT *lq1 = ws->lq1;

	REAL tmp;

	struct CORE_QP_IPM_WORKSPACE *cws = ws->core_workspace;

	ws->scale = 0;

	if(nb+ng>0)
		{
		COMPUTE_GAMMA_GAMMA_QP(qp->d->pa, qp->m->pa, cws);
		}

	if(ne>0)
		{

		// XXX needed ???
		GESE(nv, nv+nv+ng, 0.0, lq1, 0, 0); // TODO not the first part for HP and RF

		if(ws->use_hess_fact==0)
			{
			POTRF_L(nv, Hg, 0, 0, Lv+1, 0, 0);
			ws->use_hess_fact=1;
			}

		VECCP(nv, res_g, 0, lv, 0);

		if(ns>0)
			{
			COND_SLACKS_FACT_SOLVE(qp, ws);
			}
		else if(nb+ng>0)
			{
			AXPY(nb+ng,  1.0, Gamma, nb+ng, Gamma, 0, tmp_nbg+0, 0);
			AXPY(nb+ng, -1.0, gamma, nb+ng, gamma, 0, tmp_nbg+1, 0);
			}
		if(nb>0)
			{
			for(ii=0; ii<nb; ii++)
				{
				tmp = BLASFEO_DVECEL(tmp_nbg+0, ii);
				tmp = tmp>=0.0 ? tmp : 0.0;
				tmp = sqrt( tmp );
				BLASFEO_DMATEL(lq1, idxb[ii], nv+idxb[ii]) = tmp>0.0 ? tmp : 0.0;
				}
			VECAD_SP(nb, 1.0, tmp_nbg+1, 0, idxb, lv, 0);
			}
		if(ng>0)
			{
			for(ii=0; ii<ng; ii++)
				{
				tmp = BLASFEO_DVECEL(tmp_nbg+0, nb+ii);
				tmp = tmp>=0.0 ? tmp : 0.0;
				tmp = sqrt( tmp );
				BLASFEO_DVECEL(tmp_nbg+0, nb+ii) = tmp;
				}
			GEMM_R_DIAG(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+0, nb, 0.0, lq1, 0, nv+nv, lq1, 0, nv+nv);
			GEMV_N(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+1, nb, 1.0, lv, 0, lv, 0);
			}

		DIARE(nv, arg->reg_prim, lq1, 0, nv);

//blasfeo_print_dmat(nv, nv, lq1, 0, 0);
//blasfeo_print_dmat(nv, nv, lq1, 0, nv);
//blasfeo_print_dmat(nv, ng, lq1, 0, nv+ng);

#if defined(LA_HIGH_PERFORMANCE) | defined(LA_REFERENCE)
//		TRCP_L(nv, Lv+1, 0, 0, lq1, 0, 0);
//		GELQF_PD(nv, nv+nv+ng, lq1, 0, 0, lq1, 0, 0, lq_work1);
//		GELQF_PD_LA(nv, nv+ng, lq1, 0, 0, lq1, 0, nv, lq_work1);
//		GELQF_PD_LLA(nv, ng, lq1, 0, 0, lq1, 0, nv, lq1, 0, 2*nv, lq_work1);
//		TRCP_L(nv, lq1, 0, 0, Lv, 0, 0);
		TRCP_L(nv, Lv+1, 0, 0, Lv, 0, 0);
		GELQF_PD_LLA(nv, ng, Lv, 0, 0, lq1, 0, nv, lq1, 0, 2*nv, lq_work1); // TODO reduce lq1 size !!!
#else // LA_BLAS_WRAPPER
		TRCP_L(nv, Lv+1, 0, 0, lq1, 0, 0);
		GELQF(nv, nv+nv+ng, lq1, 0, 0, lq1, 0, 0, lq_work1);
		TRCP_L(nv, lq1, 0, 0, Lv, 0, 0);
		for(ii=0; ii<nv; ii++)
			if(BLASFEO_DMATEL(Lv, ii, ii) < 0)
				COLSC(nv-ii, -1.0, Lv, ii, ii);
#endif

//blasfeo_print_dmat(nv, nv, Lv, 0, 0);

		VECCP(nv, lv, 0, dv, 0);

		TRSM_RLTN(ne, nv, 1.0, Lv, 0, 0, A, 0, 0, AL, 0, 0);

		TRSV_LNN(nv, Lv, 0, 0, lv, 0, lv, 0);

		GEMV_N(ne, nv, 1.0, AL, 0, 0, lv, 0, 1.0, res_b, 0, dpi, 0);

		GECP(ne, nv, AL, 0, 0, lq0, 0, ne);

#if defined(LA_HIGH_PERFORMANCE)
//		GESE(ne, ne, 0.0, lq0, 0, 0);
//		DIARE(ne, arg->reg_dual, lq0, 0, 0);
//		GELQF_PD(ne, ne+nv, lq0, 0, 0, lq0, 0, 0, lq_work0);
//		GELQF_PD_LA(ne, nv, lq0, 0, 0, lq0, 0, ne, lq_work0);
//		TRCP_L(ne, lq0, 0, 0, Le, 0, 0);
		GESE(ne, ne, 0.0, Le, 0, 0);
		DIARE(ne, arg->reg_dual, Le, 0, 0);
		GELQF_PD_LA(ne, nv, Le, 0, 0, lq0, 0, ne, lq_work0); // TODO reduce lq0 size !!!
#elif defined(LA_REFERENCE)
//		GESE(ne, ne, 0.0, lq0, 0, 0);
//		DIARE(ne, arg->reg_dual, lq0, 0, 0);
//		GELQF_PD(ne, ne+nv, lq0, 0, 0, lq0, 0, 0, lq_work0);
//		GELQF_PD_LA(ne, nv, lq0, 0, 0, lq0, 0, ne, lq_work0);
//		TRCP_L(ne, lq0, 0, 0, Le, 0, 0);
		GESE(ne, ne, 0.0, Le, 0, 0);
		DIARE(ne, arg->reg_dual, Le, 0, 0);
		GELQF_PD_LA(ne, nv, Le, 0, 0, lq0, 0, ne, lq_work0); // TODO reduce lq0 size !!!
#else // LA_BLAS_WRAPPER
		GESE(ne, ne, 0.0, lq0, 0, 0);
		DIARE(ne, arg->reg_dual, lq0, 0, 0);
		GELQF(ne, ne+nv, lq0, 0, 0, lq0, 0, 0, lq_work0);
		TRCP_L(ne, lq0, 0, 0, Le, 0, 0);
		for(ii=0; ii<ne; ii++)
			if(BLASFEO_DMATEL(Le, ii, ii) < 0)
				COLSC(ne-ii, -1.0, Le, ii, ii);
#endif

//		blasfeo_print_dmat(ne, ne, Le, 0, 0);

		TRSV_LNN(ne, Le, 0, 0, dpi, 0, dpi, 0);
		TRSV_LTN(ne, Le, 0, 0, dpi, 0, dpi, 0);

		GEMV_T(ne, nv, 1.0, A, 0, 0, dpi, 0, -1.0, dv, 0, dv, 0);

		TRSV_LNN(nv, Lv, 0, 0, dv, 0, dv, 0);
		TRSV_LTN(nv, Lv, 0, 0, dv, 0, dv, 0);

		}
	else // ne==0
		{

		// XXX needed ???
		GESE(nv, nv+nv+ng, 0.0, lq1, 0, 0); // TODO not the first part for HP and RF

		if(ws->use_hess_fact==0)
			{
			POTRF_L(nv, Hg, 0, 0, Lv+1, 0, 0);
			ws->use_hess_fact=1;
			}

		VECCP(nv, res_g, 0, lv, 0);

		if(ns>0)
			{
			COND_SLACKS_FACT_SOLVE(qp, ws);
			}
		else if(nb+ng>0)
			{
			AXPY(nb+ng,  1.0, Gamma, nb+ng, Gamma, 0, tmp_nbg+0, 0);
			AXPY(nb+ng, -1.0, gamma, nb+ng, gamma, 0, tmp_nbg+1, 0);
			}
		if(nb>0)
			{
			for(ii=0; ii<nb; ii++)
				{
				tmp = BLASFEO_DVECEL(tmp_nbg+0, ii);
				tmp = tmp>=0.0 ? tmp : 0.0;
				tmp = sqrt( tmp );
				BLASFEO_DMATEL(lq1, idxb[ii], nv+idxb[ii]) = tmp>0.0 ? tmp : 0.0;
				}
			VECAD_SP(nb, 1.0, tmp_nbg+1, 0, idxb, lv, 0);
			}
		if(ng>0)
			{
			for(ii=0; ii<ng; ii++)
				{
				tmp = BLASFEO_DVECEL(tmp_nbg+0, nb+ii);
				tmp = tmp>=0.0 ? tmp : 0.0;
				tmp = sqrt( tmp );
				BLASFEO_DVECEL(tmp_nbg+0, nb+ii) = tmp;
				}
			GEMM_R_DIAG(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+0, nb, 0.0, lq1, 0, nv+nv, lq1, 0, nv+nv);
			GEMV_N(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+1, nb, 1.0, lv, 0, lv, 0);
			}

		DIARE(nv, arg->reg_prim, lq1, 0, nv);

#if defined(LA_HIGH_PERFORMANCE)
//		TRCP_L(nv, Lv+1, 0, 0, lq1, 0, 0);
//		GELQF_PD(nv, nv+nv+ng, lq1, 0, 0, lq1, 0, 0, lq_work1);
//		GELQF_PD_LA(nv, nv+ng, lq1, 0, 0, lq1, 0, nv, lq_work1);
//		GELQF_PD_LLA(nv, ng, lq1, 0, 0, lq1, 0, nv, lq1, 0, 2*nv, lq_work1);
//		TRCP_L(nv, lq1, 0, 0, Lv, 0, 0);
		TRCP_L(nv, Lv+1, 0, 0, Lv, 0, 0);
		GELQF_PD_LLA(nv, ng, Lv, 0, 0, lq1, 0, nv, lq1, 0, 2*nv, lq_work1); // TODO reduce lq1 size !!!
#elif defined(LA_REFERENCE)
//		TRCP_L(nv, Lv+1, 0, 0, lq1, 0, 0);
//		GELQF_PD(nv, nv+nv+ng, lq1, 0, 0, lq1, 0, 0, lq_work1);
//		GELQF_PD_LA(nv, nv+ng, lq1, 0, 0, lq1, 0, nv, lq_work1);
//		GELQF_PD_LLA(nv, ng, lq1, 0, 0, lq1, 0, nv, lq1, 0, 2*nv, lq_work1);
//		TRCP_L(nv, lq1, 0, 0, Lv, 0, 0);
		TRCP_L(nv, Lv+1, 0, 0, Lv, 0, 0);
		GELQF_PD_LLA(nv, ng, Lv, 0, 0, lq1, 0, nv, lq1, 0, 2*nv, lq_work1); // TODO reduce lq1 size !!!
#else // LA_BLAS_WRAPPER
		TRCP_L(nv, Lv+1, 0, 0, lq1, 0, 0);
		GELQF(nv, nv+nv+ng, lq1, 0, 0, lq1, 0, 0, lq_work1);
		TRCP_L(nv, lq1, 0, 0, Lv, 0, 0);
		for(ii=0; ii<nv; ii++)
			if(BLASFEO_DMATEL(Lv, ii, ii) < 0)
				COLSC(nv-ii, -1.0, Lv, ii, ii);
#endif

#if 0
if(nv<30)
{
//	blasfeo_print_dmat(nv, nv+nb+ng, lq1, 0, 0);
blasfeo_print_dmat(nv, nv, Lv, 0, 0);
exit(1);
}
#endif

		VECCP(nv, lv, 0, dv, 0);
		VECSC(nv, -1.0, dv, 0);

		TRSV_LNN(nv, Lv, 0, 0, dv, 0, dv, 0);
		TRSV_LTN(nv, Lv, 0, 0, dv, 0, dv, 0);

		} // ne>0
	
	if(nb+ng>0)
		{
		if(nb>0)
			VECEX_SP(nb, 1.0, idxb, dv, 0, dt, 0);

		VECSE(ng, 0.0, dt, nb);
		if(ng>0)
			GEMV_T(nv, ng, 1.0, Ct, 0, 0, dv, 0, 0.0, dt, nb, dt, nb);

		VECCP(nb+ng, dt, 0, dt, nb+ng);
		VECSC(nb+ng, -1.0, dt, nb+ng);

		if(ns>0)
			EXPAND_SLACKS(qp, qp_sol, ws);

		COMPUTE_LAM_T_QP(qp->d->pa, qp->m->pa, qp_sol->lam->pa, qp_sol->t->pa, cws);
		}

	return;

	}



#if 0
void FACT_SOLVE_LU_KKT_STEP_DENSE_QP(struct DENSE_QP *qp, struct DENSE_QP_SOL *qp_sol, struct DENSE_QP_IPM_ARG *arg, struct DENSE_QP_IPM_WORKSPACE *ws)
	{

	int ii;

	int nv = qp->dim->nv;
	int ne = qp->dim->ne;
	int nb = qp->dim->nb;
	int ng = qp->dim->ng;
	int ns = qp->dim->ns;

	struct STRMAT *Hg = qp->Hv;
	struct STRMAT *A = qp->A;
	struct STRMAT *Ct = qp->Ct;
	int *idxb = qp->idxb;

	struct STRVEC *res_g = qp->gz;
	struct STRVEC *res_b = qp->b;

	struct STRVEC *dv = qp_sol->v;
	struct STRVEC *dpi = qp_sol->pi;
	struct STRVEC *dt = qp_sol->t;

	struct STRMAT *Lv = ws->Lv;
	struct STRMAT *Le = ws->Le;
	struct STRMAT *Ctx = ws->Ctx;
	struct STRMAT *AL = ws->AL;
	struct STRVEC *lv = ws->lv;
	struct STRVEC *sv = ws->sv;
	struct STRVEC *se = ws->se;
	struct STRVEC *Gamma = ws->Gamma;
	struct STRVEC *gamma = ws->gamma;
	struct STRVEC *tmp_nbg = ws->tmp_nbg;
	int *ipiv_v = ws->ipiv_v;
	int *ipiv_e = ws->ipiv_e;

	REAL tmp;

	struct CORE_QP_IPM_WORKSPACE *cws = ws->core_workspace;

	if(nb+ng>0)
		{
		COMPUTE_GAMMA_GAMMA_QP(qp->d->pa, qp->m->pa, cws);
		}

	if(ne>0)
		{

		if(arg->scale)
			{

//			TRCP_L(nv, Hg, 0, 0, Lv, 0, 0);
			GECP(nv, nv, Hg, 0, 0, Lv, 0, 0);

			VECCP(nv, res_g, 0, lv, 0);

			if(ns>0)
				{
				COND_SLACKS_FACT_SOLVE(qp, ws);
				}
			else if(nb+ng>0)
				{
				AXPY(nb+ng,  1.0, Gamma, nb+ng, Gamma, 0, tmp_nbg+0, 0);
				AXPY(nb+ng, -1.0, gamma, nb+ng, gamma, 0, tmp_nbg+1, 0);
				}
			if(nb>0)
				{
				DIAAD_SP(nb, 1.0, tmp_nbg+0, 0, idxb, Lv, 0, 0);
				VECAD_SP(nb, 1.0, tmp_nbg+1, 0, idxb, lv, 0);
				}
			if(ng>0)
				{
				GEMV_N(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+1, nb, 1.0, lv, 0, lv, 0);
				GEMM_R_DIAG(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+0, nb, 0.0, Ctx, 0, 0, Ctx, 0, 0);
				SYRK_LN(nv, ng, 1.0, Ctx, 0, 0, Ct, 0, 0, 1.0, Lv, 0, 0, Lv, 0, 0);
				}

			DIAEX(nv, 1.0, Lv, 0, 0, sv, 0);
			for(ii=0; ii<nv; ii++)
				{
				tmp = sqrt(sv->pa[ii]);
//				tmp = sqrt(tmp);
//				tmp = sqrt(sv->pa[ii]+tmp);
//				tmp = 1.0;
				sv->pa[ii] = tmp==0 ? 1.0 : 1.0/tmp;
				}

			GEMM_L_DIAG(nv, nv, 1.0, sv, 0, Lv, 0, 0, 0.0, Lv, 0, 0, Lv, 0, 0);
			GEMM_R_DIAG(nv, nv, 1.0, Lv, 0, 0, sv, 0, 0.0, Lv, 0, 0, Lv, 0, 0);
			DIARE(nv, arg->reg_prim, Lv, 0, 0);
			POTRF_L(nv, Lv, 0, 0, Lv, 0, 0);

			GEMV_DIAG(nv, 1.0, sv, 0, lv, 0, 0.0, lv, 0, lv, 0);
			VECCP(nv, lv, 0, dv, 0);

			GECP(ne, nv, A, 0, 0, AL, 0, 0);
			GEMM_R_DIAG(ne, nv, 1.0, AL, 0, 0, sv, 0, 0.0, AL, 0, 0, AL, 0, 0);
			TRSM_RLTN(ne, nv, 1.0, Lv, 0, 0, AL, 0, 0, AL, 0, 0);

			TRSV_LNN(nv, Lv, 0, 0, lv, 0, lv, 0);

			GESE(ne, ne, 0.0, Le, 0, 0);
			SYRK_LN(ne, nv, 1.0, AL, 0, 0, AL, 0, 0, 1.0, Le, 0, 0, Le, 0, 0);

			DIAEX(ne, 1.0, Le, 0, 0, se, 0);
			for(ii=0; ii<ne; ii++)
				{
				tmp = sqrt(se->pa[ii]);
//				tmp = sqrt(tmp);
//				tmp = sqrt(se->pa[ii]+tmp);
//				tmp = 1.0;
				se->pa[ii] = tmp==0 ? 1.0 : 1.0/tmp;
				}

			GEMM_L_DIAG(ne, ne, 1.0, se, 0, Le, 0, 0, 0.0, Le, 0, 0, Le, 0, 0);
			GEMM_R_DIAG(ne, ne, 1.0, Le, 0, 0, se, 0, 0.0, Le, 0, 0, Le, 0, 0);
			DIARE(ne, arg->reg_prim, Le, 0, 0);
			POTRF_L(ne, Le, 0, 0, Le, 0, 0);

			GEMV_N(ne, nv, 1.0, AL, 0, 0, lv, 0, 1.0, res_b, 0, dpi, 0);

			GEMV_DIAG(ne, 1.0, se, 0, dpi, 0, 0.0, dpi, 0, dpi, 0);
			TRSV_LNN(ne, Le, 0, 0, dpi, 0, dpi, 0);
			TRSV_LTN(ne, Le, 0, 0, dpi, 0, dpi, 0);
			GEMV_DIAG(ne, 1.0, se, 0, dpi, 0, 0.0, dpi, 0, dpi, 0);

			GEMV_T(ne, nv, 1.0, A, 0, 0, dpi, 0, 0.0, lv, 0, lv, 0);
			GEMV_DIAG(nv, 1.0, sv, 0, lv, 0, -1.0, dv, 0, dv, 0);

			TRSV_LNN(nv, Lv, 0, 0, dv, 0, dv, 0);
			TRSV_LTN(nv, Lv, 0, 0, dv, 0, dv, 0);
			GEMV_DIAG(nv, 1.0, sv, 0, dv, 0, 0.0, dv, 0, dv, 0);

			}
		else // no scale
			{

//			TRCP_L(nv, Hg, 0, 0, Lv, 0, 0);
			GECP(nv, nv, Hg, 0, 0, Lv, 0, 0);

			VECCP(nv, res_g, 0, lv, 0);

			if(ns>0)
				{
				COND_SLACKS_FACT_SOLVE(qp, ws);
				}
			else if(nb+ng>0)
				{
				AXPY(nb+ng,  1.0, Gamma, nb+ng, Gamma, 0, tmp_nbg+0, 0);
				AXPY(nb+ng, -1.0, gamma, nb+ng, gamma, 0, tmp_nbg+1, 0);
				}
			if(nb>0)
				{
				DIAAD_SP(nb, 1.0, tmp_nbg+0, 0, idxb, Lv, 0, 0);
				VECAD_SP(nb, 1.0, tmp_nbg+1, 0, idxb, lv, 0);
				}
			if(ng>0)
				{
				GEMV_N(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+1, nb, 1.0, lv, 0, lv, 0);
				GEMM_R_DIAG(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+0, nb, 0.0, Ctx, 0, 0, Ctx, 0, 0);
				SYRK_LN(nv, ng, 1.0, Ctx, 0, 0, Ct, 0, 0, 1.0, Lv, 0, 0, Lv, 0, 0);
				}

			TRTR_L(nv, Lv, 0, 0, Lv, 0, 0);
			GETRF(nv, nv, Lv, 0, 0, Lv, 0, 0, ipiv_v);

			VECCP(nv, lv, 0, dv, 0);
			VECPE(nv, ipiv_v, lv, 0);
			TRSV_LNU(nv, Lv, 0, 0, lv, 0, lv, 0);

			GECP(ne, nv, A, 0, 0, AL+1, 0, 0);
			COLPE(nv, ipiv_v, AL+1);
			TRSM_RLTU(ne, nv, 1.0, Lv, 0, 0, AL+1, 0, 0, AL+1, 0, 0);
//			TRSM_RUNN(ne, nv, 1.0, Lv, 0, 0, A, 0, 0, AL+0, 0, 0);
			TRTR_U(nv, Lv, 0, 0, Lv+1, 0, 0);
			TRSM_RLTN(ne, nv, 1.0, Lv+1, 0, 0, A, 0, 0, AL+0, 0, 0);

			GEMV_N(ne, nv, 1.0, AL+0, 0, 0, lv, 0, 1.0, res_b, 0, dpi, 0);

			SYRK_LN(ne, nv, 1.0, AL+0, 0, 0, AL+1, 0, 0, 0.0, Le, 0, 0, Le, 0, 0);

#if 0
			POTRF_L(ne, Le, 0, 0, Le, 0, 0);

			TRSV_LNN(ne, Le, 0, 0, dpi, 0, dpi, 0);
			TRSV_LTN(ne, Le, 0, 0, dpi, 0, dpi, 0);
#else
			TRTR_L(ne, Le, 0, 0, Le, 0, 0);
			GETRF(ne, ne, Le, 0, 0, Le, 0, 0, ipiv_e);

			VECPE(ne, ipiv_e, dpi, 0);
			TRSV_LNU(ne, Le, 0, 0, dpi, 0, dpi, 0);
			TRSV_UNN(ne, Le, 0, 0, dpi, 0, dpi, 0);
#endif

			GEMV_T(ne, nv, 1.0, A, 0, 0, dpi, 0, -1.0, dv, 0, dv, 0);

			VECPE(nv, ipiv_v, dv, 0);
			TRSV_LNU(nv, Lv, 0, 0, dv, 0, dv, 0);
			TRSV_UNN(nv, Lv, 0, 0, dv, 0, dv, 0);

			} // scale

		}
	else // ne==0
		{

		if(arg->scale)
			{

			TRCP_L(nv, Hg, 0, 0, Lv, 0, 0);
			VECCP(nv, res_g, 0, lv, 0);

			if(ns>0)
				{
				COND_SLACKS_FACT_SOLVE(qp, ws);
				}
			else if(nb+ng>0)
				{
				AXPY(nb+ng,  1.0, Gamma, nb+ng, Gamma, 0, tmp_nbg+0, 0);
				AXPY(nb+ng, -1.0, gamma, nb+ng, gamma, 0, tmp_nbg+1, 0);
				}
			if(nb>0)
				{
				DIAAD_SP(nb, 1.0, tmp_nbg+0, 0, idxb, Lv, 0, 0);
				VECAD_SP(nb, 1.0, tmp_nbg+1, 0, idxb, lv, 0);
				}
			if(ng>0)
				{
				GEMM_R_DIAG(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+0, nb, 0.0, Ctx, 0, 0, Ctx, 0, 0);
				GEMV_N(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+1, nb, 1.0, lv, 0, lv, 0);
				SYRK_LN(nv, ng, 1.0, Ctx, 0, 0, Ct, 0, 0, 1.0, Lv, 0, 0, Lv, 0, 0);
				}

			DIAEX(nv, 1.0, Lv, 0, 0, sv, 0);
			for(ii=0; ii<nv; ii++)
				{
				tmp = sqrt(sv->pa[ii]);
//				tmp = sqrt(tmp);
//				tmp = sqrt(sv->pa[ii]+tmp);
//				tmp = 1.0;
				sv->pa[ii] = tmp==0 ? 1.0 : 1.0/tmp;
				}

			GEMM_L_DIAG(nv, nv, 1.0, sv, 0, Lv, 0, 0, 0.0, Lv, 0, 0, Lv, 0, 0);
			GEMM_R_DIAG(nv, nv, 1.0, Lv, 0, 0, sv, 0, 0.0, Lv, 0, 0, Lv, 0, 0);
			DIARE(nv, arg->reg_prim, Lv, 0, 0);
			POTRF_L_MN(nv, nv, Lv, 0, 0, Lv, 0, 0);

			VECCP(nv, lv, 0, dv, 0);
			VECSC(nv, -1.0, dv, 0);

			GEMV_DIAG(nv, 1.0, sv, 0, dv, 0, 0.0, dv, 0, dv, 0);
			TRSV_LNN(nv, Lv, 0, 0, dv, 0, dv, 0);
			TRSV_LTN(nv, Lv, 0, 0, dv, 0, dv, 0);
			GEMV_DIAG(nv, 1.0, sv, 0, dv, 0, 0.0, dv, 0, dv, 0);
			}
		else // no scale
			{
	//		TRCP_L(nv, Hg, 0, 0, Lv, 0, 0);
			GECP(nv, nv, Hg, 0, 0, Lv, 0, 0);
			ROWIN(nv, 1.0, res_g, 0, Lv, nv, 0);

			if(ns>0)
				{
				COND_SLACKS_FACT_SOLVE(qp, ws);
				}
			else if(nb+ng>0)
				{
				AXPY(nb+ng,  1.0, Gamma, nb+ng, Gamma, 0, tmp_nbg+0, 0);
				AXPY(nb+ng, -1.0, gamma, nb+ng, gamma, 0, tmp_nbg+1, 0);
				}
			if(nb>0)
				{
				DIAAD_SP(nb, 1.0, tmp_nbg+0, 0, idxb, Lv, 0, 0);
				ROWAD_SP(nb, 1.0, tmp_nbg+1, 0, idxb, Lv, nv, 0);
				}
			if(ng>0)
				{
				GEMM_R_DIAG(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+0, nb, 0.0, Ctx, 0, 0, Ctx, 0, 0);
				ROWIN(ng, 1.0, tmp_nbg+1, nb, Ctx, nv, 0);
				SYRK_LN_MN(nv+1, nv, ng, 1.0, Ctx, 0, 0, Ct, 0, 0, 1.0, Lv, 0, 0, Lv, 0, 0);
				}

			ROWEX(nv, -1.0, Lv, nv, 0, dv, 0);

			TRTR_L(nv, Lv, 0, 0, Lv, 0, 0);
			GETRF(nv, nv, Lv, 0, 0, Lv, 0, 0, ipiv_v);

			VECPE(nv, ipiv_v, dv, 0);
			TRSV_LNU(nv, Lv, 0, 0, dv, 0, dv, 0);
			TRSV_UNN(nv, Lv, 0, 0, dv, 0, dv, 0);

			} // scale

		} // ne>0

	if(nb+ng>0)
		{
		if(nb>0)
			VECEX_SP(nb, 1.0, idxb, dv, 0, dt, 0);

		if(ng>0)
			GEMV_T(nv, ng, 1.0, Ct, 0, 0, dv, 0, 0.0, dt, nb, dt, nb);

		VECCP(nb+ng, dt, 0, dt, nb+ng);
		VECSC(nb+ng, -1.0, dt, nb+ng);

		if(ns>0)
			EXPAND_SLACKS(qp, qp_sol, ws);

		COMPUTE_LAM_T_QP(qp->d->pa, qp->m->pa, qp_sol->lam->pa, qp_sol->t->pa, cws);
		}

	return;

	}
#endif



// range-space (Schur complement) method
void SOLVE_KKT_STEP_DENSE_QP(struct DENSE_QP *qp, struct DENSE_QP_SOL *qp_sol, struct DENSE_QP_IPM_ARG *arg, struct DENSE_QP_IPM_WORKSPACE *ws)
	{

	int nv = qp->dim->nv;
	int ne = qp->dim->ne;
	int nb = qp->dim->nb;
	int ng = qp->dim->ng;
	int ns = qp->dim->ns;

	struct STRMAT *A = qp->A;
	struct STRMAT *Ct = qp->Ct;
	int *idxb = qp->idxb;

//	struct STRVEC *res_g = ws->res->res_g;
//	struct STRVEC *res_b = ws->res->res_b;
	struct STRVEC *res_g = qp->gz;
	struct STRVEC *res_b = qp->b;

	struct STRVEC *dv = qp_sol->v;
	struct STRVEC *dpi = qp_sol->pi;
	struct STRVEC *dt = qp_sol->t;

	struct STRMAT *Lv = ws->Lv;
	struct STRMAT *Le = ws->Le;
	struct STRMAT *Ctx = ws->Ctx;
	struct STRMAT *AL = ws->AL;
	struct STRVEC *lv = ws->lv;
	struct STRVEC *sv = ws->sv;
	struct STRVEC *se = ws->se;
	struct STRVEC *gamma = ws->gamma;
	struct STRVEC *tmp_nbg = ws->tmp_nbg;

	struct CORE_QP_IPM_WORKSPACE *cws = ws->core_workspace;

	if(nb>0 | ng>0)
		{
		COMPUTE_GAMMA_QP(qp->d->pa, qp->m->pa, cws);
		}

	if(ne>0)
		{

		if(ws->scale)
			{

			VECCP(nv, res_g, 0, lv, 0);

			if(ns>0)
				{
				COND_SLACKS_SOLVE(qp, ws);
				}
			else if(nb+ng>0)
				{
				AXPY(nb+ng, -1.0, gamma, nb+ng, gamma, 0, tmp_nbg+1, 0);
				}
			if(nb>0)
				{
				VECAD_SP(nb, 1.0, tmp_nbg+1, 0, idxb, lv, 0);
				}
			if(ng>0)
				{
				GEMV_N(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+1, nb, 1.0, lv, 0, lv, 0);
				}

			GEMV_DIAG(nv, 1.0, sv, 0, lv, 0, 0.0, lv, 0, lv, 0);
			VECCP(nv, lv, 0, dv, 0);

			TRSV_LNN(nv, Lv, 0, 0, lv, 0, lv, 0);

			GEMV_N(ne, nv, 1.0, AL, 0, 0, lv, 0, 1.0, res_b, 0, dpi, 0);

			GEMV_DIAG(ne, 1.0, se, 0, dpi, 0, 0.0, dpi, 0, dpi, 0);
			TRSV_LNN(ne, Le, 0, 0, dpi, 0, dpi, 0);
			TRSV_LTN(ne, Le, 0, 0, dpi, 0, dpi, 0);
			GEMV_DIAG(ne, 1.0, se, 0, dpi, 0, 0.0, dpi, 0, dpi, 0);

			GEMV_T(ne, nv, 1.0, A, 0, 0, dpi, 0, 0.0, lv, 0, lv, 0);
			GEMV_DIAG(nv, 1.0, sv, 0, lv, 0, -1.0, dv, 0, dv, 0);

			TRSV_LNN(nv, Lv, 0, 0, dv, 0, dv, 0);
			TRSV_LTN(nv, Lv, 0, 0, dv, 0, dv, 0);
			GEMV_DIAG(nv, 1.0, sv, 0, dv, 0, 0.0, dv, 0, dv, 0);

			}
		else // no scale
			{

			VECCP(nv, res_g, 0, lv, 0);

			if(ns>0)
				{
				COND_SLACKS_SOLVE(qp, ws);
				}
			else if(nb+ng>0)
				{
				AXPY(nb+ng, -1.0, gamma, nb+ng, gamma, 0, tmp_nbg+1, 0);
				}
			if(nb>0)
				{
				VECAD_SP(nb, 1.0, tmp_nbg+1, 0, idxb, lv, 0);
				}
			if(ng>0)
				{
				GEMV_N(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+1, nb, 1.0, lv, 0, lv, 0);
				}

			VECCP(nv, lv, 0, dv, 0);

			TRSV_LNN(nv, Lv, 0, 0, lv, 0, lv, 0);

			GEMV_N(ne, nv, 1.0, AL, 0, 0, lv, 0, 1.0, res_b, 0, dpi, 0);

			TRSV_LNN(ne, Le, 0, 0, dpi, 0, dpi, 0);
			TRSV_LTN(ne, Le, 0, 0, dpi, 0, dpi, 0);

			GEMV_T(ne, nv, 1.0, A, 0, 0, dpi, 0, -1.0, dv, 0, dv, 0);

			TRSV_LNN(nv, Lv, 0, 0, dv, 0, dv, 0);
			TRSV_LTN(nv, Lv, 0, 0, dv, 0, dv, 0);

			} // scale

		}
	else // ne==0
		{

		if(ws->scale)
			{

			VECCP(nv, res_g, 0, lv, 0);

			if(ns>0)
				{
				COND_SLACKS_SOLVE(qp, ws);
				}
			else if(nb+ng>0)
				{
				AXPY(nb+ng, -1.0, gamma, nb+ng, gamma, 0, tmp_nbg+1, 0);
				}
			if(nb>0)
				{
				VECAD_SP(nb, 1.0, tmp_nbg+1, 0, idxb, lv, 0);
				}
			if(ng>0)
				{
				GEMV_N(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+1, nb, 1.0, lv, 0, lv, 0);
				}

			VECCP(nv, lv, 0, dv, 0);
			VECSC(nv, -1.0, dv, 0);

			GEMV_DIAG(nv, 1.0, sv, 0, dv, 0, 0.0, dv, 0, dv, 0);
			TRSV_LNN(nv, Lv, 0, 0, dv, 0, dv, 0);
			TRSV_LTN(nv, Lv, 0, 0, dv, 0, dv, 0);
			GEMV_DIAG(nv, 1.0, sv, 0, dv, 0, 0.0, dv, 0, dv, 0);

			}
		else // no scale
			{

			VECCP(nv, res_g, 0, lv, 0);

			if(ns>0)
				{
				COND_SLACKS_SOLVE(qp, ws);
				}
			else if(nb+ng>0)
				{
				AXPY(nb+ng, -1.0, gamma, nb+ng, gamma, 0, tmp_nbg+1, 0);
				}
			if(nb>0)
				{
				VECAD_SP(nb, 1.0, tmp_nbg+1, 0, idxb, lv, 0);
				}
			if(ng>0)
				{
				GEMV_N(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+1, nb, 1.0, lv, 0, lv, 0);
				}

			VECCP(nv, lv, 0, dv, 0);
			VECSC(nv, -1.0, dv, 0);

			TRSV_LNN(nv, Lv, 0, 0, dv, 0, dv, 0);
			TRSV_LTN(nv, Lv, 0, 0, dv, 0, dv, 0);

			} // scale

		} // ne>0

	if(nb+ng>0)
		{
		if(nb>0)
			VECEX_SP(nb, 1.0, idxb, dv, 0, dt, 0);

		if(ng>0)
			GEMV_T(nv, ng, 1.0, Ct, 0, 0, dv, 0, 0.0, dt, nb, dt, nb);

		VECCP(nb+ng, dt, 0, dt, nb+ng);
		VECSC(nb+ng, -1.0, dt, nb+ng);

		if(ns>0)
			EXPAND_SLACKS(qp, qp_sol, ws);

		COMPUTE_LAM_T_QP(qp->d->pa, qp->m->pa, qp_sol->lam->pa, qp_sol->t->pa, cws);
		}

	return;

	}


