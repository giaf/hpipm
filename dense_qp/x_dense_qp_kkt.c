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



void INIT_VAR_DENSE_QP(struct DENSE_QP *qp, struct DENSE_QP_SOL *qp_sol, struct IPM_DENSE_QP_WORKSPACE *ws)
	{

	struct IPM_CORE_QP_WORKSPACE *cws = ws->core_workspace;

	// extract cws members
	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;
	int ns = qp->ns;

	REAL *d = qp->d->pa;
	int *idxb = qp->idxb;

	REAL *v = qp_sol->v->pa;
	REAL *pi = qp_sol->pi->pa;
	REAL *lam = qp_sol->lam->pa;
	REAL *t = qp_sol->t->pa;

	REAL mu0 = cws->mu0;

	// local variables
	int ii;
	int idxb0;
	REAL thr0 = 0.1;

	// warm start TODO


	// cold start

	// primal variables
	for(ii=0; ii<nv+2*ns; ii++)
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
		t[0+ii]     = - d[0+ii]     + v[idxb0];
		t[nb+ng+ii] =   d[nb+ng+ii] - v[idxb0];
		if(t[0+ii]<thr0)
			{
			if(t[nb+ng+ii]<thr0)
				{
				v[idxb0] = 0.5*(d[0+ii] - d[nb+ng+ii]);
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
			v[idxb0] = d[nb+ng+ii] - thr0;
			}
		lam[0+ii]     = mu0/t[0+ii];
		lam[nb+ng+ii] = mu0/t[nb+ng+ii];
		}
	
	// inequality constraints
	GEMV_T_LIBSTR(nv, ng, 1.0, qp->Ct, 0, 0, qp_sol->v, 0, 0.0, qp_sol->t, nb, qp_sol->t, nb);
	for(ii=0; ii<ng; ii++)
		{
		t[2*nb+ng+ii] = t[nb+ii];
		t[nb+ii]      -= d[nb+ii];
		t[2*nb+ng+ii] += d[2*nb+ng+ii];
//		t[nb+ii]      = fmax( thr0, t[nb+ii] );
//		t[2*nb+ng+ii] = fmax( thr0, t[2*nb+ng+ii] );
		t[nb+ii]      = thr0>t[nb+ii]      ? thr0 : t[nb+ii];
		t[2*nb+ng+ii] = thr0>t[2*nb+ng+ii] ? thr0 : t[2*nb+ng+ii];
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



void COMPUTE_RES_DENSE_QP(struct DENSE_QP *qp, struct DENSE_QP_SOL *qp_sol, struct IPM_DENSE_QP_WORKSPACE *ws)
	{

	struct IPM_CORE_QP_WORKSPACE *cws = ws->core_workspace;

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;
	int ns = qp->ns;

	int nct = ws->core_workspace->nc;

	struct STRMAT *Hg = qp->Hg;
	struct STRVEC *g = qp->g;
	struct STRMAT *A = qp->A;
	struct STRVEC *b = qp->b;
	struct STRMAT *Ct = qp->Ct;
	struct STRVEC *d = qp->d;
	int *idxb = qp->idxb;
	struct STRVEC *Z = qp->Z;
	struct STRVEC *z = qp->z;
	int *idxs = qp->idxs;

	struct STRVEC *v = qp_sol->v;
	struct STRVEC *pi = qp_sol->pi;
	struct STRVEC *lam = qp_sol->lam;
	struct STRVEC *t = qp_sol->t;

	struct STRVEC *res_g = ws->res_g;
	struct STRVEC *res_b = ws->res_b;
	struct STRVEC *res_d = ws->res_d;
	struct STRVEC *res_m = ws->res_m;
	struct STRVEC *tmp_nbg = ws->tmp_nbg;
	struct STRVEC *tmp_ns = ws->tmp_ns;

	REAL mu;

	// res g
	SYMV_L_LIBSTR(nv, nv, 1.0, Hg, 0, 0, v, 0, 1.0, g, 0, res_g, 0);

	if(nb+ng>0)
		{
		AXPY_LIBSTR(nb+ng, -1.0, lam, 0, lam, nb+ng, tmp_nbg+0, 0);
		AXPY_LIBSTR(nb+ng,  1.0, d, 0, t, 0, res_d, 0);
		AXPY_LIBSTR(nb+ng, -1.0, d, nb+ng, t, nb+ng, res_d, nb+ng);
		// box
		if(nb>0)
			{
			VECAD_SP_LIBSTR(nb, 1.0, tmp_nbg+0, 0, idxb, res_g, 0);
			VECEX_SP_LIBSTR(nb, 1.0, idxb, v, 0, tmp_nbg+1, 0);
			}
		// general
		if(ng>0)
			{
			GEMV_NT_LIBSTR(nv, ng, 1.0, 1.0, Ct, 0, 0, tmp_nbg+0, 0, v, 0, 1.0, 0.0, res_g, 0, tmp_nbg+1, nb, res_g, 0, tmp_nbg+1, nb);
			}
		AXPY_LIBSTR(nb+ng, -1.0, tmp_nbg+1, 0, res_d, 0, res_d, 0);
		AXPY_LIBSTR(nb+ng,  1.0, tmp_nbg+1, 0, res_d, nb+ng, res_d, nb+ng);
		}
	if(ns>0)
		{
		// res_g
		GEMV_DIAG_LIBSTR(2*ns, 1.0, Z, 0, v, nv, 1.0, z, 0, res_g, nv);
		AXPY_LIBSTR(2*ns, -1.0, lam, 2*nb+2*ng, res_g, nv, res_g, nv);
		VECEX_SP_LIBSTR(ns, 1.0, idxs, lam, 0, tmp_ns, 0);
		AXPY_LIBSTR(ns, -1.0, tmp_ns, 0, res_g, nv, res_g, nv);
		VECEX_SP_LIBSTR(ns, 1.0, idxs, lam, nb+ng, tmp_ns, 0);
		AXPY_LIBSTR(ns, -1.0, tmp_ns, 0, res_g, nv, res_g, nv+ns);
		// res_d
		VECAD_SP_LIBSTR(ns, -1.0, v, nv, idxs, res_d, 0);
		VECAD_SP_LIBSTR(ns, -1.0, v, nv+ns, idxs, res_d, nb+ng);
		AXPY_LIBSTR(2*ns, -1.0, v, nv, t, 2*nb+2*ng, res_d, 2*nb+2*ng);
		}
	
	// res b, res g
	GEMV_NT_LIBSTR(ne, nv, -1.0, -1.0, A, 0, 0, v, 0, pi, 0, 1.0, 1.0, b, 0, res_g, 0, res_b, 0, res_g, 0);

	// res_mu
	mu = VECMULDOT_LIBSTR(2*nct, lam, 0, t, 0, res_m, 0); // XXX

	if(nct>0)
		ws->res_mu = mu*cws->nt_inv;
	else
		ws->res_mu = 0.0;


	return;

	}


// range-space (Schur complement) method
void FACT_SOLVE_KKT_UNCONSTR_DENSE_QP(struct DENSE_QP *qp, struct DENSE_QP_SOL *qp_sol, struct IPM_DENSE_QP_WORKSPACE *ws)
	{

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;
	struct STRMAT *Hg = qp->Hg;
	struct STRMAT *A = qp->A;
	struct STRVEC *g = qp->g;
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
		POTRF_L_LIBSTR(nv, Hg, 0, 0, Lv, 0, 0);

		GECP_LIBSTR(ne, nv, A, 0, 0, AL, 0, 0);
		TRSM_RLTN_LIBSTR(ne, nv, 1.0, Lv, 0, 0, A, 0, 0, AL, 0, 0);

		GESE_LIBSTR(ne, ne, 0.0, Le, 0, 0);
		SYRK_POTRF_LN_LIBSTR(ne, ne, nv, AL, 0, 0, AL, 0, 0, Le, 0, 0, Le, 0, 0);

		TRSV_LNN_LIBSTR(nv, Lv, 0, 0, g, 0, lv, 0);

		GEMV_N_LIBSTR(ne, nv, 1.0, AL, 0, 0, lv, 0, 1.0, b, 0, pi, 0);

		TRSV_LNN_LIBSTR(ne, Le, 0, 0, pi, 0, pi, 0);
		TRSV_LTN_LIBSTR(ne, Le, 0, 0, pi, 0, pi, 0);

		GEMV_T_LIBSTR(ne, nv, 1.0, A, 0, 0, pi, 0, -1.0, g, 0, v, 0);

		TRSV_LNN_LIBSTR(nv, Lv, 0, 0, v, 0, v, 0);
		TRSV_LTN_LIBSTR(nv, Lv, 0, 0, v, 0, v, 0);
		}
	else
		{
#if 0
		POTRF_L_LIBSTR(nv, Hg, 0, 0, Lv, 0, 0);

		VECCP_LIBSTR(nv, g, 0, v, 0);
		VECSC_LIBSTR(nv, -1.0, v, 0);

		TRSV_LNN_LIBSTR(nv, Lv, 0, 0, v, 0, v, 0);
		TRSV_LTN_LIBSTR(nv, Lv, 0, 0, v, 0, v, 0);
#else
		POTRF_L_MN_LIBSTR(nv+1, nv, Hg, 0, 0, Lv, 0, 0);

		ROWEX_LIBSTR(nv, -1.0, Lv, nv, 0, v, 0);
		TRSV_LTN_LIBSTR(nv, Lv, 0, 0, v, 0, v, 0);
#endif
		}

	return;

	}



static void COND_SLACKS_FACT_SOLVE(struct DENSE_QP *qp, struct IPM_DENSE_QP_WORKSPACE *ws)
	{

	int ii, idx;

	int nv = qp->nv;
	int nb = qp->nb;
	int ng = qp->ng;
	int ns = qp->ns;

	struct STRVEC *Z = qp->Z;
	int *idxs = qp->idxs;

	struct STRVEC *dv = ws->dv;
	struct STRVEC *res_g = ws->res_g;
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

	VECCP_LIBSTR(nb+ng, Gamma, 0, tmp_nbg+0, 0);
	VECCP_LIBSTR(nb+ng, Gamma, nb+ng, tmp_nbg+1, 0);
	VECCP_LIBSTR(nb+ng, gamma, 0, tmp_nbg+2, 0);
	VECCP_LIBSTR(nb+ng, gamma, nb+ng, tmp_nbg+3, 0);

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
	
	AXPY_LIBSTR(nb+ng,  1.0, tmp_nbg+1, 0, tmp_nbg+0, 0, tmp_nbg+0, 0);
	AXPY_LIBSTR(nb+ng, -1.0, tmp_nbg+3, 0, tmp_nbg+2, 0, tmp_nbg+1, 0);

	return;

	}



static void COND_SLACKS_SOLVE(struct DENSE_QP *qp, struct IPM_DENSE_QP_WORKSPACE *ws)
	{

	int ii, idx;

	int nv = qp->nv;
	int nb = qp->nb;
	int ng = qp->ng;
	int ns = qp->ns;

	int *idxs = qp->idxs;

	struct STRVEC *dv = ws->dv;
	struct STRVEC *res_g = ws->res_g;
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

	VECCP_LIBSTR(nb+ng, gamma, 0, tmp_nbg+2, 0);
	VECCP_LIBSTR(nb+ng, gamma, nb+ng, tmp_nbg+3, 0);

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
	
	AXPY_LIBSTR(nb+ng, -1.0, tmp_nbg+3, 0, tmp_nbg+2, 0, tmp_nbg+1, 0);

	return;

	}



static void EXPAND_SLACKS(struct DENSE_QP *qp, struct IPM_DENSE_QP_WORKSPACE *ws)
	{

	int ii, idx;

	int nv = qp->nv;
	int nb = qp->nb;
	int ng = qp->ng;
	int ns = qp->ns;

	int *idxs = qp->idxs;

	struct STRVEC *dv = ws->dv;
	struct STRVEC *dt = ws->dt;
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
void FACT_SOLVE_KKT_STEP_DENSE_QP(struct DENSE_QP *qp, struct IPM_DENSE_QP_WORKSPACE *ws)
	{

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;
	int ns = qp->ns;
	struct STRMAT *Hg = qp->Hg;
	struct STRMAT *A = qp->A;
	struct STRMAT *Ct = qp->Ct;
	int *idxb = qp->idxb;

	struct STRMAT *Lv = ws->Lv;
	struct STRMAT *Le = ws->Le;
	struct STRMAT *Ctx = ws->Ctx;
	struct STRMAT *AL = ws->AL;
	struct STRVEC *lv = ws->lv;
	struct STRVEC *dv = ws->dv;
	struct STRVEC *dpi = ws->dpi;
	struct STRVEC *dt = ws->dt;
	struct STRVEC *res_g = ws->res_g;
	struct STRVEC *res_b = ws->res_b;
	struct STRVEC *Gamma = ws->Gamma;
	struct STRVEC *gamma = ws->gamma;
	struct STRVEC *tmp_nbg = ws->tmp_nbg;

	struct IPM_CORE_QP_WORKSPACE *cws = ws->core_workspace;

	if(nb+ng>0)
		{
		COMPUTE_QX_QX_QP(cws);
		}

	if(ne>0)
		{
//		TRCP_L_LIBSTR(nv, Hg, 0, 0, Lv, 0, 0);
		GECP_LIBSTR(nv, nv, Hg, 0, 0, Lv, 0, 0);

		VECCP_LIBSTR(nv, res_g, 0, lv, 0);

		if(ns>0)
			{
			COND_SLACKS_FACT_SOLVE(qp, ws);
			}
		else if(nb+ng>0)
			{
			AXPY_LIBSTR(nb+ng,  1.0, Gamma, nb+ng, Gamma, 0, tmp_nbg+0, 0);
			AXPY_LIBSTR(nb+ng, -1.0, gamma, nb+ng, gamma, 0, tmp_nbg+1, 0);
			}
		if(nb>0)
			{
			DIAAD_SP_LIBSTR(nb, 1.0, tmp_nbg+0, 0, idxb, Lv, 0, 0);
			VECAD_SP_LIBSTR(nb, 1.0, tmp_nbg+1, 0, idxb, lv, 0);
			}
		if(ng>0)
			{
			GEMV_N_LIBSTR(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+1, nb, 1.0, lv, 0, lv, 0);
			GEMM_R_DIAG_LIBSTR(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+0, nb, 0.0, Ctx, 0, 0, Ctx, 0, 0);
			SYRK_POTRF_LN_LIBSTR(nv, nv, ng, Ctx, 0, 0, Ct, 0, 0, Lv, 0, 0, Lv, 0, 0);
			}
		else
			{
			POTRF_L_LIBSTR(nv, Lv, 0, 0, Lv, 0, 0);
			}

		VECCP_LIBSTR(nv, lv, 0, dv, 0);

		GECP_LIBSTR(ne, nv, A, 0, 0, AL, 0, 0);
		TRSM_RLTN_LIBSTR(ne, nv, 1.0, Lv, 0, 0, A, 0, 0, AL, 0, 0);

		GESE_LIBSTR(ne, ne, 0.0, Le, 0, 0);
		SYRK_POTRF_LN_LIBSTR(ne, ne, nv, AL, 0, 0, AL, 0, 0, Le, 0, 0, Le, 0, 0);

		TRSV_LNN_LIBSTR(nv, Lv, 0, 0, lv, 0, lv, 0);

		GEMV_N_LIBSTR(ne, nv, 1.0, AL, 0, 0, lv, 0, 1.0, res_b, 0, dpi, 0);

		TRSV_LNN_LIBSTR(ne, Le, 0, 0, dpi, 0, dpi, 0);
		TRSV_LTN_LIBSTR(ne, Le, 0, 0, dpi, 0, dpi, 0);

		GEMV_T_LIBSTR(ne, nv, 1.0, A, 0, 0, dpi, 0, -1.0, dv, 0, dv, 0);

		TRSV_LNN_LIBSTR(nv, Lv, 0, 0, dv, 0, dv, 0);
		TRSV_LTN_LIBSTR(nv, Lv, 0, 0, dv, 0, dv, 0);
		}
	else
		{
#if 0
		TRCP_L_LIBSTR(nv, Hg, 0, 0, Lv, 0, 0);
		VECCP_LIBSTR(nv, res_g, 0, lv, 0);

		if(nb>0)
			{
			DIAAD_SP_LIBSTR(nb, 1.0, Qx, 0, idxb, Lv, 0, 0);
			VECAD_SP_LIBSTR(nb, 1.0, qx, 0, idxb, lv, 0);
			}

		if(ng>0)
			{
			GEMM_R_DIAG_LIBSTR(nv, ng, 1.0, Ct, 0, 0, Qx, nb, 0.0, Ctx, 0, 0, Ctx, 0, 0);
			GEMV_N_LIBSTR(nv, ng, 1.0, Ct, 0, 0, qx, nb, 1.0, lv, 0, lv, 0);
			SYRK_POTRF_LN_LIBSTR(nv, nv, ng, Ctx, 0, 0, Ct, 0, 0, Lv, 0, 0, Lv, 0, 0); // TODO _mn_ routine in BLASFEO !!!
			}
		else
			{
			POTRF_L_LIBSTR(nv, Lv, 0, 0, Lv, 0, 0);
			}

		VECCP_LIBSTR(nv, lv, 0, dv, 0);
		VECSC_LIBSTR(nv, -1.0, dv, 0);

		TRSV_LNN_LIBSTR(nv, Lv, 0, 0, dv, 0, dv, 0);
		TRSV_LTN_LIBSTR(nv, Lv, 0, 0, dv, 0, dv, 0);
#else
//		TRCP_L_LIBSTR(nv, Hg, 0, 0, Lv, 0, 0);
		GECP_LIBSTR(nv, nv, Hg, 0, 0, Lv, 0, 0);
		ROWIN_LIBSTR(nv, 1.0, res_g, 0, Lv, nv, 0);

		if(ns>0)
			{
			COND_SLACKS_FACT_SOLVE(qp, ws);
			}
		else if(nb+ng>0)
			{
			AXPY_LIBSTR(nb+ng,  1.0, Gamma, nb+ng, Gamma, 0, tmp_nbg+0, 0);
			AXPY_LIBSTR(nb+ng, -1.0, gamma, nb+ng, gamma, 0, tmp_nbg+1, 0);
			}
		if(nb>0)
			{
			DIAAD_SP_LIBSTR(nb, 1.0, tmp_nbg+0, 0, idxb, Lv, 0, 0);
			ROWAD_SP_LIBSTR(nb, 1.0, tmp_nbg+1, 0, idxb, Lv, nv, 0);
			}
		if(ng>0)
			{
			GEMM_R_DIAG_LIBSTR(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+0, nb, 0.0, Ctx, 0, 0, Ctx, 0, 0);
			ROWIN_LIBSTR(ng, 1.0, tmp_nbg+1, nb, Ctx, nv, 0);
			SYRK_POTRF_LN_LIBSTR(nv+1, nv, ng, Ctx, 0, 0, Ct, 0, 0, Lv, 0, 0, Lv, 0, 0); // TODO _mn_ routine in BLASFEO !!!
			}
		else
			{
			POTRF_L_MN_LIBSTR(nv+1, nv, Lv, 0, 0, Lv, 0, 0);
			}

		ROWEX_LIBSTR(nv, -1.0, Lv, nv, 0, dv, 0);
		TRSV_LTN_LIBSTR(nv, Lv, 0, 0, dv, 0, dv, 0);
#endif
		}

	if(nb+ng>0)
		{
		if(nb>0)
			VECEX_SP_LIBSTR(nb, 1.0, idxb, dv, 0, dt, 0);

		if(ng>0)
			GEMV_T_LIBSTR(nv, ng, 1.0, Ct, 0, 0, dv, 0, 0.0, dt, nb, dt, nb);

		VECCP_LIBSTR(nb+ng, dt, 0, dt, nb+ng);
		VECSC_LIBSTR(nb+ng, -1.0, dt, nb+ng);

		if(ns>0)
			EXPAND_SLACKS(qp, ws);

		COMPUTE_LAM_T_QP(cws);
		}

	return;

	}



// range-space (Schur complement) method
void SOLVE_KKT_STEP_DENSE_QP(struct DENSE_QP *qp, struct IPM_DENSE_QP_WORKSPACE *ws)
	{

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;
	int ns = qp->ns;
	struct STRMAT *A = qp->A;
	struct STRMAT *Ct = qp->Ct;
	int *idxb = qp->idxb;

	struct STRMAT *Lv = ws->Lv;
	struct STRMAT *Le = ws->Le;
	struct STRMAT *Ctx = ws->Ctx;
	struct STRMAT *AL = ws->AL;
	struct STRVEC *lv = ws->lv;
	struct STRVEC *dv = ws->dv;
	struct STRVEC *dpi = ws->dpi;
	struct STRVEC *dt = ws->dt;
	struct STRVEC *res_g = ws->res_g;
	struct STRVEC *res_b = ws->res_b;
	struct STRVEC *gamma = ws->gamma;
	struct STRVEC *tmp_nbg = ws->tmp_nbg;

	struct IPM_CORE_QP_WORKSPACE *cws = ws->core_workspace;

	if(nb>0 | ng>0)
		{
		COMPUTE_QX_QP(cws);
		}

	if(ne>0)
		{
		VECCP_LIBSTR(nv, res_g, 0, lv, 0);

		if(ns>0)
			{
			COND_SLACKS_SOLVE(qp, ws);
			}
		else if(nb+ng>0)
			{
			AXPY_LIBSTR(nb+ng, -1.0, gamma, nb+ng, gamma, 0, tmp_nbg+1, 0);
			}
		if(nb>0)
			{
			VECAD_SP_LIBSTR(nb, 1.0, tmp_nbg+1, 0, idxb, lv, 0);
			}
		if(ng>0)
			{
			GEMV_N_LIBSTR(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+1, nb, 1.0, lv, 0, lv, 0);
			}

		VECCP_LIBSTR(nv, lv, 0, dv, 0);

		TRSV_LNN_LIBSTR(nv, Lv, 0, 0, lv, 0, lv, 0);

		GEMV_N_LIBSTR(ne, nv, 1.0, AL, 0, 0, lv, 0, 1.0, res_b, 0, dpi, 0);

		TRSV_LNN_LIBSTR(ne, Le, 0, 0, dpi, 0, dpi, 0);
		TRSV_LTN_LIBSTR(ne, Le, 0, 0, dpi, 0, dpi, 0);

		GEMV_T_LIBSTR(ne, nv, 1.0, A, 0, 0, dpi, 0, -1.0, dv, 0, dv, 0);

		TRSV_LNN_LIBSTR(nv, Lv, 0, 0, dv, 0, dv, 0);
		TRSV_LTN_LIBSTR(nv, Lv, 0, 0, dv, 0, dv, 0);
		}
	else
		{
		VECCP_LIBSTR(nv, res_g, 0, lv, 0);

		if(ns>0)
			{
			COND_SLACKS_SOLVE(qp, ws);
			}
		else if(nb+ng>0)
			{
			AXPY_LIBSTR(nb+ng, -1.0, gamma, nb+ng, gamma, 0, tmp_nbg+1, 0);
			}
		if(nb>0)
			{
			VECAD_SP_LIBSTR(nb, 1.0, tmp_nbg+1, 0, idxb, lv, 0);
			}
		if(ng>0)
			{
			GEMV_N_LIBSTR(nv, ng, 1.0, Ct, 0, 0, tmp_nbg+1, nb, 1.0, lv, 0, lv, 0);
			}

		VECCP_LIBSTR(nv, lv, 0, dv, 0);
		VECSC_LIBSTR(nv, -1.0, dv, 0);

		TRSV_LNN_LIBSTR(nv, Lv, 0, 0, dv, 0, dv, 0);
		TRSV_LTN_LIBSTR(nv, Lv, 0, 0, dv, 0, dv, 0);
		}

	if(nb+ng>0)
		{
		if(nb>0)
			VECEX_SP_LIBSTR(nb, 1.0, idxb, dv, 0, dt, 0);

		if(ng>0)
			GEMV_T_LIBSTR(nv, ng, 1.0, Ct, 0, 0, dv, 0, 0.0, dt, nb, dt, nb);

		VECCP_LIBSTR(nb+ng, dt, 0, dt, nb+ng);
		VECSC_LIBSTR(nb+ng, -1.0, dt, nb+ng);

		if(ns>0)
			EXPAND_SLACKS(qp, ws);

		COMPUTE_LAM_T_QP(cws);
		}

	return;

	}


