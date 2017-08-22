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



void INIT_VAR_HARD_DENSE_QP(struct DENSE_QP *qp, struct DENSE_QP_SOL *qp_sol, struct IPM_HARD_DENSE_QP_WORKSPACE *ws)
	{

	struct IPM_CORE_QP_WORKSPACE *rws = ws->core_workspace;

	// extract rws members
	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;

	REAL *d = qp->d->pa;
	int *idxb = qp->idxb;

	REAL *v = qp_sol->v->pa;
	REAL *pi = qp_sol->pi->pa;
	REAL *lam = qp_sol->lam->pa;
	REAL *t = qp_sol->t->pa;

	REAL mu0 = rws->mu0;

	// local variables
	int ii;
	int idxb0;
	REAL thr0 = 0.1;

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

	return;

	}



void COMPUTE_RES_HARD_DENSE_QP(struct DENSE_QP *qp, struct DENSE_QP_SOL *qp_sol, struct IPM_HARD_DENSE_QP_WORKSPACE *ws)
	{

	struct IPM_CORE_QP_WORKSPACE *cws = ws->core_workspace;

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;

	// TODO extract qp arguments !!!!!

	// TODO extract ws arguments !!!!!

	struct STRVEC *v = qp_sol->v;
	struct STRVEC *pi = qp_sol->pi;
	struct STRVEC *lam = qp_sol->lam;
	struct STRVEC *t = qp_sol->t;

	REAL mu;

	// res g
	SYMV_L_LIBSTR(nv, nv, 1.0, qp->Hg, 0, 0, v, 0, 1.0, qp->g, 0, ws->res_g, 0);

	if(nb>0)
		{
		// res_g
		AXPY_LIBSTR(nb, -1.0, lam, 0, lam, nb+ng, ws->tmp_nb, 0);
		VECAD_SP_LIBSTR(nb, 1.0, ws->tmp_nb, 0, qp->idxb, ws->res_g, 0);
		// res_d
		VECEX_SP_LIBSTR(nb, -1.0, qp->idxb, v, 0, ws->res_d, 0);
		VECCP_LIBSTR(nb, ws->res_d, 0, ws->res_d, nb+ng);
		AXPY_LIBSTR(nb, 1.0, qp->d, 0, ws->res_d, 0, ws->res_d, 0);
		AXPY_LIBSTR(nb, 1.0, qp->d, nb+ng, ws->res_d, nb+ng, ws->res_d, nb+ng);
		AXPY_LIBSTR(nb, 1.0, t, 0, ws->res_d, 0, ws->res_d, 0);
		AXPY_LIBSTR(nb, -1.0, t, nb+ng, ws->res_d, nb+ng, ws->res_d, nb+ng);
		VECSC_LIBSTR(nb, -1.0, ws->res_d, nb+ng); // TODO embed with above
		}

	if(ng>0)
		{
		AXPY_LIBSTR(ng, -1.0, lam, nb, lam, 2*nb+ng, ws->tmp_ng0, 0);
		AXPY_LIBSTR(ng, 1.0, t, nb, qp->d, nb, ws->res_d, nb);
		AXPY_LIBSTR(ng, -1.0, t, 2*nb+ng, qp->d, 2*nb+ng, ws->res_d, 2*nb+ng);
		GEMV_NT_LIBSTR(nv, ng, 1.0, 1.0, qp->Ct, 0, 0, ws->tmp_ng0, 0, v, 0, 1.0, 0.0, ws->res_g, 0, ws->tmp_ng1, 0, ws->res_g, 0, ws->tmp_ng1, 0);
		AXPY_LIBSTR(ng, -1.0, ws->tmp_ng1, 0, ws->res_d, nb, ws->res_d, nb);
		AXPY_LIBSTR(ng, -1.0, ws->tmp_ng1, 0, ws->res_d, 2*nb+ng, ws->res_d, 2*nb+ng);
		VECSC_LIBSTR(ng, -1.0, ws->res_d, 2*nb+ng); // TODO embed with above
		}
	
	// res b, res g
	GEMV_NT_LIBSTR(ne, nv, -1.0, -1.0, qp->A, 0, 0, v, 0, pi, 0, 1.0, 1.0, qp->b, 0, ws->res_g, 0, ws->res_b, 0, ws->res_g, 0);

	// res_mu
	mu = VECMULDOT_LIBSTR(2*nb+2*ng, lam, 0, t, 0, ws->res_m, 0);

	if(cws->nc>0)
		ws->res_mu = mu*cws->nt_inv;
	else
		ws->res_mu = 0.0;


	return;

	}


// range-space (Schur complement) method
void FACT_SOLVE_KKT_UNCONSTR_DENSE_QP(struct DENSE_QP *qp, struct DENSE_QP_SOL *qp_sol, struct IPM_HARD_DENSE_QP_WORKSPACE *ws)
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



// range-space (Schur complement) method
void FACT_SOLVE_KKT_STEP_HARD_DENSE_QP(struct DENSE_QP *qp, struct IPM_HARD_DENSE_QP_WORKSPACE *ws)
	{

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;
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
	struct STRVEC *Qx = ws->Qx;
	struct STRVEC *qx = ws->qx;

	struct IPM_CORE_QP_WORKSPACE *rws = ws->core_workspace;

	if(nb>0 | ng>0)
		{
		COMPUTE_QX_QX_QP(rws);
		}

	if(ne>0)
		{
//		TRCP_L_LIBSTR(nv, Hg, 0, 0, Lv, 0, 0);
		GECP_LIBSTR(nv, nv, Hg, 0, 0, Lv, 0, 0);

		VECCP_LIBSTR(nv, res_g, 0, lv, 0);

		if(nb>0)
			{
			DIAAD_SP_LIBSTR(nb, 1.0, Qx, 0, idxb, Lv, 0, 0);
			VECAD_SP_LIBSTR(nb, 1.0, qx, 0, idxb, lv, 0);
			}

		if(ng>0)
			{
			GEMV_N_LIBSTR(nv, ng, 1.0, Ct, 0, 0, qx, nb, 1.0, lv, 0, lv, 0);
			GEMM_R_DIAG_LIBSTR(nv, ng, 1.0, Ct, 0, 0, Qx, nb, 0.0, Ctx, 0, 0, Ctx, 0, 0);
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

		if(nb>0)
			{
			DIAAD_SP_LIBSTR(nb, 1.0, Qx, 0, idxb, Lv, 0, 0);
			ROWAD_SP_LIBSTR(nb, 1.0, qx, 0, idxb, Lv, nv, 0);
			}

		if(ng>0)
			{
			GEMM_R_DIAG_LIBSTR(nv, ng, 1.0, Ct, 0, 0, Qx, nb, 0.0, Ctx, 0, 0, Ctx, 0, 0);
			ROWIN_LIBSTR(ng, 1.0, qx, nb, Ctx, nv, 0);
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

	if(nb>0)
		{
		VECEX_SP_LIBSTR(nb, 1.0, idxb, dv, 0, dt, 0);
		}

	if(ng>0)
		{
		GEMV_T_LIBSTR(nv, ng, 1.0, Ct, 0, 0, dv, 0, 0.0, dt, nb, dt, nb);
		}

	if(nb>0 | ng>0)
		{
		VECCP_LIBSTR(nb+ng, dt, 0, dt, nb+ng);
		VECSC_LIBSTR(nb+ng, -1.0, dt, nb+ng);
		COMPUTE_LAM_T_QP(rws);
		}

	return;

	}



// range-space (Schur complement) method
void SOLVE_KKT_STEP_HARD_DENSE_QP(struct DENSE_QP *qp, struct IPM_HARD_DENSE_QP_WORKSPACE *ws)
	{

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;
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
	struct STRVEC *qx = ws->qx;

	struct IPM_CORE_QP_WORKSPACE *rws = ws->core_workspace;

	if(nb>0 | ng>0)
		{
		COMPUTE_QX_QP(rws);
		}

	if(ne>0)
		{
		VECCP_LIBSTR(nv, res_g, 0, lv, 0);

		if(nb>0)
			{
			VECAD_SP_LIBSTR(nb, 1.0, qx, 0, idxb, lv, 0);
			}

		if(ng>0)
			{
			GEMV_N_LIBSTR(nv, ng, 1.0, Ct, 0, 0, qx, nb, 1.0, lv, 0, lv, 0);
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

		if(nb>0)
			{
			VECAD_SP_LIBSTR(nb, 1.0, qx, 0, idxb, lv, 0);
			}

		if(ng>0)
			{
			GEMV_N_LIBSTR(nv, ng, 1.0, Ct, 0, 0, qx, nb, 1.0, lv, 0, lv, 0);
			}

		VECCP_LIBSTR(nv, lv, 0, dv, 0);
		VECSC_LIBSTR(nv, -1.0, dv, 0);

		TRSV_LNN_LIBSTR(nv, Lv, 0, 0, dv, 0, dv, 0);
		TRSV_LTN_LIBSTR(nv, Lv, 0, 0, dv, 0, dv, 0);
		}

	if(nb>0)
		{
		VECEX_SP_LIBSTR(nb, 1.0, idxb, dv, 0, dt, 0);
		}

	if(ng>0)
		{
		GEMV_T_LIBSTR(nv, ng, 1.0, Ct, 0, 0, dv, 0, 0.0, dt, nb, dt, nb);
		}

	if(nb>0 | ng>0)
		{
		VECCP_LIBSTR(nb+ng, dt, 0, dt, nb+ng);
		VECSC_LIBSTR(nb+ng, -1.0, dt, nb+ng);
		COMPUTE_LAM_T_QP(rws);
		}

	return;

	}


