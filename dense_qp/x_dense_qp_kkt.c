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



void INIT_VAR_HARD_DENSE_QP(struct DENSE_QP *qp, struct DENSE_QP_SOL *qp_sol, struct IPM_HARD_DENSE_QP_WORKSPACE *ws)
	{

	struct IPM_HARD_CORE_QP_WORKSPACE *rws = ws->core_workspace;

	// extract rws members
	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;

	REAL *d_lb = qp->d_lb->pa;
	REAL *d_ub = qp->d_ub->pa;
	REAL *d_lg = qp->d_lg->pa;
	REAL *d_ug = qp->d_ug->pa;
	int *idxb = qp->idxb;

	REAL *v = qp_sol->v->pa;
	REAL *pi = qp_sol->pi->pa;
	REAL *lam_lb = qp_sol->lam_lb->pa;
	REAL *lam_ub = qp_sol->lam_ub->pa;
	REAL *lam_lg = qp_sol->lam_lg->pa;
	REAL *lam_ug = qp_sol->lam_ug->pa;
	REAL *t_lb = qp_sol->t_lb->pa;
	REAL *t_ub = qp_sol->t_ub->pa;
	REAL *t_lg = qp_sol->t_lg->pa;
	REAL *t_ug = qp_sol->t_ug->pa;

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
	GEMV_T_LIBSTR(nv, ng, 1.0, qp->Ct, 0, 0, qp_sol->v, 0, 0.0, qp_sol->t_lg, 0, qp_sol->t_lg, 0);
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



void COMPUTE_RES_HARD_DENSE_QP(struct DENSE_QP *qp, struct DENSE_QP_SOL *qp_sol, struct IPM_HARD_DENSE_QP_WORKSPACE *ws)
	{

	struct IPM_HARD_CORE_QP_WORKSPACE *cws = ws->core_workspace;

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;

	// TODO extract qp arguments !!!!!

	struct STRVEC *v = qp_sol->v;
	struct STRVEC *pi = qp_sol->pi;
	struct STRVEC *lam_lb = qp_sol->lam_lb;
	struct STRVEC *lam_ub = qp_sol->lam_ub;
	struct STRVEC *lam_lg = qp_sol->lam_lg;
	struct STRVEC *lam_ug = qp_sol->lam_ug;
	struct STRVEC *t_lb = qp_sol->t_lb;
	struct STRVEC *t_ub = qp_sol->t_ub;
	struct STRVEC *t_lg = qp_sol->t_lg;
	struct STRVEC *t_ug = qp_sol->t_ug;

	REAL mu;

	// res g
	SYMV_L_LIBSTR(nv, nv, 1.0, qp->Hg, 0, 0, v, 0, 1.0, qp->g, 0, ws->res_g, 0);

	if(nb>0)
		{
		// res_g
		AXPY_LIBSTR(nb, -1.0, lam_lb, 0, lam_ub, 0, ws->tmp_nb, 0);
		VECAD_SP_LIBSTR(nb, 1.0, ws->tmp_nb, 0, qp->idxb, ws->res_g, 0);
		// res_d
		VECEX_SP_LIBSTR(nb, -1.0, qp->idxb, v, 0, ws->res_d_lb, 0);
		VECCP_LIBSTR(nb, ws->res_d_lb, 0, ws->res_d_ub, 0);
		AXPY_LIBSTR(nb, 1.0, qp->d_lb, 0, ws->res_d_lb, 0, ws->res_d_lb, 0);
		AXPY_LIBSTR(nb, 1.0, qp->d_ub, 0, ws->res_d_ub, 0, ws->res_d_ub, 0);
		AXPY_LIBSTR(nb, 1.0, t_lb, 0, ws->res_d_lb, 0, ws->res_d_lb, 0);
		AXPY_LIBSTR(nb, -1.0, t_ub, 0, ws->res_d_ub, 0, ws->res_d_ub, 0);
		}

	if(ng>0)
		{
		AXPY_LIBSTR(ng, -1.0, lam_lg, 0, lam_ug, 0, ws->tmp_ng0, 0);
		AXPY_LIBSTR(ng, 1.0, t_lg, 0, qp->d_lg, 0, ws->res_d_lg, 0);
		AXPY_LIBSTR(ng, -1.0, t_ug, 0, qp->d_ug, 0, ws->res_d_ug, 0);
		GEMV_NT_LIBSTR(nv, ng, 1.0, 1.0, qp->Ct, 0, 0, ws->tmp_ng0, 0, v, 0, 1.0, 0.0, ws->res_g, 0, ws->tmp_ng1, 0, ws->res_g, 0, ws->tmp_ng1, 0);
		AXPY_LIBSTR(ng, -1.0, ws->tmp_ng1, 0, ws->res_d_lg, 0, ws->res_d_lg, 0);
		AXPY_LIBSTR(ng, -1.0, ws->tmp_ng1, 0, ws->res_d_ug, 0, ws->res_d_ug, 0);
		}
	
	// res b, res g
	GEMV_NT_LIBSTR(ne, nv, -1.0, -1.0, qp->A, 0, 0, v, 0, pi, 0, 1.0, 1.0, qp->b, 0, ws->res_g, 0, ws->res_b, 0, ws->res_g, 0);

	// res_mu
	mu = VECMULDOT_LIBSTR(2*nb+2*ng, lam_lb, 0, t_lb, 0, ws->res_m, 0);

	if(cws->nb+cws->ng>0)
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
	struct STRVEC *dt_lb = ws->dt_lb;
	struct STRVEC *dt_lg = ws->dt_lg;
	struct STRVEC *res_g = ws->res_g;
	struct STRVEC *res_b = ws->res_b;
	struct STRVEC *Qx = ws->Qx;
	struct STRVEC *qx = ws->qx;

	struct IPM_HARD_CORE_QP_WORKSPACE *rws = ws->core_workspace;

	if(nb>0 | ng>0)
		{
		COMPUTE_QX_QX_HARD_QP(rws);
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
		VECEX_SP_LIBSTR(nb, 1.0, idxb, dv, 0, dt_lb, 0);
		}

	if(ng>0)
		{
		GEMV_T_LIBSTR(nv, ng, 1.0, Ct, 0, 0, dv, 0, 0.0, dt_lg, 0, dt_lg, 0);
		}

	if(nb>0 | ng>0)
		{
		COMPUTE_LAM_T_HARD_QP(rws);
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
	struct STRVEC *dt_lb = ws->dt_lb;
	struct STRVEC *dt_lg = ws->dt_lg;
	struct STRVEC *res_g = ws->res_g;
	struct STRVEC *res_b = ws->res_b;
	struct STRVEC *qx = ws->qx;

	struct IPM_HARD_CORE_QP_WORKSPACE *rws = ws->core_workspace;

	if(nb>0 | ng>0)
		{
		COMPUTE_QX_HARD_QP(rws);
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
		VECEX_SP_LIBSTR(nb, 1.0, idxb, dv, 0, dt_lb, 0);
		}

	if(ng>0)
		{
		GEMV_T_LIBSTR(nv, ng, 1.0, Ct, 0, 0, dv, 0, 0.0, dt_lg, 0, dt_lg, 0);
		}

	if(nb>0 | ng>0)
		{
		COMPUTE_LAM_T_HARD_QP(rws);
		}

	return;

	}


