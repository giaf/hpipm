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



void INIT_VAR_HARD_OCP_QP(struct OCP_QP *qp, struct OCP_QP_SOL *qp_sol, struct IPM_HARD_OCP_QP_WORKSPACE *ws)
	{

	struct IPM_HARD_CORE_QP_WORKSPACE *cws = ws->core_workspace;
	
	// loop index
	int ii, jj;

	//
	int N = qp->N;
	int *nx = qp->nx;
	int *nu = qp->nu;
	int *nb = qp->nb;
	int *ng = qp->ng;

	REAL mu0 = cws->mu0;

	//
	REAL *ux, *pi, *d_lb, *d_ub, *d_lg, *d_ug, *lam_lb, *lam_ub, *lam_lg, *lam_ug, *t_lb, *t_ub, *t_lg, *t_ug;
	int *idxb;

	REAL thr0 = 0.1;

	// warm start TODO

	// cold start

	// ux
	for(ii=0; ii<=N; ii++)
		{
		ux = qp_sol->ux[ii].pa;
		for(jj=0; jj<nu[ii]+nx[ii]; jj++)
			{
			ux[jj] = 0.0;
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
		d_lb = qp->d_lb[ii].pa;
		d_ub = qp->d_ub[ii].pa;
		lam_lb = qp_sol->lam_lb[ii].pa;
		lam_ub = qp_sol->lam_ub[ii].pa;
		t_lb = qp_sol->t_lb[ii].pa;
		t_ub = qp_sol->t_ub[ii].pa;
		idxb = qp->idxb[ii];
		for(jj=0; jj<nb[ii]; jj++)
			{
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
			lam_lb[jj] = mu0/t_lb[jj];
			lam_ub[jj] = mu0/t_ub[jj];
			}
		}
	
	// general constraints
	for(ii=0; ii<=N; ii++)
		{
		t_lg = qp_sol->t_lg[ii].pa;
		t_ug = qp_sol->t_ug[ii].pa;
		lam_lg = qp_sol->lam_lg[ii].pa;
		lam_ug = qp_sol->lam_ug[ii].pa;
		d_lg = qp->d_lg[ii].pa;
		d_ug = qp->d_ug[ii].pa;
		ux = qp_sol->ux[ii].pa;
		GEMV_T_LIBSTR(nu[ii]+nx[ii], ng[ii], 1.0, qp->DCt+ii, 0, 0, qp_sol->ux+ii, 0, 0.0, qp_sol->t_lg+ii, 0, qp_sol->t_lg+ii, 0);
		for(jj=0; jj<ng[ii]; jj++)
			{
			t_ug[jj] = - t_lg[jj];
			t_lg[jj] -= d_lg[jj];
			t_ug[jj] += d_ug[jj];
//			t_lg[jj] = fmax(thr0, t_lg[jj]);
//			t_ug[jj] = fmax(thr0, t_ug[jj]);
			t_lg[jj] = thr0>t_lg[jj] ? thr0 : t_lg[jj];
			t_ug[jj] = thr0>t_ug[jj] ? thr0 : t_ug[jj];
			lam_lg[jj] = mu0/t_lg[jj];
			lam_ug[jj] = mu0/t_ug[jj];
			}
		}


	return;

	}



void COMPUTE_RES_HARD_OCP_QP(struct OCP_QP *qp, struct OCP_QP_SOL *qp_sol, struct IPM_HARD_OCP_QP_WORKSPACE *ws)
	{

	struct IPM_HARD_CORE_QP_WORKSPACE *cws = ws->core_workspace;
	
	// loop index
	int ii;

	//
	int N = qp->N;
	int *nx = qp->nx;
	int *nu = qp->nu;
	int *nb = qp->nb;
	int *ng = qp->ng;

	int nbt = ws->core_workspace->nb;
	int ngt = ws->core_workspace->ng;

	struct STRMAT *BAbt = qp->BAbt;
	struct STRMAT *RSQrq = qp->RSQrq;
	struct STRMAT *DCt = qp->DCt;
	struct STRVEC *b = qp->b;
	struct STRVEC *rq = qp->rq;
	struct STRVEC *d_lb = qp->d_lb;
	struct STRVEC *d_ub = qp->d_ub;
	struct STRVEC *d_lg = qp->d_lg;
	struct STRVEC *d_ug = qp->d_ug;
	int **idxb = qp->idxb;

	struct STRVEC *ux = qp_sol->ux;
	struct STRVEC *pi = qp_sol->pi;
	struct STRVEC *lam_lb = qp_sol->lam_lb;
	struct STRVEC *lam_ub = qp_sol->lam_ub;
	struct STRVEC *lam_lg = qp_sol->lam_lg;
	struct STRVEC *lam_ug = qp_sol->lam_ug;
	struct STRVEC *t_lb = qp_sol->t_lb;
	struct STRVEC *t_ub = qp_sol->t_ub;
	struct STRVEC *t_lg = qp_sol->t_lg;
	struct STRVEC *t_ug = qp_sol->t_ug;

	struct STRVEC *res_g = ws->res_g;
	struct STRVEC *res_b = ws->res_b;
	struct STRVEC *res_d_lb = ws->res_d_lb;
	struct STRVEC *res_d_ub = ws->res_d_ub;
	struct STRVEC *res_d_lg = ws->res_d_lg;
	struct STRVEC *res_d_ug = ws->res_d_ug;
	struct STRVEC *tmp_ngM = ws->tmp_ngM;
	struct STRVEC *tmp_nbM = ws->tmp_nbM;

	int nx0, nx1, nu0, nu1, nb0, ng0;

	//
	REAL mu = 0.0;

	// loop over stages
	for(ii=0; ii<=N; ii++)
		{

		nx0 = nx[ii];
		nu0 = nu[ii];
		nb0 = nb[ii];
		ng0 = ng[ii];

		VECCP_LIBSTR(nu0+nx0, rq+ii, 0, res_g+ii, 0);

		if(ii>0)
			AXPY_LIBSTR(nx0, -1.0, pi+(ii-1), 0, res_g+ii, nu0, res_g+ii, nu0);

		SYMV_L_LIBSTR(nu0+nx0, nu0+nx0, 1.0, RSQrq+ii, 0, 0, ux+ii, 0, 1.0, res_g+ii, 0, res_g+ii, 0);

		if(nb0>0)
			{

			AXPY_LIBSTR(nb0, -1.0, lam_lb+ii, 0, lam_ub+ii, 0, tmp_nbM, 0);
			VECAD_SP_LIBSTR(nb0, 1.0, tmp_nbM, 0, idxb[ii], res_g+ii, 0);

			VECEX_SP_LIBSTR(nb0, -1.0, idxb[ii], ux+ii, 0, res_d_lb+ii, 0);
			VECCP_LIBSTR(nb0, res_d_lb+ii, 0, res_d_ub+ii, 0);
			AXPY_LIBSTR(nb0, 1.0, d_lb+ii, 0, res_d_lb+ii, 0, res_d_lb+ii, 0);
			AXPY_LIBSTR(nb0, 1.0, d_ub+ii, 0, res_d_ub+ii, 0, res_d_ub+ii, 0);
			AXPY_LIBSTR(nb0, 1.0, t_lb+ii, 0, res_d_lb+ii, 0, res_d_lb+ii, 0);
			AXPY_LIBSTR(nb0, -1.0, t_ub+ii, 0, res_d_ub+ii, 0, res_d_ub+ii, 0);

			}

		if(ng0>0) // TODO merge with bounds as much as possible
			{

			AXPY_LIBSTR(ng0, -1.0, lam_lg+ii, 0, lam_ug+ii, 0, tmp_ngM+0, 0);

			AXPY_LIBSTR(ng0,  1.0, t_lg+ii, 0, d_lg+ii, 0, res_d_lg+ii, 0);
			AXPY_LIBSTR(ng0, -1.0, t_ug+ii, 0, d_ug+ii, 0, res_d_ug+ii, 0);

			GEMV_NT_LIBSTR(nu0+nx0, ng0, 1.0, 1.0, DCt+ii, 0, 0, tmp_ngM+0, 0, ux+ii, 0, 1.0, 0.0, res_g+ii, 0, tmp_ngM+1, 0, res_g+ii, 0, tmp_ngM+1, 0);

			AXPY_LIBSTR(ng0, -1.0, tmp_ngM+1, 0, res_d_lg+ii, 0, res_d_lg+ii, 0);
			AXPY_LIBSTR(ng0, -1.0, tmp_ngM+1, 0, res_d_ug+ii, 0, res_d_ug+ii, 0);

			}


		if(ii<N)
			{

			nu1 = nu[ii+1];
			nx1 = nx[ii+1];

			AXPY_LIBSTR(nx1, -1.0, ux+(ii+1), nu1, b+ii, 0, res_b+ii, 0);

			GEMV_NT_LIBSTR(nu0+nx0, nx1, 1.0, 1.0, BAbt+ii, 0, 0, pi+ii, 0, ux+ii, 0, 1.0, 1.0, res_g+ii, 0, res_b+ii, 0, res_g+ii, 0, res_b+ii, 0);

			}

		}

	mu += VECMULDOT_LIBSTR(2*nbt+2*ngt, lam_lb, 0, t_lb, 0, ws->res_m, 0);

	if(cws->nb+cws->ng>0)
		ws->res_mu = mu*cws->nt_inv;
	else
		ws->res_mu = 0.0;

	return;

	}



// backward Riccati recursion
void FACT_SOLVE_KKT_UNCONSTR_OCP_QP(struct OCP_QP *qp, struct OCP_QP_SOL *qp_sol, struct IPM_HARD_OCP_QP_WORKSPACE *ws)
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
		GECP_LIBSTR(nu[N-ii-1]+nx[N-ii-1]+1, nx[N-ii], BAbt+(N-ii-1), 0, 0, AL, 0, 0);
		TRMM_RLNN_LIBSTR(nu[N-ii-1]+nx[N-ii-1]+1, nx[N-ii], 1.0, L+(N-ii), nu[N-ii], nu[N-ii], AL, 0, 0, AL, 0, 0);
		ROWEX_LIBSTR(nx[N-ii], 1.0, AL, nu[N-ii-1]+nx[N-ii-1], 0, tmp_nxM, 0);
		GEAD_LIBSTR(1, nx[N-ii], 1.0, L+(N-ii), nu[N-ii]+nx[N-ii], nu[N-ii], AL, nu[N-ii-1]+nx[N-ii-1], 0);

		SYRK_POTRF_LN_LIBSTR(nu[N-ii-1]+nx[N-ii-1]+1, nu[N-ii-1]+nx[N-ii-1], nx[N-ii], AL, 0, 0, AL, 0, 0, RSQrq+(N-ii-1), 0, 0, L+(N-ii-1), 0, 0);

//		d_print_strmat(nu[N-ii-1]+nx[N-ii-1]+1, nu[N-ii-1]+nx[N-ii-1], L+(N-ii-1), 0, 0);
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

//	d_print_tran_strvec(nu[ii]+nx[ii], ux+ii, 0);

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

//		d_print_tran_strvec(nu[ii]+nx[ii], ux+ii, 0);
		}

	return;

	}



// backward Riccati recursion
void FACT_SOLVE_KKT_STEP_HARD_OCP_QP(struct OCP_QP *qp, struct IPM_HARD_OCP_QP_WORKSPACE *ws)
	{

	int N = qp->N;
	int *nx = qp->nx;
	int *nu = qp->nu;
	int *nb = qp->nb;
	int *ng = qp->ng;

	struct STRMAT *BAbt = qp->BAbt;
	struct STRMAT *RSQrq = qp->RSQrq;
	struct STRMAT *DCt = qp->DCt;
	int **idxb = qp->idxb;

	struct STRMAT *L = ws->L;
	struct STRMAT *AL = ws->AL;
	struct STRVEC *res_b = ws->res_b;
	struct STRVEC *res_g = ws->res_g;
	struct STRVEC *dux = ws->dux;
	struct STRVEC *dpi = ws->dpi;
	struct STRVEC *dt_lb = ws->dt_lb;
	struct STRVEC *dt_lg = ws->dt_lg;
	struct STRVEC *Qx_lg = ws->Qx_lg;
	struct STRVEC *Qx_lb = ws->Qx_lb;
	struct STRVEC *qx_lg = ws->qx_lg;
	struct STRVEC *qx_lb = ws->qx_lb;
	struct STRVEC *Pb = ws->Pb;
	struct STRVEC *tmp_nxM = ws->tmp_nxM;

	//
	int ii;

	struct IPM_HARD_CORE_QP_WORKSPACE *cws = ws->core_workspace;

//	if(nb>0 | ng>0)
//		{
		COMPUTE_QX_QX_HARD_QP(cws);
//		}

	// factorization and backward substitution

	// last stage
#if defined(DOUBLE_PRECISION)
	TRCP_L_LIBSTR(nu[N]+nx[N], RSQrq+N, 0, 0, L+N, 0, 0); // TODO dtrcp_l_libstr with m and n, for m>=n
#else
	GECP_LIBSTR(nu[N]+nx[N], nu[N]+nx[N], RSQrq+N, 0, 0, L+N, 0, 0); // TODO dtrcp_l_libstr with m and n, for m>=n
#endif
	ROWIN_LIBSTR(nu[N]+nx[N], 1.0, res_g+N, 0, L+N, nu[N]+nx[N], 0);
	if(nb[N]>0)
		{
		DIAAD_SP_LIBSTR(nb[N], 1.0, Qx_lb+N, 0, idxb[N], L+N, 0, 0);
		ROWAD_SP_LIBSTR(nb[N], 1.0, qx_lb+N, 0, idxb[N], L+N, nu[N]+nx[N], 0);
		}
	if(ng[N]>0)
		{
		GEMM_R_DIAG_LIBSTR(nu[N]+nx[N], ng[N], 1.0, DCt+N, 0, 0, Qx_lg+N, 0, 0.0, AL+0, 0, 0, AL+0, 0, 0);
		ROWIN_LIBSTR(ng[N], 1.0, qx_lg+N, 0, AL+0, nu[N]+nx[N], 0);
		SYRK_POTRF_LN_LIBSTR(nu[N]+nx[N]+1, nu[N]+nx[N], ng[N], AL+0, 0, 0, DCt+N, 0, 0, L+N, 0, 0, L+N, 0, 0);
		}
	else
		{
		POTRF_L_MN_LIBSTR(nu[N]+nx[N]+1, nu[N]+nx[N], L+N, 0, 0, L+N, 0, 0);
		}

	// middle stages
	for(ii=0; ii<N; ii++)
		{
		GECP_LIBSTR(nu[N-ii-1]+nx[N-ii-1], nx[N-ii], BAbt+(N-ii-1), 0, 0, AL, 0, 0);
		ROWIN_LIBSTR(nx[N-ii], 1.0, res_b+(N-ii-1), 0, AL, nu[N-ii-1]+nx[N-ii-1], 0);
		TRMM_RLNN_LIBSTR(nu[N-ii-1]+nx[N-ii-1]+1, nx[N-ii], 1.0, L+(N-ii), nu[N-ii], nu[N-ii], AL, 0, 0, AL, 0, 0);
		ROWEX_LIBSTR(nx[N-ii], 1.0, AL, nu[N-ii-1]+nx[N-ii-1], 0, tmp_nxM, 0);
		TRMV_LNN_LIBSTR(nx[N-ii], nx[N-ii], L+(N-ii), nu[N-ii], nu[N-ii], tmp_nxM, 0, Pb+(N-ii-1), 0);
		GEAD_LIBSTR(1, nx[N-ii], 1.0, L+(N-ii), nu[N-ii]+nx[N-ii], nu[N-ii], AL, nu[N-ii-1]+nx[N-ii-1], 0);

#if defined(DOUBLE_PRECISION)
		TRCP_L_LIBSTR(nu[N-ii-1]+nx[N-ii-1], RSQrq+(N-ii-1), 0, 0, L+(N-ii-1), 0, 0);
#else
		GECP_LIBSTR(nu[N-ii-1]+nx[N-ii-1], nu[N-ii-1]+nx[N-ii-1], RSQrq+(N-ii-1), 0, 0, L+(N-ii-1), 0, 0);
#endif
		ROWIN_LIBSTR(nu[N-ii-1]+nx[N-ii-1], 1.0, res_g+(N-ii-1), 0, L+(N-ii-1), nu[N-ii-1]+nx[N-ii-1], 0);

		if(nb[N-ii-1]>0)
			{
			DIAAD_SP_LIBSTR(nb[N-ii-1], 1.0, Qx_lb+(N-ii-1), 0, idxb[N-ii-1], L+(N-ii-1), 0, 0);
			ROWAD_SP_LIBSTR(nb[N-ii-1], 1.0, qx_lb+(N-ii-1), 0, idxb[N-ii-1], L+(N-ii-1), nu[N-ii-1]+nx[N-ii-1], 0);
			}

		if(ng[N-ii-1]>0)
			{
			GEMM_R_DIAG_LIBSTR(nu[N-ii-1]+nx[N-ii-1], ng[N-ii-1], 1.0, DCt+N-ii-1, 0, 0, Qx_lg+N-ii-1, 0, 0.0, AL+0, 0, nx[N-ii], AL+0, 0, nx[N-ii]);
			ROWIN_LIBSTR(ng[N-ii-1], 1.0, qx_lg+N-ii-1, 0, AL+0, nu[N-ii-1]+nx[N-ii-1], nx[N-ii]);
			GECP_LIBSTR(nu[N-ii-1]+nx[N-ii-1], nx[N-ii], AL+0, 0, 0, AL+1, 0, 0);
			GECP_LIBSTR(nu[N-ii-1]+nx[N-ii-1], ng[N-ii-1], DCt+N-ii-1, 0, 0, AL+1, 0, nx[N-ii]);
			SYRK_POTRF_LN_LIBSTR(nu[N-ii-1]+nx[N-ii-1]+1, nu[N-ii-1]+nx[N-ii-1], nx[N-ii]+ng[N-ii-1], AL+0, 0, 0, AL+1, 0, 0, L+N-ii-1, 0, 0, L+N-ii-1, 0, 0);
			}
		else
			{
			SYRK_POTRF_LN_LIBSTR(nu[N-ii-1]+nx[N-ii-1]+1, nu[N-ii-1]+nx[N-ii-1], nx[N-ii], AL, 0, 0, AL, 0, 0, L+(N-ii-1), 0, 0, L+(N-ii-1), 0, 0);
			}

//		d_print_strmat(nu[N-ii-1]+nx[N-ii-1]+1, nu[N-ii-1]+nx[N-ii-1], L+(N-ii-1), 0, 0);
		}

	// forward substitution

	// first stage
	ii = 0;
	ROWEX_LIBSTR(nu[ii]+nx[ii], -1.0, L+(ii), nu[ii]+nx[ii], 0, dux+ii, 0);
	TRSV_LTN_LIBSTR(nu[ii]+nx[ii], L+ii, 0, 0, dux+ii, 0, dux+ii, 0);
	GEMV_T_LIBSTR(nu[ii]+nx[ii], nx[ii+1], 1.0, BAbt+ii, 0, 0, dux+ii, 0, 1.0, res_b+ii, 0, dux+(ii+1), nu[ii+1]);
	ROWEX_LIBSTR(nx[ii+1], 1.0, L+(ii+1), nu[ii+1]+nx[ii+1], nu[ii+1], tmp_nxM, 0);
	TRMV_LTN_LIBSTR(nx[ii+1], nx[ii+1], L+(ii+1), nu[ii+1], nu[ii+1], dux+(ii+1), nu[ii+1], dpi+ii, 0);
	AXPY_LIBSTR(nx[ii+1], 1.0, tmp_nxM, 0, dpi+ii, 0, dpi+ii, 0);
	TRMV_LNN_LIBSTR(nx[ii+1], nx[ii+1], L+(ii+1), nu[ii+1], nu[ii+1], dpi+ii, 0, dpi+ii, 0);

//	d_print_tran_strvec(nu[ii]+nx[ii], dux+ii, 0);

	// middle stages
	for(ii=1; ii<N; ii++)
		{
		ROWEX_LIBSTR(nu[ii], -1.0, L+(ii), nu[ii]+nx[ii], 0, dux+ii, 0);
		TRSV_LTN_MN_LIBSTR(nu[ii]+nx[ii], nu[ii], L+ii, 0, 0, dux+ii, 0, dux+ii, 0);
		GEMV_T_LIBSTR(nu[ii]+nx[ii], nx[ii+1], 1.0, BAbt+ii, 0, 0, dux+ii, 0, 1.0, res_b+ii, 0, dux+(ii+1), nu[ii+1]);
		ROWEX_LIBSTR(nx[ii+1], 1.0, L+(ii+1), nu[ii+1]+nx[ii+1], nu[ii+1], tmp_nxM, 0);
		TRMV_LTN_LIBSTR(nx[ii+1], nx[ii+1], L+(ii+1), nu[ii+1], nu[ii+1], dux+(ii+1), nu[ii+1], dpi+ii, 0);
		AXPY_LIBSTR(nx[ii+1], 1.0, tmp_nxM, 0, dpi+ii, 0, dpi+ii, 0);
		TRMV_LNN_LIBSTR(nx[ii+1], nx[ii+1], L+(ii+1), nu[ii+1], nu[ii+1], dpi+ii, 0, dpi+ii, 0);

//		d_print_tran_strvec(nu[ii]+nx[ii], dux+ii, 0);
		}



//	if(nb>0)
//		{
		for(ii=0; ii<=N; ii++)
			VECEX_SP_LIBSTR(nb[ii], 1.0, idxb[ii], dux+ii, 0, dt_lb+ii, 0);
//		}

//	if(ng>0)
//		{
		for(ii=0; ii<=N; ii++)
			GEMV_T_LIBSTR(nu[ii]+nx[ii], ng[ii], 1.0, DCt+ii, 0, 0, dux+ii, 0, 0.0, dt_lg+ii, 0, dt_lg+ii, 0);
//		}

//	if(nb>0 | ng>0)
//		{
		COMPUTE_LAM_T_HARD_QP(cws);
//		}

	return;

	}



// backward Riccati recursion
void SOLVE_KKT_STEP_HARD_OCP_QP(struct OCP_QP *qp, struct IPM_HARD_OCP_QP_WORKSPACE *ws)
	{

	int N = qp->N;
	int *nx = qp->nx;
	int *nu = qp->nu;
	int *nb = qp->nb;
	int *ng = qp->ng;

	struct STRMAT *BAbt = qp->BAbt;
	struct STRMAT *RSQrq = qp->RSQrq;
	struct STRMAT *DCt = qp->DCt;
	int **idxb = qp->idxb;

	struct STRMAT *L = ws->L;
	struct STRMAT *AL = ws->AL;
	struct STRVEC *res_b = ws->res_b;
	struct STRVEC *res_g = ws->res_g;
	struct STRVEC *dux = ws->dux;
	struct STRVEC *dpi = ws->dpi;
	struct STRVEC *dt_lb = ws->dt_lb;
	struct STRVEC *dt_lg = ws->dt_lg;
	struct STRVEC *qx_lg = ws->qx_lg;
	struct STRVEC *qx_lb = ws->qx_lb;
	struct STRVEC *Pb = ws->Pb;
	struct STRVEC *tmp_nxM = ws->tmp_nxM;

	//
	int ii;

	struct IPM_HARD_CORE_QP_WORKSPACE *cws = ws->core_workspace;

//	if(nb>0 | ng>0)
//		{
		COMPUTE_QX_HARD_QP(cws);
//		}

	// backward substitution

	// last stage
	VECCP_LIBSTR(nu[N]+nx[N], res_g+N, 0, dux+N, 0);
	if(nb[N]>0)
		{
		VECAD_SP_LIBSTR(nb[N], 1.0, qx_lb+N, 0, idxb[N], dux+N, 0);
		}
	// general constraints
	if(ng[N]>0)
		{
		GEMV_N_LIBSTR(nu[N]+nx[N], ng[N], 1.0, DCt+N, 0, 0, qx_lg+N, 0, 1.0, dux+N, 0, dux+N, 0);
		}

	// middle stages
	for(ii=0; ii<N-1; ii++)
		{
		VECCP_LIBSTR(nu[N-ii-1]+nx[N-ii-1], res_g+N-ii-1, 0, dux+N-ii-1, 0);
		if(nb[N-ii-1]>0)
			{
			VECAD_SP_LIBSTR(nb[N-ii-1], 1.0, qx_lb+N-ii-1, 0, idxb[N-ii-1], dux+N-ii-1, 0);
			}
		if(ng[N-ii-1]>0)
			{
			GEMV_N_LIBSTR(nu[N-ii-1]+nx[N-ii-1], ng[N-ii-1], 1.0, DCt+N-ii-1, 0, 0, qx_lg+N-ii-1, 0, 1.0, dux+N-ii-1, 0, dux+N-ii-1, 0);
			}
		AXPY_LIBSTR(nx[N-ii], 1.0, dux+N-ii, nu[N-ii], Pb+N-ii-1, 0, tmp_nxM, 0);
		GEMV_N_LIBSTR(nu[N-ii-1]+nx[N-ii-1], nx[N-ii], 1.0, BAbt+N-ii-1, 0, 0, tmp_nxM, 0, 1.0, dux+N-ii-1, 0, dux+N-ii-1, 0);
		TRSV_LNN_MN_LIBSTR(nu[N-ii-1]+nx[N-ii-1], nu[N-ii-1], L+N-ii-1, 0, 0, dux+N-ii-1, 0, dux+N-ii-1, 0);
		}

	// first stage
	ii = N-1;
	VECCP_LIBSTR(nu[N-ii-1]+nx[N-ii-1], res_g+N-ii-1, 0, dux+N-ii-1, 0);
	if(nb[N-ii-1]>0)
		{
		VECAD_SP_LIBSTR(nb[N-ii-1], 1.0, qx_lb+N-ii-1, 0, idxb[N-ii-1], dux+N-ii-1, 0);
		}
	if(ng[N-ii-1]>0)
		{
		GEMV_N_LIBSTR(nu[N-ii-1]+nx[N-ii-1], ng[N-ii-1], 1.0, DCt+N-ii-1, 0, 0, qx_lg+N-ii-1, 0, 1.0, dux+N-ii-1, 0, dux+N-ii-1, 0);
		}
	AXPY_LIBSTR(nx[N-ii], 1.0, dux+N-ii, nu[N-ii], Pb+N-ii-1, 0, tmp_nxM, 0);
	GEMV_N_LIBSTR(nu[N-ii-1]+nx[N-ii-1], nx[N-ii], 1.0, BAbt+N-ii-1, 0, 0, tmp_nxM, 0, 1.0, dux+N-ii-1, 0, dux+N-ii-1, 0);
	TRSV_LNN_LIBSTR(nu[N-ii-1]+nx[N-ii-1], L+N-ii-1, 0, 0, dux+N-ii-1, 0, dux+N-ii-1, 0);

	// first stage
	ii = 0;
	VECCP_LIBSTR(nx[ii+1], dux+ii+1, nu[ii+1], dpi+ii, 0);
	VECSC_LIBSTR(nu[ii]+nx[ii], -1.0, dux+ii, 0);
	TRSV_LTN_LIBSTR(nu[ii]+nx[ii], L+ii, 0, 0, dux+ii, 0, dux+ii, 0);
	GEMV_T_LIBSTR(nu[ii]+nx[ii], nx[ii+1], 1.0, BAbt+ii, 0, 0, dux+ii, 0, 1.0, res_b+ii, 0, dux+ii+1, nu[ii+1]);
	VECCP_LIBSTR(nx[ii+1], dux+ii+1, nu[ii+1], tmp_nxM, 0);
	TRMV_LTN_LIBSTR(nx[ii+1], nx[ii+1], L+ii+1, nu[ii+1], nu[ii+1], tmp_nxM, 0, tmp_nxM, 0);
	TRMV_LNN_LIBSTR(nx[ii+1], nx[ii+1], L+ii+1, nu[ii+1], nu[ii+1], tmp_nxM, 0, tmp_nxM, 0);
	AXPY_LIBSTR(nx[ii+1], 1.0, tmp_nxM, 0, dpi+ii, 0, dpi+ii, 0);

	// middle stages
	for(ii=1; ii<N; ii++)
		{
		VECCP_LIBSTR(nx[ii+1], dux+ii+1, nu[ii+1], dpi+ii, 0);
		VECSC_LIBSTR(nu[ii], -1.0, dux+ii, 0);
		TRSV_LTN_MN_LIBSTR(nu[ii]+nx[ii], nu[ii], L+ii, 0, 0, dux+ii, 0, dux+ii, 0);
		GEMV_T_LIBSTR(nu[ii]+nx[ii], nx[ii+1], 1.0, BAbt+ii, 0, 0, dux+ii, 0, 1.0, res_b+ii, 0, dux+ii+1, nu[ii+1]);
		VECCP_LIBSTR(nx[ii+1], dux+ii+1, nu[ii+1], tmp_nxM, 0);
		TRMV_LTN_LIBSTR(nx[ii+1], nx[ii+1], L+ii+1, nu[ii+1], nu[ii+1], tmp_nxM, 0, tmp_nxM, 0);
		TRMV_LNN_LIBSTR(nx[ii+1], nx[ii+1], L+ii+1, nu[ii+1], nu[ii+1], tmp_nxM, 0, tmp_nxM, 0);
		AXPY_LIBSTR(nx[ii+1], 1.0, tmp_nxM, 0, dpi+ii, 0, dpi+ii, 0);
		}



//	if(nb>0)
//		{
		for(ii=0; ii<=N; ii++)
			VECEX_SP_LIBSTR(nb[ii], 1.0, idxb[ii], dux+ii, 0, dt_lb+ii, 0);
//		}

//	if(ng>0)
//		{
		for(ii=0; ii<=N; ii++)
			GEMV_T_LIBSTR(nu[ii]+nx[ii], ng[ii], 1.0, DCt+ii, 0, 0, dux+ii, 0, 0.0, dt_lg+ii, 0, dt_lg+ii, 0);
//		}

//	if(nb>0 | ng>0)
//		{
		COMPUTE_LAM_T_HARD_QP(cws);
//		}

	return;

	}


