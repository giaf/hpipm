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



void INIT_VAR_HARD_OCP_QP(struct OCP_QP *qp, struct OCP_QP_SOL *qp_sol, struct IPM_HARD_OCP_QP_WORKSPACE *ws)
	{

	struct IPM_CORE_QP_WORKSPACE *cws = ws->core_workspace;
	
	// loop index
	int ii, jj;

	//
	int N = qp->N;
	int *nx = qp->nx;
	int *nu = qp->nu;
	int *nb = qp->nb;
	int *ng = qp->ng;
	int *ns = qp->ns;

	REAL mu0 = cws->mu0;

	//
	REAL *ux, *s_lb, *s_ub, *pi, *d_lb, *d_ub, *d_lg, *d_ug, *lam_lb, *lam_ub, *lam_lg, *lam_ug, *t_lb, *t_ub, *t_lg, *t_ug;
	int *idxb;

	REAL thr0 = 0.1;

	// warm start TODO

	// cold start

	// ux
	for(ii=0; ii<=N; ii++)
		{
		ux = qp_sol->ux[ii].pa;
		for(jj=0; jj<nu[ii]+nx[ii]+2*ns[ii]; jj++)
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
		d_lb = qp->d[ii].pa;
		d_ub = qp->d[ii].pa+nb[ii]+ng[ii];
		lam_lb = qp_sol->lam[ii].pa;
		lam_ub = qp_sol->lam[ii].pa+nb[ii]+ng[ii];
		t_lb = qp_sol->t[ii].pa;
		t_ub = qp_sol->t[ii].pa+nb[ii]+ng[ii];
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
		t_lg = qp_sol->t[ii].pa+nb[ii];
		t_ug = qp_sol->t[ii].pa+2*nb[ii]+ng[ii];
		lam_lg = qp_sol->lam[ii].pa+nb[ii];
		lam_ug = qp_sol->lam[ii].pa+2*nb[ii]+ng[ii];
		d_lg = qp->d[ii].pa+nb[ii];
		d_ug = qp->d[ii].pa+2*nb[ii]+ng[ii];
		ux = qp_sol->ux[ii].pa;
//		GEMV_T_LIBSTR(nu[ii]+nx[ii], ng[ii], 1.0, qp->DCt+ii, 0, 0, qp_sol->ux+ii, 0, 0.0, qp_sol->t+ii, nb[ii], qp_sol->t+ii, nb[ii]);
		GEMV_T_LIBSTR(nu[ii]+nx[ii], ng[ii], 0.0, qp->DCt+ii, 0, 0, qp_sol->ux+ii, 0, 0.0, qp_sol->t+ii, nb[ii], qp_sol->t+ii, nb[ii]);
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

	// soft constraints
	for(ii=0; ii<=N; ii++)
		{
		s_lb = qp_sol->ux[ii].pa+nu[ii]+nx[ii];
		s_ub = qp_sol->ux[ii].pa+nu[ii]+nx[ii]+ns[ii];
		lam_lb = qp_sol->lam[ii].pa+2*nb[ii]+2*ng[ii];
		lam_ub = qp_sol->lam[ii].pa+2*nb[ii]+2*ng[ii]+ns[ii];
		t_lb = qp_sol->t[ii].pa+2*nb[ii]+2*ng[ii];
		t_ub = qp_sol->t[ii].pa+2*nb[ii]+2*ng[ii]+ns[ii];
//		idxs = qp->idxs[ii];
		for(jj=0; jj<ns[ii]; jj++)
			{
			s_lb[jj] = thr0;
			s_ub[jj] = thr0;
			t_lb[jj] = thr0;//1.0;
			t_ub[jj] = thr0;//1.0;
			lam_lb[jj] = mu0/t_lb[jj];
			lam_ub[jj] = mu0/t_ub[jj];
			}
		}

	return;

	}



void COMPUTE_RES_HARD_OCP_QP(struct OCP_QP *qp, struct OCP_QP_SOL *qp_sol, struct IPM_HARD_OCP_QP_WORKSPACE *ws)
	{

	struct IPM_CORE_QP_WORKSPACE *cws = ws->core_workspace;
	
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
	struct STRVEC *tmp_ngM = ws->tmp_ngM;
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

		VECCP_LIBSTR(nu0+nx0, rq+ii, 0, res_g+ii, 0);

		if(ii>0)
			AXPY_LIBSTR(nx0, -1.0, pi+(ii-1), 0, res_g+ii, nu0, res_g+ii, nu0);

		SYMV_L_LIBSTR(nu0+nx0, nu0+nx0, 1.0, RSQrq+ii, 0, 0, ux+ii, 0, 1.0, res_g+ii, 0, res_g+ii, 0);

		if(nb0>0)
			{

			AXPY_LIBSTR(nb0, -1.0, lam+ii, 0, lam+ii, nb[ii]+ng[ii], tmp_nbgM, 0);
			VECAD_SP_LIBSTR(nb0, 1.0, tmp_nbgM, 0, idxb[ii], res_g+ii, 0);

			VECEX_SP_LIBSTR(nb0, -1.0, idxb[ii], ux+ii, 0, res_d+ii, 0);
			VECCP_LIBSTR(nb0, res_d+ii, 0, res_d+ii, nb0+ng0);
			AXPY_LIBSTR(nb0, 1.0, d+ii, 0, res_d+ii, 0, res_d+ii, 0);
			AXPY_LIBSTR(nb0, 1.0, d+ii, nb[ii]+ng[ii], res_d+ii, nb0+ng0, res_d+ii, nb0+ng0);
			AXPY_LIBSTR(nb0, 1.0, t+ii, 0, res_d+ii, 0, res_d+ii, 0);
			AXPY_LIBSTR(nb0, -1.0, t+ii, nb[ii]+ng[ii], res_d+ii, nb0+ng0, res_d+ii, nb0+ng0);
			VECSC_LIBSTR(nb0, -1.0, res_d+ii, nb0+ng0); // TODO embed with above

			}

		if(ng0>0) // TODO merge with bounds as much as possible
			{

//			d_print_e_tran_strvec(ng0, lam+ii, nb0);
//			d_print_e_tran_strvec(ng0, lam+ii, 2*nb0+ng0);
			AXPY_LIBSTR(ng0, -1.0, lam+ii, nb[ii], lam+ii, 2*nb[ii]+ng[ii], tmp_nbgM, nb[ii]);
//			d_print_e_tran_strvec(ng0, tmp_nbgM, nb0);

			AXPY_LIBSTR(ng0,  1.0, t+ii, nb[ii], d+ii, nb[ii], res_d+ii, nb0);
			AXPY_LIBSTR(ng0, -1.0, t+ii, 2*nb[ii]+ng[ii], d+ii, 2*nb[ii]+ng[ii], res_d+ii, 2*nb0+ng0);

//			d_print_strmat(nu0+nx0, ng0, DCt+ii, 0, 0);
			GEMV_NT_LIBSTR(nu0+nx0, ng0, 1.0, 1.0, DCt+ii, 0, 0, tmp_nbgM, nb[ii], ux+ii, 0, 1.0, 0.0, res_g+ii, 0, tmp_ngM, 0, res_g+ii, 0, tmp_ngM, 0);
//			d_print_e_tran_strvec(nu0+nx0, res_g+ii, 0);

			AXPY_LIBSTR(ng0, -1.0, tmp_ngM, 0, res_d+ii, nb0, res_d+ii, nb0);
			AXPY_LIBSTR(ng0, -1.0, tmp_ngM, 0, res_d+ii, 2*nb0+ng0, res_d+ii, 2*nb0+ng0);
			VECSC_LIBSTR(ng0, -1.0, res_d+ii, 2*nb0+ng0); // TODO embed with above

			}

		if(ns0>0)
			{

			// res_g
//			printf("\nii = %d\n", ii);
//			d_print_e_tran_strvec(2*ns0, Z+ii, 0);
//			d_print_e_tran_strvec(2*ns0, ux+ii, nu0+nx0);
//			d_print_e_tran_strvec(2*ns0, z+ii, 0);
			GEMV_DIAG_LIBSTR(2*ns0, 1.0, Z+ii, 0, ux+ii, nu0+nx0, 1.0, z+ii, 0, res_g+ii, nu0+nx0);
//			d_print_e_tran_strvec(2*ns0, res_g+ii, nu0+nx0);

			AXPY_LIBSTR(2*ns0, -1.0, lam+ii, 2*nb0+2*ng0, res_g+ii, nu0+nx0, res_g+ii, nu0+nx0);
//			d_print_e_tran_strvec(2*ns0, res_g+ii, nu0+nx0);

			VECEX_SP_LIBSTR(ns0, 1.0, idxs[ii], lam+ii, 0, tmp_nsM, 0);
			AXPY_LIBSTR(ns0, -1.0, tmp_nsM, 0, res_g+ii, nu0+nx0, res_g+ii, nu0+nx0);
//			d_print_e_tran_strvec(2*ns0, res_g+ii, nu0+nx0);

			VECEX_SP_LIBSTR(ns0, 1.0, idxs[ii], lam+ii, nb0+ng0, tmp_nsM, 0);
			AXPY_LIBSTR(ns0, -1.0, tmp_nsM, 0, res_g+ii, nu0+nx0+ns0, res_g+ii, nu0+nx0+ns0);
//			d_print_e_tran_strvec(2*ns0, res_g+ii, nu0+nx0);

			// res_d
//			VECAD_SP_LIBSTR(ns0, 1.0, t+ii, 2*nb0+2*ng0, idxs[ii], res_d+ii, 0);
//			VECAD_SP_LIBSTR(ns0, 1.0, t+ii, 2*nb0+2*ng0+ns0, idxs[ii], res_d+ii, nb0+ng0);
			VECAD_SP_LIBSTR(ns0, -1.0, ux+ii, nu0+nx0, idxs[ii], res_d+ii, 0);
			VECAD_SP_LIBSTR(ns0, -1.0, ux+ii, nu0+nx0+ns0, idxs[ii], res_d+ii, nb0+ng0);
//			VECSC_LIBSTR(2*ns0, 0.0, res_d+ii, 2*nb0+2*ng0); // TODO vecse ???
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

	mu += VECMULDOT_LIBSTR(2*nct, lam, 0, t, 0, ws->res_m, 0); // XXX

	if(nct>0)
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
		TRMM_RLNN_LIBSTR(nu[N-ii-1]+nx[N-ii-1]+1, nx[N-ii], 1.0, L+(N-ii), nu[N-ii], nu[N-ii], BAbt+(N-ii-1), 0, 0, AL, 0, 0);
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
	
	ii = N;
	ROWEX_LIBSTR(nu[ii], -1.0, L+(ii), nu[ii]+nx[ii], 0, ux+ii, 0);
	TRSV_LTN_MN_LIBSTR(nu[ii]+nx[ii], nu[ii], L+ii, 0, 0, ux+ii, 0, ux+ii, 0);

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
	struct STRVEC *tmp_ngM = ws->tmp_ngM;
	struct STRVEC *tmp_nsM = ws->tmp_nsM;

	REAL *ptr0, *ptr1, *ptr2, *ptr3;

	//
	int ii, nn, ss, idx;

	struct IPM_CORE_QP_WORKSPACE *cws = ws->core_workspace;

	COMPUTE_QX_QX_QP(cws);

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
		// Gamma lower
		VECEX_SP_LIBSTR(ns[ss], 1.0, idxs[ss], Gamma+ss, 0, Zs_inv+ss, 0);
		AXPY_LIBSTR(ns[ss], 1.0, Z+ss, 0, Zs_inv+ss, 0, Zs_inv+ss, 0);
		AXPY_LIBSTR(ns[ss], 1.0, Gamma+ss, 2*nb[ss]+2*ng[ss], Zs_inv+ss, 0, Zs_inv+ss, 0);
		ptr0 = (Zs_inv+ss)->pa;
		for(ii=0; ii<ns[ss]; ii++)
			ptr0[ii] = 1.0/ptr0[ii];
		VECCP_LIBSTR(nb[ss]+ng[ss], Gamma+ss, 0, tmp_nbgM+0, 0);
		ptr0 = (Zs_inv+ss)->pa;
		ptr1 = (tmp_nbgM+0)->pa;
		for(ii=0; ii<ns[ss]; ii++)
			{
			idx = idxs[ss][ii];
			ptr1[idx] = ptr1[idx] - ptr1[idx]*ptr0[ii]*ptr1[idx];
			}
		// Gamma upper
		VECEX_SP_LIBSTR(ns[ss], 1.0, idxs[ss], Gamma+ss, nb[ss]+ng[ss], Zs_inv+ss, ns[ss]);
		AXPY_LIBSTR(ns[ss], 1.0, Z+ss, ns[ss], Zs_inv+ss, ns[ss], Zs_inv+ss, ns[ss]);
		AXPY_LIBSTR(ns[ss], 1.0, Gamma+ss, 2*nb[ss]+2*ng[ss]+ns[ss], Zs_inv+ss, ns[ss], Zs_inv+ss, ns[ss]);
		ptr0 = (Zs_inv+ss)->pa+ns[ss];
		for(ii=0; ii<ns[ss]; ii++)
			ptr0[ii] = 1.0/ptr0[ii];
		VECCP_LIBSTR(nb[ss]+ng[ss], Gamma+ss, nb[ss]+ng[ss], tmp_nbgM+2, 0);
		ptr0 = (Zs_inv+ss)->pa+ns[ss];
		ptr1 = (tmp_nbgM+2)->pa;
		for(ii=0; ii<ns[ss]; ii++)
			{
			idx = idxs[ss][ii];
			ptr1[idx] = ptr1[idx] - ptr1[idx]*ptr0[ii]*ptr1[idx];
			}
		// gamma lower
//d_print_e_tran_strvec(ns[ss], gamma+ss, 0);
		VECEX_SP_LIBSTR(ns[ss], 1.0, idxs[ss], gamma+ss, 0, dux+ss, nu[ss]+nx[ss]);
//d_print_e_tran_strvec(ns[ss], res_g+ss, nu[ss]+nx[ss]);
		AXPY_LIBSTR(ns[ss], 1.0, res_g+ss, nu[ss]+nx[ss], dux+ss, nu[ss]+nx[ss], dux+ss, nu[ss]+nx[ss]);
//d_print_e_tran_strvec(ns[ss], gamma+ss, 2*nb[ss]+2*ng[ss]);
		AXPY_LIBSTR(ns[ss], 1.0, gamma+ss, 2*nb[ss]+2*ng[ss], dux+ss, nu[ss]+nx[ss], dux+ss, nu[ss]+nx[ss]);
//d_print_e_tran_strvec(ns[ss], dux+ss, nu[ss]+nx[ss]);
		ptr0 = (Zs_inv+ss)->pa;
		ptr1 = (dux+ss)->pa+nu[ss]+nx[ss];
		ptr2 = (tmp_nsM)->pa;
		for(ii=0; ii<ns[ss]; ii++)
			ptr2[ii] = ptr0[ii]*ptr1[ii];
		VECCP_LIBSTR(nb[ss]+ng[ss], gamma+ss, 0, tmp_nbgM+1, 0);
		ptr0 = (tmp_nsM)->pa;
		ptr1 = (Gamma+ss)->pa;
		ptr2 = (tmp_nbgM+1)->pa;
		for(ii=0; ii<ns[ss]; ii++)
			{
			idx = idxs[ss][ii];
			ptr2[idx] = ptr2[idx] - ptr1[idx]*ptr0[ii];
			}
		// gamma upper
		VECEX_SP_LIBSTR(ns[ss], 1.0, idxs[ss], gamma+ss, nb[ss]+ng[ss], dux+ss, nu[ss]+nx[ss]+ns[ss]);
		AXPY_LIBSTR(ns[ss], 1.0, res_g+ss, nu[ss]+nx[ss]+ns[ss], dux+ss, nu[ss]+nx[ss]+ns[ss], dux+ss, nu[ss]+nx[ss]+ns[ss]);
		AXPY_LIBSTR(ns[ss], 1.0, gamma+ss, 2*nb[ss]+2*ng[ss]+ns[ss], dux+ss, nu[ss]+nx[ss]+ns[ss], dux+ss, nu[ss]+nx[ss]+ns[ss]);
		ptr0 = (Zs_inv+ss)->pa+ns[ss];
		ptr1 = (dux+ss)->pa+nu[ss]+nx[ss]+ns[ss];
		ptr2 = (tmp_nsM)->pa;
		for(ii=0; ii<ns[ss]; ii++)
			ptr2[ii] = ptr0[ii]*ptr1[ii];
		VECCP_LIBSTR(nb[ss]+ng[ss], gamma+ss, nb[ss]+ng[ss], tmp_nbgM+3, 0);
		ptr0 = (tmp_nsM)->pa;
		ptr1 = (Gamma+ss)->pa+nb[ss]+ng[ss];
		ptr2 = (tmp_nbgM+3)->pa;
		for(ii=0; ii<ns[ss]; ii++)
			{
			idx = idxs[ss][ii];
			ptr2[idx] = ptr2[idx] - ptr1[idx]*ptr0[ii];
			}
//		d_print_e_tran_strvec(nb[ss]+ng[ss], tmp_nbgM+0, 0);
//		d_print_e_tran_strvec(nb[ss]+ng[ss], tmp_nbgM+2, 0);
//		d_print_e_tran_strvec(nb[ss]+ng[ss], tmp_nbgM+1, 0);
//		d_print_e_tran_strvec(nb[ss]+ng[ss], tmp_nbgM+3, 0);
//		d_print_e_tran_strvec(2*ns[ss], Zs_inv+ss, 0);
//		d_print_e_tran_strvec(2*nb[ss]+2*ng[ss], Gamma+ss, 0);
//		d_print_e_tran_strvec(2*ns[ss], dux+ss, nu[ss]+nx[ss]);
//		exit(1);
		// Gamma update
		AXPY_LIBSTR(nb[ss]+ng[ss],  1.0, tmp_nbgM+2, 0, tmp_nbgM+0, 0, tmp_nbgM+0, 0);
		AXPY_LIBSTR(nb[ss]+ng[ss], -1.0, tmp_nbgM+3, 0, tmp_nbgM+1, 0, tmp_nbgM+1, 0);
		}
	else if(nb[ss]>0 | ng[ss]>0)
		{
//d_print_e_tran_strvec(2*nb[ss]+2*ng[ss], gamma+ss, 0);
		AXPY_LIBSTR(nb[ss]+ng[ss],  1.0, Gamma+ss, nb[ss]+ng[ss], Gamma+ss, 0, tmp_nbgM+0, 0);
		AXPY_LIBSTR(nb[ss]+ng[ss], -1.0, gamma+ss, nb[ss]+ng[ss], gamma+ss, 0, tmp_nbgM+1, 0);
		}
//d_print_strmat(nu[N]+nx[N]+1, nu[N]+nx[N], L+N, 0, 0);
	if(nb[ss]>0)
		{
//		AXPY_LIBSTR(nb[ss], 1.0, Gamma+ss, nb[ss]+ng[ss], Gamma+ss, 0, tmp_nbgM+0, 0);
		DIAAD_SP_LIBSTR(nb[ss], 1.0, tmp_nbgM+0, 0, idxb[ss], L+ss, 0, 0);
//		AXPY_LIBSTR(nb[ss], -1.0, gamma+ss, nb[ss]+ng[ss], gamma+ss, 0, tmp_nbgM+1, 0);
		ROWAD_SP_LIBSTR(nb[ss], 1.0, tmp_nbgM+1, 0, idxb[ss], L+ss, nu[ss]+nx[ss], 0);
		}
	if(ng[ss]>0)
		{
//		AXPY_LIBSTR(ng[ss], 1.0, Gamma+ss, 2*nb[ss]+ng[ss], Gamma+ss, nb[ss], tmp_nbgM+0, nb[ss]);
		GEMM_R_DIAG_LIBSTR(nu[ss]+nx[ss], ng[ss], 1.0, DCt+ss, 0, 0, tmp_nbgM+0, nb[ss], 0.0, AL+0, 0, 0, AL+0, 0, 0);
//		AXPY_LIBSTR(ng[ss], -1.0, gamma+ss, 2*nb[ss]+ng[ss], gamma+ss, nb[ss], tmp_nbgM+1, nb[ss]);
		ROWIN_LIBSTR(ng[ss], 1.0, tmp_nbgM+1, nb[ss], AL+0, nu[ss]+nx[ss], 0);
		SYRK_POTRF_LN_LIBSTR(nu[ss]+nx[ss]+1, nu[ss]+nx[ss], ng[ss], AL+0, 0, 0, DCt+ss, 0, 0, L+ss, 0, 0, L+ss, 0, 0);
//		dsyrk_ln_mn_libstr(nu[ss]+nx[ss]+1, nu[ss]+nx[ss], ng[ss], 1.0, AL+0, 0, 0, DCt+ss, 0, 0, 1.0, L+ss, 0, 0, L+ss, 0, 0);
//d_print_strmat(nu[N]+nx[N]+1, nu[N]+nx[N], L+N, 0, 0);
//		dpotrf_l_mn_libstr(nu[ss]+nx[ss]+1, nu[ss]+nx[ss], L+ss, 0, 0, L+ss, 0, 0);
		}
	else
		{
		POTRF_L_MN_LIBSTR(nu[ss]+nx[ss]+1, nu[ss]+nx[ss], L+ss, 0, 0, L+ss, 0, 0);
		}
//d_print_strmat(nu[N]+nx[N]+1, nu[N]+nx[N], L+N, 0, 0);
//exit(1);
	
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
		TRCP_L_LIBSTR(nu[ss]+nx[ss], RSQrq+(ss), 0, 0, L+(ss), 0, 0);
#else
		GECP_LIBSTR(nu[ss]+nx[ss], nu[ss]+nx[ss], RSQrq+(ss), 0, 0, L+(ss), 0, 0);
#endif
		ROWIN_LIBSTR(nu[ss]+nx[ss], 1.0, res_g+(ss), 0, L+(ss), nu[ss]+nx[ss], 0);

		if(ns[ss]>0)
			{
			}
		AXPY_LIBSTR(nb[ss]+ng[ss],  1.0, Gamma+ss, nb[ss]+ng[ss], Gamma+ss, 0, tmp_nbgM+0, 0);
		AXPY_LIBSTR(nb[ss]+ng[ss], -1.0, gamma+ss, nb[ss]+ng[ss], gamma+ss, 0, tmp_nbgM+1, 0);
		if(nb[ss]>0)
			{
//			AXPY_LIBSTR(nb[ss], 1.0, Gamma+ss, nb[ss]+ng[ss], Gamma+ss, 0, tmp_nbgM+0, 0);
			DIAAD_SP_LIBSTR(nb[ss], 1.0, tmp_nbgM+0, 0, idxb[ss], L+(ss), 0, 0);
//			AXPY_LIBSTR(nb[ss], -1.0, gamma+ss, nb[ss]+ng[ss], gamma+ss, 0, tmp_nbgM+1, 0);
			ROWAD_SP_LIBSTR(nb[ss], 1.0, tmp_nbgM+1, 0, idxb[ss], L+(ss), nu[ss]+nx[ss], 0);
			}

		if(ng[ss]>0)
			{
//			AXPY_LIBSTR(ng[ss], 1.0, Gamma+ss, 2*nb[ss]+nb[ss], Gamma+ss, nb[ss], tmp_nbgM+0, nb[ss]);
			GEMM_R_DIAG_LIBSTR(nu[ss]+nx[ss], ng[ss], 1.0, DCt+ss, 0, 0, tmp_nbgM+0, nb[ss], 0.0, AL+0, 0, nx[ss+1], AL+0, 0, nx[ss+1]);
//			AXPY_LIBSTR(ng[ss], -1.0, gamma+ss, 2*nb[ss]+ng[ss], gamma+ss, nb[ss], tmp_nbgM+1, nb[ss]);
			ROWIN_LIBSTR(ng[ss], 1.0, tmp_nbgM+1, nb[ss], AL+0, nu[ss]+nx[ss], nx[ss+1]);
			GECP_LIBSTR(nu[ss]+nx[ss], nx[ss+1], AL+0, 0, 0, AL+1, 0, 0);
			GECP_LIBSTR(nu[ss]+nx[ss], ng[ss], DCt+ss, 0, 0, AL+1, 0, nx[ss+1]);
			SYRK_POTRF_LN_LIBSTR(nu[ss]+nx[ss]+1, nu[ss]+nx[ss], nx[ss+1]+ng[ss], AL+0, 0, 0, AL+1, 0, 0, L+ss, 0, 0, L+ss, 0, 0);
			}
		else
			{
			SYRK_POTRF_LN_LIBSTR(nu[ss]+nx[ss]+1, nu[ss]+nx[ss], nx[ss+1], AL, 0, 0, AL, 0, 0, L+(ss), 0, 0, L+(ss), 0, 0);
			}

//		d_print_strmat(nu[ss]+nx[ss]+1, nu[ss]+nx[ss], L+(ss), 0, 0);
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

//	d_print_tran_strvec(nu[ss]+nx[ss], dux+ss, 0);

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

//		d_print_tran_strvec(nu[ss]+nx[ss], dux+ss, 0);
		}

	ss = N;
	ROWEX_LIBSTR(nu[ss], -1.0, L+ss, nu[ss]+nx[ss], 0, dux+ss, 0);
//d_print_tran_strvec(nu[N]+nx[N], dux+N, 0);
	TRSV_LTN_MN_LIBSTR(nu[ss]+nx[ss], nu[ss], L+ss, 0, 0, dux+ss, 0, dux+ss, 0);

//d_print_strmat(nu[N]+nx[N]+1, nu[N]+nx[N], L+N, 0, 0);
//d_print_tran_strvec(nu[N]+nx[N], dux+N, 0);
//exit(1);


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
		VECEX_SP_LIBSTR(ns[ss], 1.0, idxs[ss], dt+ss, 0, tmp_nsM+0, 0); // XXX ???
		VECEX_SP_LIBSTR(ns[ss], 1.0, idxs[ss], Gamma+ss, 0, tmp_nsM+1, 0);
		ptr0 = (tmp_nsM+0)->pa;
		ptr1 = (tmp_nsM+1)->pa;
		for(ii=0; ii<ns[ss]; ii++)
			ptr1[ii] = ptr1[ii]*ptr0[ii];
		AXPY_LIBSTR(ns[ss], 1.0, tmp_nsM+1, 0, dux+ss, nu[ss]+nx[ss], dux+ss, nu[ss]+nx[ss]);

		VECEX_SP_LIBSTR(ns[ss], 1.0, idxs[ss], dt+ss, 0, tmp_nsM+0, 0); // XXX ???
		VECEX_SP_LIBSTR(ns[ss], 1.0, idxs[ss], Gamma+ss, nb[ss]+ng[ss], tmp_nsM+1, 0);
		ptr0 = (tmp_nsM+0)->pa;
		ptr1 = (tmp_nsM+1)->pa;
		for(ii=0; ii<ns[ss]; ii++)
			ptr1[ii] = ptr1[ii]*ptr0[ii];
		AXPY_LIBSTR(ns[ss], -1.0, tmp_nsM+1, 0, dux+ss, nu[ss]+nx[ss]+ns[ss], dux+ss, nu[ss]+nx[ss]+ns[ss]);

		ptr0 = (Zs_inv+ss)->pa;
		ptr1 = (dux+ss)->pa+nu[ss]+nx[ss];
		for(ii=0; ii<2*ns[ss]; ii++)
			ptr1[ii] = - ptr1[ii]*ptr0[ii];

		VECAD_SP_LIBSTR(ns[ss], 1.0, dux+ss, nu[ss]+nx[ss], idxs[ss], dt+ss, 0);
		VECAD_SP_LIBSTR(ns[ss], 1.0, dux+ss, nu[ss]+nx[ss]+ns[ss], idxs[ss], dt+ss, nb[ss]+ng[ss]);
		}
	
	for(ss=0; ss<=N; ss++)
		{
		VECCP_LIBSTR(2*ns[ss], dux+ss, nu[ss]+nx[ss], dt+ss, 2*nb[ss]+2*ng[ss]);
//		VECEX_SP_LIBSTR(ns[ss], 1.0, idxs[ss], dt+ss, 0, tmp_nsM, 0);
//		AXPY_LIBSTR(ns[ss],  1.0, tmp_nsM, 0, dt+ss, 2*nb[ss]+2*ng[ss], dt+ss, 2*nb[ss]+2*ng[ss]);
//		AXPY_LIBSTR(ns[ss], -1.0, tmp_nsM, 0, dt+ss, 2*nb[ss]+2*ng[ss]+ns[ss], dt+ss, 2*nb[ss]+2*ng[ss]+ns[ss]);
//		ptr0 = (Zs_inv+ss)->pa;
//		ptr1 = (dt+ss)->pa+2*nb[ss]+2*ng[ss];
//		for(ii=0; ii<2*ns[ss]; ii++)
//			ptr1[ii] = - ptr1[ii]*ptr0[ii];
		}

	COMPUTE_LAM_T_QP(cws);

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
	struct STRVEC *dt = ws->dt;
	struct STRVEC *gamma = ws->gamma;
	struct STRVEC *Pb = ws->Pb;
	struct STRVEC *tmp_nxM = ws->tmp_nxM;
	struct STRVEC *tmp_nbgM = ws->tmp_nbgM;
	struct STRVEC *tmp_ngM = ws->tmp_ngM;

	//
	int ii;

	struct IPM_CORE_QP_WORKSPACE *cws = ws->core_workspace;

//	if(nb>0 | ng>0)
//		{
		COMPUTE_QX_QP(cws);
//		}

	// backward substitution

	// last stage
	VECCP_LIBSTR(nu[N]+nx[N], res_g+N, 0, dux+N, 0);
	if(nb[N]>0)
		{
		AXPY_LIBSTR(nb[N], -1.0, gamma+N, nb[N]+ng[N], gamma+N, 0, tmp_nbgM, 0);
		VECAD_SP_LIBSTR(nb[N], 1.0, tmp_nbgM, 0, idxb[N], dux+N, 0);
		}
	// general constraints
	if(ng[N]>0)
		{
		AXPY_LIBSTR(ng[N], -1.0, gamma+N, 2*nb[N]+ng[N], gamma+N, nb[N], tmp_nbgM, nb[N]);
		GEMV_N_LIBSTR(nu[N]+nx[N], ng[N], 1.0, DCt+N, 0, 0, tmp_nbgM, nb[N], 1.0, dux+N, 0, dux+N, 0);
		}
	TRSV_LNN_MN_LIBSTR(nu[N]+nx[N], nu[N], L+N, 0, 0, dux+N, 0, dux+N, 0);

	// middle stages
	for(ii=0; ii<N-1; ii++)
		{
		VECCP_LIBSTR(nu[N-ii-1]+nx[N-ii-1], res_g+N-ii-1, 0, dux+N-ii-1, 0);
		if(nb[N-ii-1]>0)
			{
			AXPY_LIBSTR(nb[N-ii-1], -1.0, gamma+N-ii-1, nb[N-ii-1]+ng[N-ii-1], gamma+N-ii-1, 0, tmp_nbgM, 0);
			VECAD_SP_LIBSTR(nb[N-ii-1], 1.0, tmp_nbgM, 0, idxb[N-ii-1], dux+N-ii-1, 0);
			}
		if(ng[N-ii-1]>0)
			{
			AXPY_LIBSTR(ng[N-ii-1], -1.0, gamma+N-ii-1, 2*nb[N-ii-1]+ng[N-ii-1], gamma+N-ii-1, nb[N-ii-1], tmp_nbgM, nb[N]);
			GEMV_N_LIBSTR(nu[N-ii-1]+nx[N-ii-1], ng[N-ii-1], 1.0, DCt+N-ii-1, 0, 0, tmp_nbgM, nb[N], 1.0, dux+N-ii-1, 0, dux+N-ii-1, 0);
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
		AXPY_LIBSTR(nb[N-ii-1], -1.0, gamma+N-ii-1, nb[N-ii-1]+ng[N-ii-1], gamma+N-ii-1, 0, tmp_nbgM, 0);
		VECAD_SP_LIBSTR(nb[N-ii-1], 1.0, tmp_nbgM, 0, idxb[N-ii-1], dux+N-ii-1, 0);
		}
	if(ng[N-ii-1]>0)
		{
		AXPY_LIBSTR(ng[N-ii-1], -1.0, gamma+N-ii-1, 2*nb[N-ii-1]+ng[N-ii-1], gamma+N-ii-1, nb[N-ii-1], tmp_nbgM, nb[N-ii-1]);
		GEMV_N_LIBSTR(nu[N-ii-1]+nx[N-ii-1], ng[N-ii-1], 1.0, DCt+N-ii-1, 0, 0, tmp_nbgM, nb[N-ii-1], 1.0, dux+N-ii-1, 0, dux+N-ii-1, 0);
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

	ii = N;
	VECSC_LIBSTR(nu[ii], -1.0, dux+ii, 0);
	TRSV_LTN_MN_LIBSTR(nu[ii]+nx[ii], nu[ii], L+ii, 0, 0, dux+ii, 0, dux+ii, 0);



//	if(nb>0)
//		{
		for(ii=0; ii<=N; ii++)
			VECEX_SP_LIBSTR(nb[ii], 1.0, idxb[ii], dux+ii, 0, dt+ii, 0);
//		}

//	if(ng>0)
//		{
		for(ii=0; ii<=N; ii++)
			GEMV_T_LIBSTR(nu[ii]+nx[ii], ng[ii], 1.0, DCt+ii, 0, 0, dux+ii, 0, 0.0, dt+ii, nb[ii], dt+ii, nb[ii]);
//		}

		for(ii=0; ii<=N; ii++)
			{
			VECCP_LIBSTR(nb[ii]+ng[ii], dt+ii, 0, dt+ii, nb[ii]+ng[ii]);
			VECSC_LIBSTR(nb[ii]+ng[ii], -1.0, dt+ii, nb[ii]+ng[ii]);
			}

//		exit(2);

//	if(nb>0 | ng>0)
//		{
		COMPUTE_LAM_T_QP(cws);
//		}

	return;

	}


