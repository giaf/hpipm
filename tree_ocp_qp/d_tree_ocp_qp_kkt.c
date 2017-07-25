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

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_d_aux.h>
#include <blasfeo_d_blas.h>

#include "../include/hpipm_tree.h"
#include "../include/hpipm_d_tree_ocp_qp.h"
#include "../include/hpipm_d_tree_ocp_qp_sol.h"
#include "../include/hpipm_d_tree_ocp_qp_ipm_hard.h"
#include "../include/hpipm_d_core_qp_ipm_hard.h"
#include "../include/hpipm_d_core_qp_ipm_hard_aux.h"

#define AXPY_LIBSTR daxpy_libstr
#define COMPUTE_LAM_T_HARD_QP d_compute_lam_t_hard_qp
#define COMPUTE_QX_HARD_QP d_compute_qx_hard_qp
#define COMPUTE_QX_QX_HARD_QP d_compute_Qx_qx_hard_qp
#define DIAAD_SP_LIBSTR ddiaad_sp_libstr
#define GEAD_LIBSTR dgead_libstr
#define GECP_LIBSTR dgecp_libstr
#define GEMM_R_DIAG_LIBSTR dgemm_r_diag_libstr
#define GEMV_N_LIBSTR dgemv_n_libstr
#define GEMV_NT_LIBSTR dgemv_nt_libstr
#define GEMV_T_LIBSTR dgemv_t_libstr
#define IPM_HARD_CORE_QP_WORKSPACE d_ipm_hard_core_qp_workspace
#define IPM_HARD_TREE_OCP_QP_WORKSPACE d_ipm_hard_tree_ocp_qp_workspace
#define POTRF_L_MN_LIBSTR dpotrf_l_mn_libstr
#define REAL double
#define ROWAD_SP_LIBSTR drowad_sp_libstr
#define ROWIN_LIBSTR drowin_libstr
#define ROWEX_LIBSTR drowex_libstr
#define STRMAT d_strmat
#define STRVEC d_strvec
#define SYMV_L_LIBSTR dsymv_l_libstr
#define SYRK_LN_MN_LIBSTR dsyrk_ln_mn_libstr
#define SYRK_POTRF_LN_LIBSTR dsyrk_dpotrf_ln_libstr
#define TREE_OCP_QP d_tree_ocp_qp
#define TREE_OCP_QP_SOL d_tree_ocp_qp_sol
#define TRMM_RLNN_LIBSTR dtrmm_rlnn_libstr
#define TRMV_LNN_LIBSTR dtrmv_lnn_libstr
#define TRMV_LTN_LIBSTR dtrmv_ltn_libstr
#define TRSV_LNN_LIBSTR dtrsv_lnn_libstr
#define TRSV_LNN_MN_LIBSTR dtrsv_lnn_mn_libstr
#define TRSV_LTN_LIBSTR dtrsv_ltn_libstr
#define TRSV_LTN_MN_LIBSTR dtrsv_ltn_mn_libstr
#define VECAD_SP_LIBSTR dvecad_sp_libstr
#define VECCP_LIBSTR dveccp_libstr
#define VECEX_SP_LIBSTR dvecex_sp_libstr
#define VECMULDOT_LIBSTR dvecmuldot_libstr
#define VECSC_LIBSTR dvecsc_libstr

#define INIT_VAR_HARD_TREE_OCP_QP d_init_var_hard_tree_ocp_qp
#define COMPUTE_RES_HARD_TREE_OCP_QP d_compute_res_hard_tree_ocp_qp
#define FACT_SOLVE_KKT_UNCONSTR_TREE_OCP_QP d_fact_solve_kkt_unconstr_tree_ocp_qp
#define FACT_SOLVE_KKT_STEP_HARD_TREE_OCP_QP d_fact_solve_kkt_step_hard_tree_ocp_qp
#define SOLVE_KKT_STEP_HARD_TREE_OCP_QP d_solve_kkt_step_hard_tree_ocp_qp


void INIT_VAR_HARD_TREE_OCP_QP(struct TREE_OCP_QP *qp, struct TREE_OCP_QP_SOL *qp_sol, struct IPM_HARD_TREE_OCP_QP_WORKSPACE *ws)
	{

	struct IPM_HARD_CORE_QP_WORKSPACE *cws = ws->core_workspace;
	
	// loop index
	int ii, jj;

	//
	int Nn = qp->Nn;
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
	for(ii=0; ii<Nn; ii++)
		{
		ux = qp_sol->ux[ii].pa;
		for(jj=0; jj<nu[ii]+nx[ii]; jj++)
			{
			ux[jj] = 0.0;
			}
		}
	
	// pi
	for(ii=0; ii<Nn-1; ii++)
		{
		pi = qp_sol->pi[ii].pa;
		for(jj=0; jj<nx[ii+1]; jj++)
			{
			pi[jj] = 0.0;
			}
		}
	
	// box constraints
	for(ii=0; ii<Nn; ii++)
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
	for(ii=0; ii<Nn; ii++)
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



void COMPUTE_RES_HARD_TREE_OCP_QP(struct TREE_OCP_QP *qp, struct TREE_OCP_QP_SOL *qp_sol, struct IPM_HARD_TREE_OCP_QP_WORKSPACE *ws)
	{

	struct IPM_HARD_CORE_QP_WORKSPACE *cws = ws->core_workspace;

	struct tree *ttree = qp->ttree;
	
	// loop index
	int ii, jj;

	int nkids, idxkid;

	//
	int Nn = qp->Nn;
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

	// loop over nodes
	for(ii=0; ii<Nn; ii++)
		{

		nx0 = nx[ii];
		nu0 = nu[ii];
		nb0 = nb[ii];
		ng0 = ng[ii];

		VECCP_LIBSTR(nu0+nx0, rq+ii, 0, res_g+ii, 0);

		// if not root
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

		// work on kids
		nkids = (ttree->root+ii)->nkids;
		for(jj=0; jj<nkids; jj++)
			{

			idxkid = (ttree->root+ii)->kids[jj];

			nu1 = nu[idxkid];
			nx1 = nx[idxkid];

			AXPY_LIBSTR(nx1, -1.0, ux+idxkid, nu1, b+idxkid-1, 0, res_b+idxkid-1, 0);

			GEMV_NT_LIBSTR(nu0+nx0, nx1, 1.0, 1.0, BAbt+idxkid-1, 0, 0, pi+idxkid-1, 0, ux+ii, 0, 1.0, 1.0, res_g+ii, 0, res_b+idxkid-1, 0, res_g+ii, 0, res_b+idxkid-1, 0);

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
void FACT_SOLVE_KKT_UNCONSTR_TREE_OCP_QP(struct TREE_OCP_QP *qp, struct TREE_OCP_QP_SOL *qp_sol, struct IPM_HARD_TREE_OCP_QP_WORKSPACE *ws)
	{

	int Nn = qp->Nn;
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
	int ii, jj;

	int idx, nkids, idxkid;

	struct IPM_HARD_CORE_QP_WORKSPACE *cws = ws->core_workspace;

	// backward factorization and substitution

	// loop over nodes, starting from the end

	for(ii=0; ii<Nn; ii++)
		{

		idx = Nn-ii-1;

		nkids = (qp->ttree->root+idx)->nkids;

#if defined(DOUBLE_PRECISION)
		TRCP_L_LIBSTR(nu[idx]+nx[idx], RSQrq+idx, 0, 0, L+idx, 0, 0); // TODO dtrcp_l_libstr with m and n, for m>=n
		GECP_LIBSTR(1, nu[idx]+nx[idx], RSQrq+idx, nu[idx]+nx[idx], 0, L+idx, nu[idx]+nx[idx], 0);
#else
		GECP_LIBSTR(nu[idx]+nx[idx]+1, nu[idx]+nx[idx], RSQrq+idx, 0, 0, L+idx, 0, 0); // TODO dtrcp_l_libstr with m and n, for m>=n
#endif

		for(jj=0; jj<nkids; jj++)
			{

			idxkid = (qp->ttree->root+idx)->kids[jj];

			TRMM_RLNN_LIBSTR(nu[idx]+nx[idx]+1, nx[idxkid], 1.0, L+idxkid, nu[idxkid], nu[idxkid], BAbt+idxkid-1, 0, 0, AL, 0, 0);
			GEAD_LIBSTR(1, nx[idxkid], 1.0, L+idxkid, nu[idxkid]+nx[idxkid], nu[idxkid], AL, nu[idx]+nx[idx], 0);

			SYRK_LN_MN_LIBSTR(nu[idx]+nx[idx]+1, nu[idx]+nx[idx], nx[idxkid], 1.0, AL, 0, 0, AL, 0, 0, 1.0, L+idx, 0, 0, L+idx, 0, 0);

			}

		POTRF_L_MN_LIBSTR(nu[idx]+nx[idx]+1, nu[idx]+nx[idx], L+idx, 0, 0, L+idx, 0, 0);

		}


	// forward substitution

	// loop over nodes, starting from the root

	// root
	ii = 0;

	idx = ii;
	nkids = (qp->ttree->root+idx)->nkids;

	ROWEX_LIBSTR(nu[idx]+nx[idx], -1.0, L+idx, nu[idx]+nx[idx], 0, ux+idx, 0);
	TRSV_LTN_LIBSTR(nu[idx]+nx[idx], L+idx, 0, 0, ux+idx, 0, ux+idx, 0);

	for(jj=0; jj<nkids; jj++)
		{

		idxkid = (qp->ttree->root+idx)->kids[jj];

		GEMV_T_LIBSTR(nu[idx]+nx[idx], nx[idxkid], 1.0, BAbt+idxkid-1, 0, 0, ux+idx, 0, 1.0, b+idxkid-1, 0, ux+idxkid, nu[idxkid]);
		ROWEX_LIBSTR(nx[idxkid], 1.0, L+idxkid, nu[idxkid]+nx[idxkid], nu[idxkid], tmp_nxM, 0);
		TRMV_LTN_LIBSTR(nx[idxkid], nx[idxkid], L+idxkid, nu[idxkid], nu[idxkid], ux+idxkid, nu[idxkid], pi+idxkid-1, 0);
		AXPY_LIBSTR(nx[idxkid], 1.0, tmp_nxM, 0, pi+idxkid-1, 0, pi+idxkid-1, 0);
		TRMV_LNN_LIBSTR(nx[idxkid], nx[idxkid], L+idxkid, nu[idxkid], nu[idxkid], pi+idxkid-1, 0, pi+idxkid-1, 0);

		}

	// other nodes
	for(ii=1; ii<Nn; ii++)
		{

		idx = ii;
		nkids = (qp->ttree->root+idx)->nkids;

		ROWEX_LIBSTR(nu[idx], -1.0, L+idx, nu[idx]+nx[idx], 0, ux+idx, 0);
		TRSV_LTN_MN_LIBSTR(nu[idx]+nx[idx], nu[idx], L+idx, 0, 0, ux+idx, 0, ux+idx, 0);

		for(jj=0; jj<nkids; jj++)
			{

			idxkid = (qp->ttree->root+idx)->kids[jj];

			GEMV_T_LIBSTR(nu[idx]+nx[idx], nx[idxkid], 1.0, BAbt+idxkid-1, 0, 0, ux+idx, 0, 1.0, b+idxkid-1, 0, ux+idxkid, nu[idxkid]);
			ROWEX_LIBSTR(nx[idxkid], 1.0, L+idxkid, nu[idxkid]+nx[idxkid], nu[idxkid], tmp_nxM, 0);
			TRMV_LTN_LIBSTR(nx[idxkid], nx[idxkid], L+idxkid, nu[idxkid], nu[idxkid], ux+idxkid, nu[idxkid], pi+idxkid-1, 0);
			AXPY_LIBSTR(nx[idxkid], 1.0, tmp_nxM, 0, pi+idxkid-1, 0, pi+idxkid-1, 0);
			TRMV_LNN_LIBSTR(nx[idxkid], nx[idxkid], L+idxkid, nu[idxkid], nu[idxkid], pi+idxkid-1, 0, pi+idxkid-1, 0);

			}

		}

	return;

	}



// backward Riccati recursion
void FACT_SOLVE_KKT_STEP_HARD_TREE_OCP_QP(struct TREE_OCP_QP *qp, struct IPM_HARD_TREE_OCP_QP_WORKSPACE *ws)
	{

	int Nn = qp->Nn;
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
	int ii, jj;

	int idx, nkids, idxkid;

	struct IPM_HARD_CORE_QP_WORKSPACE *cws = ws->core_workspace;


	COMPUTE_QX_QX_HARD_QP(cws);

	// backward factorization and substitution

	// loop over nodes, starting from the end

	for(ii=0; ii<Nn; ii++)
		{

		idx = Nn-ii-1;

		nkids = (qp->ttree->root+idx)->nkids;

#if defined(DOUBLE_PRECISION)
		TRCP_L_LIBSTR(nu[idx]+nx[idx], RSQrq+idx, 0, 0, L+idx, 0, 0); // TODO dtrcp_l_libstr with m and n, for m>=n
#else
		GECP_LIBSTR(nu[idx]+nx[idx], nu[idx]+nx[idx], RSQrq+idx, 0, 0, L+idx, 0, 0); // TODO dtrcp_l_libstr with m and n, for m>=n
#endif
		ROWIN_LIBSTR(nu[idx]+nx[idx], 1.0, res_g+idx, 0, L+idx, nu[idx]+nx[idx], 0);

		for(jj=0; jj<nkids; jj++)
			{

			idxkid = (qp->ttree->root+idx)->kids[jj];

			GECP_LIBSTR(nu[idx]+nx[idx], nx[idxkid], BAbt+idxkid-1, 0, 0, AL, 0, 0);
			ROWIN_LIBSTR(nx[idxkid], 1.0, res_b+idxkid-1, 0, AL, nu[idx]+nx[idx], 0);
			TRMM_RLNN_LIBSTR(nu[idx]+nx[idx]+1, nx[idxkid], 1.0, L+idxkid, nu[idxkid], nu[idxkid], AL, 0, 0, AL, 0, 0);
			ROWEX_LIBSTR(nx[idxkid], 1.0, AL, nu[idx]+nx[idx], 0, tmp_nxM, 0);
			TRMV_LNN_LIBSTR(nx[idxkid], nx[idxkid], L+idxkid, nu[idxkid], nu[idxkid], tmp_nxM, 0, Pb+idxkid-1, 0);
			GEAD_LIBSTR(1, nx[idxkid], 1.0, L+idxkid, nu[idxkid]+nx[idxkid], nu[idxkid], AL, nu[idx]+nx[idx], 0);

			SYRK_LN_MN_LIBSTR(nu[idx]+nx[idx]+1, nu[idx]+nx[idx], nx[idxkid], 1.0, AL, 0, 0, AL, 0, 0, 1.0, L+idx, 0, 0, L+idx, 0, 0);

			}

		if(nb[idx]>0)
			{
			DIAAD_SP_LIBSTR(nb[idx], 1.0, Qx_lb+idx, 0, idxb[idx], L+idx, 0, 0);
			ROWAD_SP_LIBSTR(nb[idx], 1.0, qx_lb+idx, 0, idxb[idx], L+idx, nu[idx]+nx[idx], 0);
			}
		if(ng[idx]>0)
			{
			GEMM_R_DIAG_LIBSTR(nu[idx]+nx[idx], ng[idx], 1.0, DCt+idx, 0, 0, Qx_lg+idx, 0, 0.0, AL+0, 0, 0, AL+0, 0, 0);
			ROWIN_LIBSTR(ng[idx], 1.0, qx_lg+idx, 0, AL+0, nu[idx]+nx[idx], 0);
			SYRK_POTRF_LN_LIBSTR(nu[idx]+nx[idx]+1, nu[idx]+nx[idx], ng[idx], AL+0, 0, 0, DCt+idx, 0, 0, L+idx, 0, 0, L+idx, 0, 0);
			}
		else
			{
			POTRF_L_MN_LIBSTR(nu[idx]+nx[idx]+1, nu[idx]+nx[idx], L+idx, 0, 0, L+idx, 0, 0);
			}

		}

	// forward substitution

	// loop over nodes, starting from the root

	// root
	ii = 0;

	idx = ii;
	nkids = (qp->ttree->root+idx)->nkids;

	ROWEX_LIBSTR(nu[idx]+nx[idx], -1.0, L+idx, nu[idx]+nx[idx], 0, dux+idx, 0);
	TRSV_LTN_LIBSTR(nu[idx]+nx[idx], L+idx, 0, 0, dux+idx, 0, dux+idx, 0);

	for(jj=0; jj<nkids; jj++)
		{

		idxkid = (qp->ttree->root+idx)->kids[jj];

		GEMV_T_LIBSTR(nu[idx]+nx[idx], nx[idxkid], 1.0, BAbt+idxkid-1, 0, 0, dux+idx, 0, 1.0, res_b+idxkid-1, 0, dux+idxkid, nu[idxkid]);
		ROWEX_LIBSTR(nx[idxkid], 1.0, L+idxkid, nu[idxkid]+nx[idxkid], nu[idxkid], tmp_nxM, 0);
		TRMV_LTN_LIBSTR(nx[idxkid], nx[idxkid], L+idxkid, nu[idxkid], nu[idxkid], dux+idxkid, nu[idxkid], dpi+idxkid-1, 0);
		AXPY_LIBSTR(nx[idxkid], 1.0, tmp_nxM, 0, dpi+idxkid-1, 0, dpi+idxkid-1, 0);
		TRMV_LNN_LIBSTR(nx[idxkid], nx[idxkid], L+idxkid, nu[idxkid], nu[idxkid], dpi+idxkid-1, 0, dpi+idxkid-1, 0);

		}

	// other nodes
	for(ii=1; ii<Nn; ii++)
		{

		idx = ii;
		nkids = (qp->ttree->root+idx)->nkids;

		ROWEX_LIBSTR(nu[idx], -1.0, L+idx, nu[idx]+nx[idx], 0, dux+idx, 0);
		TRSV_LTN_MN_LIBSTR(nu[idx]+nx[idx], nu[idx], L+idx, 0, 0, dux+idx, 0, dux+idx, 0);

		for(jj=0; jj<nkids; jj++)
			{

			idxkid = (qp->ttree->root+idx)->kids[jj];

			GEMV_T_LIBSTR(nu[idx]+nx[idx], nx[idxkid], 1.0, BAbt+idxkid-1, 0, 0, dux+idx, 0, 1.0, res_b+idxkid-1, 0, dux+idxkid, nu[idxkid]);
			ROWEX_LIBSTR(nx[idxkid], 1.0, L+idxkid, nu[idxkid]+nx[idxkid], nu[idxkid], tmp_nxM, 0);
			TRMV_LTN_LIBSTR(nx[idxkid], nx[idxkid], L+idxkid, nu[idxkid], nu[idxkid], dux+idxkid, nu[idxkid], dpi+idxkid-1, 0);
			AXPY_LIBSTR(nx[idxkid], 1.0, tmp_nxM, 0, dpi+idxkid-1, 0, dpi+idxkid-1, 0);
			TRMV_LNN_LIBSTR(nx[idxkid], nx[idxkid], L+idxkid, nu[idxkid], nu[idxkid], dpi+idxkid-1, 0, dpi+idxkid-1, 0);

			}

		}



	for(ii=0; ii<Nn; ii++)
		VECEX_SP_LIBSTR(nb[ii], 1.0, idxb[ii], dux+ii, 0, dt_lb+ii, 0);

	for(ii=0; ii<Nn; ii++)
		GEMV_T_LIBSTR(nu[ii]+nx[ii], ng[ii], 1.0, DCt+ii, 0, 0, dux+ii, 0, 0.0, dt_lg+ii, 0, dt_lg+ii, 0);

	COMPUTE_LAM_T_HARD_QP(cws);

	return;

	}



// backward Riccati recursion
void SOLVE_KKT_STEP_HARD_TREE_OCP_QP(struct TREE_OCP_QP *qp, struct IPM_HARD_TREE_OCP_QP_WORKSPACE *ws)
	{

	int Nn = qp->Nn;
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
	int ii, jj;

	int idx, nkids, idxkid;

	struct IPM_HARD_CORE_QP_WORKSPACE *cws = ws->core_workspace;

	COMPUTE_QX_HARD_QP(cws);


	// backward substitution

	// loop over nodes, starting from the end

	// middle stages
	for(ii=0; ii<Nn-1; ii++)
		{

		idx = Nn-ii-1;

		nkids = (qp->ttree->root+idx)->nkids;

		VECCP_LIBSTR(nu[idx]+nx[idx], res_g+idx, 0, dux+idx, 0);

		for(jj=0; jj<nkids; jj++)
			{

			idxkid = (qp->ttree->root+idx)->kids[jj];

			AXPY_LIBSTR(nx[idxkid], 1.0, dux+idxkid, nu[idxkid], Pb+idxkid-1, 0, tmp_nxM, 0);
			GEMV_N_LIBSTR(nu[idx]+nx[idx], nx[idxkid], 1.0, BAbt+idxkid-1, 0, 0, tmp_nxM, 0, 1.0, dux+idx, 0, dux+idx, 0);

			}

		if(nb[idx]>0)
			{
			VECAD_SP_LIBSTR(nb[idx], 1.0, qx_lb+idx, 0, idxb[idx], dux+idx, 0);
			}
		if(ng[idx]>0)
			{
			GEMV_N_LIBSTR(nu[idx]+nx[idx], ng[idx], 1.0, DCt+idx, 0, 0, qx_lg+idx, 0, 1.0, dux+idx, 0, dux+idx, 0);
			}
		TRSV_LNN_MN_LIBSTR(nu[idx]+nx[idx], nu[idx], L+idx, 0, 0, dux+idx, 0, dux+idx, 0);
		}

	// root
	ii = Nn-1;

	idx = Nn-ii-1;

	nkids = (qp->ttree->root+idx)->nkids;

	VECCP_LIBSTR(nu[idx]+nx[idx], res_g+idx, 0, dux+idx, 0);

	for(jj=0; jj<nkids; jj++)
		{

		idxkid = (qp->ttree->root+idx)->kids[jj];

		AXPY_LIBSTR(nx[idxkid], 1.0, dux+idxkid, nu[idxkid], Pb+idxkid-1, 0, tmp_nxM, 0);
		GEMV_N_LIBSTR(nu[idx]+nx[idx], nx[idxkid], 1.0, BAbt+idxkid-1, 0, 0, tmp_nxM, 0, 1.0, dux+idx, 0, dux+idx, 0);

		}

	if(nb[idx]>0)
		{
		VECAD_SP_LIBSTR(nb[idx], 1.0, qx_lb+idx, 0, idxb[idx], dux+idx, 0);
		}
	if(ng[idx]>0)
		{
		GEMV_N_LIBSTR(nu[idx]+nx[idx], ng[idx], 1.0, DCt+idx, 0, 0, qx_lg+idx, 0, 1.0, dux+idx, 0, dux+idx, 0);
		}
	TRSV_LNN_LIBSTR(nu[idx]+nx[idx], L+idx, 0, 0, dux+idx, 0, dux+idx, 0);


	// forward substitution

	// root
	ii = 0;

	idx = ii;
	nkids = (qp->ttree->root+idx)->nkids;

	VECSC_LIBSTR(nu[idx]+nx[idx], -1.0, dux+idx, 0);
	TRSV_LTN_LIBSTR(nu[idx]+nx[idx], L+idx, 0, 0, dux+idx, 0, dux+idx, 0);

	for(jj=0; jj<nkids; jj++)
		{

		idxkid = (qp->ttree->root+idx)->kids[jj];

		VECCP_LIBSTR(nx[idxkid], dux+idxkid, nu[idxkid], dpi+idxkid-1, 0);
		GEMV_T_LIBSTR(nu[idx]+nx[idx], nx[idxkid], 1.0, BAbt+idxkid-1, 0, 0, dux+idx, 0, 1.0, res_b+idxkid-1, 0, dux+idxkid, nu[idxkid]);
		VECCP_LIBSTR(nx[idxkid], dux+idxkid, nu[idxkid], tmp_nxM, 0);
		TRMV_LTN_LIBSTR(nx[idxkid], nx[idxkid], L+idxkid, nu[idxkid], nu[idxkid], tmp_nxM, 0, tmp_nxM, 0);
		TRMV_LNN_LIBSTR(nx[idxkid], nx[idxkid], L+idxkid, nu[idxkid], nu[idxkid], tmp_nxM, 0, tmp_nxM, 0);
		AXPY_LIBSTR(nx[idxkid], 1.0, tmp_nxM, 0, dpi+idxkid-1, 0, dpi+idxkid-1, 0);

		}

	// other nodes
	for(ii=1; ii<Nn; ii++)
		{

		idx = ii;
		nkids = (qp->ttree->root+idx)->nkids;

		VECSC_LIBSTR(nu[idx], -1.0, dux+idx, 0);
		TRSV_LTN_MN_LIBSTR(nu[idx]+nx[idx], nu[idx], L+idx, 0, 0, dux+idx, 0, dux+idx, 0);

		for(jj=0; jj<nkids; jj++)
			{

			idxkid = (qp->ttree->root+idx)->kids[jj];

			VECCP_LIBSTR(nx[idxkid], dux+idxkid, nu[idxkid], dpi+idxkid-1, 0);
			GEMV_T_LIBSTR(nu[idx]+nx[idx], nx[idxkid], 1.0, BAbt+idxkid-1, 0, 0, dux+idx, 0, 1.0, res_b+idxkid-1, 0, dux+idxkid, nu[idxkid]);
			VECCP_LIBSTR(nx[idxkid], dux+idxkid, nu[idxkid], tmp_nxM, 0);
			TRMV_LTN_LIBSTR(nx[idxkid], nx[idxkid], L+idxkid, nu[idxkid], nu[idxkid], tmp_nxM, 0, tmp_nxM, 0);
			TRMV_LNN_LIBSTR(nx[idxkid], nx[idxkid], L+idxkid, nu[idxkid], nu[idxkid], tmp_nxM, 0, tmp_nxM, 0);
			AXPY_LIBSTR(nx[idxkid], 1.0, tmp_nxM, 0, dpi+idxkid-1, 0, dpi+idxkid-1, 0);

			}

		}



	for(ii=0; ii<Nn; ii++)
		VECEX_SP_LIBSTR(nb[ii], 1.0, idxb[ii], dux+ii, 0, dt_lb+ii, 0);

	for(ii=0; ii<Nn; ii++)
		GEMV_T_LIBSTR(nu[ii]+nx[ii], ng[ii], 1.0, DCt+ii, 0, 0, dux+ii, 0, 0.0, dt_lg+ii, 0, dt_lg+ii, 0);

	COMPUTE_LAM_T_HARD_QP(cws);

	return;

	}


