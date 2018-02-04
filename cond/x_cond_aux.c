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

void COND_BABT(struct OCP_QP *ocp_qp, struct STRMAT *BAbt2, struct STRVEC *b2, struct COND_QP_OCP2DENSE_WORKSPACE *cond_ws)
	{

	int N = ocp_qp->dim->N;

	// early return
	if(N<0)
		return;

	// extract input members
	int *nx = ocp_qp->dim->nx;
	int *nu = ocp_qp->dim->nu;
	struct STRMAT *BAbt = ocp_qp->BAbt;

	// extract memory members
	struct STRMAT *Gamma = cond_ws->Gamma;
	struct STRVEC *Gammab = cond_ws->Gammab;

	int ii, jj;

	int nu_tmp;

	nu_tmp = 0;
	ii = 0;
	// B & A & b
	GECP_LIBSTR(nu[0]+nx[0]+1, nx[1], &BAbt[0], 0, 0, &Gamma[0], 0, 0);
	// b
	ROWEX_LIBSTR(nx[1], 1.0, &Gamma[0], nu[0]+nx[0], 0, &Gammab[0], 0);

	nu_tmp += nu[0];
	ii++;

	for(ii=1; ii<N; ii++)
		{
		// TODO check for equal pointers and avoid copy

		// Gamma * A^T
		GEMM_NN_LIBSTR(nu_tmp+nx[0]+1, nx[ii+1], nx[ii], 1.0, &Gamma[ii-1], 0, 0, &BAbt[ii], nu[ii], 0, 0.0, &Gamma[ii], nu[ii], 0, &Gamma[ii], nu[ii], 0); // Gamma * A^T

		GECP_LIBSTR(nu[ii], nx[ii+1], &BAbt[ii], 0, 0, &Gamma[ii], 0, 0);

		nu_tmp += nu[ii];

		GEAD_LIBSTR(1, nx[ii+1], 1.0, &BAbt[ii], nu[ii]+nx[ii], 0, &Gamma[ii], nu_tmp+nx[0], 0);

		ROWEX_LIBSTR(nx[ii+1], 1.0, &Gamma[ii], nu_tmp+nx[0], 0, &Gammab[ii], 0);
		}
	
	if(cond_ws->cond_last_stage==0)
		{
		// B & A & b
		GECP_LIBSTR(nu_tmp+nx[0]+1, nx[N], &Gamma[N-1], 0, 0, &BAbt2[0], 0, 0);
		// b
		ROWEX_LIBSTR(nx[N], 1.0, &BAbt2[0], nu_tmp+nx[0], 0, &b2[0], 0);
		}

	return;

	}



void COND_B(struct OCP_QP *ocp_qp, struct STRVEC *b2, struct COND_QP_OCP2DENSE_WORKSPACE *cond_ws)
	{

	int N = ocp_qp->dim->N;

	// early return
	if(N<0)
		return;

	// extract input members
	int *nx = ocp_qp->dim->nx;
	int *nu = ocp_qp->dim->nu;
	struct STRMAT *BAbt = ocp_qp->BAbt;
	struct STRVEC *b = ocp_qp->b;

	// extract memory members
	struct STRVEC *Gammab = cond_ws->Gammab;

	int ii, jj;

	int nu_tmp;

	ii = 0;
	// b
	VECCP_LIBSTR(nx[1], b+0, 0, Gammab+0, 0);

	ii++;

	for(ii=1; ii<N; ii++)
		{
		// TODO check for equal pointers and avoid copy

		// A * Gammab
		GEMV_T_LIBSTR(nx[ii], nx[ii+1], 1.0, BAbt+ii, nu[ii], 0, Gammab+(ii-1), 0, 1.0, b+ii, 0, Gammab+ii, 0);
		}
	
	if(cond_ws->cond_last_stage==0)
		{
		// b
		VECCP_LIBSTR(nx[N], Gammab+(N-1), 0, b2+0, 0);
		}

	return;

	}



void COND_RSQRQ_N2NX3(struct OCP_QP *ocp_qp, struct STRMAT *RSQrq2, struct STRVEC *rq2, struct COND_QP_OCP2DENSE_WORKSPACE *cond_ws)
	{

	int N = ocp_qp->dim->N;
	if(cond_ws->cond_last_stage==0)
		N -= 1;

	// early return
	if(N<0)
		return;

	// extract input members
	int *nx = ocp_qp->dim->nx;
	int *nu = ocp_qp->dim->nu;

	struct STRMAT *BAbt = ocp_qp->BAbt;
	struct STRMAT *RSQrq = ocp_qp->RSQrq;

	// extract memory members
	struct STRMAT *Gamma = cond_ws->Gamma;
	struct STRMAT *L = cond_ws->L;
	struct STRMAT *Lx = cond_ws->Lx;
	struct STRMAT *AL = cond_ws->AL;

	// early return
	if(N==0)
		{
		GECP_LIBSTR(nu[0]+nx[0]+1, nu[0]+nx[0], &RSQrq[0], 0, 0, &RSQrq2[0], 0, 0); // XXX when removing rq from RSQ !!!
		ROWEX_LIBSTR(nu[0]+nx[0], 1.0, &RSQrq[0], nu[0]+nx[0], 0, &rq2[0], 0); // XXX when removing rq from RSQ !!!
		return;
		}

	int nn;

	int nu2 = 0; // sum of all nu
	for(nn=0; nn<=N; nn++)
		nu2 += nu[nn];
	
	int nub = nu2; // backward partial sum
	int nuf = 0; // forward partial sum

	// final stage 
	nub -= nu[N];

	// TODO g is not part of H !!!!!
	GECP_LIBSTR(nu[N]+nx[N]+1, nu[N]+nx[N], &RSQrq[N], 0, 0, &L[N], 0, 0);

	// D
	TRCP_L_LIBSTR(nu[N], &L[N], 0, 0, &RSQrq2[0], nuf, nuf);

	GEMM_NN_LIBSTR(nub+nx[0]+1, nu[N], nx[N], 1.0, &Gamma[N-1], 0, 0, &L[N], nu[N], 0, 0.0, &RSQrq2[0], nuf+nu[N], nuf, &RSQrq2[0], nuf+nu[N], nuf);

	// m
	GEAD_LIBSTR(1, nu[N], 1.0, &L[N], nu[N]+nx[N], 0, &RSQrq2[0], nu2+nx[0], nuf);

	nuf += nu[N];



	// middle stages 
	for(nn=0; nn<N-1; nn++)
		{	
		nub -= nu[N-nn-1];

#if defined(LA_HIGH_PERFORMANCE)
		GECP_LIBSTR(nx[N-nn]+1, nx[N-nn], &L[N-nn], nu[N-nn], nu[N-nn], Lx, 0, 0);

		POTRF_L_MN_LIBSTR(nx[N-nn]+1, nx[N-nn], Lx, 0, 0, Lx, 0, 0);
#else
		POTRF_L_MN_LIBSTR(nx[N-nn]+1, nx[N-nn], &L[N-nn], nu[N-nn], nu[N-nn], Lx, 0, 0);
#endif
		TRMM_RLNN_LIBSTR(nu[N-nn-1]+nx[N-nn-1]+1, nx[N-nn], 1.0, Lx, 0, 0, &BAbt[N-nn-1], 0, 0, AL, 0, 0);
		GEAD_LIBSTR(1, nx[N-nn], 1.0, Lx, nx[N-nn], 0, AL, nu[N-nn-1]+nx[N-nn-1], 0);

		SYRK_LN_MN_LIBSTR(nu[N-nn-1]+nx[N-nn-1]+1, nu[N-nn-1]+nx[N-nn-1], nx[N-nn], 1.0, AL, 0, 0, AL, 0, 0, 1.0, &RSQrq[N-nn-1], 0, 0, &L[N-nn-1], 0, 0);

		// D
		TRCP_L_LIBSTR(nu[N-nn-1], &L[N-nn-1], 0, 0, &RSQrq2[0], nuf, nuf);

		GEMM_NN_LIBSTR(nub+nx[0]+1, nu[N-nn-1], nx[N-nn-1], 1.0, &Gamma[N-nn-2], 0, 0, &L[N-nn-1], nu[N-nn-1], 0, 0.0, &RSQrq2[0], nuf+nu[N-nn-1], nuf, &RSQrq2[0], nuf+nu[N-nn-1], nuf);

		// m
		GEAD_LIBSTR(1, nu[N-nn-1], 1.0, &L[N-nn-1], nu[N-nn-1]+nx[N-nn-1], 0, &RSQrq2[0], nu2+nx[0], nuf);

		nuf += nu[N-nn-1];

		}

	// first stage
	nn = N-1;

#if defined(LA_HIGH_PERFORMANCE)
	GECP_LIBSTR(nx[N-nn]+1, nx[N-nn], &L[N-nn], nu[N-nn], nu[N-nn], Lx, 0, 0);

	POTRF_L_MN_LIBSTR(nx[N-nn]+1, nx[N-nn], Lx, 0, 0, Lx, 0, 0);
#else
	POTRF_L_MN_LIBSTR(nx[N-nn]+1, nx[N-nn], &L[N-nn], nu[N-nn], nu[N-nn], Lx, 0, 0);
#endif
	TRMM_RLNN_LIBSTR(nu[N-nn-1]+nx[N-nn-1]+1, nx[N-nn], 1.0, Lx, 0, 0, &BAbt[N-nn-1], 0, 0, AL, 0, 0);
	GEAD_LIBSTR(1, nx[N-nn], 1.0, Lx, nx[N-nn], 0, AL, nu[N-nn-1]+nx[N-nn-1], 0);

	SYRK_LN_MN_LIBSTR(nu[N-nn-1]+nx[N-nn-1]+1, nu[N-nn-1]+nx[N-nn-1], nx[N-nn], 1.0, AL, 0, 0, AL, 0, 0, 1.0, &RSQrq[N-nn-1], 0, 0, &L[N-nn-1], 0, 0);

	// D, M, m, P, p
//	GECP_LIBSTR(nu[0]+nx[0]+1, nu[0]+nx[0], &L[N-nn-1], 0, 0, &RSQrq2[0], nuf, nuf); // TODO dtrcp for 'rectangular' matrices
	TRCP_L_LIBSTR(nu[0]+nx[0], &L[N-nn-1], 0, 0, &RSQrq2[0], nuf, nuf); // TODO dtrcp for 'rectangular' matrices
	GECP_LIBSTR(1, nu[0]+nx[0], &L[N-nn-1], nu[0]+nx[0], 0, &RSQrq2[0], nuf+nu[0]+nx[0], nuf); // TODO dtrcp for 'rectangular' matrices
	// m p
	ROWEX_LIBSTR(nu2+nx[0], 1.0, &RSQrq2[0], nu2+nx[0], 0, &rq2[0], 0);

	return;

	}



void COND_RQ_N2NX3(struct OCP_QP *ocp_qp, struct STRVEC *rq2, struct COND_QP_OCP2DENSE_WORKSPACE *cond_ws)
	{

	int N = ocp_qp->dim->N;
	if(cond_ws->cond_last_stage==0)
		N -= 1;

	// early return
	if(N<0)
		return;

	// extract input members
	int *nx = ocp_qp->dim->nx;
	int *nu = ocp_qp->dim->nu;

	struct STRMAT *BAbt = ocp_qp->BAbt;
	struct STRVEC *b = ocp_qp->b;
	struct STRVEC *rq = ocp_qp->rq;

	// extract memory members
	struct STRMAT *L = cond_ws->L;
	struct STRVEC *Gammab = cond_ws->Gammab;
	struct STRVEC *l = cond_ws->l;
	struct STRVEC *tmp_nuxM = cond_ws->tmp_nuxM;

	// early return
	if(N==0)
		{
		VECCP_LIBSTR(nu[0]+nx[0], rq+0, 0, rq2+0, 0);
		return;
		}

	int nn;

	int nu2 = 0; // sum of all nu
	for(nn=0; nn<=N; nn++)
		nu2 += nu[nn];
	
	int nub = nu2; // backward partial sum
	int nuf = 0; // forward partial sum

	// final stage 
	nub -= nu[N];

	VECCP_LIBSTR(nu[N]+nx[N], rq+N, 0, l+N, 0);

	GEMV_T_LIBSTR(nx[N], nu[N], 1.0, L+N, nu[N], 0, Gammab+(N-1), 0, 1.0, l+N, 0, rq2+0, nuf);

	nuf += nu[N];



	// middle stages 
	for(nn=0; nn<N-1; nn++)
		{	
		nub -= nu[N-nn-1];

		SYMV_L_LIBSTR(nx[N-nn], nx[N-nn], 1.0, L+(N-nn), nu[N-nn], nu[N-nn], b+(N-nn-1), 0, 1.0, l+(N-nn), nu[N-nn], tmp_nuxM, 0);

		GEMV_N_LIBSTR(nu[N-nn-1]+nx[N-nn-1], nx[N-nn], 1.0, BAbt+(N-nn-1), 0, 0, tmp_nuxM, 0, 1.0, rq+(N-nn-1), 0, l+(N-nn-1), 0);

		GEMV_T_LIBSTR(nx[N-nn-1], nu[N-nn-1], 1.0, L+(N-nn-1), nu[N-nn-1], 0, Gammab+(N-nn-2), 0, 1.0, l+(N-nn-1), 0, rq2+0, nuf);

		nuf += nu[N-nn-1];

		}

	// first stage
	nn = N-1;

	SYMV_L_LIBSTR(nx[N-nn], nx[N-nn], 1.0, L+(N-nn), nu[N-nn], nu[N-nn], b+(N-nn-1), 0, 1.0, l+(N-nn), nu[N-nn], tmp_nuxM, 0);

	GEMV_N_LIBSTR(nu[N-nn-1]+nx[N-nn-1], nx[N-nn], 1.0, BAbt+(N-nn-1), 0, 0, tmp_nuxM, 0, 1.0, rq+(N-nn-1), 0, l+(N-nn-1), 0);

	// m p
	VECCP_LIBSTR(nu[0]+nx[0], l+(N-nn-1), 0, rq2+0, nuf);

	return;

	}



void COND_DCTD(struct OCP_QP *ocp_qp, int *idxb2, struct STRMAT *DCt2, struct STRVEC *d2, int *idxs2, struct STRVEC *Z2, struct STRVEC *z2, struct COND_QP_OCP2DENSE_WORKSPACE *cond_ws)
	{

	int N = ocp_qp->dim->N;
	if(cond_ws->cond_last_stage==0)
		N -= 1;

	// early return
	if(N<0)
		return;

	// extract input members
	int *nx = ocp_qp->dim->nx;
	int *nu = ocp_qp->dim->nu;
	int *nb = ocp_qp->dim->nb;
	int *ng = ocp_qp->dim->ng;
	int *ns = ocp_qp->dim->ns;

	int **idxb = ocp_qp->idxb;
	struct STRVEC *d = ocp_qp->d;
	struct STRMAT *DCt = ocp_qp->DCt;
	int **idxs = ocp_qp->idxs;
	struct STRVEC *Z = ocp_qp->Z;
	struct STRVEC *z = ocp_qp->z;

	// extract memory members
	struct STRMAT *Gamma = cond_ws->Gamma;
	struct STRVEC *Gammab = cond_ws->Gammab;
	struct STRVEC *tmp_ngM = cond_ws->tmp_ngM;
	int *idxs_rev = cond_ws->idxs_rev;


	REAL *ptr_d_lb;
	REAL *ptr_d_ub;
	
	int nu_tmp, ng_tmp;

	int ii, jj;

	int nu0, nx0, nb0, ng0, ns0;

	// problem size

	int nbb = nb[0]; // box that remain box constraints
	int nbg = 0; // box that becomes general constraints
	for(ii=1; ii<=N; ii++)
		for(jj=0; jj<nb[ii]; jj++)
			if(idxb[ii][jj]<nu[ii])
				nbb++;
			else
				nbg++;
	
	int nx2 = nx[0];
	int nu2 = nu[0];
	int ns2 = ns[0];
	int ngg = ng[0];
	for(ii=1; ii<=N; ii++)
		{
		nu2 += nu[ii];
		ns2 += ns[ii];
		ngg += ng[ii];
		}
	int ng2 = nbg + ngg;
	int nb2 = nbb;
	int nt2 = nb2 + ng2;

	REAL *d_lb3 = d2->pa+0;
	REAL *d_ub3 = d2->pa+nb2+ng2;
	REAL *d_lg3 = d2->pa+nb2;
	REAL *d_ug3 = d2->pa+2*nb2+ng2;

	// set constraint matrix to zero (it's 2 lower triangular matrices atm)
	GESE_LIBSTR(nu2+nx2, ng2, 0.0, &DCt2[0], 0, 0);

	// box constraints

	int idx_gammab = nx[0];
	for(ii=0; ii<N; ii++)
		idx_gammab += nu[ii];

	int ib = 0;
	int ig = 0;
	int is = 0; // XXX

	int idxb0, idxg0;

	REAL *ptr_Z;
	REAL *ptr_z;
	REAL *ptr_Z2 = Z2->pa;
	REAL *ptr_z2 = z2->pa;

	REAL tmp;
	int idx_g;

	// middle stages
	nu_tmp = 0;
	for(ii=0; ii<N; ii++)
		{
		nu0 = nu[N-ii];
		nb0 = nb[N-ii];
		ng0 = ng[N-ii];
		ns0 = ns[N-ii];
		for(jj=0; jj<nb0+ng0; jj++)
			idxs_rev[jj] = -1;
		if(ns0>0)
			{
			for(jj=0; jj<ns0; jj++)
				idxs_rev[idxs[N-ii][jj]] = jj;
			ptr_Z = Z[N-ii].pa;
			ptr_z = z[N-ii].pa;
			}
//		int_print_mat(1, nb0+ng0, idxs_rev, 1);
		nu_tmp += nu0;
		ptr_d_lb = d[N-ii].pa+0;
		ptr_d_ub = d[N-ii].pa+nb0+ng0;
		for(jj=0; jj<nb0; jj++)
			{
			idxb0 = idxb[N-ii][jj];
			if(idxb0<nu0) // input: box constraint
				{
				d_lb3[ib] = ptr_d_lb[jj];
				d_ub3[ib] = ptr_d_ub[jj];
				idxb2[ib] = nu_tmp - nu0 + idxb[N-ii][jj];
				if(idxs_rev[jj]>=0)
					{
					idxs2[is] = ib;
					ptr_Z2[0+is]   = ptr_Z[0+idxs_rev[jj]];
					ptr_Z2[ns2+is] = ptr_Z[ns[N-ii]+idxs_rev[jj]];
					ptr_z2[0+is]   = ptr_z[0+idxs_rev[jj]];
					ptr_z2[ns2+is] = ptr_z[ns[N-ii]+idxs_rev[jj]];
					is++;
					}
				ib++;
				}
			else // state: general constraint
				{
				idx_g = idxb0-nu0;
				tmp = GEEX1_LIBSTR(&Gamma[N-1-ii], idx_gammab, idx_g);
				d_lg3[ig] = ptr_d_lb[jj] - tmp;
//				d_ug3[ig] = ptr_d_ub[jj] - tmp;
				d_ug3[ig] = ptr_d_ub[jj] + tmp; // XXX
				GECP_LIBSTR(idx_gammab, 1, &Gamma[N-ii-1], 0, idx_g, &DCt2[0], nu_tmp, ig);
				if(idxs_rev[jj]>=0)
					{
					idxs2[is] = nb2+ig;
					ptr_Z2[0+is]   = ptr_Z[0+idxs_rev[jj]];
					ptr_Z2[ns2+is] = ptr_Z[ns0+idxs_rev[jj]];
					ptr_z2[0+is]   = ptr_z[0+idxs_rev[jj]];
					ptr_z2[ns2+is] = ptr_z[ns0+idxs_rev[jj]];
					is++;
					}
				ig++;
				}
			}
		idx_gammab -= nu[N-1-ii];
		}

	// initial stage: both inputs and states as box constraints
	nu0 = nu[0];
	nb0 = nb[0];
	ng0 = ng[0];
	ns0 = ns[0];
	for(jj=0; jj<nb0+ng0; jj++)
		idxs_rev[jj] = -1;
	if(ns0>0)
		{
		for(jj=0; jj<ns0; jj++)
			idxs_rev[idxs[0][jj]] = jj;
		ptr_Z = Z[0].pa;
		ptr_z = z[0].pa;
		}
//	int_print_mat(1, nb0+ng0, idxs_rev, 1);
	nu_tmp += nu0;
	ptr_d_lb = d[0].pa+0;
	ptr_d_ub = d[0].pa+nb0+ng0;
	for(jj=0; jj<nb0; jj++)
		{
		idxb0 = idxb[0][jj];
		d_lb3[ib] = ptr_d_lb[jj];
		d_ub3[ib] = ptr_d_ub[jj];
		idxb2[ib] = nu_tmp - nu0 + idxb0;
		if(idxs_rev[jj]>=0) // XXX
			{
			idxs2[is] = ib;
			ptr_Z2[0+is]   = ptr_Z[0+idxs_rev[jj]];
			ptr_Z2[ns2+is] = ptr_Z[ns0+idxs_rev[jj]];
			ptr_z2[0+is]   = ptr_z[0+idxs_rev[jj]];
			ptr_z2[ns2+is] = ptr_z[ns0+idxs_rev[jj]];
			is++;
			}
		ib++;
		}

	// XXX for now, just shift after box-to-general constraints
	// better interleave them, to keep the block lower trianlgular structure !!!

	// general constraints

	char *c_ptr;

	nu_tmp = 0;
	ng_tmp = 0;
	for(ii=0; ii<N; ii++)
		{

		nx0 = nx[N-ii];
		nu0 = nu[N-ii];
		nb0 = nb[N-ii];
		ng0 = ng[N-ii];
		ns0 = ns[N-ii];
		for(jj=0; jj<nb0+ng0; jj++)
			idxs_rev[jj] = -1;
		if(ns0>0)
			{
			for(jj=0; jj<ns0; jj++)
				idxs_rev[idxs[N-ii][jj]] = jj;
			ptr_Z = Z[N-ii].pa;
			ptr_z = z[N-ii].pa;
			}
//		int_print_mat(1, nb0+ng0, idxs_rev, 1);

		if(ng0>0)
			{

			for(ig=0; ig<ng0; ig++)
				{

				idxg0 = nb0+ig;

				if(idxs_rev[idxg0]>=0)
					{
					idxs2[is] = nb2+nbg+ng_tmp+ig;
					ptr_Z2[0+is]   = ptr_Z[0+idxs_rev[idxg0]];
					ptr_Z2[ns2+is] = ptr_Z[ns0+idxs_rev[idxg0]];
					ptr_z2[0+is]   = ptr_z[0+idxs_rev[idxg0]];
					ptr_z2[ns2+is] = ptr_z[ns0+idxs_rev[idxg0]];
					is++;
					}

				}

			GECP_LIBSTR(nu0, ng0, &DCt[N-ii], 0, 0, DCt2, nu_tmp, nbg+ng_tmp);

			nu_tmp += nu0;

			GEMM_NN_LIBSTR(nu2+nx[0]-nu_tmp, ng0, nx0, 1.0, &Gamma[N-1-ii], 0, 0, &DCt[N-ii], nu0, 0, 0.0, DCt2, nu_tmp, nbg+ng_tmp, DCt2, nu_tmp, nbg+ng_tmp);

			VECCP_LIBSTR(ng0, &d[N-ii], nb0, d2, nb2+nbg+ng_tmp);
			VECCP_LIBSTR(ng0, &d[N-ii], 2*nb0+ng0, d2, 2*nb2+ng2+nbg+ng_tmp);

			GEMV_T_LIBSTR(nx0, ng0, 1.0, &DCt[N-ii], nu0, 0, &Gammab[N-ii-1], 0, 0.0, tmp_ngM, 0, tmp_ngM, 0);

			AXPY_LIBSTR(ng0, -1.0, tmp_ngM, 0, d2, nb2+nbg+ng_tmp, d2, nb2+nbg+ng_tmp);
//			AXPY_LIBSTR(ng0, -1.0, tmp_ngM, 0, d2, 2*nb2+ng2+nbg+ng_tmp, d2, 2*nb2+ng2+nbg+ng_tmp);
			AXPY_LIBSTR(ng0,  1.0, tmp_ngM, 0, d2, 2*nb2+ng2+nbg+ng_tmp, d2, 2*nb2+ng2+nbg+ng_tmp); // XXX

			ng_tmp += ng0;
			
			}
		else
			{

			nu_tmp += nu0;

			}

		}

	ii = N;

	nx0 = nx[0];
	nu0 = nu[0];
	nb0 = nb[0];
	ng0 = ng[0];
	ns0 = ns[0];
	for(jj=0; jj<nb0+ng0; jj++)
		idxs_rev[jj] = -1;
	if(ns0>0)
		{
		for(jj=0; jj<ns0; jj++)
			idxs_rev[idxs[N-ii][jj]] = jj;
		ptr_Z = Z[N-ii].pa;
		ptr_z = z[N-ii].pa;
		}
//	int_print_mat(1, nb0+ng0, idxs_rev, 1);

	if(ng0>0)
		{

		for(ig=0; ig<ng0; ig++)
			{

			idxg0 = nb0+ig;

			if(idxs_rev[idxg0]>=0)
				{
				idxs2[is] = nb2+nbg+ng_tmp+ig;
				ptr_Z2[0+is]   = ptr_Z[0+idxs_rev[idxg0]];
				ptr_Z2[ns2+is] = ptr_Z[ns0+idxs_rev[idxg0]];
				ptr_z2[0+is]   = ptr_z[0+idxs_rev[idxg0]];
				ptr_z2[ns2+is] = ptr_z[ns0+idxs_rev[idxg0]];
				is++;
				}

			}

		GECP_LIBSTR(nu0+nx0, ng0, &DCt[0], 0, 0, DCt2, nu_tmp, nbg+ng_tmp);

		VECCP_LIBSTR(ng0, &d[0], nb0, d2, nb2+nbg+ng_tmp);
		VECCP_LIBSTR(ng0, &d[0], 2*nb0+ng0, d2, 2*nb2+ng2+nbg+ng_tmp);

//		ng_tmp += ng[N-ii];

		}

	return;

	}



void COND_D(struct OCP_QP *ocp_qp, struct STRVEC *d2, struct STRVEC *z2, struct COND_QP_OCP2DENSE_WORKSPACE *cond_ws)
	{

	int N = ocp_qp->dim->N;
	if(cond_ws->cond_last_stage==0)
		N -= 1;

	// early return
	if(N<0)
		return;

	// extract input members
	int *nx = ocp_qp->dim->nx;
	int *nu = ocp_qp->dim->nu;
	int *nb = ocp_qp->dim->nb;
	int *ng = ocp_qp->dim->ng;
	int *ns = ocp_qp->dim->ns;

	int **idxb = ocp_qp->idxb;
	struct STRVEC *d = ocp_qp->d;
	struct STRMAT *DCt = ocp_qp->DCt;
	int **idxs = ocp_qp->idxs;
	struct STRVEC *Z = ocp_qp->Z;
	struct STRVEC *z = ocp_qp->z;

	// extract memory members
	struct STRMAT *Gamma = cond_ws->Gamma;
	struct STRVEC *Gammab = cond_ws->Gammab;
	struct STRVEC *tmp_ngM = cond_ws->tmp_ngM;
	int *idxs_rev = cond_ws->idxs_rev;


	REAL *ptr_d_lb;
	REAL *ptr_d_ub;
	
	int nu_tmp, ng_tmp;

	int ii, jj;

	int nu0, nx0, nb0, ng0, ns0;

	// problem size

	int nbb = nb[0]; // box that remain box constraints
	int nbg = 0; // box that becomes general constraints
	for(ii=1; ii<=N; ii++)
		for(jj=0; jj<nb[ii]; jj++)
			if(idxb[ii][jj]<nu[ii])
				nbb++;
			else
				nbg++;
	
	int nx2 = nx[0];
	int nu2 = nu[0];
	int ns2 = ns[0];
	int ngg = ng[0];
	for(ii=1; ii<=N; ii++)
		{
		nu2 += nu[ii];
		ns2 += ns[ii];
		ngg += ng[ii];
		}
	int ng2 = nbg + ngg;
	int nb2 = nbb;
	int nt2 = nb2 + ng2;

	REAL *d_lb3 = d2->pa+0;
	REAL *d_ub3 = d2->pa+nb2+ng2;
	REAL *d_lg3 = d2->pa+nb2;
	REAL *d_ug3 = d2->pa+2*nb2+ng2;

	// box constraints

	int idx_gammab = nx[0];
	for(ii=0; ii<N; ii++)
		idx_gammab += nu[ii];

	int ib = 0;
	int ig = 0;
	int is = 0; // XXX

	int idxb0, idxg0;

	REAL *ptr_z;
	REAL *ptr_z2 = z2->pa;

	REAL tmp;
	int idx_g;

	// middle stages
	nu_tmp = 0;
	for(ii=0; ii<N; ii++)
		{
		nu0 = nu[N-ii];
		nb0 = nb[N-ii];
		ng0 = ng[N-ii];
		ns0 = ns[N-ii];
		for(jj=0; jj<nb0+ng0; jj++)
			idxs_rev[jj] = -1;
		if(ns0>0)
			{
			for(jj=0; jj<ns0; jj++)
				idxs_rev[idxs[N-ii][jj]] = jj;
			ptr_z = z[N-ii].pa;
			}
		nu_tmp += nu0;
		ptr_d_lb = d[N-ii].pa+0;
		ptr_d_ub = d[N-ii].pa+nb0+ng0;
		for(jj=0; jj<nb0; jj++)
			{
			idxb0 = idxb[N-ii][jj];
			if(idxb0<nu0) // input: box constraint
				{
				d_lb3[ib] = ptr_d_lb[jj];
				d_ub3[ib] = ptr_d_ub[jj];
				if(idxs_rev[jj]>=0)
					{
					ptr_z2[0+is]   = ptr_z[0+idxs_rev[jj]];
					ptr_z2[ns2+is] = ptr_z[ns[N-ii]+idxs_rev[jj]];
					is++;
					}
				ib++;
				}
			else // state: general constraint
				{
				idx_g = idxb0-nu0;
				tmp = GEEX1_LIBSTR(&Gamma[N-1-ii], idx_gammab, idx_g);
				d_lg3[ig] = ptr_d_lb[jj] - tmp;
				d_ug3[ig] = ptr_d_ub[jj] - tmp;
				if(idxs_rev[jj]>=0)
					{
					ptr_z2[0+is]   = ptr_z[0+idxs_rev[jj]];
					ptr_z2[ns2+is] = ptr_z[ns0+idxs_rev[jj]];
					is++;
					}
				ig++;
				}
			}
		idx_gammab -= nu[N-1-ii];
		}

	// initial stage: both inputs and states as box constraints
	nu0 = nu[0];
	nb0 = nb[0];
	ng0 = ng[0];
	ns0 = ns[0];
	for(jj=0; jj<nb0+ng0; jj++)
		idxs_rev[jj] = -1;
	if(ns0>0)
		{
		for(jj=0; jj<ns0; jj++)
			idxs_rev[idxs[0][jj]] = jj;
		ptr_z = z[0].pa;
		}
	nu_tmp += nu0;
	ptr_d_lb = d[0].pa+0;
	ptr_d_ub = d[0].pa+nb0+ng0;
	for(jj=0; jj<nb0; jj++)
		{
		idxb0 = idxb[0][jj];
		d_lb3[ib] = ptr_d_lb[jj];
		d_ub3[ib] = ptr_d_ub[jj];
		if(idxs_rev[jj]>=0) // XXX
			{
			ptr_z2[0+is]   = ptr_z[0+idxs_rev[jj]];
			ptr_z2[ns2+is] = ptr_z[ns0+idxs_rev[jj]];
			is++;
			}
		ib++;
		}

	// XXX for now, just shift after box-to-general constraints
	// better interleave them, to keep the block lower trianlgular structure !!!

	// general constraints

	char *c_ptr;

	nu_tmp = 0;
	ng_tmp = 0;
	for(ii=0; ii<N; ii++)
		{

		nx0 = nx[N-ii];
		nu0 = nu[N-ii];
		nb0 = nb[N-ii];
		ng0 = ng[N-ii];
		ns0 = ns[N-ii];
		for(jj=0; jj<nb0+ng0; jj++)
			idxs_rev[jj] = -1;
		if(ns0>0)
			{
			for(jj=0; jj<ns0; jj++)
				idxs_rev[idxs[N-ii][jj]] = jj;
			ptr_z = z[N-ii].pa;
			}

		if(ng0>0)
			{

			for(ig=0; ig<ng0; ig++)
				{

				idxg0 = nb0+ig;

				if(idxs_rev[idxg0]>=0)
					{
					ptr_z2[0+is]   = ptr_z[0+idxs_rev[idxg0]];
					ptr_z2[ns2+is] = ptr_z[ns0+idxs_rev[idxg0]];
					is++;
					}

				}

			nu_tmp += nu0;

			VECCP_LIBSTR(ng0, &d[N-ii], nb0, d2, nb2+nbg+ng_tmp);
			VECCP_LIBSTR(ng0, &d[N-ii], 2*nb0+ng0, d2, 2*nb2+ng2+nbg+ng_tmp);

			GEMV_T_LIBSTR(nx0, ng0, 1.0, &DCt[N-ii], nu0, 0, &Gammab[N-ii-1], 0, 0.0, tmp_ngM, 0, tmp_ngM, 0);

			AXPY_LIBSTR(ng0, -1.0, tmp_ngM, 0, d2, nb2+nbg+ng_tmp, d2, nb2+nbg+ng_tmp);
			AXPY_LIBSTR(ng0, -1.0, tmp_ngM, 0, d2, 2*nb2+ng2+nbg+ng_tmp, d2, 2*nb2+ng2+nbg+ng_tmp);

			ng_tmp += ng0;
			
			}
		else
			{

			nu_tmp += nu0;

			}

		}

	ii = N;

	nx0 = nx[0];
	nu0 = nu[0];
	nb0 = nb[0];
	ng0 = ng[0];
	ns0 = ns[0];
	for(jj=0; jj<nb0+ng0; jj++)
		idxs_rev[jj] = -1;
	if(ns0>0)
		{
		for(jj=0; jj<ns0; jj++)
			idxs_rev[idxs[N-ii][jj]] = jj;
		ptr_z = z[N-ii].pa;
		}

	if(ng0>0)
		{

		for(ig=0; ig<ng0; ig++)
			{

			idxg0 = nb0+ig;

			if(idxs_rev[idxg0]>=0)
				{
				ptr_z2[0+is]   = ptr_z[0+idxs_rev[idxg0]];
				ptr_z2[ns2+is] = ptr_z[ns0+idxs_rev[idxg0]];
				is++;
				}

			}

		VECCP_LIBSTR(ng0, &d[0], nb0, d2, nb2+nbg+ng_tmp);
		VECCP_LIBSTR(ng0, &d[0], 2*nb0+ng0, d2, 2*nb2+ng2+nbg+ng_tmp);

//		ng_tmp += ng[N-ii];

		}

	return;

	}



void EXPAND_SOL(struct OCP_QP *ocp_qp, struct DENSE_QP_SOL *dense_qp_sol, struct OCP_QP_SOL *ocp_qp_sol, struct COND_QP_OCP2DENSE_WORKSPACE *cond_ws)
	{

	int N = ocp_qp->dim->N;
	int Np = N;
	if(cond_ws->cond_last_stage==0)
		N -= 1;

	int ii, jj;

	int *nu = ocp_qp->dim->nu;
	int *nx = ocp_qp->dim->nx;
	int *nb = ocp_qp->dim->nb;
	int *ng = ocp_qp->dim->ng;
	int *ns = ocp_qp->dim->ns;

	struct STRMAT *BAbt = ocp_qp->BAbt;
	struct STRVEC *b = ocp_qp->b;
	int **idxb = ocp_qp->idxb;
	struct STRMAT *RSQrq = ocp_qp->RSQrq;
	struct STRVEC *rq = ocp_qp->rq;
	struct STRMAT *DCt = ocp_qp->DCt;
	int **idxs = ocp_qp->idxs;

	struct STRVEC *vc = dense_qp_sol->v;
	struct STRVEC *pic = dense_qp_sol->pi;
	struct STRVEC *lamc = dense_qp_sol->lam;
	struct STRVEC *tc = dense_qp_sol->t;

	struct STRVEC *ux = ocp_qp_sol->ux;
	struct STRVEC *pi = ocp_qp_sol->pi;
	struct STRVEC *lam = ocp_qp_sol->lam;
	struct STRVEC *t = ocp_qp_sol->t;

	struct STRVEC *tmp_nuxM = cond_ws->tmp_nuxM;
	struct STRVEC *tmp_ngM = cond_ws->tmp_ngM;
	int *idxs_rev = cond_ws->idxs_rev;

	// problem size

	int nbb = nb[0]; // box that remain box constraints
	int nbg = 0; // box that becomes general constraints
	for(ii=1; ii<=N; ii++)
		for(jj=0; jj<nb[ii]; jj++)
			if(idxb[ii][jj]<nu[ii])
				nbb++;
			else
				nbg++;
	
	int nx2 = nx[0];
	int nu2 = nu[0];
	int ngg = ng[0];
	int ns2 = ns[0];
	for(ii=1; ii<=N; ii++)
		{
		nu2 += nu[ii];
		ngg += ng[ii];
		ns2 += ns[ii];
		}
	int ng2 = nbg + ngg;
	int nb2 = nbb;
	int nt2 = nb2 + ng2;

	// inputs & initial states
	int nu_tmp = 0;
	// final stages: copy only input
	for(ii=0; ii<N; ii++)
		{
		VECCP_LIBSTR(nu[N-ii], vc, nu_tmp, ux+(N-ii), 0);
		nu_tmp += nu[N-ii];
		}
	// first stage: copy input and state
	VECCP_LIBSTR(nu[0]+nx[0], vc, nu_tmp, ux+0, 0);

	// compute missing states by simulation within each block
	for(ii=0; ii<N; ii++)
		{
		GEMV_T_LIBSTR(nu[ii]+nx[ii], nx[ii+1], 1.0, BAbt+ii, 0, 0, ux+ii, 0, 1.0, b+ii, 0, ux+(ii+1), nu[ii+1]);
		}

	// slack variables and ineq lagrange multipliers
	nbb = 0;
	nbg = 0;
	ngg = 0;
	REAL *ptr_lam;
	REAL *ptr_lam_lb;
	REAL *ptr_lam_ub;
	REAL *ptr_lam_lg;
	REAL *ptr_lam_ug;
	REAL *ptr_t;
	REAL *ptr_t_lb;
	REAL *ptr_t_ub;
	REAL *ptr_t_lg;
	REAL *ptr_t_ug;
	REAL *ptr_lamc = lamc->pa;
	REAL *ptr_lam_lbc = lamc->pa+0;
	REAL *ptr_lam_ubc = lamc->pa+nb2+ng2;
	REAL *ptr_lam_lgc = lamc->pa+nb2;
	REAL *ptr_lam_ugc = lamc->pa+2*nb2+ng2;
	REAL *ptr_tc = tc->pa;
	REAL *ptr_t_lbc = tc->pa+0;
	REAL *ptr_t_ubc = tc->pa+nb2+ng2;
	REAL *ptr_t_lgc = tc->pa+nb2;
	REAL *ptr_t_ugc = tc->pa+2*nb2+ng2;
	REAL *ptr_ux;
	REAL *ptr_vc = vc->pa;
	int is = 0;
	int nx0, nu0, nb0, ng0, ns0;
	int idxb0, idxg0;

	// box constraints
	// final stages
	for(ii=0; ii<N; ii++)
		{
		nx0 = nx[N-ii];
		nu0 = nu[N-ii];
		nb0 = nb[N-ii];
		ng0 = ng[N-ii];
		ns0 = ns[N-ii];
		for(jj=0; jj<nb0+ng0; jj++)
			idxs_rev[jj] = -1;
		if(ns0>0)
			{
			for(jj=0; jj<ns0; jj++)
				idxs_rev[idxs[N-ii][jj]] = jj;
			}
		ptr_ux = (ux+N-ii)->pa;
		ptr_lam = (lam+N-ii)->pa;
		ptr_t = (t+N-ii)->pa;
		ptr_lam_lb = (lam+N-ii)->pa+0;
		ptr_lam_ub = (lam+N-ii)->pa+nb0+ng0;
		ptr_t_lb = (t+N-ii)->pa+0;
		ptr_t_ub = (t+N-ii)->pa+nb0+ng0;
		for(jj=0; jj<nb0; jj++)
			{
			idxb0 = idxb[N-ii][jj];
			if(idxb0<nu0)
				{
				// box as box
				ptr_lam_lb[jj] = ptr_lam_lbc[nbb];
				ptr_lam_ub[jj] = ptr_lam_ubc[nbb];
				ptr_t_lb[jj] = ptr_t_lbc[nbb];
				ptr_t_ub[jj] = ptr_t_ubc[nbb];
				nbb++;
				}
			else
				{
				// box as general XXX change when decide where nbg are placed wrt ng
				ptr_lam_lb[jj] = ptr_lam_lgc[nbg];
				ptr_lam_ub[jj] = ptr_lam_ugc[nbg];
				ptr_t_lb[jj] = ptr_t_lgc[nbg];
				ptr_t_ub[jj] = ptr_t_ugc[nbg];
				nbg++;
				}
			if(idxs_rev[jj]>=0)
				{
				ptr_lam[2*nb0+2*ng0+idxs_rev[jj]]     = ptr_lamc[2*nb2+2*ng2+is];
				ptr_lam[2*nb0+2*ng0+ns0+idxs_rev[jj]] = ptr_lamc[2*nb2+2*ng2+ns2+is];
				ptr_t[2*nb0+2*ng0+idxs_rev[jj]]     = ptr_tc[2*nb2+2*ng2+is];
				ptr_t[2*nb0+2*ng0+ns0+idxs_rev[jj]] = ptr_tc[2*nb2+2*ng2+ns2+is];
				ptr_ux[nu0+nx0+idxs_rev[jj]] = ptr_vc[nu2+nx2+is];
				ptr_ux[nu0+nx0+ns0+idxs_rev[jj]] = ptr_vc[nu2+nx2+ns2+is];
				is++;
				}
			}
		}
	// first stage
	// all box as box
	nx0 = nx[0];
	nu0 = nu[0];
	nb0 = nb[0];
	ng0 = ng[0];
	ns0 = ns[0];
	VECCP_LIBSTR(nb0, lamc, 0+nbb, lam+0, 0);
	VECCP_LIBSTR(nb0, lamc, nb2+ng2+nbb, lam+0, nb0+ng0);
	VECCP_LIBSTR(nb0, tc, 0+nbb, t+0, 0);
	VECCP_LIBSTR(nb0, tc, nb2+ng2+nbb, t+0, nb0+ng0);
	nbb += nb0;
	for(jj=0; jj<nb0+ng0; jj++)
		idxs_rev[jj] = -1;
	if(ns0>0)
		{
		for(jj=0; jj<ns0; jj++)
			idxs_rev[idxs[0][jj]] = jj;
		}
	ptr_ux = (ux+0)->pa;
	ptr_lam = (lam+0)->pa;
	ptr_t = (t+0)->pa;
	for(jj=0; jj<nb0; jj++)
		{
//		idxb0 = idxb[0][jj];
		if(idxs_rev[jj]>=0)
			{
			ptr_lam[2*nb0+2*ng0+idxs_rev[jj]]     = ptr_lamc[2*nb2+2*ng2+is];
			ptr_lam[2*nb0+2*ng0+ns0+idxs_rev[jj]] = ptr_lamc[2*nb2+2*ng2+ns2+is];
			ptr_t[2*nb0+2*ng0+idxs_rev[jj]]     = ptr_tc[2*nb2+2*ng2+is];
			ptr_t[2*nb0+2*ng0+ns0+idxs_rev[jj]] = ptr_tc[2*nb2+2*ng2+ns2+is];
			ptr_ux[nu0+nx0+idxs_rev[jj]] = ptr_vc[nu2+nx2+is];
			ptr_ux[nu0+nx0+ns0+idxs_rev[jj]] = ptr_vc[nu2+nx2+ns2+is];
			is++;
			}
		}

	// general constraints
	// process as vectors ???
	for(ii=0; ii<N; ii++)
		{
		nx0 = nx[N-ii];
		nu0 = nu[N-ii];
		nb0 = nb[N-ii];
		ng0 = ng[N-ii];
		ns0 = ns[N-ii];
		for(jj=0; jj<nb0+ng0; jj++)
			idxs_rev[jj] = -1;
		if(ns0>0)
			{
			for(jj=0; jj<ns0; jj++)
				idxs_rev[idxs[N-ii][jj]] = jj;
			}
		ptr_ux = (ux+(N-ii))->pa;
		ptr_lam = (lam+(N-ii))->pa;
		ptr_t = (t+(N-ii))->pa;
		ptr_lam_lg = (lam+(N-ii))->pa+nb0;
		ptr_lam_ug = (lam+(N-ii))->pa+2*nb0+ng0;
		ptr_t_lg = (t+(N-ii))->pa+nb0;
		ptr_t_ug = (t+(N-ii))->pa+2*nb0+ng0;
		for(jj=0; jj<ng0; jj++)
			{
			// genenral as general
			ptr_lam_lg[jj] = ptr_lam_lgc[nbg+ngg];
			ptr_lam_ug[jj] = ptr_lam_ugc[nbg+ngg];
			ptr_t_lg[jj] = ptr_t_lgc[nbg+ngg];
			ptr_t_ug[jj] = ptr_t_ugc[nbg+ngg];
			ngg++;

			idxg0 = nb0+jj;
			if(idxs_rev[idxg0]>=0)
				{
				ptr_lam[2*nb0+2*ng0+idxs_rev[idxg0]]     = ptr_lamc[2*nb2+2*ng2+is];
				ptr_lam[2*nb0+2*ng0+ns0+idxs_rev[idxg0]] = ptr_lamc[2*nb2+2*ng2+ns2+is];
				ptr_t[2*nb0+2*ng0+idxs_rev[idxg0]]     = ptr_tc[2*nb2+2*ng2+is];
				ptr_t[2*nb0+2*ng0+ns0+idxs_rev[idxg0]] = ptr_tc[2*nb2+2*ng2+ns2+is];
				ptr_ux[nu0+nx0+idxs_rev[idxg0]] = ptr_vc[nu2+nx2+is];
				ptr_ux[nu0+nx0+ns0+idxs_rev[idxg0]] = ptr_vc[nu2+nx2+ns2+is];
				is++;
				}
			}
		}
	// first stage
	// all general as general
	nx0 = nx[0];
	nu0 = nu[0];
	nb0 = nb[0];
	ng0 = ng[0];
	ns0 = ns[0];
	VECCP_LIBSTR(ng[0], lamc, nb2+nbg+ngg, lam+0, nb0);
	VECCP_LIBSTR(ng[0], lamc, 2*nb2+ng2+nbg+ngg, lam+0, 2*nb0+ng0);
	VECCP_LIBSTR(ng[0], tc, nb2+nbg+ngg, t+0, nb0);
	VECCP_LIBSTR(ng[0], tc, 2*nb2+ng2+nbg+ngg, t+0, 2*nb0+ng0);
	ngg += ng0;
	for(jj=0; jj<nb0+ng0; jj++)
		idxs_rev[jj] = -1;
	if(ns0>0)
		{
		for(jj=0; jj<ns0; jj++)
			idxs_rev[idxs[0][jj]] = jj;
		}
	ptr_ux = (ux+(0))->pa;
	ptr_lam = (lam+(0))->pa;
	ptr_t = (t+(0))->pa;
	for(jj=0; jj<ng0; jj++)
		{
		idxg0 = nb0+jj;
		if(idxs_rev[idxg0]>=0)
			{
			ptr_lam[2*nb0+2*ng0+idxs_rev[idxg0]]     = ptr_lamc[2*nb2+2*ng2+is];
			ptr_lam[2*nb0+2*ng0+ns0+idxs_rev[idxg0]] = ptr_lamc[2*nb2+2*ng2+ns2+is];
			ptr_t[2*nb0+2*ng0+idxs_rev[idxg0]]     = ptr_tc[2*nb2+2*ng2+is];
			ptr_t[2*nb0+2*ng0+ns0+idxs_rev[idxg0]] = ptr_tc[2*nb2+2*ng2+ns2+is];
			ptr_ux[nu0+nx0+idxs_rev[idxg0]] = ptr_vc[nu2+nx2+is];
			ptr_ux[nu0+nx0+ns0+idxs_rev[idxg0]] = ptr_vc[nu2+nx2+ns2+is];
			is++;
			}
		}

	// lagrange multipliers of equality constraints
	// TODO avoid to multiply by R and B (i.e. the u part)
	REAL *ptr_nuxM = tmp_nuxM->pa;
	REAL *ptr_ngM = tmp_ngM->pa;
	// last stage
	if(cond_ws->cond_last_stage==0)
		VECCP_LIBSTR(nx[Np], pic, 0, pi+Np-1, 0);
	else
		{
//		SYMV_L_LIBSTR(nx[Np], nx[Np], 1.0, RSQrq+Np, nu[Np], nu[Np], ux+Np, nu[Np], 1.0, rq+Np, nu[Np], pi+Np-1, 0);
		ptr_lam_lb = (lam+Np)->pa+0;
		ptr_lam_ub = (lam+Np)->pa+nb[Np]+ng[Np];
		ptr_lam_lg = (lam+Np)->pa+nb[Np];
		ptr_lam_ug = (lam+Np)->pa+2*nb[Np]+ng[Np];
		VECCP_LIBSTR(nu[Np]+nx[Np], rq+(Np), 0, tmp_nuxM, 0);
		for(jj=0; jj<nb[Np]; jj++)
			ptr_nuxM[idxb[Np][jj]] += ptr_lam_ub[jj] - ptr_lam_lb[jj];
		SYMV_L_LIBSTR(nu[Np]+nx[Np], nu[Np]+nx[Np], 1.0, RSQrq+(Np), 0, 0, ux+(Np), 0, 1.0, tmp_nuxM, 0, tmp_nuxM, 0);
		for(jj=0; jj<ng[Np]; jj++)
			ptr_ngM[jj] = ptr_lam_ug[jj] - ptr_lam_lg[jj];
		GEMV_N_LIBSTR(nu[Np]+nx[Np], ng[Np], 1.0, DCt+(Np), 0, 0, tmp_ngM, 0, 1.0, tmp_nuxM, 0, tmp_nuxM, 0);

		VECCP_LIBSTR(nx[Np], tmp_nuxM, nu[Np], pi+(Np-1), 0);
		}

	for(ii=0; ii<Np-1; ii++)
		{
		ptr_lam_lb = (lam+Np-1-ii)->pa+0;
		ptr_lam_ub = (lam+Np-1-ii)->pa+nb[Np-1-ii]+ng[Np-1-ii];
		ptr_lam_lg = (lam+Np-1-ii)->pa+nb[Np-1-ii];
		ptr_lam_ug = (lam+Np-1-ii)->pa+2*nb[Np-1-ii]+ng[Np-1-ii];
		VECCP_LIBSTR(nu[Np-1-ii]+nx[Np-1-ii], rq+(Np-1-ii), 0, tmp_nuxM, 0);
		for(jj=0; jj<nb[Np-1-ii]; jj++)
			ptr_nuxM[idxb[Np-1-ii][jj]] += ptr_lam_ub[jj] - ptr_lam_lb[jj];
		SYMV_L_LIBSTR(nu[Np-1-ii]+nx[Np-1-ii], nu[Np-1-ii]+nx[Np-1-ii], 1.0, RSQrq+(Np-1-ii), 0, 0, ux+(Np-1-ii), 0, 1.0, tmp_nuxM, 0, tmp_nuxM, 0);
		GEMV_N_LIBSTR(nu[Np-1-ii]+nx[Np-1-ii], nx[Np-ii], 1.0, BAbt+(Np-1-ii), 0, 0, pi+(Np-1-ii), 0, 1.0, tmp_nuxM, 0, tmp_nuxM, 0);
		for(jj=0; jj<ng[Np-1-ii]; jj++)
			ptr_ngM[jj] = ptr_lam_ug[jj] - ptr_lam_lg[jj];
		GEMV_N_LIBSTR(nu[Np-1-ii]+nx[Np-1-ii], ng[Np-1-ii], 1.0, DCt+(Np-1-ii), 0, 0, tmp_ngM, 0, 1.0, tmp_nuxM, 0, tmp_nuxM, 0);

		VECCP_LIBSTR(nx[Np-1-ii], tmp_nuxM, nu[Np-1-ii], pi+(Np-2-ii), 0);
		}
	

	return;

	}



void EXPAND_PRIMAL_SOL(struct OCP_QP *ocp_qp, struct DENSE_QP_SOL *dense_qp_sol, struct OCP_QP_SOL *ocp_qp_sol, struct COND_QP_OCP2DENSE_WORKSPACE *cond_ws)
	{

	int N = ocp_qp->dim->N;
	int Np = N;
	if(cond_ws->cond_last_stage==0)
		N -= 1;

	int ii, jj;

	int *nu = ocp_qp->dim->nu;
	int *nx = ocp_qp->dim->nx;
	int *nb = ocp_qp->dim->nb;
	int *ng = ocp_qp->dim->ng;
	int *ns = ocp_qp->dim->ns;

	struct STRMAT *BAbt = ocp_qp->BAbt;
	struct STRVEC *b = ocp_qp->b;
	int **idxb = ocp_qp->idxb;
	struct STRMAT *RSQrq = ocp_qp->RSQrq;
	struct STRVEC *rq = ocp_qp->rq;
	struct STRMAT *DCt = ocp_qp->DCt;
	int **idxs = ocp_qp->idxs;

	struct STRVEC *vc = dense_qp_sol->v;
	struct STRVEC *pic = dense_qp_sol->pi;
	struct STRVEC *lamc = dense_qp_sol->lam;
	struct STRVEC *tc = dense_qp_sol->t;

	struct STRVEC *ux = ocp_qp_sol->ux;
	struct STRVEC *pi = ocp_qp_sol->pi;
	struct STRVEC *lam = ocp_qp_sol->lam;
	struct STRVEC *t = ocp_qp_sol->t;

	struct STRVEC *tmp_nuxM = cond_ws->tmp_nuxM;
	struct STRVEC *tmp_ngM = cond_ws->tmp_ngM;
	int *idxs_rev = cond_ws->idxs_rev;

	// problem size

	int nbb = nb[0]; // box that remain box constraints
	int nbg = 0; // box that becomes general constraints
	for(ii=1; ii<=N; ii++)
		for(jj=0; jj<nb[ii]; jj++)
			if(idxb[ii][jj]<nu[ii])
				nbb++;
			else
				nbg++;
	
	int nx2 = nx[0];
	int nu2 = nu[0];
	int ngg = ng[0];
	int ns2 = ns[0];
	for(ii=1; ii<=N; ii++)
		{
		nu2 += nu[ii];
		ngg += ng[ii];
		ns2 += ns[ii];
		}
	int ng2 = nbg + ngg;
	int nb2 = nbb;
	int nt2 = nb2 + ng2;

	// inputs & initial states
	int nu_tmp = 0;
	// final stages: copy only input
	for(ii=0; ii<N; ii++)
		{
		VECCP_LIBSTR(nu[N-ii], vc, nu_tmp, ux+(N-ii), 0);
		nu_tmp += nu[N-ii];
		}
	// first stage: copy input and state
	VECCP_LIBSTR(nu[0]+nx[0], vc, nu_tmp, ux+0, 0);

	// compute missing states by simulation within each block
	for(ii=0; ii<N; ii++)
		{
		GEMV_T_LIBSTR(nu[ii]+nx[ii], nx[ii+1], 1.0, BAbt+ii, 0, 0, ux+ii, 0, 1.0, b+ii, 0, ux+(ii+1), nu[ii+1]);
		}

	return;

	}
