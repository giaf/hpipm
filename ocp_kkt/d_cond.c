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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_d_blas.h>
#include <blasfeo_d_aux.h>

#include "../include/hpipm_ocp_kkt.h"



void d_cond_BAbt(int N, struct d_ocp_qp *str_in, int idx_in, struct d_ocp_qp *str_out, int idx_out, struct d_cond_mem *str_mem, int idx_mem)
	{

	// early return
	if(N<0)
		return;

	// extract input members
	int *nx = str_in->nx + idx_in;
	int *nu = str_in->nu + idx_in;
	struct d_strmat *sBAbt = str_in->sBAbt + idx_in;

	// extract output members
	struct d_strmat *sBAbt2 = str_out->sBAbt + idx_out;
	struct d_strvec *sb2 = str_out->sb + idx_out;

	// extract memory members
	struct d_strmat *sGamma = str_mem->sGamma + idx_mem;

	int ii, jj;

	int nu_tmp;

	nu_tmp = 0;
	ii = 0;
	// B & A & b
	dgecp_libstr(nu[0]+nx[0]+1, nx[1], 1.0, &sBAbt[0], 0, 0, &sGamma[0], 0, 0);
	//
	nu_tmp += nu[0];
	ii++;

	for(ii=1; ii<N; ii++)
		{
		// TODO check for equal pointers and avoid copy

		// Gamma * A^T
		dgemm_nn_libstr(nu_tmp+nx[0]+1, nx[ii+1], nx[ii], 1.0, &sGamma[ii-1], 0, 0, &sBAbt[ii], nu[ii], 0, 0.0, &sGamma[ii], nu[ii], 0, &sGamma[ii], nu[ii], 0); // Gamma * A^T

		dgecp_libstr(nu[ii], nx[ii+1], 1.0, &sBAbt[ii], 0, 0, &sGamma[ii], 0, 0);

		nu_tmp += nu[ii];

		dgead_libstr(1, nx[ii+1], 1.0, &sBAbt[ii], nu[ii]+nx[ii], 0, &sGamma[ii], nu_tmp+nx[0], 0);
		}
	
	// B & A & b
	dgecp_libstr(nu_tmp+nx[0]+1, nx[N], 1.0, &sGamma[N-1], 0, 0, &sBAbt2[0], 0, 0);
	// b
	drowex_libstr(nx[N], 1.0, &sBAbt2[0], 0, 0, &sb2[0], 0);

	return;

	}



// workspace_sizes[0] : Lx
// workspace_sizes[1] : BAbtL
int d_cond_RSQrq_N2nx3_workspace_size(int N, struct d_ocp_qp *str_in, int idx_in, struct d_cond_work *str_work, int idx_work)
	{

	// early return
	if(N<0)
		return 0;

	// extract input members
	int *nx = str_in->nx + idx_in;
	int *nu = str_in->nu + idx_in;

	// extract workspace members
	int *sizes = str_work->cond_RSQrq_N2nx3_sizes + idx_work;

	int ii;
	int tmp;

	int size = 0;

	// Lx
	sizes[0] = 0;
	for(ii=1; ii<=N; ii++)
		{
		tmp = d_size_strmat(nx[ii]+1, nx[ii]);
		sizes[0] = tmp>sizes[0] ? tmp : sizes[0];
		}

	// BAbtL
	sizes[1] = 0;
	for(ii=0; ii<N; ii++)
		{
		tmp = d_size_strmat(nu[ii]+nx[ii]+1, nx[ii+1]);
		sizes[1] = tmp>sizes[1] ? tmp : sizes[1];
		}
	
	size += sizes[0];
	size += sizes[1];

	return size;
	}



void d_cond_RSQrq_N2nx3(int N, struct d_ocp_qp *str_in, int idx_in, struct d_ocp_qp *str_out, int idx_out, struct d_cond_mem *str_mem, int idx_mem, struct d_cond_work *str_work, int idx_work, void *workspace)
	{

	// early return
	if(N<0)
		return;

	// extract input members
	int *nx = str_in->nx + idx_in;
	int *nu = str_in->nu + idx_in;
	struct d_strmat *sBAbt = str_in->sBAbt + idx_in;
	struct d_strmat *sRSQrq = str_in->sRSQrq + idx_in;

	// extract output members
	struct d_strmat *sRSQrq2 = str_out->sRSQrq + idx_out;
	struct d_strvec *srq2 = str_out->srq + idx_out;

	// extract memory members
	struct d_strmat *sGamma = str_mem->sGamma + idx_mem;
	struct d_strmat *sL = str_mem->sL + idx_mem;

	// extract workspace members
	int *sizes = str_work->cond_RSQrq_N2nx3_sizes;

	// declare workspace matrices
	struct d_strmat sLx;
	struct d_strmat sBAbtL;

	// early return
	if(N==0)
		{
		dgecp_libstr(nu[0]+nx[0]+1, nu[0]+nx[0], 1.0, &sRSQrq[0], 0, 0, &sRSQrq2[0], 0, 0);
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

	dgecp_libstr(nu[N]+nx[N]+1, nu[N]+nx[N], 1.0, &sRSQrq[N], 0, 0, &sL[N], 0, 0);

	// D
	dtrcp_l_libstr(nu[N], 1.0, &sL[N], 0, 0, &sRSQrq2[0], nuf, nuf);

	dgemm_nn_libstr(nub+nx[0]+1, nu[N], nx[N], 1.0, &sGamma[N-1], 0, 0, &sL[N], nu[N], 0, 0.0, &sRSQrq2[0], nuf+nu[N], nuf, &sRSQrq2[0], nuf+nu[N], nuf);

	// m
	dgead_libstr(1, nu[N], 1.0, &sL[N], nu[N]+nx[N], 0, &sRSQrq2[0], nu2+nx[0], nuf);

	nuf += nu[N];



	// middle stages 
	for(nn=0; nn<N-1; nn++)
		{	
		nub -= nu[N-nn-1];

		d_create_strmat(nx[N-nn]+1, nx[N-nn], &sLx, workspace+0);
		d_create_strmat(nu[N-nn-1]+nx[N-nn-1]+1, nx[N-nn], &sBAbtL, workspace+sizes[0]);

#if defined(LA_HIGH_PERFORMANCE)
		dgecp_libstr(nx[N-nn]+1, nx[N-nn], 1.0, &sL[N-nn], nu[N-nn], nu[N-nn], &sLx, 0, 0);

		dpotrf_l_mn_libstr(nx[N-nn]+1, nx[N-nn], &sLx, 0, 0, &sLx, 0, 0);

		dtrmm_rlnn_libstr(nu[N-nn-1]+nx[N-nn-1]+1, nx[N-nn], 1.0, &sLx, 0, 0, &sBAbt[N-nn-1], 0, 0, &sBAbtL, 0, 0);
#else
		dpotrf_l_mn_libstr(nx[N-nn]+1, nx[N-nn], &sL[N-nn], nu[N-nn], nu[N-nn], &sLx, 0, 0);

		dtrmm_rlnn_libstr(nu[N-nn-1]+nx[N-nn-1]+1, nx[N-nn], 1.0, &sLx, 0, 0, &sBAbt[N-nn-1], 0, 0, &sBAbtL, 0, 0);
#endif
		dgead_libstr(1, nx[N-nn], 1.0, &sLx, nx[N-nn], 0, &sBAbtL, nu[N-nn-1]+nx[N-nn-1], 0);

		dsyrk_ln_mn_libstr(nu[N-nn-1]+nx[N-nn-1]+1, nu[N-nn-1]+nx[N-nn-1], nx[N-nn], 1.0, &sBAbtL, 0, 0, &sBAbtL, 0, 0, 1.0, &sRSQrq[N-nn-1], 0, 0, &sL[N-nn-1], 0, 0);

		// D
		dtrcp_l_libstr(nu[N-nn-1], 1.0, &sL[N-nn-1], 0, 0, &sRSQrq2[0], nuf, nuf);

		dgemm_nn_libstr(nub+nx[0]+1, nu[N-nn-1], nx[N-nn-1], 1.0, &sGamma[N-nn-2], 0, 0, &sL[N-nn-1], nu[N-nn-1], 0, 0.0, &sRSQrq2[0], nuf+nu[N-nn-1], nuf, &sRSQrq2[0], nuf+nu[N-nn-1], nuf);

		// m
		dgead_libstr(1, nu[N-nn-1], 1.0, &sL[N-nn-1], nu[N-nn-1]+nx[N-nn-1], 0, &sRSQrq2[0], nu2+nx[0], nuf);

		nuf += nu[N-nn-1];

		}

	// first stage
	nn = N-1;

	d_create_strmat(nx[N-nn]+1, nx[N-nn], &sLx, workspace+0);
	d_create_strmat(nu[N-nn-1]+nx[N-nn-1]+1, nx[N-nn], &sBAbtL, workspace+sizes[0]);
	
#if defined(LA_HIGH_PERFORMANCE)
	dgecp_libstr(nx[N-nn]+1, nx[N-nn], 1.0, &sL[N-nn], nu[N-nn], nu[N-nn], &sLx, 0, 0);

	dpotrf_l_mn_libstr(nx[N-nn]+1, nx[N-nn], &sLx, 0, 0, &sLx, 0, 0);

	dtrmm_rlnn_libstr(nu[N-nn-1]+nx[N-nn-1]+1, nx[N-nn], 1.0, &sLx, 0, 0, &sBAbt[N-nn-1], 0, 0, &sBAbtL, 0, 0);
#else
	dpotrf_l_mn_libstr(nx[N-nn]+1, nx[N-nn], &sL[N-nn], nu[N-nn], nu[N-nn], &sLx, 0, 0);

	dtrmm_rlnn_libstr(nu[N-nn-1]+nx[N-nn-1]+1, nx[N-nn], 1.0, &sLx, 0, 0, &sBAbt[N-nn-1], 0, 0, &sBAbtL, 0, 0);
#endif
	dgead_libstr(1, nx[N-nn], 1.0, &sLx, nx[N-nn], 0, &sBAbtL, nu[N-nn-1]+nx[N-nn-1], 0);

	dsyrk_ln_mn_libstr(nu[N-nn-1]+nx[N-nn-1]+1, nu[N-nn-1]+nx[N-nn-1], nx[N-nn], 1.0, &sBAbtL, 0, 0, &sBAbtL, 0, 0, 1.0, &sRSQrq[N-nn-1], 0, 0, &sL[N-nn-1], 0, 0);

	// D, M, m, P, p
//	dgecp_libstr(nu[0]+nx[0]+1, nu[0]+nx[0], 1.0, &sL[N-nn-1], 0, 0, &sRSQrq2[0], nuf, nuf); // TODO dtrcp for 'rectangular' matrices
	dtrcp_l_libstr(nu[0]+nx[0], 1.0, &sL[N-nn-1], 0, 0, &sRSQrq2[0], nuf, nuf); // TODO dtrcp for 'rectangular' matrices
	dgecp_libstr(1, nu[0]+nx[0], 1.0, &sL[N-nn-1], nu[0]+nx[0], 0, &sRSQrq2[0], nuf+nu[0]+nx[0], nuf); // TODO dtrcp for 'rectangular' matrices
	// m p
	drowex_libstr(nu[0]+nx[0], 1.0, &sRSQrq2[0], 0, 0, &srq2[0], 0);

	return;

	}



int d_cond_DCtd_workspace_size(int N, struct d_ocp_qp *str_in, int idx_in)
	{

	// early return
	if(N<0)
		return 0;

	// extract input members
	int *nx = str_in->nx + idx_in;
	int *ng = str_in->ng + idx_in;

	int ii;
	int tmp;

	int size = 0;

	for(ii=0; ii<=N; ii++)
		{
		tmp = d_size_strmat(ng[ii], nx[ii]);
		tmp += d_size_strvec(nx[ii]);
		tmp += d_size_strvec(ng[ii]);
		size = tmp > size ? tmp : size;
		}

	return size;

	}



void d_cond_DCtd(int N, struct d_ocp_qp *str_in, int idx_in, struct d_ocp_qp *str_out, int idx_out, struct d_cond_mem *str_mem, int idx_mem, void *workspace)
	{

	// early return
	if(N<0)
		return;

	// extract input members
	int *nx = str_in->nx + idx_in;
	int *nu = str_in->nu + idx_in;
	int *nb = str_in->nb + idx_in;
	int *ng = str_in->ng + idx_in;
	int **idxb = str_in->idxb + idx_in;
	struct d_strmat *sDCt = str_in->sDCt + idx_in;
	struct d_strvec *sd = str_in->sd + idx_in;

	// extract output members
	struct d_strmat *sDCt2 = str_out->sDCt + idx_out;
	struct d_strvec *sd2 = str_out->sd + idx_out;
	int **idxb2 = str_out->idxb + idx_out;

	// extract memory members
	struct d_strmat *sGamma = str_mem->sGamma + idx_mem;


	double *d2 = sd2->pa; // TODO use VECEL instead !!!
	double *ptr_d;
	
	int nu_tmp, ng_tmp;

	int ii, jj;

	int nu0, nx0, nb0, ng0, nt0;

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
	for(ii=1; ii<=N; ii++)
		{
		nu2 += nu[ii];
		ngg += ng[ii];
		}
	int ng2 = nbg + ngg;
	int nb2 = nbb;
	int nt2 = nb2 + ng2;

	// set constraint matrix to zero (it's 2 lower triangular matrices atm)
	dgese_libstr(nu2+nx2, ng2, 0.0, &sDCt2[0], 0, 0);

	// box constraints

	int idx_gammab = nx[0];
	for(ii=0; ii<N; ii++)
		idx_gammab += nu[ii];

	int ib = 0;
	int ig = 0;

	double tmp;
	int idx_g;

	// middle stages
	nu_tmp = 0;
	for(ii=0; ii<N; ii++)
		{
		nu0 = nu[N-ii];
		nb0 = nb[N-ii];
		ng0 = ng[N-ii];
		nt0 = nb0 + ng0;
		nu_tmp += nu0;
		ptr_d = sd[N-ii].pa;
		for(jj=0; jj<nb0; jj++)
			{
			if(idxb[N-ii][jj]<nu0) // input: box constraint
				{
				d2[0*nt2+ib] = ptr_d[0*nt0+jj];
				d2[1*nt2+ib] = ptr_d[1*nt0+jj];
				idxb2[0][ib] = nu_tmp - nu0 + idxb[N-ii][jj];
				ib++;
				}
			else // state: general constraint
				{
				idx_g = idxb[N-ii][jj]-nu0;
				tmp = dgeex1_libstr(&sGamma[N-1-ii], idx_gammab, idx_g);
				d2[nb2+0*nt2+ig] = ptr_d[0*nt0+jj] - tmp;
				d2[nb2+1*nt2+ig] = ptr_d[1*nt0+jj] - tmp;
				dgecp_libstr(idx_gammab, 1, 1.0, &sGamma[N-ii-1], 0, idx_g, &sDCt2[0], nu_tmp, ig);
				ig++;
				}
			}
		idx_gammab -= nu[N-1-ii];
		}

	// initial stage: both inputs and states as box constraints
	nu0 = nu[0];
	nb0 = nb[0];
	ng0 = ng[0];
	nt0 = nb0 + ng0;
	nu_tmp += nu0;
	ptr_d = sd[0].pa;
	for(jj=0; jj<nb0; jj++)
		{
		d2[0*nt2+ib] = ptr_d[0*nt0+jj];
		d2[1*nt2+ib] = ptr_d[1*nt0+jj];
		idxb2[0][ib] = nu_tmp - nu0 + idxb[0][jj];
		ib++;
		}

	// XXX for now, just shift after box-to-general constraints
	// better interleave them, to keep the block lower trianlgular structure !!!

	// general constraints

	struct d_strmat sC;
	struct d_strvec sGammab;
	struct d_strvec sCGammab;

	char *c_ptr;

	nu_tmp = 0;
	ng_tmp = 0;
	for(ii=0; ii<N; ii++)
		{

		nx0 = nx[N-ii];
		nu0 = nu[N-ii];
		ng0 = ng[N-ii];
		nt0 = nb0 + ng0;

		if(ng0>0)
			{

			c_ptr = (char *) workspace;

			dgecp_libstr(nu0, ng0, 1.0, &sDCt[N-ii], 0, 0, &sDCt2[0], nu_tmp, nbg+ng_tmp);

			nu_tmp += nu0;

			d_create_strmat(ng0, nx0, &sC, (void *) c_ptr);
			c_ptr += sC.memory_size;

			dgetr_libstr(nx0, ng0, 1.0, &sDCt[N-ii], nu0, 0, &sC, 0, 0);

			dgemm_nt_libstr(nu2+nx[0]-nu_tmp, ng0, nx0, 1.0, &sGamma[N-1-ii], 0, 0, &sC, 0, 0, 0.0, &sDCt2[0], nu_tmp, nbg+ng_tmp, &sDCt2[0], nu_tmp, nbg+ng_tmp);

			dveccp_libstr(ng0, 1.0, &sd[N-ii], nb0+0*nt0, sd2, nb2+0*nt2+nbg+ng_tmp);
			dveccp_libstr(ng0, 1.0, &sd[N-ii], nb0+1*nt0, sd2, nb2+1*nt2+nbg+ng_tmp);

			d_create_strvec(nx0, &sGammab, (void *) c_ptr);
			c_ptr += sGammab.memory_size;
			d_create_strvec(ng0, &sCGammab, (void *) c_ptr);
			c_ptr += sCGammab.memory_size;

			drowex_libstr(nx0, 1.0, &sGamma[N-1-ii], nu2+nx[0]-nu_tmp, 0, &sGammab, 0);

			dgemv_n_libstr(ng0, nx0, 1.0, &sC, 0, 0, &sGammab, 0, 0.0, &sCGammab, 0, &sCGammab, 0);

			daxpy_libstr(ng0, -1.0, &sCGammab, 0, sd2, nb2+0*nt2+nbg+ng_tmp, sd2, nb2+0*nt2+nbg+ng_tmp);
			daxpy_libstr(ng0, -1.0, &sCGammab, 0, sd2, nb2+1*nt2+nbg+ng_tmp, sd2, nb2+1*nt2+nbg+ng_tmp);

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
	ng0 = ng[0];
	nt0 = nb0 + ng0;

	if(ng0>0)
		{

		dgecp_libstr(nu0+nx0, ng0, 1.0, &sDCt[0], 0, 0, &sDCt2[0], nu_tmp, nbg+ng_tmp);

		dveccp_libstr(ng0, 1.0, &sd[0], nb0+0*nt0, sd2, nb2+0*nt2+nbg+ng_tmp);
		dveccp_libstr(ng0, 1.0, &sd[0], nb0+1*nt0, sd2, nb2+1*nt2+nbg+ng_tmp);

//		ng_tmp += ng[N-ii];

		}

	return;

	}




