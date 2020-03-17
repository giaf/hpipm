/**************************************************************************************************
*                                                                                                 *
* This file is part of HPIPM.                                                                     *
*                                                                                                 *
* HPIPM -- High-Performance Interior Point Method.                                                *
* Copyright (C) 2019 by Gianluca Frison.                                                          *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
*                                                                                                 *
* The 2-Clause BSD License                                                                        *
*                                                                                                 *
* Redistribution and use in source and binary forms, with or without                              *
* modification, are permitted provided that the following conditions are met:                     *
*                                                                                                 *
* 1. Redistributions of source code must retain the above copyright notice, this                  *
*    list of conditions and the following disclaimer.                                             *
* 2. Redistributions in binary form must reproduce the above copyright notice,                    *
*    this list of conditions and the following disclaimer in the documentation                    *
*    and/or other materials provided with the distribution.                                       *
*                                                                                                 *
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND                 *
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED                   *
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE                          *
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR                 *
* ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES                  *
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;                    *
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND                     *
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT                      *
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS                   *
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                                    *
*                                                                                                 *
* Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             *
*                                                                                                 *
**************************************************************************************************/



int OCP_QP_REDUCE_EQ_DOF_WORK_MEMSIZE(struct OCP_QP_DIM *dim)
	{

	int ii;

	// extract dim
	int N = dim->N;
	int *nu = dim->nu;
	int *nx = dim->nx;
	int *nb = dim->nb;
	int *ng = dim->ng;

	int nuxM = nu[0]+nx[0];
	int nbgM = nb[0]+ng[0];
	for(ii=1; ii<=N; ii++)
		{
		nuxM = nu[ii]+nx[ii]>nuxM ? nu[ii]+nx[ii] : nuxM;
		nbgM = nb[ii]+ng[ii]>nbgM ? nb[ii]+ng[ii] : nbgM;
		}

	int size = 0;

	size += nuxM*sizeof(int); // e_imask_ux
	size += nbgM*sizeof(int); // e_imask_d

	size += 1*sizeof(struct STRVEC); // tmp_nuxM

	size += 1*SIZE_STRVEC(nuxM); // tmp_nuxM

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 64; // align to typical cache line size

	return size;

	}



void OCP_QP_REDUCE_EQ_DOF_WORK_CREATE(struct OCP_QP_DIM *dim, struct OCP_QP_REDUCE_EQ_DOF_WORK *work, void *mem)
	{

	int ii;

	// zero memory (to avoid corrupted memory like e.g. NaN)
	int memsize = OCP_QP_REDUCE_EQ_DOF_WORK_MEMSIZE(dim);
	hpipm_zero_memset(memsize, mem);

	// extract dim
	int N = dim->N;
	int *nu = dim->nu;
	int *nx = dim->nx;
	int *nb = dim->nb;
	int *ng = dim->ng;

	int nuxM = nu[0]+nx[0];
	int nbgM = nb[0]+ng[0];
	for(ii=1; ii<=N; ii++)
		{
		nuxM = nu[ii]+nx[ii]>nuxM ? nu[ii]+nx[ii] : nuxM;
		nbgM = nb[ii]+ng[ii]>nbgM ? nb[ii]+ng[ii] : nbgM;
		}


	// vector struct stuff
	struct STRVEC *sv_ptr = (struct STRVEC *) mem;

	// tmp_nuxM
	work->tmp_nuxM = sv_ptr;
	sv_ptr += 1;


	// integer stuff
	int *i_ptr;
	i_ptr = (int *) sv_ptr;

	// e_imask_ux
	work->e_imask_ux = i_ptr;
	i_ptr += nuxM;

	// e_imask_d
	work->e_imask_d = i_ptr;
	i_ptr += nbgM;


	// align to typical cache line size
	long long l_ptr = (long long) i_ptr;
	l_ptr = (l_ptr+63)/64*64;


	// floating point stuff
	char *c_ptr;
	c_ptr = (char *) l_ptr;


	// tmp_nuxM
	CREATE_STRVEC(nuxM, work->tmp_nuxM, c_ptr);
	c_ptr += (work->tmp_nuxM)->memsize;


	return;

	}



void OCP_QP_REDUCE_EQ_DOF(struct OCP_QP *qp, struct OCP_QP *qp_red, struct OCP_QP_REDUCE_EQ_DOF_WORK *work)
	{

	int ii, jj, idx;

	struct OCP_QP_DIM *dim = qp->dim;
	int N = dim->N;
	int *nx = dim->nx;
	int *nu = dim->nu;
	int *nb = dim->nb;
	int *nbx = dim->nbx;
	int *nbu = dim->nbu;
	int *ng = dim->ng;
	int *ns = dim->ns;
	int *nbue = dim->nbue;
	int *nbxe = dim->nbxe;
	int *nge = dim->nge;

	struct OCP_QP_DIM *dim_red = qp_red->dim;
	int *nb_red = dim_red->nb;
	int *ng_red = dim_red->ng;
	int *ns_red = dim_red->ns;

	// TODO handle case of softed equalities !!!!!!!!!!!!!!!!

	// initial stage
	ii = 0;
	if(nbue[ii]+nbxe[ii]>0) // reduce inputs and states
		{
		VECSE(nu[ii]+nx[ii], 0.0, work->tmp_nuxM, 0);
		for(jj=0; jj<nu[ii]+nx[ii]; jj++)
			work->e_imask_ux[jj] = 0;
		for(jj=0; jj<nb[ii]; jj++)
			work->e_imask_d[jj] = 0;
		// TODO e_mask_d !!!!!!!!!!!!!!!!!!!!!
		for(jj=0; jj<nbue[ii]+nbxe[ii]; jj++)
			{
#ifdef DOUBLE_PRECISION
			BLASFEO_DVECEL(work->tmp_nuxM, qp->idxb[ii][qp->idxe[ii][jj]]) = BLASFEO_DVECEL(qp->d, qp->idxe[ii][jj]);
#else
			BLASFEO_SVECEL(work->tmp_nuxM, qp->idxb[ii][qp->idxe[ii][jj]]) = BLASFEO_DVECEL(qp->d, qp->idxe[ii][jj]);
#endif
			work->e_imask_ux[qp->idxb[ii][qp->idxe[ii][jj]]] = 1;
			work->e_imask_d[qp->idxe[ii][jj]] = 1;
			}
		// XXX if ii<N !!!
		int_print_mat(1, nu[ii]+nx[ii], work->e_imask_ux, 1);
		int_print_mat(1, nb[ii]+ng[ii], work->e_imask_d, 1);
		// TODO check first and last non-zero in e_mask and only multiply between them
		// BAt
		idx = 0;
		for(jj=0; jj<nu[ii]+nx[ii]; jj++)
			{
			if(work->e_imask_ux[jj]==0)
				{
				GECP(1, nx[ii+1], qp->BAbt+ii, jj, 0, qp_red->BAbt+ii, idx, 0);
				idx++;
				}
			}
		// b
		GEMV_T(nu[ii]+nx[ii], nx[ii+1], 1.0, qp->BAbt+ii, 0, 0, work->tmp_nuxM, 0, 1.0, qp->b+ii, 0, qp_red->b+ii, 0);
		// RSQ
		// rq
		// d
		idx = 0;
		for(jj=0; jj<nb[ii]; jj++)
			{
			if(work->e_imask_d[jj]==0)
				{
#ifdef DOUBLE_PRECISION
				BLASFEO_DVECEL(qp_red->d+ii, idx) = BLASFEO_DVECEL(qp->d+ii, jj);
				BLASFEO_DVECEL(qp_red->d+ii, nb_red[ii]+ng_red[ii]+idx) = BLASFEO_DVECEL(qp->d+ii, nb[ii]+ng[ii]+jj);
#else
				BLASFEO_SVECEL(qp_red->d+ii, idx) = BLASFEO_SVECEL(qp->d+ii, jj);
				BLASFEO_SVECEL(qp_red->d+ii, nb_red[ii]+ng_red[ii]+idx) = BLASFEO_SVECEL(qp->d+ii, nb[ii]+ng[ii]+jj);
#endif
				idx++;
				}
			}
		VECCP(ng[ii], qp->d+ii, nb[ii], qp_red->d+ii, nb_red[ii]);
		VECCP(ng[ii]+2*ns[ii], qp->d+ii, 2*nb[ii]+ng[ii], qp_red->d+ii, 2*nb_red[ii]+2*ng_red[ii]);
		// idxb
		// DCt
		}
	else // copy everything
		{
		GECP(nu[ii]+nx[ii]+1, nx[ii+1], qp->BAbt+ii, 0, 0, qp_red->BAbt+ii, 0, 0);
		VECCP(nx[ii+1], qp->b+ii, 0, qp_red->b+ii, 0);
		GECP(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], qp->RSQrq+ii, 0, 0, qp_red->RSQrq+ii, 0, 0);
		VECCP(2*ns[ii], qp->Z+ii, 0, qp_red->Z+ii, 0);
		VECCP(nu[ii]+nx[ii]+2*ns[ii], qp->rqz+ii, 0, qp_red->rqz+ii, 0);
		for(jj=0; jj<nb[ii]; jj++)
			qp_red->idxb[ii][jj] = qp->idxb[ii][jj];
		GECP(nu[ii]+nx[ii], ng[ii], qp->DCt+ii, 0, 0, qp_red->DCt+ii, 0, 0);
		VECCP(2*nb[ii]+2*ng[ii]+2*ns[ii], qp->d+ii, 0, qp_red->d+ii, 0);
		VECCP(2*nb[ii]+2*ng[ii]+2*ns[ii], qp->d_mask+ii, 0, qp_red->d_mask+ii, 0);
		VECCP(2*nb[ii]+2*ng[ii]+2*ns[ii], qp->m+ii, 0, qp_red->m+ii, 0);
		for(jj=0; jj<ns[ii]; jj++)
			qp_red->idxs[ii][jj] = qp->idxs[ii][jj];
		for(jj=0; jj<nb[ii]+ng[ii]; jj++)
			qp_red->idxs_rev[ii][jj] = qp->idxs_rev[ii][jj];
		for(jj=0; jj<nbue[ii]+nbxe[ii]+nge[ii]; jj++)
			qp_red->idxe[ii][jj] = qp->idxe[ii][jj];
		}
	// other stages
	for(ii=1; ii<=N; ii++)
		{
		if(nbue[ii]>0) // reduce inputs and states
			{
			// TODO
			}
		else // copy everything
			{
			if(ii<N)
				{
				GECP(nu[ii]+nx[ii]+1, nx[ii+1], qp->BAbt+ii, 0, 0, qp_red->BAbt+ii, 0, 0);
				VECCP(nx[ii+1], qp->b+ii, 0, qp_red->b+ii, 0);
				}
			GECP(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], qp->RSQrq+ii, 0, 0, qp_red->RSQrq+ii, 0, 0);
			VECCP(2*ns[ii], qp->Z+ii, 0, qp_red->Z+ii, 0);
			VECCP(nu[ii]+nx[ii]+2*ns[ii], qp->rqz+ii, 0, qp_red->rqz+ii, 0);
			for(jj=0; jj<nb[ii]; jj++)
				qp_red->idxb[ii][jj] = qp->idxb[ii][jj];
			GECP(nu[ii]+nx[ii], ng[ii], qp->DCt+ii, 0, 0, qp_red->DCt+ii, 0, 0);
			VECCP(2*nb[ii]+2*ng[ii]+2*ns[ii], qp->d+ii, 0, qp_red->d+ii, 0);
			VECCP(2*nb[ii]+2*ng[ii]+2*ns[ii], qp->d_mask+ii, 0, qp_red->d_mask+ii, 0);
			VECCP(2*nb[ii]+2*ng[ii]+2*ns[ii], qp->m+ii, 0, qp_red->m+ii, 0);
			for(jj=0; jj<ns[ii]; jj++)
				qp_red->idxs[ii][jj] = qp->idxs[ii][jj];
			for(jj=0; jj<nb[ii]+ng[ii]; jj++)
				qp_red->idxs_rev[ii][jj] = qp->idxs_rev[ii][jj];
			for(jj=0; jj<nbue[ii]+nbxe[ii]+nge[ii]; jj++)
				qp_red->idxe[ii][jj] = qp->idxe[ii][jj];
			}
		}

	return;

	}
