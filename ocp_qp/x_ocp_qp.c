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



int OCP_QP_STRSIZE()
	{
	return sizeof(struct OCP_QP);
	}



int OCP_QP_MEMSIZE(struct OCP_QP_DIM *dim)
	{

	// extract dim
	int N = dim->N;
	int *nx = dim->nx;
	int *nu = dim->nu;
	int *nb = dim->nb;
	int *ng = dim->ng;
	int *ns = dim->ns;
	int *nbue = dim->nbue;
	int *nbxe = dim->nbxe;
	int *nge = dim->nge;

	// loop index
	int ii;

	// compute core qp size
	int nvt = 0;
	int net = 0;
	int nct = 0;
	for(ii=0; ii<N; ii++)
		{
		nvt += nx[ii]+nu[ii]+2*ns[ii];
		net += nx[ii+1];
		nct += 2*nb[ii]+2*ng[ii]+2*ns[ii];
		}
	nvt += nx[ii]+nu[ii]+2*ns[ii];
	nct += 2*nb[ii]+2*ng[ii]+2*ns[ii];

	int size = 0;

	size += 5*(N+1)*sizeof(int); // nx nu nb ng ns
	size += 3*(N+1)*sizeof(int *); // idxb idxs_rev idxe
	size += 2*(N+1)*sizeof(struct STRMAT); // RSqrq DCt
	size += 1*N*sizeof(struct STRMAT); // BAbt
	size += 5*(N+1)*sizeof(struct STRVEC); // rqz d m Z d_mask
	size += 1*N*sizeof(struct STRVEC); // b

	for(ii=0; ii<N; ii++)
		{
		size += nb[ii]*sizeof(int); // idxb
		size += (nb[ii]+ng[ii])*sizeof(int); // idxs_rev
		size += (nbue[ii]+nbxe[ii]+nge[ii])*sizeof(int); // idxe
		size += SIZE_STRMAT(nu[ii]+nx[ii]+1, nx[ii+1]); // BAbt
		size += SIZE_STRMAT(nu[ii]+nx[ii]+1, nu[ii]+nx[ii]); // RSQrq
		size += SIZE_STRMAT(nu[ii]+nx[ii], ng[ii]); // DCt
		size += SIZE_STRVEC(2*ns[ii]); // Z
		}
	ii = N;
	size += nb[ii]*sizeof(int); // idxb
	size += (nb[ii]+ng[ii])*sizeof(int); // idxs_rev
	size += (nbue[ii]+nbxe[ii]+nge[ii])*sizeof(int); // idxe
	size += SIZE_STRMAT(nu[ii]+nx[ii]+1, nu[ii]+nx[ii]); // RSQrq
	size += SIZE_STRMAT(nu[ii]+nx[ii], ng[ii]); // DCt
	size += SIZE_STRVEC(2*ns[ii]); // Z

	size += 1*SIZE_STRVEC(nvt); // rqz
	size += 1*SIZE_STRVEC(net); // b
	size += 3*SIZE_STRVEC(nct); // d m d_mask

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 64; // align to typical cache line size

	return size;

	}



void OCP_QP_CREATE(struct OCP_QP_DIM *dim, struct OCP_QP *qp, void *mem)
	{

	// loop index
	int ii, jj;

	// zero memory (to avoid corrupted memory like e.g. NaN)
	int memsize = OCP_QP_MEMSIZE(dim);
	hpipm_zero_memset(memsize, mem);

	// extract dim
	int N = dim->N;
	int *nx = dim->nx;
	int *nu = dim->nu;
	int *nb = dim->nb;
	int *ng = dim->ng;
	int *ns = dim->ns;
	int *nbue = dim->nbue;
	int *nbxe = dim->nbxe;
	int *nge = dim->nge;

	// compute core qp size
	int nvt = 0;
	int net = 0;
	int nct = 0;
	for(ii=0; ii<N; ii++)
		{
		nvt += nx[ii]+nu[ii]+2*ns[ii];
		net += nx[ii+1];
		nct += 2*nb[ii]+2*ng[ii]+2*ns[ii];
		}
	nvt += nx[ii]+nu[ii]+2*ns[ii];
	nct += 2*nb[ii]+2*ng[ii]+2*ns[ii];



	// int pointer stuff
	int **ip_ptr;
	ip_ptr = (int **) mem;

	// idxb
	qp->idxb = ip_ptr;
	ip_ptr += N+1;

	// idxs_rev
	qp->idxs_rev = ip_ptr;
	ip_ptr += N+1;

	// idxe
	qp->idxe = ip_ptr;
	ip_ptr += N+1;


	// matrix struct stuff
	struct STRMAT *sm_ptr = (struct STRMAT *) ip_ptr;

	// BAbt
	qp->BAbt = sm_ptr;
	sm_ptr += N;

	// RSQrq
	qp->RSQrq = sm_ptr;
	sm_ptr += N+1;

	// DCt
	qp->DCt = sm_ptr;
	sm_ptr += N+1;


	// vector struct stuff
	struct STRVEC *sv_ptr = (struct STRVEC *) sm_ptr;

	// b
	qp->b = sv_ptr;
	sv_ptr += N;

	// rqz
	qp->rqz = sv_ptr;
	sv_ptr += N+1;

	// d
	qp->d = sv_ptr;
	sv_ptr += N+1;

	// d_mask
	qp->d_mask = sv_ptr;
	sv_ptr += N+1;

	// m
	qp->m = sv_ptr;
	sv_ptr += N+1;

	// Z
	qp->Z = sv_ptr;
	sv_ptr += N+1;


	// integer stuff
	int *i_ptr;
	i_ptr = (int *) sv_ptr;

	// idxb
	for(ii=0; ii<=N; ii++)
		{
		(qp->idxb)[ii] = i_ptr;
		i_ptr += nb[ii];
		}

	// idxs_rev
	for(ii=0; ii<=N; ii++)
		{
		(qp->idxs_rev)[ii] = i_ptr;
		i_ptr += nb[ii]+ng[ii];
		// default value: -1
		for(jj=0; jj<nb[ii]+ng[ii]; jj++)
			qp->idxs_rev[ii][jj] = -1;
		}

	// idxe
	for(ii=0; ii<=N; ii++)
		{
		(qp->idxe)[ii] = i_ptr;
		i_ptr += nbue[ii]+nbxe[ii]+nge[ii];
		}


	// align to typical cache line size
	long long l_ptr = (long long) i_ptr;
	l_ptr = (l_ptr+63)/64*64;


	// floating point stuff
	char *c_ptr;
	c_ptr = (char *) l_ptr;

	char *tmp_ptr;

	// BAbt
	for(ii=0; ii<N; ii++)
		{
		CREATE_STRMAT(nu[ii]+nx[ii]+1, nx[ii+1], qp->BAbt+ii, c_ptr);
		c_ptr += (qp->BAbt+ii)->memsize;
		GESE(nu[ii]+nx[ii]+1, nx[ii+1], 0.0, qp->BAbt+ii, 0, 0);
		}

	// RSQrq
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRMAT(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], qp->RSQrq+ii, c_ptr);
		c_ptr += (qp->RSQrq+ii)->memsize;
		GESE(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], 0.0, qp->RSQrq+ii, 0, 0);
		}

	// DCt
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRMAT(nu[ii]+nx[ii], ng[ii], qp->DCt+ii, c_ptr);
		c_ptr += (qp->DCt+ii)->memsize;
		GESE(nu[ii]+nx[ii], ng[ii], 0.0, qp->DCt+ii, 0, 0);
		}

	// Z
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(2*ns[ii], qp->Z+ii, c_ptr);
		c_ptr += (qp->Z+ii)->memsize;
		VECSE(2*ns[ii], 0.0, qp->Z+ii, 0);
		}

	// g
	tmp_ptr = c_ptr;
	c_ptr += SIZE_STRVEC(nvt);
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(nu[ii]+nx[ii]+2*ns[ii], qp->rqz+ii, tmp_ptr);
		tmp_ptr += nu[ii]*sizeof(REAL);
		tmp_ptr += nx[ii]*sizeof(REAL);
		tmp_ptr += ns[ii]*sizeof(REAL);
		tmp_ptr += ns[ii]*sizeof(REAL);
		VECSE(nu[ii]+nx[ii]+2*ns[ii], 0.0, qp->rqz+ii, 0);
		}

	// b
	tmp_ptr = c_ptr;
	c_ptr += SIZE_STRVEC(net);
	for(ii=0; ii<N; ii++)
		{
		CREATE_STRVEC(nx[ii+1], qp->b+ii, tmp_ptr);
		tmp_ptr += nx[ii+1]*sizeof(REAL);
		VECSE(nx[ii+1], 0.0, qp->b+ii, 0);
		}

	// d
	tmp_ptr = c_ptr;
	c_ptr += SIZE_STRVEC(nct);
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(2*nb[ii]+2*ng[ii]+2*ns[ii], qp->d+ii, tmp_ptr);
		tmp_ptr += nb[ii]*sizeof(REAL); // lb
		tmp_ptr += ng[ii]*sizeof(REAL); // lg
		tmp_ptr += nb[ii]*sizeof(REAL); // ub
		tmp_ptr += ng[ii]*sizeof(REAL); // ug
		tmp_ptr += ns[ii]*sizeof(REAL); // ls
		tmp_ptr += ns[ii]*sizeof(REAL); // us
		VECSE(2*nb[ii]+2*ng[ii]+2*ns[ii], 0.0, qp->d+ii, 0);
		}

	// d_mask
	tmp_ptr = c_ptr;
	c_ptr += SIZE_STRVEC(nct);
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(2*nb[ii]+2*ng[ii]+2*ns[ii], qp->d_mask+ii, tmp_ptr);
		tmp_ptr += nb[ii]*sizeof(REAL); // lb
		tmp_ptr += ng[ii]*sizeof(REAL); // lg
		tmp_ptr += nb[ii]*sizeof(REAL); // ub
		tmp_ptr += ng[ii]*sizeof(REAL); // ug
		tmp_ptr += ns[ii]*sizeof(REAL); // ls
		tmp_ptr += ns[ii]*sizeof(REAL); // us
		VECSE(2*nb[ii]+2*ng[ii]+2*ns[ii], 1.0, qp->d_mask+ii, 0);
		}

	// m
	tmp_ptr = c_ptr;
	c_ptr += SIZE_STRVEC(nct);
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(2*nb[ii]+2*ng[ii]+2*ns[ii], qp->m+ii, tmp_ptr);
		tmp_ptr += nb[ii]*sizeof(REAL); // lb
		tmp_ptr += ng[ii]*sizeof(REAL); // lg
		tmp_ptr += nb[ii]*sizeof(REAL); // ub
		tmp_ptr += ng[ii]*sizeof(REAL); // ug
		tmp_ptr += ns[ii]*sizeof(REAL); // ls
		tmp_ptr += ns[ii]*sizeof(REAL); // us
		VECSE(2*nb[ii]+2*ng[ii]+2*ns[ii], 0.0, qp->m+ii, 0);
		}

	qp->dim = dim;

	qp->memsize = OCP_QP_MEMSIZE(dim);


#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + qp->memsize)
		{
		printf("\nerror: OCP_QP_CREATE: outside memory bounds!\n\n");
		exit(1);
		}
#endif


	return;

	}



void OCP_QP_COPY_ALL(struct OCP_QP *qp_orig, struct OCP_QP *qp_dest)
	{

	// extract dim
	int N = qp_orig->dim->N;
	int *nx = qp_orig->dim->nx;
	int *nu = qp_orig->dim->nu;
	int *nb = qp_orig->dim->nb;
	int *nbx = qp_orig->dim->nbx;
	int *nbu = qp_orig->dim->nbu;
	int *ng = qp_orig->dim->ng;
	int *ns = qp_orig->dim->ns;
	int *nbxe = qp_orig->dim->nbxe;
	int *nbue = qp_orig->dim->nbue;
	int *nge = qp_orig->dim->nge;

	int ii, jj;

	// copy dim pointer
//	qp_dest->dim = qp_orig->dim;

	for(ii=0; ii<N; ii++)
		{
		GECP(nu[ii]+nx[ii]+1, nx[ii+1], qp_orig->BAbt+ii, 0, 0, qp_dest->BAbt+ii, 0, 0);
		VECCP(nx[ii+1], qp_orig->b+ii, 0, qp_dest->b+ii, 0);
		}

	for(ii=0; ii<=N; ii++)
		{
		GECP(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], qp_orig->RSQrq+ii, 0, 0, qp_dest->RSQrq+ii, 0, 0);
		VECCP(2*ns[ii], qp_orig->Z+ii, 0, qp_dest->Z+ii, 0);
		VECCP(nu[ii]+nx[ii]+2*ns[ii], qp_orig->rqz+ii, 0, qp_dest->rqz+ii, 0);
		for(jj=0; jj<nb[ii]; jj++)
			qp_dest->idxb[ii][jj] = qp_orig->idxb[ii][jj];
		GECP(nu[ii]+nx[ii], ng[ii], qp_orig->DCt+ii, 0, 0, qp_dest->DCt+ii, 0, 0);
		VECCP(2*nb[ii]+2*ng[ii]+2*ns[ii], qp_orig->d+ii, 0, qp_dest->d+ii, 0);
		VECCP(2*nb[ii]+2*ng[ii]+2*ns[ii], qp_orig->d_mask+ii, 0, qp_dest->d_mask+ii, 0);
		VECCP(2*nb[ii]+2*ng[ii]+2*ns[ii], qp_orig->m+ii, 0, qp_dest->m+ii, 0);
		for(jj=0; jj<nb[ii]+ng[ii]; jj++)
			qp_dest->idxs_rev[ii][jj] = qp_orig->idxs_rev[ii][jj];
		for(jj=0; jj<nbue[ii]+nbxe[ii]+nge[ii]; jj++)
			qp_dest->idxe[ii][jj] = qp_orig->idxe[ii][jj];
		}

	return;

	}



void OCP_QP_SET_ALL_ZERO(struct OCP_QP *qp)
	{

	// extract dim
	int N = qp->dim->N;
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;
	int *nb = qp->dim->nb;
	int *nbx = qp->dim->nbx;
	int *nbu = qp->dim->nbu;
	int *ng = qp->dim->ng;
	int *ns = qp->dim->ns;
	int *nbxe = qp->dim->nbxe;
	int *nbue = qp->dim->nbue;
	int *nge = qp->dim->nge;

	int ii, jj;

	for(ii=0; ii<N; ii++)
		{
		GESE(nu[ii]+nx[ii]+1, nx[ii+1], 0.0, qp->BAbt+ii, 0, 0);
		VECSE(nx[ii+1], 0.0, qp->b+ii, 0);
		}

	for(ii=0; ii<=N; ii++)
		{
		GESE(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], 0.0, qp->RSQrq+ii, 0, 0);
		VECSE(2*ns[ii], 0.0, qp->Z+ii, 0);
		VECSE(nu[ii]+nx[ii]+2*ns[ii], 0.0, qp->rqz+ii, 0);
		for(jj=0; jj<nb[ii]; jj++)
			qp->idxb[ii][jj] = 0;
		GESE(nu[ii]+nx[ii], ng[ii], 0.0, qp->DCt+ii, 0, 0);
		VECSE(2*nb[ii]+2*ng[ii]+2*ns[ii], 0.0, qp->d+ii, 0);
		VECSE(2*nb[ii]+2*ng[ii]+2*ns[ii], 1.0, qp->d_mask+ii, 0);
		VECSE(2*nb[ii]+2*ng[ii]+2*ns[ii], 0.0, qp->m+ii, 0);
		for(jj=0; jj<ns[ii]; jj++)
			qp->idxs_rev[ii][jj] = -1;
		for(jj=0; jj<nbue[ii]+nbxe[ii]+nge[ii]; jj++)
			qp->idxe[ii][jj] = 0;
		}

	return;

	}



void OCP_QP_SET_RHS_ZERO(struct OCP_QP *qp)
	{

	// extract dim
	int N = qp->dim->N;
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;
	int *nb = qp->dim->nb;
	int *nbx = qp->dim->nbx;
	int *nbu = qp->dim->nbu;
	int *ng = qp->dim->ng;
	int *ns = qp->dim->ns;

	int ii, jj;

	for(ii=0; ii<N; ii++)
		{
		VECSE(nx[ii+1], 0.0, qp->b+ii, 0);
		}

	for(ii=0; ii<=N; ii++)
		{
		VECSE(2*ns[ii], 0.0, qp->Z+ii, 0);
		VECSE(nu[ii]+nx[ii]+2*ns[ii], 0.0, qp->rqz+ii, 0);
		VECSE(2*nb[ii]+2*ng[ii]+2*ns[ii], 0.0, qp->d+ii, 0);
		VECSE(2*nb[ii]+2*ng[ii]+2*ns[ii], 1.0, qp->d_mask+ii, 0);
		VECSE(2*nb[ii]+2*ng[ii]+2*ns[ii], 0.0, qp->m+ii, 0);
		}

	return;

	}



// TODO deprecate and remove ???
void OCP_QP_SET_ALL(REAL **A, REAL **B, REAL **b, REAL **Q, REAL **S, REAL **R, REAL **q, REAL **r, int **idxbx, REAL **d_lbx, REAL **d_ubx, int **idxbu, REAL **d_lbu, REAL **d_ubu, REAL **C, REAL **D, REAL **d_lg, REAL **d_ug, REAL **Zl, REAL **Zu, REAL **zl, REAL **zu, int **idxs, REAL **d_ls, REAL **d_us, struct OCP_QP *qp)
	{

	// extract dim
	int N = qp->dim->N;
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;
	int *nb = qp->dim->nb;
	int *nbx = qp->dim->nbx;
	int *nbu = qp->dim->nbu;
	int *ng = qp->dim->ng;
	int *ns = qp->dim->ns;

	int ii, jj;

	for(ii=0; ii<N; ii++)
		{
		CVT_TRAN_MAT2STRMAT(nx[ii+1], nu[ii], B[ii], nx[ii+1], qp->BAbt+ii, 0, 0);
		CVT_TRAN_MAT2STRMAT(nx[ii+1], nx[ii], A[ii], nx[ii+1], qp->BAbt+ii, nu[ii], 0);
		CVT_TRAN_MAT2STRMAT(nx[ii+1], 1, b[ii], nx[ii+1], qp->BAbt+ii, nu[ii]+nx[ii], 0); // XXX remove ???
		CVT_VEC2STRVEC(nx[ii+1], b[ii], qp->b+ii, 0);
		}

	for(ii=0; ii<=N; ii++)
		{
		CVT_MAT2STRMAT(nu[ii], nu[ii], R[ii], nu[ii], qp->RSQrq+ii, 0, 0);
		CVT_TRAN_MAT2STRMAT(nu[ii], nx[ii], S[ii], nu[ii], qp->RSQrq+ii, nu[ii], 0);
		CVT_MAT2STRMAT(nx[ii], nx[ii], Q[ii], nx[ii], qp->RSQrq+ii, nu[ii], nu[ii]);
		CVT_TRAN_MAT2STRMAT(nu[ii], 1, r[ii], nu[ii], qp->RSQrq+ii, nu[ii]+nx[ii], 0); // XXX remove ???
		CVT_TRAN_MAT2STRMAT(nx[ii], 1, q[ii], nx[ii], qp->RSQrq+ii, nu[ii]+nx[ii], nu[ii]); // XXX remove ???
		CVT_VEC2STRVEC(nu[ii], r[ii], qp->rqz+ii, 0);
		CVT_VEC2STRVEC(nx[ii], q[ii], qp->rqz+ii, nu[ii]);
		}

	for(ii=0; ii<=N; ii++)
		{
		if(nbu[ii]>0)
			{
			for(jj=0; jj<nbu[ii]; jj++)
				qp->idxb[ii][jj] = idxbu[ii][jj];
			CVT_VEC2STRVEC(nbu[ii], d_lbu[ii], qp->d+ii, 0);
			CVT_VEC2STRVEC(nbu[ii], d_ubu[ii], qp->d+ii, nb[ii]+ng[ii]);
			}
		if(nbx[ii>0])
			{
			for(jj=0; jj<nbx[ii]; jj++)
				qp->idxb[ii][nbu[ii]+jj] = nu[ii]+idxbx[ii][jj];
			CVT_VEC2STRVEC(nbx[ii], d_lbx[ii], qp->d+ii, nbu[ii]);
			CVT_VEC2STRVEC(nbx[ii], d_ubx[ii], qp->d+ii, nb[ii]+ng[ii]+nbu[ii]);
			}
		if(nb[ii>0])
			{
			VECSC(nb[ii], -1.0, qp->d+ii, nb[ii]+ng[ii]);
			VECSE(nb[ii], 0.0, qp->m+ii, 0);
			VECSE(nb[ii], 0.0, qp->m+ii, nb[ii]+ng[ii]);
			}
		}

	for(ii=0; ii<=N; ii++)
		{
		if(ng[ii]>0)
			{
			CVT_TRAN_MAT2STRMAT(ng[ii], nu[ii], D[ii], ng[ii], qp->DCt+ii, 0, 0);
			CVT_TRAN_MAT2STRMAT(ng[ii], nx[ii], C[ii], ng[ii], qp->DCt+ii, nu[ii], 0);
			CVT_VEC2STRVEC(ng[ii], d_lg[ii], qp->d+ii, nb[ii]);
			CVT_VEC2STRVEC(ng[ii], d_ug[ii], qp->d+ii, 2*nb[ii]+ng[ii]);
			VECSC(ng[ii], -1.0, qp->d+ii, 2*nb[ii]+ng[ii]);
			VECSE(ng[ii], 0.0, qp->m+ii, nb[ii]);
			VECSE(ng[ii], 0.0, qp->m+ii, 2*nb[ii]+ng[ii]);
			}
		}

	for(ii=0; ii<=N; ii++)
		{
		if(ns[ii]>0)
			{
			for(jj=0; jj<ns[ii]; jj++)
				{
				qp->idxs_rev[ii][idxs[ii][jj]] = jj;
				}
			CVT_VEC2STRVEC(ns[ii], Zl[ii], qp->Z+ii, 0);
			CVT_VEC2STRVEC(ns[ii], Zu[ii], qp->Z+ii, ns[ii]);
			CVT_VEC2STRVEC(ns[ii], zl[ii], qp->rqz+ii, nu[ii]+nx[ii]);
			CVT_VEC2STRVEC(ns[ii], zu[ii], qp->rqz+ii, nu[ii]+nx[ii]+ns[ii]);
			CVT_VEC2STRVEC(ns[ii], d_ls[ii], qp->d+ii, 2*nb[ii]+2*ng[ii]);
			CVT_VEC2STRVEC(ns[ii], d_us[ii], qp->d+ii, 2*nb[ii]+2*ng[ii]+ns[ii]);
			VECSE(ns[ii], 0.0, qp->m+ii, 2*nb[ii]+2*ng[ii]);
			VECSE(ns[ii], 0.0, qp->m+ii, 2*nb[ii]+2*ng[ii]+ns[ii]);
			}
		}

	return;

	}



void OCP_QP_SET_ALL_ROWMAJ(REAL **A, REAL **B, REAL **b, REAL **Q, REAL **S, REAL **R, REAL **q, REAL **r, int **idxbx, REAL **d_lbx, REAL **d_ubx, int **idxbu, REAL **d_lbu, REAL **d_ubu, REAL **C, REAL **D, REAL **d_lg, REAL **d_ug, REAL **Zl, REAL **Zu, REAL **zl, REAL **zu, int **idxs, REAL **d_ls, REAL **d_us, struct OCP_QP *qp)
	{

	// extract dim
	int N = qp->dim->N;
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;
	int *nb = qp->dim->nb;
	int *nbx = qp->dim->nbx;
	int *nbu = qp->dim->nbu;
	int *ng = qp->dim->ng;
	int *ns = qp->dim->ns;

	int ii, jj;

	for(ii=0; ii<N; ii++)
		{
		CVT_MAT2STRMAT(nu[ii], nx[ii+1], B[ii], nu[ii], qp->BAbt+ii, 0, 0);
		CVT_MAT2STRMAT(nx[ii], nx[ii+1], A[ii], nx[ii], qp->BAbt+ii, nu[ii], 0);
		CVT_TRAN_MAT2STRMAT(nx[ii+1], 1, b[ii], nx[ii+1], qp->BAbt+ii, nu[ii]+nx[ii], 0); // XXX remove ???
		CVT_VEC2STRVEC(nx[ii+1], b[ii], qp->b+ii, 0);
		}

	for(ii=0; ii<=N; ii++)
		{
		CVT_MAT2STRMAT(nu[ii], nu[ii], R[ii], nu[ii], qp->RSQrq+ii, 0, 0);
		CVT_MAT2STRMAT(nx[ii], nu[ii], S[ii], nx[ii], qp->RSQrq+ii, nu[ii], 0);
		CVT_MAT2STRMAT(nx[ii], nx[ii], Q[ii], nx[ii], qp->RSQrq+ii, nu[ii], nu[ii]);
		CVT_TRAN_MAT2STRMAT(nu[ii], 1, r[ii], nu[ii], qp->RSQrq+ii, nu[ii]+nx[ii], 0); // XXX remove ???
		CVT_TRAN_MAT2STRMAT(nx[ii], 1, q[ii], nx[ii], qp->RSQrq+ii, nu[ii]+nx[ii], nu[ii]); // XXX remove ???
		CVT_VEC2STRVEC(nu[ii], r[ii], qp->rqz+ii, 0);
		CVT_VEC2STRVEC(nx[ii], q[ii], qp->rqz+ii, nu[ii]);
		}

	for(ii=0; ii<=N; ii++)
		{
		if(nbu[ii]>0)
			{
			for(jj=0; jj<nbu[ii]; jj++)
				qp->idxb[ii][jj] = idxbu[ii][jj];
			CVT_VEC2STRVEC(nbu[ii], d_lbu[ii], qp->d+ii, 0);
			CVT_VEC2STRVEC(nbu[ii], d_ubu[ii], qp->d+ii, nb[ii]+ng[ii]);
			}
		if(nbx[ii>0])
			{
			for(jj=0; jj<nbx[ii]; jj++)
				qp->idxb[ii][nbu[ii]+jj] = nu[ii]+idxbx[ii][jj];
			CVT_VEC2STRVEC(nbx[ii], d_lbx[ii], qp->d+ii, nbu[ii]);
			CVT_VEC2STRVEC(nbx[ii], d_ubx[ii], qp->d+ii, nb[ii]+ng[ii]+nbu[ii]);
			}
		if(nb[ii>0])
			{
			VECSC(nb[ii], -1.0, qp->d+ii, nb[ii]+ng[ii]);
			VECSE(nb[ii], 0.0, qp->m+ii, 0);
			VECSE(nb[ii], 0.0, qp->m+ii, nb[ii]+ng[ii]);
			}
		}

	for(ii=0; ii<=N; ii++)
		{
		if(ng[ii]>0)
			{
			CVT_MAT2STRMAT(nu[ii], ng[ii], D[ii], nu[ii], qp->DCt+ii, 0, 0);
			CVT_MAT2STRMAT(nx[ii], ng[ii], C[ii], nx[ii], qp->DCt+ii, nu[ii], 0);
			CVT_VEC2STRVEC(ng[ii], d_lg[ii], qp->d+ii, nb[ii]);
			CVT_VEC2STRVEC(ng[ii], d_ug[ii], qp->d+ii, 2*nb[ii]+ng[ii]);
			VECSC(ng[ii], -1.0, qp->d+ii, 2*nb[ii]+ng[ii]);
			VECSE(ng[ii], 0.0, qp->m+ii, nb[ii]);
			VECSE(ng[ii], 0.0, qp->m+ii, 2*nb[ii]+ng[ii]);
			}
		}

	for(ii=0; ii<=N; ii++)
		{
		if(ns[ii]>0)
			{
			for(jj=0; jj<ns[ii]; jj++)
				{
				qp->idxs_rev[ii][idxs[ii][jj]] = jj;
				}
			CVT_VEC2STRVEC(ns[ii], Zl[ii], qp->Z+ii, 0);
			CVT_VEC2STRVEC(ns[ii], Zu[ii], qp->Z+ii, ns[ii]);
			CVT_VEC2STRVEC(ns[ii], zl[ii], qp->rqz+ii, nu[ii]+nx[ii]);
			CVT_VEC2STRVEC(ns[ii], zu[ii], qp->rqz+ii, nu[ii]+nx[ii]+ns[ii]);
			CVT_VEC2STRVEC(ns[ii], d_ls[ii], qp->d+ii, 2*nb[ii]+2*ng[ii]);
			CVT_VEC2STRVEC(ns[ii], d_us[ii], qp->d+ii, 2*nb[ii]+2*ng[ii]+ns[ii]);
			VECSE(ns[ii], 0.0, qp->m+ii, 2*nb[ii]+2*ng[ii]);
			VECSE(ns[ii], 0.0, qp->m+ii, 2*nb[ii]+2*ng[ii]+ns[ii]);
			}
		}

	return;

	}



void OCP_QP_SET(char *field, int stage, void *value, struct OCP_QP *qp)
	{
	REAL *r_ptr;
	int *i_ptr;
    
	// matrices
	if(hpipm_strcmp(field, "A")) 
		{
		OCP_QP_SET_A(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "B")) 
		{
		OCP_QP_SET_B(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "Q")) 
		{
		OCP_QP_SET_Q(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "S")) 
		{
		OCP_QP_SET_S(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "R")) 
		{
		OCP_QP_SET_R(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "C")) 
		{
		OCP_QP_SET_C(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "D")) 
		{
		OCP_QP_SET_D(stage, value, qp);
		}
	// vectors
	else if(hpipm_strcmp(field, "b"))
		{ 
		OCP_QP_SET_BVEC(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "q"))
		{ 
		OCP_QP_SET_QVEC(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "r"))
		{ 
		OCP_QP_SET_RVEC(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "lb"))
		{ 
		OCP_QP_SET_LB(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "lb_mask"))
		{ 
		OCP_QP_SET_LB_MASK(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "lbu") | hpipm_strcmp(field, "lu"))
		{ 
		OCP_QP_SET_LBU(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "lbu_mask"))
		{ 
		OCP_QP_SET_LBU_MASK(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "lbx") | hpipm_strcmp(field, "lx"))
		{ 
		OCP_QP_SET_LBX(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "lbx_mask"))
		{ 
		OCP_QP_SET_LBX_MASK(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "ub"))
		{ 
		OCP_QP_SET_UB(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "ub_mask"))
		{ 
		OCP_QP_SET_UB_MASK(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "ubu") | hpipm_strcmp(field, "uu"))
		{ 
		OCP_QP_SET_UBU(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "ubu_mask"))
		{ 
		OCP_QP_SET_UBU_MASK(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "ubx") | hpipm_strcmp(field, "ux"))
		{ 
		OCP_QP_SET_UBX(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "ubx_mask"))
		{ 
		OCP_QP_SET_UBX_MASK(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "lg"))
		{ 
		OCP_QP_SET_LG(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "lg_mask"))
		{ 
		OCP_QP_SET_LG_MASK(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "ug"))
		{ 
		OCP_QP_SET_UG(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "ug_mask"))
		{ 
		OCP_QP_SET_UG_MASK(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "Zl"))
		{ 
		OCP_QP_SET_ZL(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "Zu"))
		{ 
		OCP_QP_SET_ZU(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "zl"))
		{ 
		OCP_QP_SET_ZLVEC(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "zu"))
		{ 
		OCP_QP_SET_ZUVEC(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "lls"))
		{ 
		OCP_QP_SET_LLS(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "lls_mask"))
		{ 
		OCP_QP_SET_LLS_MASK(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "lus"))
		{ 
		OCP_QP_SET_LUS(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "lus_mask"))
		{ 
		OCP_QP_SET_LUS_MASK(stage, value, qp);
		}
	// int
	else if(hpipm_strcmp(field, "idxb"))
		{
		OCP_QP_SET_IDXB(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "idxbx"))
		{
		OCP_QP_SET_IDXBX(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "Jbx") | hpipm_strcmp(field, "Jx"))
		{
		OCP_QP_SET_JBX(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "idxbu"))
		{
		OCP_QP_SET_IDXBU(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "Jbu") | hpipm_strcmp(field, "Ju"))
		{
		OCP_QP_SET_JBU(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "idxs"))
		{
		OCP_QP_SET_IDXS(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "idxs_rev"))
		{
		OCP_QP_SET_IDXS_REV(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "Jsbu") | hpipm_strcmp(field, "Jsu"))
		{
		OCP_QP_SET_JSBU(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "Jsbx") | hpipm_strcmp(field, "Jsx"))
		{
		OCP_QP_SET_JSBX(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "Jsg"))
		{
		OCP_QP_SET_JSG(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "idxe"))
		{
		OCP_QP_SET_IDXE(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "idxbue"))
		{
		OCP_QP_SET_IDXBUE(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "idxbxe"))
		{
		OCP_QP_SET_IDXBXE(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "idxge"))
		{
		OCP_QP_SET_IDXGE(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "Jbue"))
		{
		OCP_QP_SET_JBUE(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "Jbxe"))
		{
		OCP_QP_SET_JBXE(stage, value, qp);
		}
	else if(hpipm_strcmp(field, "Jge"))
		{
		OCP_QP_SET_JGE(stage, value, qp);
		}
	else
		{
		printf("error: OCP_QP_SET: wrong field name '%s'. Exiting.\n", field);
		exit(1);	
		}
	return;
	}



void OCP_QP_SET_EL(char *field, int stage, int index, void *elem, struct OCP_QP *qp)
	{
	REAL *r_ptr;
	int *i_ptr;
    
	// matrices
	if(hpipm_strcmp(field, "lbx") | hpipm_strcmp(field, "lx"))
		{ 
		OCP_QP_SET_EL_LBX(stage, index, elem, qp);
		}
	else if(hpipm_strcmp(field, "ubx") | hpipm_strcmp(field, "ux"))
		{ 
		OCP_QP_SET_EL_UBX(stage, index, elem, qp);
		}
	else
		{
		printf("error: OCP_QP_SET_EL: wrong field%s\n", field);
		exit(1);	
		}
	return;
	}



void OCP_QP_SET_A(int stage, REAL *A, struct OCP_QP *qp)
	{
	// extract dim
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;

	CVT_TRAN_MAT2STRMAT(nx[stage+1], nx[stage], A, nx[stage+1], qp->BAbt+stage, nu[stage], 0);

	return;
	}



void OCP_QP_SET_B(int stage, REAL *B, struct OCP_QP *qp)
	{
	// extract dim
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;

	CVT_TRAN_MAT2STRMAT(nx[stage+1], nu[stage], B, nx[stage+1], qp->BAbt+stage, 0, 0);

	return;
	}



void OCP_QP_SET_BVEC(int stage, REAL *b, struct OCP_QP *qp)
	{
	// extract dim
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;

	CVT_TRAN_MAT2STRMAT(nx[stage+1], 1, b, nx[stage+1], &(qp->BAbt[stage]), nu[stage]+nx[stage], 0); // TODO remove ???
	CVT_VEC2STRVEC(nx[stage+1], b, qp->b+stage, 0);

	return;
	}



void OCP_QP_SET_Q(int stage, REAL *Q, struct OCP_QP *qp)
	{
	// extract dim
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;

	CVT_MAT2STRMAT(nx[stage], nx[stage], Q, nx[stage], qp->RSQrq+stage, nu[stage], nu[stage]);

	return;
	}



void OCP_QP_SET_S(int stage, REAL *S, struct OCP_QP *qp)
	{
	// extract dim
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;

	CVT_TRAN_MAT2STRMAT(nu[stage], nx[stage], S, nu[stage], qp->RSQrq+stage, nu[stage], 0);

	return;
	}



void OCP_QP_SET_R(int stage, REAL *R, struct OCP_QP *qp)
	{
	// extract dim
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;

	CVT_MAT2STRMAT(nu[stage], nu[stage], R, nu[stage], qp->RSQrq+stage, 0, 0);

	return;
	}



void OCP_QP_SET_QVEC(int stage, REAL *q, struct OCP_QP *qp)
	{
	// extract dim
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;

 	CVT_TRAN_MAT2STRMAT(nx[stage], 1, q, nx[stage], &(qp->RSQrq[stage]), nu[stage]+nx[stage], nu[stage]); // TODO remove ???
	CVT_VEC2STRVEC(nx[stage], q, qp->rqz+stage, nu[stage]);

	return;
	}



void OCP_QP_SET_RVEC(int stage, REAL *r, struct OCP_QP *qp)
	{
	// extract dim
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;

	CVT_TRAN_MAT2STRMAT(nu[stage], 1, r, nu[stage], &(qp->RSQrq[stage]), nu[stage]+nx[stage], 0); // TODO remove ???
	CVT_VEC2STRVEC(nu[stage], r, qp->rqz+stage, 0);

	return;
	}



void OCP_QP_SET_LB(int stage, REAL *lb, struct OCP_QP *qp)
	{
	// extract dim
	int *nb = qp->dim->nb;

	CVT_VEC2STRVEC(nb[stage], lb, qp->d+stage, 0);

	return;
	}



void OCP_QP_SET_LB_MASK(int stage, REAL *lb_mask, struct OCP_QP *qp)
	{
	// extract dim
	int *nb = qp->dim->nb;

	CVT_VEC2STRVEC(nb[stage], lb_mask, qp->d_mask+stage, 0);

	return;
	}



void OCP_QP_SET_LBX(int stage, REAL *lbx, struct OCP_QP *qp)
	{
	// extract dim
	int *nbu = qp->dim->nbu;
	int *nbx = qp->dim->nbx;

	CVT_VEC2STRVEC(nbx[stage], lbx, qp->d+stage, nbu[stage]);

	return;
	}



void OCP_QP_SET_EL_LBX(int stage, int index, REAL *elem, struct OCP_QP *qp)
	{
	// extract dim
	int *nbu = qp->dim->nbu;

#ifdef DOUBLE_PRECISION
	BLASFEO_DVECEL(qp->d+stage, nbu[stage]+index) = *elem;
#else
	BLASFEO_SVECEL(qp->d+stage, nbu[stage]+index) = *elem;
#endif

	return;
	}



void OCP_QP_SET_LBX_MASK(int stage, REAL *lbx_mask, struct OCP_QP *qp)
	{
	// extract dim
	int *nbu = qp->dim->nbu;
	int *nbx = qp->dim->nbx;

	CVT_VEC2STRVEC(nbx[stage], lbx_mask, qp->d_mask+stage, nbu[stage]);

	return;
	}



void OCP_QP_SET_LBU(int stage, REAL *lbu, struct OCP_QP *qp)
	{
	// extract dim
	int *nbu = qp->dim->nbu;

	CVT_VEC2STRVEC(nbu[stage], lbu, qp->d+stage, 0);

	return;
	}



void OCP_QP_SET_LBU_MASK(int stage, REAL *lbu_mask, struct OCP_QP *qp)
	{
	// extract dim
	int *nbu = qp->dim->nbu;

	CVT_VEC2STRVEC(nbu[stage], lbu_mask, qp->d_mask+stage, 0);

	return;
	}



void OCP_QP_SET_UB(int stage, REAL *ub, struct OCP_QP *qp)
	{
	// extract dim
	int *nb = qp->dim->nb;
	int *ng = qp->dim->ng;

	CVT_VEC2STRVEC(nb[stage], ub, qp->d+stage, nb[stage]+ng[stage]);
	VECSC(nb[stage], -1.0, qp->d+stage, nb[stage]+ng[stage]);

	return;
	}



void OCP_QP_SET_UB_MASK(int stage, REAL *ub_mask, struct OCP_QP *qp)
	{
	// extract dim
	int *nb = qp->dim->nb;
	int *ng = qp->dim->ng;

	CVT_VEC2STRVEC(nb[stage], ub_mask, qp->d_mask+stage, nb[stage]+ng[stage]);

	return;
	}



void OCP_QP_SET_UBX(int stage, REAL *ubx, struct OCP_QP *qp)
	{
	// extract dim
	int *nb = qp->dim->nb;
	int *nbx = qp->dim->nbx;
	int *nbu = qp->dim->nbu;
	int *ng = qp->dim->ng;

	CVT_VEC2STRVEC(nbx[stage], ubx, qp->d+stage, nb[stage]+ng[stage]+nbu[stage]);
	VECSC(nbx[stage], -1.0, qp->d+stage, nb[stage]+ng[stage]+nbu[stage]);

	return;
	}



void OCP_QP_SET_EL_UBX(int stage, int index, REAL *elem, struct OCP_QP *qp)
	{
	// extract dim
	int *nb = qp->dim->nb;
	int *nbu = qp->dim->nbu;
	int *ng = qp->dim->ng;

#ifdef DOUBLE_PRECISION
	BLASFEO_DVECEL(qp->d+stage, nb[stage]+ng[stage]+nbu[stage]+index) = - *elem;
#else
	BLASFEO_SVECEL(qp->d+stage, nb[stage]+ng[stage]+nbu[stage]+index) = - *elem;
#endif

	return;
	}



void OCP_QP_SET_UBX_MASK(int stage, REAL *ubx_mask, struct OCP_QP *qp)
	{
	// extract dim
	int *nb = qp->dim->nb;
	int *nbx = qp->dim->nbx;
	int *nbu = qp->dim->nbu;
	int *ng = qp->dim->ng;

	CVT_VEC2STRVEC(nbx[stage], ubx_mask, qp->d_mask+stage, nb[stage]+ng[stage]+nbu[stage]);

	return;
	}



void OCP_QP_SET_UBU(int stage, REAL *ubu, struct OCP_QP *qp)
	{
	// extract dim
	int *nb = qp->dim->nb;
	int *nbu = qp->dim->nbu;
	int *ng = qp->dim->ng;

	CVT_VEC2STRVEC(nbu[stage], ubu, qp->d+stage, nb[stage]+ng[stage]);
	VECSC(nbu[stage], -1.0, qp->d+stage, nb[stage]+ng[stage]);

	return;
	}



void OCP_QP_SET_UBU_MASK(int stage, REAL *ubu_mask, struct OCP_QP *qp)
	{
	// extract dim
	int *nb = qp->dim->nb;
	int *nbu = qp->dim->nbu;
	int *ng = qp->dim->ng;

	CVT_VEC2STRVEC(nbu[stage], ubu_mask, qp->d_mask+stage, nb[stage]+ng[stage]);

	return;
	}



void OCP_QP_SET_IDXB(int stage, int *idxb, struct OCP_QP *qp)
	{
	// extract dim
	int *nb = qp->dim->nb;

	int ii;
	for(ii=0; ii<nb[stage]; ii++)
		qp->idxb[stage][ii] = idxb[ii];

	return;
	}



void OCP_QP_SET_IDXBX(int stage, int *idxbx, struct OCP_QP *qp)
	{
	// extract dim
	int *nu = qp->dim->nu;
	int *nbx = qp->dim->nbx;
	int *nbu = qp->dim->nbu;

	int ii;
	for(ii=0; ii<nbx[stage]; ii++)
		{
		qp->idxb[stage][nbu[stage]+ii] = nu[stage] + idxbx[ii];
		}

	return;
	}



void OCP_QP_SET_JBX(int stage, REAL *Jbx, struct OCP_QP *qp)
	{
	// extract dim
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;
	int *nbx = qp->dim->nbx;
	int *nbu = qp->dim->nbu;

	int ii, jj, jj0;
	for(ii=0; ii<nbx[stage]; ii++)
		{
		jj0 = -1;
		for(jj=0; jj<nx[stage]; jj++)
			{
			if(jj0==-1 & Jbx[ii+jj*nbx[stage]]!=0.0)
				{
				jj0 = jj;
				qp->idxb[stage][nbu[stage]+ii] = nu[stage]+jj;
				}
			}
		}
	return;
	}



void OCP_QP_SET_IDXBU(int stage, int *idxbx, struct OCP_QP *qp)
	{
	// extract dim
	int *nbu = qp->dim->nbu;

	int ii;
	for(ii=0; ii<nbu[stage]; ii++)
		{
		qp->idxb[stage][ii] = idxbx[ii];
		}

	return;
	}



void OCP_QP_SET_JBU(int stage, REAL *Jbu, struct OCP_QP *qp)
	{
	// extract dim
	int *nu = qp->dim->nu;
	int *nbu = qp->dim->nbu;

	int ii, jj, jj0;
	for(ii=0; ii<nbu[stage]; ii++)
		{
		jj0 = -1;
		for(jj=0; jj<nu[stage]; jj++)
			{
			if(jj0==-1 & Jbu[ii+jj*nbu[stage]]!=0.0)
				{
				jj0 = jj;
				qp->idxb[stage][ii] = jj;
				}
			}
		}
	return;
	}



void OCP_QP_SET_C(int stage, REAL *C, struct OCP_QP *qp)
	{
	// extract dim
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;
	int *ng = qp->dim->ng;

	CVT_TRAN_MAT2STRMAT(ng[stage], nx[stage], C, ng[stage], qp->DCt+stage, nu[stage], 0);

	return;
	}



void OCP_QP_SET_D(int stage, REAL *D, struct OCP_QP *qp)
	{
	// extract dim
	int *nu = qp->dim->nu;
	int *ng = qp->dim->ng;

	CVT_TRAN_MAT2STRMAT(ng[stage], nu[stage], D, ng[stage], qp->DCt+stage, 0, 0);

	return;
	}



void OCP_QP_SET_LG(int stage, REAL *lg, struct OCP_QP *qp)
	{
	// extract dim
	int *nb = qp->dim->nb;
	int *ng = qp->dim->ng;

	CVT_VEC2STRVEC(ng[stage], lg, qp->d+stage, nb[stage]);

	return;
	}



void OCP_QP_SET_LG_MASK(int stage, REAL *lg_mask, struct OCP_QP *qp)
	{
	// extract dim
	int *nb = qp->dim->nb;
	int *ng = qp->dim->ng;

	CVT_VEC2STRVEC(ng[stage], lg_mask, qp->d_mask+stage, nb[stage]);

	return;
	}



void OCP_QP_SET_UG(int stage, REAL *ug, struct OCP_QP *qp)
	{
	// extract dim
	int *nb = qp->dim->nb;
	int *ng = qp->dim->ng;

	CVT_VEC2STRVEC(ng[stage], ug, qp->d+stage, 2*nb[stage]+ng[stage]);
	VECSC(ng[stage], -1.0, qp->d+stage, 2*nb[stage]+ng[stage]);

	return;
	}



void OCP_QP_SET_UG_MASK(int stage, REAL *ug_mask, struct OCP_QP *qp)
	{
	// extract dim
	int *nb = qp->dim->nb;
	int *ng = qp->dim->ng;

	CVT_VEC2STRVEC(ng[stage], ug_mask, qp->d_mask+stage, 2*nb[stage]+ng[stage]);

	return;
	}



void OCP_QP_SET_ZL(int stage, REAL *Zl, struct OCP_QP *qp)
	{
	// extract dim
	int *ns = qp->dim->ns;

	CVT_VEC2STRVEC(ns[stage], Zl, qp->Z+stage, 0);

	return;
	}



void OCP_QP_SET_ZU(int stage, REAL *Zu, struct OCP_QP *qp)
	{
	// extract dim
	int *ns = qp->dim->ns;

	CVT_VEC2STRVEC(ns[stage], Zu, qp->Z+stage, ns[stage]);

	return;
	}



void OCP_QP_SET_ZLVEC(int stage, REAL *zl, struct OCP_QP *qp)
	{
	// extract dim
	int *nu = qp->dim->nu;
	int *nx = qp->dim->nx;
	int *ns = qp->dim->ns;

	CVT_VEC2STRVEC(ns[stage], zl, qp->rqz+stage, nu[stage]+nx[stage]);

	return;
	}



void OCP_QP_SET_ZUVEC(int stage, REAL *zu, struct OCP_QP *qp)
	{
	// extract dim
	int *nu = qp->dim->nu;
	int *nx = qp->dim->nx;
	int *ns = qp->dim->ns;

	CVT_VEC2STRVEC(ns[stage], zu, qp->rqz+stage, nu[stage]+nx[stage]+ns[stage]);

	return;
	}



void OCP_QP_SET_IDXS(int stage, int *idxs, struct OCP_QP *qp)
	{
	// extract dim
	int *ns = qp->dim->ns;

	int ii;
	for(ii=0; ii<ns[stage]; ii++)
		{
		qp->idxs_rev[stage][idxs[ii]] = ii;
		}

	return;
	}



void OCP_QP_SET_IDXS_REV(int stage, int *idxs_rev, struct OCP_QP *qp)
	{
	// extract dim
	int *nb = qp->dim->nb;
	int *ng = qp->dim->ng;

	int ii;
	for(ii=0; ii<nb[stage]+ng[stage]; ii++)
		{
		qp->idxs_rev[stage][ii] = idxs_rev[ii];
		}

	return;
	}



void OCP_QP_SET_JSBU(int stage, REAL *Jsbu, struct OCP_QP *qp)
	{
	// extract dim
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;
	int *nb = qp->dim->nb;
	int *nbx = qp->dim->nbx;
	int *nbu = qp->dim->nbu;
	int *ng = qp->dim->ng;
	int *ns = qp->dim->ns;

	int ii, jj, jj0, idx_tmp;
	// compute nbu part of idxs_rev
	for(ii=0; ii<nbu[stage]; ii++)
		{
		jj0 = -1;
		for(jj=0; jj<ns[stage]; jj++)
			{
			if(jj0==-1 & Jsbu[ii+jj*nbu[stage]]!=0.0)
				{
				jj0 = jj;
				qp->idxs_rev[stage][0+ii] = jj;
				}
			}
		}
	return;
	}



void OCP_QP_SET_JSBX(int stage, REAL *Jsbx, struct OCP_QP *qp)
	{
	// extract dim
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;
	int *nb = qp->dim->nb;
	int *nbx = qp->dim->nbx;
	int *nbu = qp->dim->nbu;
	int *ng = qp->dim->ng;
	int *ns = qp->dim->ns;

	int ii, jj, jj0, idx_tmp;
	// compute nbx part of idxs_rev
	for(ii=0; ii<nbx[stage]; ii++)
		{
		jj0 = -1;
		for(jj=0; jj<ns[stage]; jj++)
			{
			if(jj0==-1 & Jsbx[ii+jj*nbx[stage]]!=0.0)
				{
				jj0 = jj;
				qp->idxs_rev[stage][nbu[stage]+ii] = jj;
				}
			}
		}
	return;
	}



void OCP_QP_SET_JSG(int stage, REAL *Jsg, struct OCP_QP *qp)
	{
	// extract dim
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;
	int *nb = qp->dim->nb;
	int *nbx = qp->dim->nbx;
	int *nbu = qp->dim->nbu;
	int *ng = qp->dim->ng;
	int *ns = qp->dim->ns;

	int ii, jj, jj0, idx_tmp;
	// compute nbx part of idxs_rev
	for(ii=0; ii<ng[stage]; ii++)
		{
		jj0 = -1;
		for(jj=0; jj<ns[stage]; jj++)
			{
			if(jj0==-1 & Jsg[ii+jj*ng[stage]]!=0.0)
				{
				jj0 = jj;
				qp->idxs_rev[stage][nb[stage]+ii] = jj;
				}
			}
		}
	return;
	}



void OCP_QP_SET_LLS(int stage, REAL *ls, struct OCP_QP *qp)
	{
	// extract dim
	int *nb = qp->dim->nb;
	int *ng = qp->dim->ng;
	int *ns = qp->dim->ns;

	CVT_VEC2STRVEC(ns[stage], ls, qp->d+stage, 2*nb[stage]+2*ng[stage]);

	return;
	}



void OCP_QP_SET_LLS_MASK(int stage, REAL *ls_mask, struct OCP_QP *qp)
	{
	// extract dim
	int *nb = qp->dim->nb;
	int *ng = qp->dim->ng;
	int *ns = qp->dim->ns;

	CVT_VEC2STRVEC(ns[stage], ls_mask, qp->d_mask+stage, 2*nb[stage]+2*ng[stage]);

	return;
	}



void OCP_QP_SET_LUS(int stage, REAL *us, struct OCP_QP *qp)
	{
	// extract dim
	int *nb = qp->dim->nb;
	int *ng = qp->dim->ng;
	int *ns = qp->dim->ns;

	CVT_VEC2STRVEC(ns[stage], us, qp->d+stage, 2*nb[stage]+2*ng[stage]+ns[stage]);

	return;
	}



void OCP_QP_SET_LUS_MASK(int stage, REAL *lus_mask, struct OCP_QP *qp)
	{
	// extract dim
	int *nb = qp->dim->nb;
	int *ng = qp->dim->ng;
	int *ns = qp->dim->ns;

	CVT_VEC2STRVEC(ns[stage], lus_mask, qp->d_mask+stage, 2*nb[stage]+2*ng[stage]+ns[stage]);

	return;
	}



void OCP_QP_SET_IDXE(int stage, int *idxe, struct OCP_QP *qp)
	{
	// extract dim
	int *nbxe = qp->dim->nbxe;
	int *nbue = qp->dim->nbue;
	int *nge = qp->dim->nge;

	int ii;
	for(ii=0; ii<nbxe[stage]+nbue[stage]+nge[stage]; ii++)
		qp->idxe[stage][ii] = idxe[ii];

	return;
	}



void OCP_QP_SET_IDXBUE(int stage, int *idxbue, struct OCP_QP *qp)
	{
	// extract dim
	int *nbxe = qp->dim->nbxe;
	int *nbue = qp->dim->nbue;
	int *nge = qp->dim->nge;

	int ii;
	for(ii=0; ii<nbue[stage]; ii++)
		qp->idxe[stage][ii] = idxbue[ii];

	return;
	}



void OCP_QP_SET_IDXBXE(int stage, int *idxbxe, struct OCP_QP *qp)
	{
	// extract dim
	int *nbu = qp->dim->nbu;
	int *nbxe = qp->dim->nbxe;
	int *nbue = qp->dim->nbue;
	int *nge = qp->dim->nge;

	int ii;
	for(ii=0; ii<nbxe[stage]; ii++)
		qp->idxe[stage][nbue[stage]+ii] = nbu[stage] + idxbxe[ii];

	return;
	}



void OCP_QP_SET_IDXGE(int stage, int *idxge, struct OCP_QP *qp)
	{
	// extract dim
	int *nbu = qp->dim->nbu;
	int *nbx = qp->dim->nbx;
	int *nbxe = qp->dim->nbxe;
	int *nbue = qp->dim->nbue;
	int *nge = qp->dim->nge;

	int ii;
	for(ii=0; ii<nge[stage]; ii++)
		qp->idxe[stage][nbue[stage]+nbxe[stage]+ii] = nbu[stage] + nbx[stage] + idxge[ii];

	return;
	}



void OCP_QP_SET_JBUE(int stage, REAL *Jbue, struct OCP_QP *qp)
	{
	// extract dim
	int *nbu = qp->dim->nbu;
	int *nbue = qp->dim->nbue;

	int ii, idx;
	idx = 0;
	for(ii=0; ii<nbu[stage]; ii++)
		{
		if(idx<nbue[stage] | Jbue[ii]!=0.0)
			{
			qp->idxe[stage][idx] = ii;
			idx++;
			}
		}
	return;
	}



void OCP_QP_SET_JBXE(int stage, REAL *Jbxe, struct OCP_QP *qp)
	{
	// extract dim
	int *nbu = qp->dim->nbu;
	int *nbx = qp->dim->nbx;
	int *nbue = qp->dim->nbue;
	int *nbxe = qp->dim->nbxe;

	int ii, idx;
	idx = 0;
	for(ii=0; ii<nbx[stage]; ii++)
		{
		if(idx<nbxe[stage] | Jbxe[ii]!=0.0)
			{
			qp->idxe[stage][nbue[stage]+idx] = nbu[stage] + ii;
			idx++;
			}
		}
	return;
	}



void OCP_QP_SET_JGE(int stage, REAL *Jge, struct OCP_QP *qp)
	{
	// extract dim
	int *nbu = qp->dim->nbu;
	int *nbx = qp->dim->nbx;
	int *ng = qp->dim->ng;
	int *nbue = qp->dim->nbue;
	int *nbxe = qp->dim->nbxe;
	int *nge = qp->dim->nge;

	int ii, idx;
	idx = 0;
	for(ii=0; ii<ng[stage]; ii++)
		{
		if(idx<nge[stage] | Jge[ii]!=0.0)
			{
			qp->idxe[stage][nbue[stage]+nbxe[stage]+idx] = nbu[stage] + nbx[stage] + ii;
			idx++;
			}
		}
	return;
	}



void OCP_QP_GET(char *field, int stage, struct OCP_QP *qp, void *value)
	{
	// matrices
	if(hpipm_strcmp(field, "A")) 
		{
		OCP_QP_GET_A(stage, qp, value);
		}
	// vectors
	else if(hpipm_strcmp(field, "lbx") | hpipm_strcmp(field, "lx"))
		{ 
		OCP_QP_GET_LBX(stage, qp, value);
		}
	else if(hpipm_strcmp(field, "ubx") | hpipm_strcmp(field, "ux"))
		{ 
		OCP_QP_GET_UBX(stage, qp, value);
		}
	// int
	else
		{
		printf("error: OCP_QP_GET: wrong field %s\n", field);
		exit(1);	
		}
	return;
	}



void OCP_QP_GET_A(int stage, struct OCP_QP *qp, REAL *A)
	{
	// extract dim
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;

	CVT_TRAN_STRMAT2MAT(nx[stage], nx[stage+1], qp->BAbt+stage, nu[stage], 0, A, nx[stage+1]);

	return;
	}



void OCP_QP_GET_B(int stage, struct OCP_QP *qp, REAL *B)
	{
	// extract dim
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;

	CVT_TRAN_STRMAT2MAT(nu[stage], nx[stage+1], qp->BAbt+stage, 0, 0, B, nx[stage+1]);

	return;
	}



void OCP_QP_GET_BVEC(int stage, struct OCP_QP *qp, REAL *b)
	{
	// extract dim
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;

	CVT_STRVEC2VEC(nx[stage+1], qp->b+stage, 0, b);

	return;
	}



void OCP_QP_GET_Q(int stage, struct OCP_QP *qp, REAL *Q)
	{
	// extract dim
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;

	CVT_STRMAT2MAT(nx[stage], nx[stage], qp->RSQrq+stage, nu[stage], nu[stage], Q, nx[stage]);

	return;
	}



void OCP_QP_GET_S(int stage, struct OCP_QP *qp, REAL *S)
	{
	// extract dim
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;

	CVT_TRAN_STRMAT2MAT(nx[stage], nu[stage], qp->RSQrq+stage, nu[stage], 0, S, nu[stage]);

	return;
	}



void OCP_QP_GET_R(int stage, struct OCP_QP *qp, REAL *R)
	{
	// extract dim
	int *nu = qp->dim->nu;

	CVT_STRMAT2MAT(nu[stage], nu[stage], qp->RSQrq+stage, 0, 0, R, nu[stage]);

	return;
	}



void OCP_QP_GET_QVEC(int stage, struct OCP_QP *qp, REAL *q)
	{
	// extract dim
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;

	CVT_STRVEC2VEC(nx[stage], qp->rqz+stage, nu[stage], q);

	return;
	}



void OCP_QP_GET_RVEC(int stage, struct OCP_QP *qp, REAL *r)
	{
	// extract dim
	int *nu = qp->dim->nu;

	CVT_STRVEC2VEC(nu[stage], qp->rqz+stage, 0, r);

	return;
	}



void OCP_QP_GET_LB(int stage, struct OCP_QP *qp, REAL *lb)
	{
	// extract dim
	int *nb = qp->dim->nb;

	int i;

	CVT_STRVEC2VEC(nb[stage], qp->d+stage, 0, lb);

	return;
	}



void OCP_QP_GET_UB(int stage, struct OCP_QP *qp, REAL *ub)
	{
	// extract dim
	int *nb = qp->dim->nb;
	int *ng = qp->dim->ng;

	int i;

	CVT_STRVEC2VEC(nb[stage], qp->d+stage, nb[stage]+ng[stage], ub);
	for(i=0; i<nb[stage]; i++)
		{
		ub[i] = -ub[i];
		}

	return;
	}



void OCP_QP_GET_LBX(int stage, struct OCP_QP *qp, REAL *lbx)
	{
	// extract dim
	int *nbu = qp->dim->nbu;
	int *nbx = qp->dim->nbx;

	CVT_STRVEC2VEC(nbx[stage], qp->d+stage, nbu[stage], lbx);

	return;
	}



void OCP_QP_GET_UBX(int stage, struct OCP_QP *qp, REAL *ubx)
	{
	// extract dim
	int *nb = qp->dim->nb;
	int *nbx = qp->dim->nbx;
	int *nbu = qp->dim->nbu;
	int *ng = qp->dim->ng;

	int i;

	CVT_STRVEC2VEC(nbx[stage], qp->d+stage, nb[stage]+ng[stage]+nbu[stage], ubx);
	for(i=0; i<nbx[stage]; i++)
		{
		ubx[i] = -ubx[i];
		}

	return;
	}



void OCP_QP_GET_LBU(int stage, struct OCP_QP *qp, REAL *lbu)
	{
	// extract dim
	int *nbu = qp->dim->nbu;

	CVT_STRVEC2VEC(nbu[stage], qp->d+stage, 0, lbu);

	return;
	}



void OCP_QP_GET_UBU(int stage, struct OCP_QP *qp, REAL *ubu)
	{
	// extract dim
	int *nb = qp->dim->nb;
	int *nbu = qp->dim->nbu;
	int *ng = qp->dim->ng;

	int i;

	CVT_STRVEC2VEC(nbu[stage], qp->d+stage, nb[stage]+ng[stage], ubu);
	for(i=0; i<nbu[stage]; i++)
		{
		ubu[i] = -ubu[i];
		}

	return;
	}



void OCP_QP_GET_IDXB(int stage, struct OCP_QP *qp, int *idxb)
	{
	// extract dim
	int *nb = qp->dim->nb;

	int ii;
	for(ii=0; ii<nb[stage]; ii++)
		idxb[ii] = qp->idxb[stage][ii];

	return;
	}



//void OCP_QP_GET_IDXBX(int stage, struct OCP_QP *qp, int *idxb)
//	{
//	TODO
//	return;
//	}



//void OCP_QP_GET_JBX(int stage, struct OCP_QP *qp, int *Jbx)
//	{
//	TODO
//	return;
//	}



//void OCP_QP_GET_IDXBU(int stage, struct OCP_QP *qp, int *idxbu)
//	{
//	TODO
//	return;
//	}



//void OCP_QP_GET_JBU(int stage, struct OCP_QP *qp, int *Jbu)
//	{
//	TODO
//	return;
//	}



void OCP_QP_GET_C(int stage, struct OCP_QP *qp, REAL *C)
	{
	// extract dim
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;
	int *ng = qp->dim->ng;

	CVT_TRAN_STRMAT2MAT(nx[stage], ng[stage], qp->DCt+stage, nu[stage], 0, C, ng[stage]);

	return;
	}



void OCP_QP_GET_D(int stage, struct OCP_QP *qp, REAL *D)
	{
	// extract dim
	int *nu = qp->dim->nu;
	int *ng = qp->dim->ng;

	CVT_TRAN_STRMAT2MAT(nu[stage], ng[stage], qp->DCt+stage, 0, 0, D, ng[stage]);

	return;
	}



void OCP_QP_GET_LG(int stage, struct OCP_QP *qp, REAL *lg)
	{
	// extract dim
	int *nb = qp->dim->nb;
	int *ng = qp->dim->ng;

	CVT_STRVEC2VEC(ng[stage], qp->d+stage, nb[stage], lg);

	return;
	}



void OCP_QP_GET_UG(int stage, struct OCP_QP *qp, REAL *ug)
	{
	// extract dim
	int *nb = qp->dim->nb;
	int *ng = qp->dim->ng;

	int i;

	CVT_STRVEC2VEC(ng[stage], qp->d+stage, 2*nb[stage]+ng[stage], ug);
	for(i=0; i<ng[stage]; i++)
		{
		ug[i] = -ug[i];
		}

	return;
	}



void OCP_QP_GET_ZL(int stage, struct OCP_QP *qp, REAL *Zl)
	{
	// extract dim
	int *ns = qp->dim->ns;

	CVT_STRVEC2VEC(ns[stage], qp->Z+stage, 0, Zl);

	return;
	}



void OCP_QP_GET_ZU(int stage, struct OCP_QP *qp, REAL *Zu)
	{
	// extract dim
	int *ns = qp->dim->ns;

	CVT_STRVEC2VEC(ns[stage], qp->Z+stage, ns[stage], Zu);

	return;
	}



void OCP_QP_GET_ZLVEC(int stage, struct OCP_QP *qp, REAL *zl)
	{
	// extract dim
	int *nu = qp->dim->nu;
	int *nx = qp->dim->nx;
	int *ns = qp->dim->ns;

	CVT_STRVEC2VEC(ns[stage], qp->rqz+stage, nu[stage]+nx[stage], zl);

	return;
	}



void OCP_QP_GET_ZUVEC(int stage, struct OCP_QP *qp, REAL *zu)
	{
	// extract dim
	int *nu = qp->dim->nu;
	int *nx = qp->dim->nx;
	int *ns = qp->dim->ns;

	CVT_STRVEC2VEC(ns[stage], qp->rqz+stage, nu[stage]+nx[stage]+ns[stage], zu);

	return;
	}



// XXX only valid if there is one slack per softed constraint !!!
void OCP_QP_GET_IDXS(int stage, struct OCP_QP *qp, int *idxs)
	{
	// extract dim
	int *nb = qp->dim->nb;
	int *ng = qp->dim->ng;
	int *ns = qp->dim->ns;

	int ii, idx_tmp;

	for(ii=0; ii<nb[stage]+ng[stage]; ii++)
		{
		idx_tmp = qp->idxs_rev[stage][ii];
		if(idx_tmp!=-1)
			{
			idxs[idx_tmp] = ii;
			}
		}

	return;
	}



void OCP_QP_GET_IDXS_REV(int stage, struct OCP_QP *qp, int *idxs_rev)
	{
	// extract dim
	int *nb = qp->dim->nb;
	int *ng = qp->dim->ng;

	int ii;
	for(ii=0; ii<nb[stage]+ng[stage]; ii++)
		idxs_rev[ii] = qp->idxs_rev[stage][ii];

	return;
	}



//void OCP_QP_GET_JSBX(int stage, struct OCP_QP *qp, int *Jsbx)
//	{
//	TODO
//	return;
//	}



//void OCP_QP_GET_JSBX(int stage, struct OCP_QP *qp, int *Jsbx)
//	{
//	TODO
//	return;
//	}



//void OCP_QP_GET_JSG(int stage, struct OCP_QP *qp, int *Jsg)
//	{
//	TODO
//	return;
//	}



void OCP_QP_GET_LLS(int stage, struct OCP_QP *qp, REAL *ls)
	{
	// extract dim
	int *nb = qp->dim->nb;
	int *ng = qp->dim->ng;
	int *ns = qp->dim->ns;

	CVT_STRVEC2VEC(ns[stage], qp->d+stage, 2*nb[stage]+2*ng[stage], ls);

	return;
	}



void OCP_QP_GET_LUS(int stage, struct OCP_QP *qp, REAL *us)
	{
	// extract dim
	int *nb = qp->dim->nb;
	int *ng = qp->dim->ng;
	int *ns = qp->dim->ns;

	int i;

	CVT_STRVEC2VEC(ns[stage], qp->d+stage, 2*nb[stage]+2*ng[stage]+ns[stage], us);

	return;
	}



