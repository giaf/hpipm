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



int MEMSIZE_OCP_QP(int N, int *nx, int *nu, int *nb, int *ng)
	{

	int ii;

	int size = 0;

	size += 4*(N+1)*sizeof(int); // nx nu nb ng
	size += (N+1)*sizeof(int *); // idxb
	size += (1*N+2*(N+1))*sizeof(struct STRMAT); // BAbt ; RSqrq DCt
	size += (1*N+5*(N+1))*sizeof(struct STRVEC); // b ; rq d_lb d_ub d_lg d_ug

	for(ii=0; ii<N; ii++)
		{
		size += nb[ii]*sizeof(int); // idxb
		size += SIZE_STRMAT(nu[ii]+nx[ii]+1, nx[ii+1]); // BAbt
		size += SIZE_STRVEC(nx[ii+1]); // b
		size += SIZE_STRMAT(nu[ii]+nx[ii]+1, nu[ii]+nx[ii]); // RSQrq
		size += SIZE_STRVEC(nu[ii]+nx[ii]); // rq
		size += SIZE_STRMAT(nu[ii]+nx[ii], ng[ii]); // DCt
		size += 2*SIZE_STRVEC(nb[ii]); // d_lb d_ub
		size += 2*SIZE_STRVEC(ng[ii]); // d_lg d_ug
		}
	ii = N;
	size += nb[ii]*sizeof(int); // idxb
	size += SIZE_STRMAT(nu[ii]+nx[ii]+1, nu[ii]+nx[ii]); // RSQrq
	size += SIZE_STRVEC(nu[ii]+nx[ii]); // rq
	size += SIZE_STRMAT(nu[ii]+nx[ii], ng[ii]); // DCt
	size += 2*SIZE_STRVEC(nb[ii]); // d_lb d_ub
	size += 2*SIZE_STRVEC(ng[ii]); // d_lg d_ug

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 64; // align to typical cache line size
	
	return size;

	}



void CREATE_OCP_QP(int N, int *nx, int *nu, int *nb, int *ng, struct OCP_QP *qp, void *memory)
	{

	int ii;


	// memsize
	qp->memsize = MEMSIZE_OCP_QP(N, nx, nu, nb, ng);


	// horizon length
	qp->N = N;


	// int pointer stuff
	int **ip_ptr;
	ip_ptr = (int **) memory;

	// idxb
	qp->idxb = ip_ptr;
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

	// rq
	qp->rq = sv_ptr;
	sv_ptr += N+1;

	// d_lb
	qp->d_lb = sv_ptr;
	sv_ptr += N+1;

	// d_ub
	qp->d_ub = sv_ptr;
	sv_ptr += N+1;

	// d_lg
	qp->d_lg = sv_ptr;
	sv_ptr += N+1;

	// d_ug
	qp->d_ug = sv_ptr;
	sv_ptr += N+1;


	// integer stuff
	int *i_ptr;
	i_ptr = (int *) sv_ptr;

	// nx
	qp->nx = i_ptr;
	for(ii=0; ii<=N; ii++)
		{
		i_ptr[ii] = nx[ii];
		}
	i_ptr += N+1;
	
	// nu
	qp->nu = i_ptr;
	for(ii=0; ii<=N; ii++)
		{
		i_ptr[ii] = nu[ii];
		}
	i_ptr += N+1;
	
	// nb
	qp->nb = i_ptr;
	for(ii=0; ii<=N; ii++)
		{
		i_ptr[ii] = nb[ii];
		}
	i_ptr += N+1;

	// ng
	qp->ng = i_ptr;
	for(ii=0; ii<=N; ii++)
		{
		i_ptr[ii] = ng[ii];
		}
	i_ptr += N+1;
	
	// idxb
	for(ii=0; ii<=N; ii++)
		{
		(qp->idxb)[ii] = i_ptr;
		i_ptr += nb[ii];
		}


	// align to typical cache line size
	long long l_ptr = (long long) i_ptr;
	l_ptr = (l_ptr+63)/64*64;


	// double stuff
	char *c_ptr;
	c_ptr = (char *) l_ptr;

	// BAbt
	for(ii=0; ii<N; ii++)
		{
		CREATE_STRMAT(nu[ii]+nx[ii]+1, nx[ii+1], qp->BAbt+ii, c_ptr);
		c_ptr += (qp->BAbt+ii)->memory_size;
		}

	// RSQrq
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRMAT(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], qp->RSQrq+ii, c_ptr);
		c_ptr += (qp->RSQrq+ii)->memory_size;
		}

	// DCt
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRMAT(nu[ii]+nx[ii], ng[ii], qp->DCt+ii, c_ptr);
		c_ptr += (qp->DCt+ii)->memory_size;
		}

	// b
	for(ii=0; ii<N; ii++)
		{
		CREATE_STRVEC(nx[ii+1], qp->b+ii, c_ptr);
		c_ptr += (qp->b+ii)->memory_size;
		}

	// rq
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(nu[ii]+nx[ii], qp->rq+ii, c_ptr);
		c_ptr += (qp->rq+ii)->memory_size;
		}

	// d_lb
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(nb[ii], qp->d_lb+ii, c_ptr);
		c_ptr += (qp->d_lb+ii)->memory_size;
		}

	// d_ub
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(nb[ii], qp->d_ub+ii, c_ptr);
		c_ptr += (qp->d_ub+ii)->memory_size;
		}

	// d_lg
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(ng[ii], qp->d_lg+ii, c_ptr);
		c_ptr += (qp->d_lg+ii)->memory_size;
		}

	// d_ug
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(ng[ii], qp->d_ug+ii, c_ptr);
		c_ptr += (qp->d_ug+ii)->memory_size;
		}

	return;

	}



void CAST_OCP_QP(int N, int *nx, int *nu, int *nb, int **idxb, int *ng, struct STRMAT *BAbt, struct STRVEC *b, struct STRMAT *RSQrq, struct STRVEC *rq, struct STRMAT *DCt, struct STRVEC *d_lb, struct STRVEC *d_ub, struct STRVEC *d_lg, struct STRVEC *d_ug, struct OCP_QP *qp)
	{

	qp->N = N;
	qp->nx = nx;
	qp->nu = nu;
	qp->nb = nb;
	qp->idxb = idxb;
	qp->ng = ng;
	qp->BAbt = BAbt;
	qp->b = b;
	qp->RSQrq = RSQrq;
	qp->rq = rq;
	qp->DCt = DCt;
	qp->d_lb = d_lb;
	qp->d_ub = d_ub;
	qp->d_lg = d_lg;
	qp->d_ug = d_ug;

	return;

	}



void CVT_COLMAJ_TO_OCP_QP(REAL **A, REAL **B, REAL **b, REAL **Q, REAL **S, REAL **R, REAL **q, REAL **r, int **idxb, REAL **d_lb, REAL **d_ub, REAL **C, REAL **D, REAL **d_lg, REAL **d_ug, struct OCP_QP *qp)
	{

	int N = qp->N;
	int *nx = qp->nx;
	int *nu = qp->nu;
	int *nb = qp->nb;
	int *ng = qp->ng;

	int ii, jj;

	for(ii=0; ii<N; ii++)
		{
		CVT_TRAN_MAT2STRMAT(nx[ii+1], nu[ii], B[ii], nx[ii+1], qp->BAbt+ii, 0, 0);
		CVT_TRAN_MAT2STRMAT(nx[ii+1], nx[ii], A[ii], nx[ii+1], qp->BAbt+ii, nu[ii], 0);
		CVT_TRAN_MAT2STRMAT(nx[ii+1], 1, b[ii], nx[ii+1], qp->BAbt+ii, nu[ii]+nx[ii], 0);
		CVT_VEC2STRVEC(nx[ii+1], b[ii], qp->b+ii, 0);
		}
	
	for(ii=0; ii<=N; ii++)
		{
		CVT_MAT2STRMAT(nu[ii], nu[ii], R[ii], nu[ii], qp->RSQrq+ii, 0, 0);
		CVT_TRAN_MAT2STRMAT(nu[ii], nx[ii], S[ii], nu[ii], qp->RSQrq+ii, nu[ii], 0);
		CVT_MAT2STRMAT(nx[ii], nx[ii], Q[ii], nx[ii], qp->RSQrq+ii, nu[ii], nu[ii]);
		CVT_TRAN_MAT2STRMAT(nu[ii], 1, r[ii], nu[ii], qp->RSQrq+ii, nu[ii]+nx[ii], 0);
		CVT_TRAN_MAT2STRMAT(nx[ii], 1, q[ii], nx[ii], qp->RSQrq+ii, nu[ii]+nx[ii], nu[ii]);
		CVT_VEC2STRVEC(nu[ii], r[ii], qp->rq+ii, 0);
		CVT_VEC2STRVEC(nx[ii], q[ii], qp->rq+ii, nu[ii]);
		}
	
	for(ii=0; ii<=N; ii++)
		{
		for(jj=0; jj<nb[ii]; jj++)
			qp->idxb[ii][jj] = idxb[ii][jj];
		CVT_VEC2STRVEC(nb[ii], d_lb[ii], qp->d_lb+ii, 0);
		CVT_VEC2STRVEC(nb[ii], d_ub[ii], qp->d_ub+ii, 0);
		}
	
	for(ii=0; ii<=N; ii++)
		{
		CVT_TRAN_MAT2STRMAT(ng[ii], nu[ii], D[ii], ng[ii], qp->DCt+ii, 0, 0);
		CVT_TRAN_MAT2STRMAT(ng[ii], nx[ii], C[ii], ng[ii], qp->DCt+ii, nu[ii], 0);
		CVT_VEC2STRVEC(ng[ii], d_lg[ii], qp->d_lg+ii, 0);
		CVT_VEC2STRVEC(ng[ii], d_ug[ii], qp->d_ug+ii, 0);
		}

	return;

	}



void CVT_ROWMAJ_TO_OCP_QP(REAL **A, REAL **B, REAL **b, REAL **Q, REAL **S, REAL **R, REAL **q, REAL **r, int **idxb, REAL **d_lb, REAL **d_ub, REAL **C, REAL **D, REAL **d_lg, REAL **d_ug, struct OCP_QP *qp)
	{

	int N = qp->N;
	int *nx = qp->nx;
	int *nu = qp->nu;
	int *nb = qp->nb;
	int *ng = qp->ng;

	int ii, jj;

	for(ii=0; ii<N; ii++)
		{
		CVT_MAT2STRMAT(nu[ii], nx[ii+1], B[ii], nu[ii], qp->BAbt+ii, 0, 0);
		CVT_MAT2STRMAT(nx[ii], nx[ii+1], A[ii], nx[ii], qp->BAbt+ii, nu[ii], 0);
		CVT_TRAN_MAT2STRMAT(nx[ii+1], 1, b[ii], nx[ii+1], qp->BAbt+ii, nu[ii]+nx[ii], 0);
		CVT_VEC2STRVEC(nx[ii+1], b[ii], qp->b+ii, 0);
		}
	
	for(ii=0; ii<=N; ii++)
		{
		CVT_TRAN_MAT2STRMAT(nu[ii], nu[ii], R[ii], nu[ii], qp->RSQrq+ii, 0, 0);
		CVT_MAT2STRMAT(nx[ii], nu[ii], S[ii], nx[ii], qp->RSQrq+ii, nu[ii], 0);
		CVT_TRAN_MAT2STRMAT(nx[ii], nx[ii], Q[ii], nx[ii], qp->RSQrq+ii, nu[ii], nu[ii]);
		CVT_TRAN_MAT2STRMAT(nu[ii], 1, r[ii], nu[ii], qp->RSQrq+ii, nu[ii]+nx[ii], 0);
		CVT_TRAN_MAT2STRMAT(nx[ii], 1, q[ii], nx[ii], qp->RSQrq+ii, nu[ii]+nx[ii], nu[ii]);
		CVT_VEC2STRVEC(nu[ii], r[ii], qp->rq+ii, 0);
		CVT_VEC2STRVEC(nx[ii], q[ii], qp->rq+ii, nu[ii]);
		}
	
	for(ii=0; ii<=N; ii++)
		{
		for(jj=0; jj<nb[ii]; jj++)
			qp->idxb[ii][jj] = idxb[ii][jj];
		CVT_VEC2STRVEC(nb[ii], d_lb[ii], qp->d_lb+ii, 0);
		CVT_VEC2STRVEC(nb[ii], d_ub[ii], qp->d_ub+ii, 0);
		}
	
	for(ii=0; ii<=N; ii++)
		{
		CVT_MAT2STRMAT(nu[ii], ng[ii], D[ii], nu[ii], qp->DCt+ii, 0, 0);
		CVT_MAT2STRMAT(nx[ii], ng[ii], C[ii], nx[ii], qp->DCt+ii, nu[ii], 0);
		CVT_VEC2STRVEC(ng[ii], d_lg[ii], qp->d_lg+ii, 0);
		CVT_VEC2STRVEC(ng[ii], d_ug[ii], qp->d_ug+ii, 0);
		}

	return;

	}



void COPY_OCP_QP(struct OCP_QP *qp_in, struct OCP_QP *qp_out)
	{

	int N = qp_in->N;
	int *nx = qp_in->nx;
	int *nu = qp_in->nu;
	int *nb = qp_in->nb;
	int *ng = qp_in->ng;

	int ii, jj;

#if defined(RUNTIME_CHECKS)
	if(qp_out->N != qp_in->N)
		{
		printf("\nError : x_copy_ocp_qp : qp_out->N != qp_out->N : %d != %d\n\n", qp_out->N, qp_in->N);
		exit(1);
		}
	for(ii=0; ii<=N; ii++)
		{
		if(qp_out->nx[ii] != qp_in->nx[ii])
			{
			printf("\nError : x_copy_ocp_qp : qp_out->nx[%d] != qp_out->nx[%d] : %d != %d\n\n", ii, ii, qp_out->nx[ii], qp_in->nx[ii]);
			exit(1);
			}
		if(qp_out->nu[ii] != qp_in->nu[ii])
			{
			printf("\nError : x_copy_ocp_qp : qp_out->nu[%d] != qp_out->nu[%d] : %d != %d\n\n", ii, ii, qp_out->nu[ii], qp_in->nu[ii]);
			exit(1);
			}
		if(qp_out->nb[ii] != qp_in->nb[ii])
			{
			printf("\nError : x_copy_ocp_qp : qp_out->nb[%d] != qp_out->nb[%d] : %d != %d\n\n", ii, ii, qp_out->nb[ii], qp_in->nb[ii]);
			exit(1);
			}
		if(qp_out->ng[ii] != qp_in->ng[ii])
			{
			printf("\nError : x_copy_ocp_qp : qp_out->ng[%d] != qp_out->ng[%d] : %d != %d\n\n", ii, ii, qp_out->ng[ii], qp_in->ng[ii]);
			exit(1);
			}
		}
#endif

	for(ii=0; ii<N; ii++)
		{
		for(jj=0; jj<qp_in->nb[ii]; jj++) qp_out->idxb[ii][jj] = qp_in->idxb[ii][jj];
		GECP_LIBSTR(nx[ii]+nu[ii]+1, nx[ii+1], qp_in->BAbt+ii, 0, 0, qp_out->BAbt+ii, 0, 0);
		VECCP_LIBSTR(nx[ii+1], qp_in->b+ii, 0, qp_out->b+ii, 0);
		GECP_LIBSTR(nx[ii]+nu[ii]+1, nu[ii]+nx[ii], qp_in->RSQrq+ii, 0, 0, qp_out->RSQrq+ii, 0, 0);
		VECCP_LIBSTR(nu[ii]+nx[ii], qp_in->rq+ii, 0, qp_out->rq+ii, 0);
		GECP_LIBSTR(nx[ii]+nu[ii], ng[ii], qp_in->DCt+ii, 0, 0, qp_out->DCt+ii, 0, 0);
		VECCP_LIBSTR(nb[ii], qp_in->d_lb+ii, 0, qp_out->d_lb+ii, 0);
		VECCP_LIBSTR(nb[ii], qp_in->d_ub+ii, 0, qp_out->d_ub+ii, 0);
		VECCP_LIBSTR(ng[ii], qp_in->d_lg+ii, 0, qp_out->d_lg+ii, 0);
		VECCP_LIBSTR(ng[ii], qp_in->d_ug+ii, 0, qp_out->d_ug+ii, 0);
		}

	return;

	}


