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



int SIZE_OCP_QP(int N, int *nx, int *nu, int *nb, int *ng)
	{

	int ii;

	int size = 0;

	size += 4*(N+1)*sizeof(int); // nx nu nb ng
	size += (N+1)*sizeof(int *); // idxb
	size += (1*N+2*(N+1))*sizeof(struct STRMAT); // BAbt ; RSqrq DCt
	size += (1*N+5*(N+1))*sizeof(struct STRVEC); // b ; rq lb ub lg ug

	for(ii=0; ii<N; ii++)
		{
		size += nb[ii]*sizeof(int); // idxb
		size += SIZE_STRMAT(nu[ii]+nx[ii]+1, nx[ii+1]); // BAbt
		size += SIZE_STRVEC(nx[ii+1]); // b
		size += SIZE_STRMAT(nu[ii]+nx[ii]+1, nu[ii]+nx[ii]); // RSQrq
		size += SIZE_STRVEC(nu[ii]+nx[ii]); // rq
		size += SIZE_STRMAT(nu[ii]+nx[ii], ng[ii]); // DCt
		size += 2*SIZE_STRVEC(nb[ii]); // lb ub
		size += 2*SIZE_STRVEC(ng[ii]); // lg ug
		}
	ii = N;
	size += nb[ii]*sizeof(int); // idxb
	size += SIZE_STRMAT(nu[ii]+nx[ii]+1, nu[ii]+nx[ii]); // RSQrq
	size += SIZE_STRVEC(nu[ii]+nx[ii]); // rq
	size += SIZE_STRMAT(nu[ii]+nx[ii], ng[ii]); // DCt
	size += 2*SIZE_STRVEC(nb[ii]); // lb ub
	size += 2*SIZE_STRVEC(ng[ii]); // lg ug

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 64; // align to typical cache line size
	
	return size;

	}



void CREATE_OCP_QP(int N, int *nx, int *nu, int *nb, int *ng, struct OCP_QP *str_out, void *memory)
	{

	int ii;


	// horizon length
	str_out->NN = N;


	// int pointer stuff
	int **ip_ptr;
	ip_ptr = (int **) memory;

	// idxb
	str_out->idxb = ip_ptr;
	ip_ptr += N+1;


	// matrix struct stuff
	struct STRMAT *sm_ptr = (struct STRMAT *) ip_ptr;

	// BAbt
	str_out->BAbt = sm_ptr;
	sm_ptr += N;

	// RSQrq
	str_out->RSQrq = sm_ptr;
	sm_ptr += N+1;

	// DCt
	str_out->DCt = sm_ptr;
	sm_ptr += N+1;


	// vector struct stuff
	struct STRVEC *sv_ptr = (struct STRVEC *) sm_ptr;

	// b
	str_out->b = sv_ptr;
	sv_ptr += N;

	// rq
	str_out->rq = sv_ptr;
	sv_ptr += N+1;

	// lb
	str_out->lb = sv_ptr;
	sv_ptr += N+1;

	// ub
	str_out->ub = sv_ptr;
	sv_ptr += N+1;

	// lg
	str_out->lg = sv_ptr;
	sv_ptr += N+1;

	// ug
	str_out->ug = sv_ptr;
	sv_ptr += N+1;


	// integer stuff
	int *i_ptr;
	i_ptr = (int *) sv_ptr;

	// nx
	str_out->nx = i_ptr;
	for(ii=0; ii<=N; ii++)
		{
		i_ptr[ii] = nx[ii];
		}
	i_ptr += N+1;
	
	// nu
	str_out->nu = i_ptr;
	for(ii=0; ii<=N; ii++)
		{
		i_ptr[ii] = nu[ii];
		}
	i_ptr += N+1;
	
	// nb
	str_out->nb = i_ptr;
	for(ii=0; ii<=N; ii++)
		{
		i_ptr[ii] = nb[ii];
		}
	i_ptr += N+1;

	// ng
	str_out->ng = i_ptr;
	for(ii=0; ii<=N; ii++)
		{
		i_ptr[ii] = ng[ii];
		}
	i_ptr += N+1;
	
	// idxb
	for(ii=0; ii<=N; ii++)
		{
		(str_out->idxb)[ii] = i_ptr;
		i_ptr += nb[ii];
		}


	// align to typical cache line size
	long long l_ptr = (long long) i_ptr;
	l_ptr = (l_ptr+63)/64*64;

	// double stuff
	void *v_ptr;
	v_ptr = (void *) l_ptr;

	// BAbt
	for(ii=0; ii<N; ii++)
		{
		CREATE_STRMAT(nu[ii]+nx[ii]+1, nx[ii+1], str_out->BAbt+ii, v_ptr);
		v_ptr += (str_out->BAbt+ii)->memory_size;
		}

	// RSQrq
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRMAT(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], str_out->RSQrq+ii, v_ptr);
		v_ptr += (str_out->RSQrq+ii)->memory_size;
		}

	// DCt
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRMAT(nu[ii]+nx[ii], ng[ii], str_out->DCt+ii, v_ptr);
		v_ptr += (str_out->DCt+ii)->memory_size;
		}

	// b
	for(ii=0; ii<N; ii++)
		{
		CREATE_STRVEC(nx[ii+1], str_out->b+ii, v_ptr);
		v_ptr += (str_out->b+ii)->memory_size;
		}

	// rq
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(nu[ii]+nx[ii], str_out->rq+ii, v_ptr);
		v_ptr += (str_out->rq+ii)->memory_size;
		}

	// lb
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(nb[ii], str_out->lb+ii, v_ptr);
		v_ptr += (str_out->lb+ii)->memory_size;
		}

	// ub
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(nb[ii], str_out->ub+ii, v_ptr);
		v_ptr += (str_out->ub+ii)->memory_size;
		}

	// lg
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(ng[ii], str_out->lg+ii, v_ptr);
		v_ptr += (str_out->lg+ii)->memory_size;
		}

	// ug
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(ng[ii], str_out->ug+ii, v_ptr);
		v_ptr += (str_out->ug+ii)->memory_size;
		}

	return;

	}



void d_cast_ocp_qp(int N, int *nx, int *nu, int *nb, int **idxb, int *ng, struct STRMAT *BAbt, struct STRVEC *b, struct STRMAT *RSQrq, struct STRVEC *rq, struct STRMAT *DCt, struct STRVEC *lb, struct STRVEC *ub, struct STRVEC *lg, struct STRVEC *ug, struct OCP_QP *qp)
	{

	qp->NN = N;
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
	qp->lb = lb;
	qp->ub = ub;
	qp->lg = lg;
	qp->ug = ug;

	return;

	}



void COPY_OCP_QP(struct OCP_QP *qp_in, struct OCP_QP *qp_out)
	{

	int N = qp_in->NN;
	int *nx = qp_in->nx;
	int *nu = qp_in->nu;
	int *nb = qp_in->nb;
	int *ng = qp_in->ng;

	int ii, jj;

#if defined(RUNTIME_CHECKS)
	if(qp_out->NN != qp_in->NN)
		{
		printf("\nError : x_copy_ocp_qp : qp_out->NN != qp_out->NN : %d != %d\n\n", qp_out->NN, qp_in->NN);
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
		VECCP_LIBSTR(nb[ii], qp_in->lb+ii, 0, qp_out->lb+ii, 0);
		VECCP_LIBSTR(nb[ii], qp_in->ub+ii, 0, qp_out->ub+ii, 0);
		VECCP_LIBSTR(ng[ii], qp_in->lg+ii, 0, qp_out->lg+ii, 0);
		VECCP_LIBSTR(ng[ii], qp_in->ug+ii, 0, qp_out->ug+ii, 0);
		}

	return;

	}


