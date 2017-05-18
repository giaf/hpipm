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



#if defined(RUNTIME_CHECKS)
#include <stdlib.h>
#endif

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_d_aux.h>

#include "../include/hpipm_ocp_kkt.h"



int d_size_ocp_qp(int N, int *nx, int *nu, int *nb, int *ng)
	{

	int ii;

	int size = 0;

	size += 4*(N+1)*sizeof(int); // nx nu nb ng
	size += (N+1)*sizeof(int *); // idxb
	size += (1*N+2*(N+1))*sizeof(struct d_strmat); // BAbt ; RSqrq DCt
	size += (1*N+5*(N+1))*sizeof(struct d_strvec); // b ; rq lb ub lg ug

	for(ii=0; ii<N; ii++)
		{
		size += nb[ii]*sizeof(int); // idxb
		size += d_size_strmat(nu[ii]+nx[ii]+1, nx[ii+1]); // BAbt
		size += d_size_strvec(nx[ii+1]); // b
		size += d_size_strmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii]); // RSQrq
		size += d_size_strvec(nu[ii]+nx[ii]); // rq
		size += d_size_strmat(nu[ii]+nx[ii], ng[ii]); // DCt
		size += 2*d_size_strvec(nb[ii]); // lb ub
		size += 2*d_size_strvec(ng[ii]); // lg ug
		}
	ii = N;
	size += nb[ii]*sizeof(int); // idxb
	size += d_size_strmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii]); // RSQrq
	size += d_size_strvec(nu[ii]+nx[ii]); // rq
	size += d_size_strmat(nu[ii]+nx[ii], ng[ii]); // DCt
	size += 2*d_size_strvec(nb[ii]); // lb ub
	size += 2*d_size_strvec(ng[ii]); // lg ug

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 64; // align to typical cache line size
	
	return size;

	}



void d_create_ocp_qp(int N, int *nx, int *nu, int *nb, int *ng, struct d_ocp_qp *str_out, void *memory)
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
	struct d_strmat *sm_ptr = (struct d_strmat *) ip_ptr;

	// BAbt
	str_out->sBAbt = sm_ptr;
	sm_ptr += N;

	// RSQrq
	str_out->sRSQrq = sm_ptr;
	sm_ptr += N+1;

	// DCt
	str_out->sDCt = sm_ptr;
	sm_ptr += N+1;


	// vector struct stuff
	struct d_strvec *sv_ptr = (struct d_strvec *) sm_ptr;

	// b
	str_out->sb = sv_ptr;
	sv_ptr += N;

	// rq
	str_out->srq = sv_ptr;
	sv_ptr += N+1;

	// lb
	str_out->slb = sv_ptr;
	sv_ptr += N+1;

	// ub
	str_out->sub = sv_ptr;
	sv_ptr += N+1;

	// lg
	str_out->slg = sv_ptr;
	sv_ptr += N+1;

	// ug
	str_out->sug = sv_ptr;
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
		d_create_strmat(nu[ii]+nx[ii]+1, nx[ii+1], str_out->sBAbt+ii, v_ptr);
		v_ptr += (str_out->sBAbt+ii)->memory_size;
		}

	// RSQrq
	for(ii=0; ii<=N; ii++)
		{
		d_create_strmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], str_out->sRSQrq+ii, v_ptr);
		v_ptr += (str_out->sRSQrq+ii)->memory_size;
		}

	// DCt
	for(ii=0; ii<=N; ii++)
		{
		d_create_strmat(nu[ii]+nx[ii], ng[ii], str_out->sDCt+ii, v_ptr);
		v_ptr += (str_out->sDCt+ii)->memory_size;
		}

	// b
	for(ii=0; ii<N; ii++)
		{
		d_create_strvec(nx[ii+1], str_out->sb+ii, v_ptr);
		v_ptr += (str_out->sb+ii)->memory_size;
		}

	// rq
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(nu[ii]+nx[ii], str_out->srq+ii, v_ptr);
		v_ptr += (str_out->srq+ii)->memory_size;
		}

	// lb
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(nb[ii], str_out->slb+ii, v_ptr);
		v_ptr += (str_out->slb+ii)->memory_size;
		}

	// ub
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(nb[ii], str_out->sub+ii, v_ptr);
		v_ptr += (str_out->sub+ii)->memory_size;
		}

	// lg
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(ng[ii], str_out->slg+ii, v_ptr);
		v_ptr += (str_out->slg+ii)->memory_size;
		}

	// ug
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(ng[ii], str_out->sug+ii, v_ptr);
		v_ptr += (str_out->sug+ii)->memory_size;
		}

	return;

	}



void d_cast_ocp_qp(int N, int *nx, int *nu, int *nb, int **idxb, int *ng, struct d_strmat *sBAbt, struct d_strvec *sb, struct d_strmat *sRSQrq, struct d_strvec *srq, struct d_strmat *sDCt, struct d_strvec *slb, struct d_strvec *sub, struct d_strvec *slg, struct d_strvec *sug, struct d_ocp_qp *str_out)
	{

	str_out->NN = N;
	str_out->nx = nx;
	str_out->nu = nu;
	str_out->nb = nb;
	str_out->idxb = idxb;
	str_out->ng = ng;
	str_out->sBAbt = sBAbt;
	str_out->sb = sb;
	str_out->sRSQrq = sRSQrq;
	str_out->srq = srq;
	str_out->sDCt = sDCt;
	str_out->slb = slb;
	str_out->sub = sub;
	str_out->slg = slg;
	str_out->sug = sug;

	return;

	}



void d_copy_ocp_qp(struct d_ocp_qp *str_in, struct d_ocp_qp *str_out)
	{

	int N = str_in->NN;
	int *nx = str_in->nx;
	int *nu = str_in->nu;
	int *nb = str_in->nb;
	int *ng = str_in->ng;

	int ii, jj;

#if defined(RUNTIME_CHECKS)
	if(str_out->NN != str_in->NN)
		{
		printf("\nError : d_copy_ocp_qp : str_out->NN != str_out->NN : %d != %d\n\n", str_out->NN, str_in->NN);
		exit(1);
		}
	for(ii=0; ii<=N; ii++)
		{
		if(str_out->nx[ii] != str_in->nx[ii])
			{
			printf("\nError : d_copy_ocp_qp : str_out->nx[%d] != str_out->nx[%d] : %d != %d\n\n", ii, ii, str_out->nx[ii], str_in->nx[ii]);
			exit(1);
			}
		if(str_out->nu[ii] != str_in->nu[ii])
			{
			printf("\nError : d_copy_ocp_qp : str_out->nu[%d] != str_out->nu[%d] : %d != %d\n\n", ii, ii, str_out->nu[ii], str_in->nu[ii]);
			exit(1);
			}
		if(str_out->nb[ii] != str_in->nb[ii])
			{
			printf("\nError : d_copy_ocp_qp : str_out->nb[%d] != str_out->nb[%d] : %d != %d\n\n", ii, ii, str_out->nb[ii], str_in->nb[ii]);
			exit(1);
			}
		if(str_out->ng[ii] != str_in->ng[ii])
			{
			printf("\nError : d_copy_ocp_qp : str_out->ng[%d] != str_out->ng[%d] : %d != %d\n\n", ii, ii, str_out->ng[ii], str_in->ng[ii]);
			exit(1);
			}
		}
#endif

	for(ii=0; ii<N; ii++)
		{
		for(jj=0; jj<str_in->nb[ii]; jj++) str_out->idxb[ii][jj] = str_in->idxb[ii][jj];
		dgecp_libstr(nx[ii]+nu[ii]+1, nx[ii+1], str_in->sBAbt+ii, 0, 0, str_out->sBAbt+ii, 0, 0);
		dveccp_libstr(nx[ii+1], str_in->sb+ii, 0, str_out->sb+ii, 0);
		dgecp_libstr(nx[ii]+nu[ii]+1, nu[ii]+nx[ii], str_in->sRSQrq+ii, 0, 0, str_out->sRSQrq+ii, 0, 0);
		dveccp_libstr(nu[ii]+nx[ii], str_in->srq+ii, 0, str_out->srq+ii, 0);
		dgecp_libstr(nx[ii]+nu[ii], ng[ii], str_in->sDCt+ii, 0, 0, str_out->sDCt+ii, 0, 0);
		dveccp_libstr(nb[ii], str_in->slb+ii, 0, str_out->slb+ii, 0);
		dveccp_libstr(nb[ii], str_in->sub+ii, 0, str_out->sub+ii, 0);
		dveccp_libstr(ng[ii], str_in->slg+ii, 0, str_out->slg+ii, 0);
		dveccp_libstr(ng[ii], str_in->sug+ii, 0, str_out->sug+ii, 0);
		}

	return;

	}

