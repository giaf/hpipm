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



hpipm_size_t OCP_QCQP_DIM_STRSIZE()
	{

	hpipm_size_t size = 0;

	size += sizeof(struct OCP_QCQP_DIM);

	return size;

	}



hpipm_size_t OCP_QCQP_DIM_MEMSIZE(int N)
	{

	hpipm_size_t size = 0;

	size += 16*(N+1)*sizeof(int);

	size += 1*sizeof(struct OCP_QP_DIM);
	size += 1*OCP_QP_DIM_MEMSIZE(N);

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

	}



void OCP_QCQP_DIM_CREATE(int N, struct OCP_QCQP_DIM *dim, void *mem)
	{

	// loop index
	int ii;

	// zero memory (to avoid corrupted memory like e.g. NaN)
	hpipm_size_t memsize = OCP_QCQP_DIM_MEMSIZE(N);
	hpipm_zero_memset(memsize, mem);

	// qp_dim struct
	struct OCP_QP_DIM *dim_ptr = mem;

	dim->qp_dim = dim_ptr;
	dim_ptr += 1;

	// align to typical cache line size
	hpipm_size_t s_ptr = (hpipm_size_t) dim_ptr;
	s_ptr = (s_ptr+63)/64*64;

	// void
	char *c_ptr = (char *) s_ptr;

	OCP_QP_DIM_CREATE(N, dim->qp_dim, c_ptr);
	c_ptr += dim->qp_dim->memsize;

	// nx
	dim->nx = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	// nu
	dim->nu = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	// nb
	dim->nb = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	// nbx
	dim->nbx = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	// nbu
	dim->nbu = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	// ng
	dim->ng = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	// nq
	dim->nq = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	// ns
	dim->ns = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	// nsbx
	dim->nsbx = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	// nsbu
	dim->nsbu = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	// nsg
	dim->nsg = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	// nsq
	dim->nsq = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	// nbxe
	dim->nbxe = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	// nbue
	dim->nbue = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	// nge
	dim->nge = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	// nqe
	dim->nqe = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);

	// N
	dim->N = N;

	// initialize dims to zero by default
	// XXX already zero-initialized while zeroing all memory

	dim->memsize = memsize; //OCP_QCQP_DIM_MEMSIZE(N);

	return;

	}


void OCP_QCQP_DIM_COPY_ALL(struct OCP_QCQP_DIM *dim_orig, struct OCP_QCQP_DIM *dim_dest)
	{

#if defined(RUNTIME_CHECKS)
	if(dim_orig->N!=dim_dest->N)
		{
#ifdef EXT_DEP
		printf("\nerror: OCP_QCQP_DIM_COPY_ALL: dim_orig->N != dim_dest->N\n");
#endif
		exit(1);
		}
#endif

	// loop index
	int ii;

	// N
	int N = dim_orig->N;

	// copy qp dim
	for(ii=0; ii<=N; ii++)
		OCP_QCQP_DIM_SET_NX(ii, dim_orig->nx[ii], dim_dest);
	for(ii=0; ii<=N; ii++)
		OCP_QCQP_DIM_SET_NU(ii, dim_orig->nu[ii], dim_dest);
	for(ii=0; ii<=N; ii++)
		OCP_QCQP_DIM_SET_NBX(ii, dim_orig->nbx[ii], dim_dest);
	for(ii=0; ii<=N; ii++)
		OCP_QCQP_DIM_SET_NBU(ii, dim_orig->nbu[ii], dim_dest);
	for(ii=0; ii<=N; ii++)
		OCP_QCQP_DIM_SET_NQ(ii, dim_orig->nq[ii], dim_dest);
	for(ii=0; ii<=N; ii++)
		OCP_QCQP_DIM_SET_NG(ii, dim_orig->ng[ii], dim_dest);
	for(ii=0; ii<=N; ii++)
		OCP_QCQP_DIM_SET_NSBX(ii, dim_orig->nsbx[ii], dim_dest);
	for(ii=0; ii<=N; ii++)
		OCP_QCQP_DIM_SET_NSBU(ii, dim_orig->nsbu[ii], dim_dest);
	for(ii=0; ii<=N; ii++)
		OCP_QCQP_DIM_SET_NSG(ii, dim_orig->nsg[ii], dim_dest);
	for(ii=0; ii<=N; ii++)
		OCP_QCQP_DIM_SET_NSQ(ii, dim_orig->nsq[ii], dim_dest);
	for(ii=0; ii<=N; ii++)
		OCP_QCQP_DIM_SET_NBXE(ii, dim_orig->nbxe[ii], dim_dest);
	for(ii=0; ii<=N; ii++)
		OCP_QCQP_DIM_SET_NBUE(ii, dim_orig->nbue[ii], dim_dest);
	for(ii=0; ii<=N; ii++)
		OCP_QCQP_DIM_SET_NGE(ii, dim_orig->nge[ii], dim_dest);
	for(ii=0; ii<=N; ii++)
		OCP_QCQP_DIM_SET_NQE(ii, dim_orig->nqe[ii], dim_dest);

	return;

	}



void OCP_QCQP_DIM_SET(char *field_name, int stage, int value, struct OCP_QCQP_DIM *dim)
	{
	if(hpipm_strcmp(field_name, "nx"))
		{ 
		OCP_QCQP_DIM_SET_NX(stage, value, dim);
		}
	else if(hpipm_strcmp(field_name, "nu"))
		{ 
		OCP_QCQP_DIM_SET_NU(stage, value, dim);
		}
	else if(hpipm_strcmp(field_name, "nbx"))
		{
		OCP_QCQP_DIM_SET_NBX(stage, value, dim);
		}
	else if(hpipm_strcmp(field_name, "nbu"))
		{
		OCP_QCQP_DIM_SET_NBU(stage, value, dim);
		}
	else if(hpipm_strcmp(field_name, "ng"))
		{
		OCP_QCQP_DIM_SET_NG(stage, value, dim);
		}
	else if(hpipm_strcmp(field_name, "nq"))
		{
		OCP_QCQP_DIM_SET_NQ(stage, value, dim);
		}
	else if(hpipm_strcmp(field_name, "ns"))
		{
		OCP_QCQP_DIM_SET_NS(stage, value, dim);
		}
	else if(hpipm_strcmp(field_name, "nsbx"))
		{
		OCP_QCQP_DIM_SET_NSBX(stage, value, dim);
		}
	else if(hpipm_strcmp(field_name, "nsbu"))
		{
		OCP_QCQP_DIM_SET_NSBU(stage, value, dim);
		}
	else if(hpipm_strcmp(field_name, "nsg"))
		{
		OCP_QCQP_DIM_SET_NSG(stage, value, dim);
		}
	else if(hpipm_strcmp(field_name, "nsq"))
		{
		OCP_QCQP_DIM_SET_NSQ(stage, value, dim);
		}
	else if(hpipm_strcmp(field_name, "nbxe"))
		{
		OCP_QCQP_DIM_SET_NBXE(stage, value, dim);
		}
	else if(hpipm_strcmp(field_name, "nbue"))
		{
		OCP_QCQP_DIM_SET_NBUE(stage, value, dim);
		}
	else if(hpipm_strcmp(field_name, "nge"))
		{
		OCP_QCQP_DIM_SET_NGE(stage, value, dim);
		}
	else if(hpipm_strcmp(field_name, "nqe"))
		{
		OCP_QCQP_DIM_SET_NQE(stage, value, dim);
		}
	else 
		{
#ifdef EXT_DEP
		printf("error: OCP_QCQP_DIM_SET: wrong field %s\n", field_name);
#endif
		exit(1);
		}
	return;
	}



void OCP_QCQP_DIM_SET_NX(int stage, int value, struct OCP_QCQP_DIM *dim)
	{
	dim->nx[stage] = value;
	OCP_QP_DIM_SET_NX(stage, dim->nx[stage], dim->qp_dim);
	return;
	}



void OCP_QCQP_DIM_SET_NU(int stage, int value, struct OCP_QCQP_DIM *dim)
	{
	dim->nu[stage] = value;
	OCP_QP_DIM_SET_NU(stage, dim->nu[stage], dim->qp_dim);
	return;
	}



void OCP_QCQP_DIM_SET_NBX(int stage, int value, struct OCP_QCQP_DIM *dim)
	{
	dim->nbx[stage] = value;
	dim->nb[stage] = dim->nbx[stage] + dim->nbu[stage];
	OCP_QP_DIM_SET_NBX(stage, dim->nbx[stage], dim->qp_dim);
	return;
	}



void OCP_QCQP_DIM_SET_NBU(int stage, int value, struct OCP_QCQP_DIM *dim)
	{
	dim->nbu[stage] = value;
	dim->nb[stage] = dim->nbx[stage] + dim->nbu[stage];
	OCP_QP_DIM_SET_NBU(stage, dim->nbu[stage], dim->qp_dim);
	return;
	}



void OCP_QCQP_DIM_SET_NG(int stage, int value, struct OCP_QCQP_DIM *dim)
	{
	dim->ng[stage] = value;
	OCP_QP_DIM_SET_NG(stage, dim->ng[stage]+dim->nq[stage], dim->qp_dim);
	return;
	}



void OCP_QCQP_DIM_SET_NQ(int stage, int value, struct OCP_QCQP_DIM *dim)
	{
	dim->nq[stage] = value;
	OCP_QP_DIM_SET_NG(stage, dim->ng[stage]+dim->nq[stage], dim->qp_dim);
	return;
	}



void OCP_QCQP_DIM_SET_NS(int stage, int value, struct OCP_QCQP_DIM *dim)
	{
	dim->ns[stage] = value;
	OCP_QP_DIM_SET_NS(stage, dim->ns[stage], dim->qp_dim);
	return;
	}



void OCP_QCQP_DIM_SET_NSBX(int stage, int value, struct OCP_QCQP_DIM *dim)
	{
	dim->nsbx[stage] = value;
	dim->ns[stage] = dim->nsbx[stage] + dim->nsbu[stage] + dim->nsg[stage] + dim->nsq[stage];
	OCP_QP_DIM_SET_NSBX(stage, dim->nsbx[stage], dim->qp_dim);
	return;
	}



void OCP_QCQP_DIM_SET_NSBU(int stage, int value, struct OCP_QCQP_DIM *dim)
	{
	dim->nsbu[stage] = value;
	dim->ns[stage] = dim->nsbx[stage] + dim->nsbu[stage] + dim->nsg[stage] + dim->nsq[stage];
	OCP_QP_DIM_SET_NSBU(stage, dim->nsbu[stage], dim->qp_dim);
	return;
	}



void OCP_QCQP_DIM_SET_NSG(int stage, int value, struct OCP_QCQP_DIM *dim)
	{
	dim->nsg[stage] = value;
	dim->ns[stage] = dim->nsbx[stage] + dim->nsbu[stage] + dim->nsg[stage] + dim->nsq[stage];
	OCP_QP_DIM_SET_NSG(stage, dim->nsg[stage]+dim->nsq[stage], dim->qp_dim);
	return;
	}



void OCP_QCQP_DIM_SET_NSQ(int stage, int value, struct OCP_QCQP_DIM *dim)
	{
	dim->nsq[stage] = value;
	dim->ns[stage] = dim->nsbx[stage] + dim->nsbu[stage] + dim->nsg[stage] + dim->nsq[stage];
	OCP_QP_DIM_SET_NSG(stage, dim->nsg[stage]+dim->nsq[stage], dim->qp_dim);
	return;
	}



void OCP_QCQP_DIM_SET_NBXE(int stage, int value, struct OCP_QCQP_DIM *dim)
	{
	dim->nbxe[stage] = value;
	OCP_QP_DIM_SET_NBXE(stage, dim->nbxe[stage], dim->qp_dim);
	return;
	}



void OCP_QCQP_DIM_SET_NBUE(int stage, int value, struct OCP_QCQP_DIM *dim)
	{
	dim->nbue[stage] = value;
	OCP_QP_DIM_SET_NBUE(stage, dim->nbue[stage], dim->qp_dim);
	return;
	}



void OCP_QCQP_DIM_SET_NGE(int stage, int value, struct OCP_QCQP_DIM *dim)
	{
	dim->nge[stage] = value;
	OCP_QP_DIM_SET_NGE(stage, dim->nge[stage]+dim->nqe[stage], dim->qp_dim);
	return;
	}



void OCP_QCQP_DIM_SET_NQE(int stage, int value, struct OCP_QCQP_DIM *dim)
	{
	dim->nqe[stage] = value;
	OCP_QP_DIM_SET_NGE(stage, dim->nge[stage]+dim->nqe[stage], dim->qp_dim);
	return;
	}



void OCP_QCQP_DIM_GET(struct OCP_QCQP_DIM *dim, char *field_name, int stage, int *value)
	{
	if(hpipm_strcmp(field_name, "nx"))
		{
		OCP_QCQP_DIM_GET_NX(dim, stage, value);
		}
	else if(hpipm_strcmp(field_name, "nu"))
		{
		OCP_QCQP_DIM_GET_NU(dim, stage, value);
		}
	else
		{
#ifdef EXT_DEP
		printf("error: OCP_QCQP_DIM_GET: wrong field %s\n", field_name);
#endif
		exit(1);
		}
	return;
	}



void OCP_QCQP_DIM_GET_N(struct OCP_QCQP_DIM *dim, int *value)
	{
	*value = dim->N;
	return;
	}



void OCP_QCQP_DIM_GET_NX(struct OCP_QCQP_DIM *dim, int stage, int *value)
	{
	*value = dim->nx[stage];
	return;
	}



void OCP_QCQP_DIM_GET_NU(struct OCP_QCQP_DIM *dim, int stage, int *value)
	{
	*value = dim->nu[stage];
	return;
	}



