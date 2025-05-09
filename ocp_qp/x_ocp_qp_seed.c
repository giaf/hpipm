/**************************************************************************************************
*                                                                                                 *
* This file is part of HPIPM.                                                                     *
*                                                                                                 *
* HPIPM -- High-Performance Interior Point Method.                                                *
* Copyright (C) 2019 by Gianluca Frison.                                                          *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights seederved.                                                                            *
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



hpipm_size_t OCP_QP_SEED_STRSIZE()
	{
	return sizeof(struct OCP_QP_SEED);
	}



hpipm_size_t OCP_QP_SEED_MEMSIZE(struct OCP_QP_DIM *dim)
	{

	// loop index
	int ii;

	// extract ocp qp size
	int N = dim->N;
	int *nx = dim->nx;
	int *nu = dim->nu;
	int *nb = dim->nb;
	int *ng = dim->ng;
	int *ns = dim->ns;

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

	hpipm_size_t size = 0;

	size += 3*(N+1)*sizeof(struct STRVEC); // seed_g seed_d seed_m
	size += 1*N*sizeof(struct STRVEC); // seed_b

	size += 1*SIZE_STRVEC(nvt); // seed_g
	size += 1*SIZE_STRVEC(net); // seed_b
	size += 2*SIZE_STRVEC(nct); // seed_d seed_m

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

	}



void OCP_QP_SEED_CREATE(struct OCP_QP_DIM *dim, struct OCP_QP_SEED *seed, void *mem)
	{

	// loop index
	int ii;

	// zero memory (to avoid corrupted memory like e.g. NaN)
	hpipm_size_t memsize = OCP_QP_SEED_MEMSIZE(dim);
	hpipm_zero_memset(memsize, mem);

	// extract ocp qp size
	int N = dim->N;
	int *nx = dim->nx;
	int *nu = dim->nu;
	int *nb = dim->nb;
	int *ng = dim->ng;
	int *ns = dim->ns;

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


	// vector struct
	struct STRVEC *sv_ptr = (struct STRVEC *) mem;

	seed->seed_g = sv_ptr;
	sv_ptr += N+1;
	seed->seed_b = sv_ptr;
	sv_ptr += N;
	seed->seed_d = sv_ptr;
	sv_ptr += N+1;
	seed->seed_m = sv_ptr;
	sv_ptr += N+1;


	// align to typical cache line size
	hpipm_size_t s_ptr = (hpipm_size_t) sv_ptr;
	s_ptr = (s_ptr+63)/64*64;


	// void stuf
	char *c_ptr = (char *) s_ptr;

	CREATE_STRVEC(nvt, seed->seed_g, c_ptr);
	c_ptr += SIZE_STRVEC(nvt);

	CREATE_STRVEC(net, seed->seed_b, c_ptr);
	c_ptr += SIZE_STRVEC(net);

	CREATE_STRVEC(nct, seed->seed_d, c_ptr);
	c_ptr += SIZE_STRVEC(nct);

	CREATE_STRVEC(nct, seed->seed_m, c_ptr);
	c_ptr += SIZE_STRVEC(nct);

	// alias
	//
	c_ptr = (char *) seed->seed_g->pa;
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(nu[ii]+nx[ii]+2*ns[ii], seed->seed_g+ii, c_ptr);
		c_ptr += nu[ii]*sizeof(REAL);
		c_ptr += nx[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		}
	//
	c_ptr = (char *) seed->seed_b->pa;
	for(ii=0; ii<N; ii++)
		{
		CREATE_STRVEC(nx[ii+1], seed->seed_b+ii, c_ptr);
		c_ptr += (nx[ii+1])*sizeof(REAL);
		}
	//
	c_ptr = (char *) seed->seed_d->pa;
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(2*nb[ii]+2*ng[ii]+2*ns[ii], seed->seed_d+ii, c_ptr);
		c_ptr += nb[ii]*sizeof(REAL);
		c_ptr += ng[ii]*sizeof(REAL);
		c_ptr += nb[ii]*sizeof(REAL);
		c_ptr += ng[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		}
	//
	c_ptr = (char *) seed->seed_m->pa;
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(2*nb[ii]+2*ng[ii]+2*ns[ii], seed->seed_m+ii, c_ptr);
		c_ptr += nb[ii]*sizeof(REAL);
		c_ptr += ng[ii]*sizeof(REAL);
		c_ptr += nb[ii]*sizeof(REAL);
		c_ptr += ng[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		}



	seed->dim = dim;

	seed->memsize = memsize; //OCP_QP_SEED_MEMSIZE(dim);


#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + seed->memsize)
		{
#ifdef EXT_DEP
		printf("\ncreate_ocp_qp_seed: outside memory bounds!\n\n");
#endif
		exit(1);
		}
#endif


	return;

	}



void OCP_QP_SEED_SET_ZERO(struct OCP_QP_SEED *seed)
	{

	int N = seed->dim->N;
	int *nx = seed->dim->nx;
	int *nu = seed->dim->nu;
	int *nb = seed->dim->nb;
	int *ng = seed->dim->ng;
	int *ns = seed->dim->ns;

	int ii;

	// seed_g
	for(ii=0; ii<=N; ii++)
		{
		VECSE(nu[ii]+nx[ii]+2*ns[ii], 0.0, seed->seed_g+ii, 0);
		}
	// seed_b
	for(ii=0; ii<N; ii++)
		{
		VECSE(nx[ii+1], 0.0, seed->seed_b+ii, 0);
		}
	// seed_d
	for(ii=0; ii<=N; ii++)
		{
		VECSE(2*nb[ii]+2*ng[ii]+2*ns[ii], 0.0, seed->seed_d+ii, 0);
		}
	// seed_m
	for(ii=0; ii<=N; ii++)
		{
		VECSE(2*nb[ii]+2*ng[ii]+2*ns[ii], 0.0, seed->seed_m+ii, 0);
		}

	return;

	}



void OCP_QP_SEED_SET(char *field, int stage, REAL *vec, struct OCP_QP_SEED *seed)
	{
	if(hpipm_strcmp(field, "seed_r"))
		{
		OCP_QP_SEED_SET_SEED_R(stage, vec, seed);
		}
	else if(hpipm_strcmp(field, "seed_q"))
		{
		OCP_QP_SEED_SET_SEED_Q(stage, vec, seed);
		}
	else if(hpipm_strcmp(field, "seed_zl"))
		{
		OCP_QP_SEED_SET_SEED_ZL(stage, vec, seed);
		}
	else if(hpipm_strcmp(field, "seed_zu"))
		{
		OCP_QP_SEED_SET_SEED_ZU(stage, vec, seed);
		}
	else if(hpipm_strcmp(field, "seed_b"))
		{
		OCP_QP_SEED_SET_SEED_B(stage, vec, seed);
		}
	else if(hpipm_strcmp(field, "seed_lb"))
		{
		OCP_QP_SEED_SET_SEED_LB(stage, vec, seed);
		}
	else if(hpipm_strcmp(field, "seed_lbu"))
		{
		OCP_QP_SEED_SET_SEED_LBU(stage, vec, seed);
		}
	else if(hpipm_strcmp(field, "seed_lbx"))
		{
		OCP_QP_SEED_SET_SEED_LBX(stage, vec, seed);
		}
	else if(hpipm_strcmp(field, "seed_ub"))
		{
		OCP_QP_SEED_SET_SEED_UB(stage, vec, seed);
		}
	else if(hpipm_strcmp(field, "seed_ubu"))
		{
		OCP_QP_SEED_SET_SEED_UBU(stage, vec, seed);
		}
	else if(hpipm_strcmp(field, "seed_ubx"))
		{
		OCP_QP_SEED_SET_SEED_UBX(stage, vec, seed);
		}
	else if(hpipm_strcmp(field, "seed_lg"))
		{
		OCP_QP_SEED_SET_SEED_LG(stage, vec, seed);
		}
	else if(hpipm_strcmp(field, "seed_ug"))
		{
		OCP_QP_SEED_SET_SEED_UG(stage, vec, seed);
		}
	else if(hpipm_strcmp(field, "seed_ls"))
		{
		OCP_QP_SEED_SET_SEED_LS(stage, vec, seed);
		}
	else if(hpipm_strcmp(field, "seed_us"))
		{
		OCP_QP_SEED_SET_SEED_US(stage, vec, seed);
		}
	else
		{
#ifdef EXT_DEP
		printf("error [OCP_QP_SEED_SET]: unknown field name '%s'. Exiting.\n", field);
#endif
		exit(1);
		}
	return;
	}



void OCP_QP_SEED_SET_SEED_R(int stage, REAL *vec, struct OCP_QP_SEED *seed)
	{
	int *nu = seed->dim->nu;
	PACK_VEC(nu[stage], vec, 1, seed->seed_g+stage, 0);
	return;
	}



void OCP_QP_SEED_SET_SEED_Q(int stage, REAL *vec, struct OCP_QP_SEED *seed)
	{
	int *nu = seed->dim->nu;
	int *nx = seed->dim->nx;
	PACK_VEC(nx[stage], vec, 1, seed->seed_g+stage, nu[stage]);
	return;
	}



void OCP_QP_SEED_SET_SEED_ZL(int stage, REAL *vec, struct OCP_QP_SEED *seed)
	{
	int *nu = seed->dim->nu;
	int *nx = seed->dim->nx;
	int *ns = seed->dim->ns;
	PACK_VEC(ns[stage], vec, 1, seed->seed_g+stage, nu[stage]+nx[stage]);
	return;
	}



void OCP_QP_SEED_SET_SEED_ZU(int stage, REAL *vec, struct OCP_QP_SEED *seed)
	{
	int *nu = seed->dim->nu;
	int *nx = seed->dim->nx;
	int *ns = seed->dim->ns;
	PACK_VEC(ns[stage], vec, 1, seed->seed_g+stage, nu[stage]+nx[stage]+ns[stage]);
	return;
	}



void OCP_QP_SEED_SET_SEED_B(int stage, REAL *vec, struct OCP_QP_SEED *seed)
	{
	int *nx = seed->dim->nx;
	PACK_VEC(nx[stage+1], vec, 1, seed->seed_b+stage, 0);
	}



void OCP_QP_SEED_SET_SEED_LB(int stage, REAL *vec, struct OCP_QP_SEED *seed)
	{
	int *nb = seed->dim->nb;
	PACK_VEC(nb[stage], vec, 1, seed->seed_d+stage, 0);
	}



void OCP_QP_SEED_SET_SEED_LBU(int stage, REAL *vec, struct OCP_QP_SEED *seed)
	{
	int *nbu = seed->dim->nbu;
	PACK_VEC(nbu[stage], vec, 1, seed->seed_d+stage, 0);
	}



void OCP_QP_SEED_SET_SEED_LBX(int stage, REAL *vec, struct OCP_QP_SEED *seed)
	{
	int *nbu = seed->dim->nbu;
	int *nbx = seed->dim->nbx;
	PACK_VEC(nbx[stage], vec, 1, seed->seed_d+stage, nbu[stage]);
	}



void OCP_QP_SEED_SET_SEED_UB(int stage, REAL *vec, struct OCP_QP_SEED *seed)
	{
	int *nb = seed->dim->nb;
	int *ng = seed->dim->ng;
	PACK_VEC(nb[stage], vec, 1, seed->seed_d+stage, nb[stage]+ng[stage]);
	VECSC(nb[stage], -1.0, seed->seed_d+stage, nb[stage]+ng[stage]);
	}



void OCP_QP_SEED_SET_SEED_UBU(int stage, REAL *vec, struct OCP_QP_SEED *seed)
	{
	int *nb = seed->dim->nb;
	int *nbu = seed->dim->nbu;
	int *ng = seed->dim->ng;
	PACK_VEC(nbu[stage], vec, 1, seed->seed_d+stage, nb[stage]+ng[stage]);
	VECSC(nbu[stage], -1.0, seed->seed_d+stage, nb[stage]+ng[stage]);
	}



void OCP_QP_SEED_SET_SEED_UBX(int stage, REAL *vec, struct OCP_QP_SEED *seed)
	{
	int *nb = seed->dim->nb;
	int *nbu = seed->dim->nbu;
	int *nbx = seed->dim->nbx;
	int *ng = seed->dim->ng;
	PACK_VEC(nbx[stage], vec, 1, seed->seed_d+stage, nb[stage]+ng[stage]+nbu[stage]);
	VECSC(nbx[stage], -1.0, seed->seed_d+stage, nb[stage]+ng[stage]+nbu[stage]);
	}



void OCP_QP_SEED_SET_SEED_LG(int stage, REAL *vec, struct OCP_QP_SEED *seed)
	{
	int *nb = seed->dim->nb;
	int *ng = seed->dim->ng;
	PACK_VEC(ng[stage], vec, 1, seed->seed_d+stage, nb[stage]);
	}



void OCP_QP_SEED_SET_SEED_UG(int stage, REAL *vec, struct OCP_QP_SEED *seed)
	{
	int *nb = seed->dim->nb;
	int *ng = seed->dim->ng;
	PACK_VEC(ng[stage], vec, 1, seed->seed_d+stage, 2*nb[stage]+ng[stage]);
	VECSC(ng[stage], -1.0, seed->seed_d+stage, 2*nb[stage]+ng[stage]);
	}



void OCP_QP_SEED_SET_SEED_LS(int stage, REAL *vec, struct OCP_QP_SEED *seed)
	{
	int *nb = seed->dim->nb;
	int *ng = seed->dim->ng;
	int *ns = seed->dim->ns;
	PACK_VEC(ns[stage], vec, 1, seed->seed_d+stage, 2*nb[stage]+2*ng[stage]);
	}



void OCP_QP_SEED_SET_SEED_US(int stage, REAL *vec, struct OCP_QP_SEED *seed)
	{
	int *nb = seed->dim->nb;
	int *ng = seed->dim->ng;
	int *ns = seed->dim->ns;
	PACK_VEC(ns[stage], vec, 1, seed->seed_d+stage, 2*nb[stage]+2*ng[stage]+ns[stage]);
	}




