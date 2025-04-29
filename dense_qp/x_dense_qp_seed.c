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



hpipm_size_t DENSE_QP_SEED_MEMSIZE(struct DENSE_QP_DIM *dim)
	{

	// loop index
	int ii;

	// extract ocp qp size
	int nv = dim->nv;
	int ne = dim->ne;
	int nb = dim->nb;
	int ng = dim->ng;
	int ns = dim->ns;

	hpipm_size_t size = 0;

	size += 4*sizeof(struct STRVEC); // seed_g seed_b seed_d seed_m

	size += 1*SIZE_STRVEC(nv+2*ns); // seed_g
	size += 1*SIZE_STRVEC(ne); // seed_b
	size += 2*SIZE_STRVEC(2*nb+2*ng+2*ns); // seed_d seed_m

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

	}



void DENSE_QP_SEED_CREATE(struct DENSE_QP_DIM *dim, struct DENSE_QP_SEED *seed, void *mem)
	{

	// loop index
	int ii;

	// zero memory (to avoid corrupted memory like e.g. NaN)
	hpipm_size_t memsize = DENSE_QP_SEED_MEMSIZE(dim);
	hpipm_zero_memset(memsize, mem);

	// extract ocp qp size
	int nv = dim->nv;
	int ne = dim->ne;
	int nb = dim->nb;
	int ng = dim->ng;
	int ns = dim->ns;


	// vector struct
	struct STRVEC *sv_ptr = (struct STRVEC *) mem;

	seed->seed_g = sv_ptr;
	sv_ptr += 1;
	seed->seed_b = sv_ptr;
	sv_ptr += 1;
	seed->seed_d = sv_ptr;
	sv_ptr += 1;
	seed->seed_m = sv_ptr;
	sv_ptr += 1;


	// align to typical cache line size
	hpipm_size_t s_ptr = (hpipm_size_t) sv_ptr;
	s_ptr = (s_ptr+63)/64*64;


	// void stuf
	char *c_ptr = (char *) s_ptr;

	CREATE_STRVEC(nv+2*ns, seed->seed_g, c_ptr);
	c_ptr += (seed->seed_g)->memsize;

	CREATE_STRVEC(ne, seed->seed_b, c_ptr);
	c_ptr += (seed->seed_b)->memsize;

	CREATE_STRVEC(2*nb+2*ng+2*ns, seed->seed_d, c_ptr);
	c_ptr += (seed->seed_d)->memsize;

	CREATE_STRVEC(2*nb+2*ng+2*ns, seed->seed_m, c_ptr);
	c_ptr += (seed->seed_m)->memsize;

	seed->dim = dim;

	seed->memsize = DENSE_QP_SEED_MEMSIZE(dim);


#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + seed->memsize)
		{
#ifdef EXT_DEP
		printf("\ncreate_dense_qp_seed: outsize memory bounds!\n\n");
#endif
		exit(1);
		}
#endif


	return;

	}



void DENSE_QP_SEED_SET(char *field, void *value, struct DENSE_QP_SEED *qp_seed)
	{
	if(hpipm_strcmp(field, "seed_g")) 
		{
		DENSE_QP_SEED_SET_SEED_G(value, qp_seed);
		}
	else if(hpipm_strcmp(field, "seed_zl"))
		{ 
		DENSE_QP_SEED_SET_SEED_ZL(value, qp_seed);
		}
	else if(hpipm_strcmp(field, "seed_zu"))
		{ 
		DENSE_QP_SEED_SET_SEED_ZU(value, qp_seed);
		}
	else if(hpipm_strcmp(field, "seed_b"))
		{
		DENSE_QP_SEED_SET_SEED_B(value, qp_seed);
		}
	else if(hpipm_strcmp(field, "seed_lb"))
		{
		DENSE_QP_SEED_SET_SEED_LB(value, qp_seed);
		}
	else if(hpipm_strcmp(field, "seed_ub"))
		{
		DENSE_QP_SEED_SET_SEED_UB(value, qp_seed);
		}
	else if(hpipm_strcmp(field, "seed_lg"))
		{
		DENSE_QP_SEED_SET_SEED_LG(value, qp_seed);
		}
	else if(hpipm_strcmp(field, "seed_ug"))
		{
		DENSE_QP_SEED_SET_SEED_UG(value, qp_seed);
		}
	else if(hpipm_strcmp(field, "seed_ls"))
		{ 
		DENSE_QP_SEED_SET_SEED_LS(value, qp_seed);
		}
	else if(hpipm_strcmp(field, "seed_us"))
		{ 
		DENSE_QP_SEED_SET_SEED_US(value, qp_seed);
		}
	else
		{
#ifdef EXT_DEP
		printf("error: DENSE_QP_SEED_SET: wrong field name '%s'. Exiting.\n", field);
#endif
		exit(1);	
		}
	return;
	}



void DENSE_QP_SEED_SET_SEED_G(REAL *g, struct DENSE_QP_SEED *qp_seed)
	{
	int nv = qp_seed->dim->nv;
	PACK_VEC(nv, g, 1, qp_seed->seed_g, 0);
	}



void DENSE_QP_SEED_SET_SEED_ZL(REAL *zl, struct DENSE_QP_SEED *qp_seed)
	{
	int nv = qp_seed->dim->nv;
	int ns = qp_seed->dim->ns;
	PACK_VEC(ns, zl, 1, qp_seed->seed_g, nv);
	}



void DENSE_QP_SEED_SET_SEED_ZU(REAL *zu, struct DENSE_QP_SEED *qp_seed)
	{
	int nv = qp_seed->dim->nv;
	int ns = qp_seed->dim->ns;
	PACK_VEC(ns, zu, 1, qp_seed->seed_g, nv+ns);
	}



void DENSE_QP_SEED_SET_SEED_B(REAL *b, struct DENSE_QP_SEED *qp_seed)
	{
	int ne = qp_seed->dim->ne;
	PACK_VEC(ne, b, 1, qp_seed->seed_b, 0);
	}



void DENSE_QP_SEED_SET_SEED_LB(REAL *lb, struct DENSE_QP_SEED *qp_seed)
	{
	int nb = qp_seed->dim->nb;
	PACK_VEC(nb, lb, 1, qp_seed->seed_d, 0);
	}



void DENSE_QP_SEED_SET_SEED_UB(REAL *ub, struct DENSE_QP_SEED *qp_seed)
	{
	int nb = qp_seed->dim->nb;
	int ng = qp_seed->dim->ng;
	PACK_VEC(nb, ub, 1, qp_seed->seed_d, nb+ng);
	VECSC(nb, -1.0, qp_seed->seed_d, nb+ng);
	}



void DENSE_QP_SEED_SET_SEED_LG(REAL *lg, struct DENSE_QP_SEED *qp_seed)
	{
	int nb = qp_seed->dim->nb;
	int ng = qp_seed->dim->ng;
	PACK_VEC(ng, lg, 1, qp_seed->seed_d, nb);
	}



void DENSE_QP_SEED_SET_SEED_UG(REAL *ug, struct DENSE_QP_SEED *qp_seed)
	{
	int nb = qp_seed->dim->nb;
	int ng = qp_seed->dim->ng;
	PACK_VEC(ng, ug, 1, qp_seed->seed_d, 2*nb+ng);
	VECSC(ng, -1.0, qp_seed->seed_d, 2*nb+ng);
	}



void DENSE_QP_SEED_SET_SEED_LS(REAL *ls, struct DENSE_QP_SEED *qp_seed)
	{
	int nb = qp_seed->dim->nb;
	int ng = qp_seed->dim->ng;
	int ns = qp_seed->dim->ns;
	PACK_VEC(ns, ls, 1, qp_seed->seed_d, 2*nb+2*ng);
	}



void DENSE_QP_SEED_SET_SEED_US(REAL *us, struct DENSE_QP_SEED *qp_seed)
	{
	int nb = qp_seed->dim->nb;
	int ng = qp_seed->dim->ng;
	int ns = qp_seed->dim->ns;
	PACK_VEC(ns, us, 1, qp_seed->seed_d, 2*nb+2*ng+ns);
	}
