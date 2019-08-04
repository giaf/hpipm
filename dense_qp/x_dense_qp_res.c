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



int MEMSIZE_DENSE_QP_RES(struct DENSE_QP_DIM *dim)
	{

	// loop index
	int ii;

	// extract ocp qp size
	int nv = dim->nv;
	int ne = dim->ne;
	int nb = dim->nb;
	int ng = dim->ng;
	int ns = dim->ns;

	int size = 0;

	size += 4*sizeof(struct STRVEC); // res_g res_b res_d res_m

	size += 1*SIZE_STRVEC(nv+2*ns); // res_g
	size += 1*SIZE_STRVEC(ne); // res_b
	size += 2*SIZE_STRVEC(2*nb+2*ng+2*ns); // res_d res_m

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

	}



void CREATE_DENSE_QP_RES(struct DENSE_QP_DIM *dim, struct DENSE_QP_RES *res, void *mem)
	{

	// loop index
	int ii;

	// extract ocp qp size
	int nv = dim->nv;
	int ne = dim->ne;
	int nb = dim->nb;
	int ng = dim->ng;
	int ns = dim->ns;


	// vector struct
	struct STRVEC *sv_ptr = (struct STRVEC *) mem;

	res->res_g = sv_ptr;
	sv_ptr += 1;
	res->res_b = sv_ptr;
	sv_ptr += 1;
	res->res_d = sv_ptr;
	sv_ptr += 1;
	res->res_m = sv_ptr;
	sv_ptr += 1;


	// align to typicl cache line size
	size_t s_ptr = (size_t) sv_ptr;
	s_ptr = (s_ptr+63)/64*64;


	// void stuf
	char *c_ptr = (char *) s_ptr;

	CREATE_STRVEC(nv+2*ns, res->res_g, c_ptr);
	c_ptr += (res->res_g)->memsize;

	CREATE_STRVEC(ne, res->res_b, c_ptr);
	c_ptr += (res->res_b)->memsize;

	CREATE_STRVEC(2*nb+2*ng+2*ns, res->res_d, c_ptr);
	c_ptr += (res->res_d)->memsize;

	CREATE_STRVEC(2*nb+2*ng+2*ns, res->res_m, c_ptr);
	c_ptr += (res->res_m)->memsize;

	res->dim = dim;

	res->memsize = MEMSIZE_DENSE_QP_RES(dim);


#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + res->memsize)
		{
		printf("\ncreate_dense_qp_res: outsize memory bounds!\n\n");
		exit(1);
		}
#endif


	return;

	}



int MEMSIZE_DENSE_QP_RES_WORKSPACE(struct DENSE_QP_DIM *dim)
	{

	// loop index
	int ii;

	// extract ocp qp size
	int nv = dim->nv;
	int ne = dim->ne;
	int nb = dim->nb;
	int ng = dim->ng;
	int ns = dim->ns;

	int size = 0;

	size += 3*sizeof(struct STRVEC); // 2*tmp_nbg tmp_ns

	size += 2*SIZE_STRVEC(nb+ng); // tmp_nbg
	size += 1*SIZE_STRVEC(ns); // tmp_ns

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

	}



void CREATE_DENSE_QP_RES_WORKSPACE(struct DENSE_QP_DIM *dim, struct DENSE_QP_RES_WORKSPACE *ws, void *mem)
	{

	// loop index
	int ii;

	// extract ocp qp size
	int nv = dim->nv;
	int ne = dim->ne;
	int nb = dim->nb;
	int ng = dim->ng;
	int ns = dim->ns;


	// vector struct
	struct STRVEC *sv_ptr = (struct STRVEC *) mem;

	ws->tmp_nbg = sv_ptr;
	sv_ptr += 2;
	ws->tmp_ns = sv_ptr;
	sv_ptr += 1;


	// align to typicl cache line size
	size_t s_ptr = (size_t) sv_ptr;
	s_ptr = (s_ptr+63)/64*64;


	// void stuf
	char *c_ptr = (char *) s_ptr;


	CREATE_STRVEC(nb+ng, ws->tmp_nbg+0, c_ptr);
	c_ptr += (ws->tmp_nbg+0)->memsize;

	CREATE_STRVEC(nb+ng, ws->tmp_nbg+1, c_ptr);
	c_ptr += (ws->tmp_nbg+1)->memsize;

	CREATE_STRVEC(ns, ws->tmp_ns+0, c_ptr);
	c_ptr += (ws->tmp_ns+0)->memsize;

	ws->memsize = MEMSIZE_DENSE_QP_RES(dim);


#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + ws->memsize)
		{
		printf("\ncreate_dense_qp_res_workspace: outsize memory bounds!\n\n");
		exit(1);
		}
#endif


	return;

	}



void CVT_DENSE_QP_RES_TO_COLMAJ(struct DENSE_QP_RES *res, REAL *res_g, REAL *res_ls, REAL *res_us, REAL *res_b, REAL *res_d_lb, REAL *res_d_ub, REAL *res_d_lg, REAL *res_d_ug, REAL *res_d_ls, REAL *res_d_us, REAL *res_m_lb, REAL *res_m_ub, REAL *res_m_lg, REAL *res_m_ug, REAL *res_m_ls, REAL *res_m_us)
	{

	int nv = res->dim->nv;
	int ne = res->dim->ne;
	int nb = res->dim->nb;
	int ng = res->dim->ng;
	int ns = res->dim->ns;

	CVT_STRVEC2VEC(nv, res->res_g, 0, res_g);
	CVT_STRVEC2VEC(ns, res->res_g, nv, res_ls);
	CVT_STRVEC2VEC(ns, res->res_g, nv+ns, res_us);

	CVT_STRVEC2VEC(ne, res->res_b, 0, res_b);
	CVT_STRVEC2VEC(nb, res->res_d, 0, res_d_lb);
	CVT_STRVEC2VEC(ng, res->res_d, nb, res_d_lg);
	CVT_STRVEC2VEC(nb, res->res_d, nb+ng, res_d_ub);
	CVT_STRVEC2VEC(ng, res->res_d, 2*nb+ng, res_d_ug);
	CVT_STRVEC2VEC(ns, res->res_d, 2*nb+2*ng, res_d_ls);
	CVT_STRVEC2VEC(ns, res->res_d, 2*nb+2*ng+ns, res_d_us);
	CVT_STRVEC2VEC(nb, res->res_m, 0, res_m_lb);
	CVT_STRVEC2VEC(ng, res->res_m, nb, res_m_lg);
	CVT_STRVEC2VEC(nb, res->res_m, nb+ng, res_m_ub);
	CVT_STRVEC2VEC(ng, res->res_m, 2*nb+ng, res_m_ug);
	CVT_STRVEC2VEC(ns, res->res_m, 2*nb+2*ng, res_m_ls);
	CVT_STRVEC2VEC(ns, res->res_m, 2*nb+2*ng+ns, res_m_us);

	return;

	}



void CVT_DENSE_QP_RES_TO_ROWMAJ(struct DENSE_QP_RES *res, REAL *res_g, REAL *res_ls, REAL *res_us, REAL *res_b, REAL *res_d_lb, REAL *res_d_ub, REAL *res_d_lg, REAL *res_d_ug, REAL *res_d_ls, REAL *res_d_us, REAL *res_m_lb, REAL *res_m_ub, REAL *res_m_lg, REAL *res_m_ug, REAL *res_m_ls, REAL *res_m_us)
	{

	int nv = res->dim->nv;
	int ne = res->dim->ne;
	int nb = res->dim->nb;
	int ng = res->dim->ng;
	int ns = res->dim->ns;

	CVT_STRVEC2VEC(nv, res->res_g, 0, res_g);
	CVT_STRVEC2VEC(ns, res->res_g, nv, res_ls);
	CVT_STRVEC2VEC(ns, res->res_g, nv+ns, res_us);

	CVT_STRVEC2VEC(ne, res->res_b, 0, res_b);
	CVT_STRVEC2VEC(nb, res->res_d, 0, res_d_lb);
	CVT_STRVEC2VEC(ng, res->res_d, nb, res_d_lg);
	CVT_STRVEC2VEC(nb, res->res_d, nb+ng, res_d_ub);
	CVT_STRVEC2VEC(ng, res->res_d, 2*nb+ng, res_d_ug);
	CVT_STRVEC2VEC(ns, res->res_d, 2*nb+2*ng, res_d_ls);
	CVT_STRVEC2VEC(ns, res->res_d, 2*nb+2*ng+ns, res_d_us);
	CVT_STRVEC2VEC(nb, res->res_m, 0, res_m_lb);
	CVT_STRVEC2VEC(ng, res->res_m, nb, res_m_lg);
	CVT_STRVEC2VEC(nb, res->res_m, nb+ng, res_m_ub);
	CVT_STRVEC2VEC(ng, res->res_m, 2*nb+ng, res_m_ug);
	CVT_STRVEC2VEC(ns, res->res_m, 2*nb+2*ng, res_m_ls);
	CVT_STRVEC2VEC(ns, res->res_m, 2*nb+2*ng+ns, res_m_us);

	return;

	}


