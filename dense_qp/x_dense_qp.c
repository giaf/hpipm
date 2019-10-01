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



int DENSE_QP_MEMSIZE(struct DENSE_QP_DIM *dim)
	{

	int nv = dim->nv;
	int ne = dim->ne;
	int nb = dim->nb;
	int ng = dim->ng;
	int ns = dim->ns;

	int size = 0;

	size += 5*sizeof(struct STRVEC); // gz b d m Z
	size += 3*sizeof(struct STRMAT); // Hv A Ct

	size += 1*SIZE_STRVEC(nv+2*ns); // g
	size += 1*SIZE_STRVEC(ne); // b
	size += 2*SIZE_STRVEC(2*nb+2*ng+2*ns); // d m
	size += 1*SIZE_STRVEC(2*ns); // Z
	size += 1*nb*sizeof(int); // idxb
	size += 1*ns*sizeof(int); // idxb

	size += 1*SIZE_STRMAT(nv+1, nv); // Hv
	size += 1*SIZE_STRMAT(ne, nv); // A
	size += 1*SIZE_STRMAT(nv, ng); // Ct

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

	}



void DENSE_QP_CREATE(struct DENSE_QP_DIM *dim, struct DENSE_QP *qp, void *mem)
	{

	// TODO set memory to zero !!!!!!!

	int ii;

	int nv = dim->nv;
	int ne = dim->ne;
	int nb = dim->nb;
	int ng = dim->ng;
	int ns = dim->ns;


	// matrix struct stuff
	struct STRMAT *sm_ptr = (struct STRMAT *) mem;

	qp->Hv = sm_ptr;
	sm_ptr += 1;

	qp->A = sm_ptr;
	sm_ptr += 1;

	qp->Ct = sm_ptr;
	sm_ptr += 1;


	// vector struct stuff
	struct STRVEC *sv_ptr = (struct STRVEC *) sm_ptr;

	qp->gz = sv_ptr;
	sv_ptr += 1;

	qp->b = sv_ptr;
	sv_ptr += 1;

	qp->d = sv_ptr;
	sv_ptr += 1;

	qp->m = sv_ptr;
	sv_ptr += 1;

	qp->Z = sv_ptr;
	sv_ptr += 1;


	// int stuff
	int *i_ptr;
	i_ptr = (int *) sv_ptr;

	// idxb
	qp->idxb = i_ptr;
	i_ptr += nb;

	// idxs
	qp->idxs = i_ptr;
	i_ptr += ns;


	// align to typical cache line size
	size_t s_ptr = (size_t) i_ptr;
	s_ptr = (s_ptr+63)/64*64;


	//  stuff
	char *c_ptr;
	c_ptr = (char *) s_ptr;

	CREATE_STRMAT(nv+1, nv, qp->Hv, c_ptr);
	c_ptr += qp->Hv->memsize;

	CREATE_STRMAT(ne, nv, qp->A, c_ptr);
	c_ptr += qp->A->memsize;

	CREATE_STRMAT(nv, ng, qp->Ct, c_ptr);
	c_ptr += qp->Ct->memsize;

	CREATE_STRVEC(nv+2*ns, qp->gz, c_ptr);
	c_ptr += qp->gz->memsize;

	CREATE_STRVEC(ne, qp->b, c_ptr);
	c_ptr += qp->b->memsize;

	CREATE_STRVEC(2*nb+2*ng+2*ns, qp->d, c_ptr);
	c_ptr += qp->d->memsize;

	CREATE_STRVEC(2*nb+2*ng+2*ns, qp->m, c_ptr);
	c_ptr += qp->m->memsize;

	CREATE_STRVEC(2*ns, qp->Z, c_ptr);
	c_ptr += qp->Z->memsize;


	qp->dim = dim;

	qp->memsize = DENSE_QP_MEMSIZE(dim);


#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + qp->memsize)
		{
		printf("\nCreate_ocp_qp: outside memory bounds!\n\n");
		exit(1);
		}
#endif


	return;

	}



void DENSE_QP_SET_ALL(REAL *H, REAL *g, REAL *A, REAL *b, int *idxb, REAL *d_lb, REAL *d_ub, REAL *C, REAL *d_lg, REAL *d_ug, REAL *Zl, REAL *Zu, REAL *zl, REAL *zu, int *idxs, REAL *d_ls, REAL *d_us, struct DENSE_QP *qp)
	{

	int ii, jj;

	int nv = qp->dim->nv;
	int ne = qp->dim->ne;
	int nb = qp->dim->nb;
	int ng = qp->dim->ng;
	int ns = qp->dim->ns;

	CVT_MAT2STRMAT(nv, nv, H, nv, qp->Hv, 0, 0);
	CVT_VEC2STRVEC(nv, g, qp->gz, 0);
	if(ne>0)
		{
		CVT_MAT2STRMAT(ne, nv, A, ne, qp->A, 0, 0);
		CVT_VEC2STRVEC(ne, b, qp->b, 0);
		}
	if(nb>0)
		{
		for(ii=0; ii<nb; ii++) qp->idxb[ii] = idxb[ii];
		CVT_VEC2STRVEC(nb, d_lb, qp->d, 0);
		CVT_VEC2STRVEC(nb, d_ub, qp->d, nb+ng);
		VECSC_LIBSTR(nb, -1.0, qp->d, nb+ng);
		VECSE_LIBSTR(nb, 0.0, qp->m, 0);
		VECSE_LIBSTR(nb, 0.0, qp->m, nb+ng);
		}
	if(ng>0)
		{
		CVT_TRAN_MAT2STRMAT(ng, nv, C, ng, qp->Ct, 0, 0);
		CVT_VEC2STRVEC(ng, d_lg, qp->d, nb);
		CVT_VEC2STRVEC(ng, d_ug, qp->d, 2*nb+ng);
		VECSC_LIBSTR(ng, -1.0, qp->d, 2*nb+ng);
		VECSE_LIBSTR(ng, 0.0, qp->m, nb);
		VECSE_LIBSTR(ng, 0.0, qp->m, 2*nb+ng);
		}
	if(ns>0)
		{
		for(ii=0; ii<ns; ii++) qp->idxs[ii] = idxs[ii];
		CVT_VEC2STRVEC(ns, Zl, qp->Z, 0);
		CVT_VEC2STRVEC(ns, Zu, qp->Z, ns);
		CVT_VEC2STRVEC(ns, zl, qp->gz, nv);
		CVT_VEC2STRVEC(ns, zu, qp->gz, nv+ns);
		CVT_VEC2STRVEC(ns, d_ls, qp->d, 2*nb+2*ng);
		CVT_VEC2STRVEC(ns, d_us, qp->d, 2*nb+2*ng+ns);
		VECSE_LIBSTR(ns, 0.0, qp->m, 2*nb+2*ng);
		VECSE_LIBSTR(ns, 0.0, qp->m, 2*nb+2*ng+ns);
		}

	return;

	}



void DENSE_QP_GET_ALL(struct DENSE_QP *qp, REAL *H, REAL *g, REAL *A, REAL *b, int *idxb, REAL *d_lb, REAL *d_ub, REAL *C, REAL *d_lg, REAL *d_ug, REAL *Zl, REAL *Zu, REAL *zl, REAL *zu, int *idxs, REAL *d_ls, REAL *d_us)
	{

	int ii;

	int nv = qp->dim->nv;
	int ne = qp->dim->ne;
	int nb = qp->dim->nb;
	int ng = qp->dim->ng;
	int ns = qp->dim->ns;

	CVT_STRMAT2MAT(nv, nv, qp->Hv, 0, 0, H, nv);
	CVT_STRVEC2VEC(nv, qp->gz, 0, g);
	if(ne>0)
		{
		CVT_STRMAT2MAT(ne, nv, qp->A, 0, 0, A, ne);
		CVT_STRVEC2VEC(ne, qp->b, 0, b);
		}
	if(nb>0)
		{
		for(ii=0; ii<nb; ii++) idxb[ii] = qp->idxb[ii];
		CVT_STRVEC2VEC(nb, qp->d, 0, d_lb);
		CVT_STRVEC2VEC(nb, qp->d, nb+ng, d_ub);
		for(ii=0; ii<nb; ii++) d_ub[ii] = - d_ub[ii];
		}
	if(ng>0)
		{
		CVT_TRAN_STRMAT2MAT(nv, ng, qp->Ct, 0, 0, C, ng);
		CVT_STRVEC2VEC(ng, qp->d, nb, d_lg);
		CVT_STRVEC2VEC(ng, qp->d, 2*nb+ng, d_ug);
		for(ii=0; ii<ng; ii++) d_ug[ii] = - d_ug[ii];
		}
	if(ns>0)
		{
		for(ii=0; ii<ns; ii++) idxs[ii] = qp->idxs[ii];
		CVT_STRVEC2VEC(ns, qp->Z, 0, Zl);
		CVT_STRVEC2VEC(ns, qp->Z, ns, Zu);
		CVT_STRVEC2VEC(ns, qp->gz, nv, zl);
		CVT_STRVEC2VEC(ns, qp->gz, nv+ns, zu);
		CVT_STRVEC2VEC(ns, qp->d, 2*nb+2*ng, d_ls);
		CVT_STRVEC2VEC(ns, qp->d, 2*nb+2*ng+ns, d_us);
		}

	return;

	}



void DENSE_QP_SET_H(REAL *H, struct DENSE_QP *qp)
	{

	int nv = qp->dim->nv;

	CVT_MAT2STRMAT(nv, nv, H, nv, qp->Hv, 0, 0);

	return;

	}



void DENSE_QP_SET_G(REAL *g, struct DENSE_QP *qp)
	{

	int nv = qp->dim->nv;

	CVT_VEC2STRVEC(nv, g, qp->gz, 0);

	return;

	}



void DENSE_QP_SET_A(REAL *A, struct DENSE_QP *qp)
	{

	int nv = qp->dim->nv;
	int ne = qp->dim->ne;

	CVT_MAT2STRMAT(ne, nv, A, ne, qp->A, 0, 0);

	return;

	}



void DENSE_QP_SET_B(REAL *b, struct DENSE_QP *qp)
	{

	int ne = qp->dim->ne;

	CVT_VEC2STRVEC(ne, b, qp->b, 0);

	return;

	}



void DENSE_QP_SET_IDXB(int *idxb, struct DENSE_QP *qp)
	{

	int ii;
	int nb = qp->dim->nb;

	for(ii=0; ii<nb; ii++) qp->idxb[ii] = idxb[ii];

	return;

	}



void DENSE_QP_SET_LB(REAL *lb, struct DENSE_QP *qp)
	{

	int nb = qp->dim->nb;

	CVT_VEC2STRVEC(nb, lb, qp->d, 0);
	VECSE_LIBSTR(nb, 0.0, qp->m, 0);

	return;

	}



void DENSE_QP_SET_UB(REAL *ub, struct DENSE_QP *qp)
	{

	int nb = qp->dim->nb;
	int ng = qp->dim->ng;

	CVT_VEC2STRVEC(nb, ub, qp->d, nb+ng);
	VECSC_LIBSTR(nb, -1.0, qp->d, nb+ng);
	VECSE_LIBSTR(nb, 0.0, qp->m, nb+ng);

	return;

	}



void DENSE_QP_SET_C(REAL *C, struct DENSE_QP *qp)
	{

	int nv = qp->dim->nv;
	int ng = qp->dim->ng;

	CVT_TRAN_MAT2STRMAT(ng, nv, C, ng, qp->Ct, 0, 0);

	return;

	}



void DENSE_QP_SET_LG(REAL *lg, struct DENSE_QP *qp)
	{

	int nb = qp->dim->nb;
	int ng = qp->dim->ng;

	CVT_VEC2STRVEC(ng, lg, qp->d, nb);
	VECSE_LIBSTR(ng, 0.0, qp->m, nb);

	return;

	}



void DENSE_QP_SET_UG(REAL *ug, struct DENSE_QP *qp)
	{

	int nb = qp->dim->nb;
	int ng = qp->dim->ng;

	CVT_VEC2STRVEC(ng, ug, qp->d, 2*nb+ng);
	VECSC_LIBSTR(ng, -1.0, qp->d, 2*nb+ng);
	VECSE_LIBSTR(ng, 0.0, qp->m, 2*nb+ng);

	return;

	}



void DENSE_QP_SET_IDXS(int *idxs, struct DENSE_QP *qp)
	{

	int ii;
	int ns = qp->dim->ns;

	for(ii=0; ii<ns; ii++) qp->idxs[ii] = idxs[ii];

	return;

	}



void DENSE_QP_SET_ZZL(REAL *Zl, struct DENSE_QP *qp)
	{

	int ns = qp->dim->ns;

	CVT_VEC2STRVEC(ns, Zl, qp->Z, 0);

	return;

	}



void DENSE_QP_SET_ZZU(REAL *Zu, struct DENSE_QP *qp)
	{

	int ns = qp->dim->ns;

	CVT_VEC2STRVEC(ns, Zu, qp->Z, ns);

	return;

	}




void DENSE_QP_SET_ZL(REAL *zl, struct DENSE_QP *qp)
	{

	int ns = qp->dim->ns;
	int nv = qp->dim->nv;

	CVT_VEC2STRVEC(ns, zl, qp->gz, nv);

	return;

	}



void DENSE_QP_SET_ZU(REAL *zu, struct DENSE_QP *qp)
	{

	int ns = qp->dim->ns;
	int nv = qp->dim->nv;

	CVT_VEC2STRVEC(ns, zu, qp->gz, nv+ns);

	return;

	}


void DENSE_QP_SET_LS(REAL *ls, struct DENSE_QP *qp)
	{

	int ns = qp->dim->ns;
	int nb = qp->dim->nb;
	int ng = qp->dim->ng;

	CVT_VEC2STRVEC(ns, ls, qp->d, 2*nb+2*ng);
	VECSE_LIBSTR(ns, 0.0, qp->m, 2*nb+2*ng);

	return;

	}



void DENSE_QP_SET_US(REAL *us, struct DENSE_QP *qp)
	{

	int ns = qp->dim->ns;
	int nb = qp->dim->nb;
	int ng = qp->dim->ng;

	CVT_VEC2STRVEC(ns, us, qp->d, 2*nb+2*ng+ns);
	VECSE_LIBSTR(ns, 0.0, qp->m, 2*nb+2*ng+ns);

	return;

	}



void DENSE_QP_GET_H(struct DENSE_QP *qp, REAL *H)
	{

	int nv = qp->dim->nv;

	CVT_STRMAT2MAT(nv, nv, qp->Hv, 0, 0, H, nv);

	return;

	}



void DENSE_QP_GET_G(struct DENSE_QP *qp, REAL *g)
	{

	int nv = qp->dim->nv;

	CVT_STRVEC2VEC(nv, qp->gz, 0, g);

	return;

	}



void DENSE_QP_GET_A(struct DENSE_QP *qp, REAL *A)
	{

	int nv = qp->dim->nv;
	int ne = qp->dim->ne;

	CVT_STRMAT2MAT(ne, nv, qp->A, 0, 0, A, ne);

	return;

	}



void DENSE_QP_GET_B(struct DENSE_QP *qp, REAL *b)
	{

	int ne = qp->dim->ne;

	CVT_STRVEC2VEC(ne, qp->b, 0, b);

	return;

	}



void DENSE_QP_GET_IDXB(struct DENSE_QP *qp, int *idxb)
	{

	int ii;
	int nb = qp->dim->nb;

	for(ii=0; ii<nb; ii++) idxb[ii] = qp->idxb[ii];

	return;

	}



void DENSE_QP_GET_LB(struct DENSE_QP *qp, REAL *lb)
	{

	int nb = qp->dim->nb;

	CVT_STRVEC2VEC(nb, qp->d, 0, lb);

	return;

	}



void DENSE_QP_GET_UB(struct DENSE_QP *qp, REAL *ub)
	{

	int ii;
	int nb = qp->dim->nb;
	int ng = qp->dim->ng;

	CVT_STRVEC2VEC(nb, qp->d, nb+ng, ub);
	for(ii=0; ii<nb; ii++) ub[ii] = - ub[ii];

	return;

	}



void DENSE_QP_GET_C(struct DENSE_QP *qp, REAL *C)
	{

	int nv = qp->dim->nv;
	int ng = qp->dim->ng;

	CVT_TRAN_STRMAT2MAT(nv, ng, qp->Ct, 0, 0, C, ng);

	return;

	}



void DENSE_QP_GET_LG(struct DENSE_QP *qp, REAL *lg)
	{

	int nb = qp->dim->nb;
	int ng = qp->dim->ng;

	CVT_STRVEC2VEC(ng, qp->d, nb, lg);

	return;

	}



void DENSE_QP_GET_UG(struct DENSE_QP *qp, REAL *ug)
	{

	int ii;
	int nb = qp->dim->nb;
	int ng = qp->dim->ng;

	CVT_STRVEC2VEC(ng, qp->d, 2*nb+ng, ug);
	for(ii=0; ii<ng; ii++) ug[ii] = - ug[ii];

	return;

	}



void DENSE_QP_GET_IDXS(struct DENSE_QP *qp, int *idxs)
	{

	int ii;
	int ns = qp->dim->ns;

	for(ii=0; ii<ns; ii++) idxs[ii] = qp->idxs[ii];

	return;

	}



void DENSE_QP_GET_ZZL(struct DENSE_QP *qp, REAL *Zl)
	{

	int ns = qp->dim->ns;

	CVT_STRVEC2VEC(ns, qp->Z, 0, Zl);

	return;

	}



void DENSE_QP_GET_ZZU(struct DENSE_QP *qp, REAL *Zu)
	{

	int ns = qp->dim->ns;

	CVT_STRVEC2VEC(ns, qp->Z, ns, Zu);

	return;

	}




void DENSE_QP_GET_ZL(struct DENSE_QP *qp, REAL *zl)
	{

	int ns = qp->dim->ns;
	int nv = qp->dim->nv;

	CVT_STRVEC2VEC(ns, qp->gz, nv, zl);

	return;

	}



void DENSE_QP_GET_ZU(struct DENSE_QP *qp, REAL *zu)
	{

	int ns = qp->dim->ns;
	int nv = qp->dim->nv;

	CVT_STRVEC2VEC(ns, qp->gz, nv+ns, zu);

	return;

	}


void DENSE_QP_GET_LS(struct DENSE_QP *qp, REAL *ls)
	{

	int ns = qp->dim->ns;
	int nb = qp->dim->nb;
	int ng = qp->dim->ng;

	CVT_STRVEC2VEC(ns, qp->d, 2*nb+2*ng, ls);

	return;

	}



void DENSE_QP_GET_US(struct DENSE_QP *qp, REAL *us)
	{

	int ns = qp->dim->ns;
	int nb = qp->dim->nb;
	int ng = qp->dim->ng;

	CVT_STRVEC2VEC(ns, qp->d, 2*nb+2*ng+ns, us);

	return;

	}



void DENSE_QP_SET_ALL_ROWMAJ(REAL *H, REAL *g, REAL *A, REAL *b, int *idxb, REAL *d_lb, REAL *d_ub, REAL *C, REAL *d_lg, REAL *d_ug, REAL *Zl, REAL *Zu, REAL *zl, REAL *zu, int *idxs, REAL *d_ls, REAL *d_us, struct DENSE_QP *qp)
	{

	int ii;

	int nv = qp->dim->nv;
	int ne = qp->dim->ne;
	int nb = qp->dim->nb;
	int ng = qp->dim->ng;
	int ns = qp->dim->ns;

	CVT_TRAN_MAT2STRMAT(nv, nv, H, nv, qp->Hv, 0, 0);
	CVT_VEC2STRVEC(nv, g, qp->gz, 0);
	if(ne>0)
		{
		CVT_TRAN_MAT2STRMAT(nv, ne, A, nv, qp->A, 0, 0);
		CVT_VEC2STRVEC(ne, b, qp->b, 0);
		}
	if(nb>0)
		{
		for(ii=0; ii<nb; ii++) qp->idxb[ii] = idxb[ii];
		CVT_VEC2STRVEC(nb, d_lb, qp->d, 0);
		CVT_VEC2STRVEC(nb, d_ub, qp->d, nb+ng);
		VECSC_LIBSTR(nb, -1.0, qp->d, nb+ng);
		VECSE_LIBSTR(nb, 0.0, qp->m, 0);
		VECSE_LIBSTR(nb, 0.0, qp->m, nb+ng);
		}
	if(ng>0)
		{
		CVT_MAT2STRMAT(nv, ng, C, nv, qp->Ct, 0, 0);
		CVT_VEC2STRVEC(ng, d_lg, qp->d, nb);
		CVT_VEC2STRVEC(ng, d_ug, qp->d, 2*nb+ng);
		VECSC_LIBSTR(ng, -1.0, qp->d, 2*nb+ng);
		VECSE_LIBSTR(ng, 0.0, qp->m, nb);
		VECSE_LIBSTR(ng, 0.0, qp->m, 2*nb+ng);
		}
	if(ns>0)
		{
		for(ii=0; ii<ns; ii++) qp->idxs[ii] = idxs[ii];
		CVT_VEC2STRVEC(ns, Zl, qp->Z, 0);
		CVT_VEC2STRVEC(ns, Zu, qp->Z, ns);
		CVT_VEC2STRVEC(ns, zl, qp->gz, nv);
		CVT_VEC2STRVEC(ns, zu, qp->gz, nv+ns);
		CVT_VEC2STRVEC(ns, d_ls, qp->d, 2*nb+2*ng);
		CVT_VEC2STRVEC(ns, d_us, qp->d, 2*nb+2*ng+ns);
		VECSE_LIBSTR(ns, 0.0, qp->m, 2*nb+2*ng);
		VECSE_LIBSTR(ns, 0.0, qp->m, 2*nb+2*ng+ns);
		}

	return;

	}



void DENSE_QP_GET_ALL_ROWMAJ(struct DENSE_QP *qp, REAL *H, REAL *g, REAL *A, REAL *b, int *idxb, REAL *d_lb, REAL *d_ub, REAL *C, REAL *d_lg, REAL *d_ug, REAL *Zl, REAL *Zu, REAL *zl, REAL *zu, int *idxs, REAL *d_ls, REAL *d_us)
	{

	int ii;

	int nv = qp->dim->nv;
	int ne = qp->dim->ne;
	int nb = qp->dim->nb;
	int ng = qp->dim->ng;
	int ns = qp->dim->ns;

	CVT_TRAN_STRMAT2MAT(nv, nv, qp->Hv, 0, 0, H, nv);
	CVT_STRVEC2VEC(nv, qp->gz, 0, g);
	if(ne>0)
		{
		CVT_TRAN_STRMAT2MAT(ne, nv, qp->A, 0, 0, A, nv);
		CVT_STRVEC2VEC(ne, qp->b, 0, b);
		}
	if(nb>0)
		{
		for(ii=0; ii<nb; ii++) idxb[ii] = qp->idxb[ii];
		CVT_STRVEC2VEC(nb, qp->d, 0, d_lb);
		CVT_STRVEC2VEC(nb, qp->d, nb+ng, d_ub);
		for(ii=0; ii<nb; ii++) d_ub[ii] = - d_ub[ii];
		}
	if(ng>0)
		{
		CVT_STRMAT2MAT(nv, ng, qp->Ct, 0, 0, C, nv);
		CVT_STRVEC2VEC(ng, qp->d, nb, d_lg);
		CVT_STRVEC2VEC(ng, qp->d, 2*nb+ng, d_ug);
		for(ii=0; ii<ng; ii++) d_ug[ii] = - d_ug[ii];
		}
	if(ns>0)
		{
		for(ii=0; ii<ns; ii++) idxs[ii] = qp->idxs[ii];
		CVT_STRVEC2VEC(ns, qp->Z, 0, Zl);
		CVT_STRVEC2VEC(ns, qp->Z, ns, Zu);
		CVT_STRVEC2VEC(ns, qp->gz, nv, zl);
		CVT_STRVEC2VEC(ns, qp->gz, nv+ns, zu);
		CVT_STRVEC2VEC(ns, qp->d, 2*nb+2*ng, d_ls);
		CVT_STRVEC2VEC(ns, qp->d, 2*nb+2*ng+ns, d_us);
		}

	return;

	}

