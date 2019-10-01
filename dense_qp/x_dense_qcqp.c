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



int DENSE_QCQP_MEMSIZE(struct DENSE_QCQP_DIM *dim)
	{

	int nv = dim->nv;
	int ne = dim->ne;
	int nb = dim->nb;
	int ng = dim->ng;
	int nq = dim->nq;
	int ns = dim->ns;

	int size = 0;

	size += (nq+5)*sizeof(struct STRVEC); // gz b d m Z gq
	size += (nq+3)*sizeof(struct STRMAT); // Hv A Ct Hq

	size += 1*SIZE_STRVEC(nv+2*ns); // g
	size += 1*SIZE_STRVEC(ne); // b
	size += 2*SIZE_STRVEC(2*nb+2*ng+2*nq+2*ns); // d m
	size += 1*SIZE_STRVEC(2*ns); // Z
	size += nq*SIZE_STRVEC(nv); // gq
	size += 1*nb*sizeof(int); // idxb
	size += 1*ns*sizeof(int); // idxb

	size += 1*SIZE_STRMAT(nv+1, nv); // Hv
	size += nq*SIZE_STRMAT(nv, nv); // Hq
	size += 1*SIZE_STRMAT(ne, nv); // A
	size += 1*SIZE_STRMAT(nv, ng); // Ct

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

	}



void DENSE_QCQP_CREATE(struct DENSE_QCQP_DIM *dim, struct DENSE_QCQP *qp, void *mem)
	{

	// TODO set memory to zero !!!!!!!

	int ii;

	int nv = dim->nv;
	int ne = dim->ne;
	int nb = dim->nb;
	int ng = dim->ng;
	int nq = dim->nq;
	int ns = dim->ns;


	// matrix struct stuff
	struct STRMAT *sm_ptr = (struct STRMAT *) mem;

	qp->Hv = sm_ptr;
	sm_ptr += 1;

	qp->A = sm_ptr;
	sm_ptr += 1;

	qp->Ct = sm_ptr;
	sm_ptr += 1;

	qp->Hq = sm_ptr;
	sm_ptr += nq;


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

	qp->gq = sv_ptr;
	sv_ptr += nq;


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

	for(ii=0; ii<nq; ii++)
		{
		CREATE_STRMAT(nv, nv, qp->Hq+ii, c_ptr);
		c_ptr += (qp->Hq+ii)->memsize;
		}

	CREATE_STRVEC(nv+2*ns, qp->gz, c_ptr);
	c_ptr += qp->gz->memsize;

	CREATE_STRVEC(ne, qp->b, c_ptr);
	c_ptr += qp->b->memsize;

	CREATE_STRVEC(2*nb+2*ng+2*nq+2*ns, qp->d, c_ptr);
	c_ptr += qp->d->memsize;

	CREATE_STRVEC(2*nb+2*ng+2*nq+2*ns, qp->m, c_ptr);
	c_ptr += qp->m->memsize;

	CREATE_STRVEC(2*ns, qp->Z, c_ptr);
	c_ptr += qp->Z->memsize;

	for(ii=0; ii<nq; ii++)
		{
		CREATE_STRVEC(nv, qp->gq+ii, c_ptr);
		c_ptr += (qp->gq+ii)->memsize;
		}
	
	// default init
	// TODO put to a larger value, and check that it doesn't make convergence worse
	// TODO compute unconstr sol and set it just below
	REAL inf = 1e3;
	VECSE(nq, -inf, qp->d, nb+ng);


	qp->dim = dim;

	qp->memsize = DENSE_QCQP_MEMSIZE(dim);


#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + qp->memsize)
		{
		printf("\ndense_qcqp_create: outside memory bounds!\n\n");
		exit(1);
		}
#endif


	return;

	}



void DENSE_QCQP_SET_H(REAL *H, struct DENSE_QCQP *qp)
	{

	int nv = qp->dim->nv;

	CVT_MAT2STRMAT(nv, nv, H, nv, qp->Hv, 0, 0);

	return;

	}



void DENSE_QCQP_SET_G(REAL *g, struct DENSE_QCQP *qp)
	{

	int nv = qp->dim->nv;

	CVT_VEC2STRVEC(nv, g, qp->gz, 0);

	return;

	}



void DENSE_QCQP_SET_A(REAL *A, struct DENSE_QCQP *qp)
	{

	int nv = qp->dim->nv;
	int ne = qp->dim->ne;

	CVT_MAT2STRMAT(ne, nv, A, ne, qp->A, 0, 0);

	return;

	}



void DENSE_QCQP_SET_B(REAL *b, struct DENSE_QCQP *qp)
	{

	int ne = qp->dim->ne;

	CVT_VEC2STRVEC(ne, b, qp->b, 0);

	return;

	}



void DENSE_QCQP_SET_IDXB(int *idxb, struct DENSE_QCQP *qp)
	{

	int ii;
	int nb = qp->dim->nb;

	for(ii=0; ii<nb; ii++) qp->idxb[ii] = idxb[ii];

	return;

	}



void DENSE_QCQP_SET_LB(REAL *lb, struct DENSE_QCQP *qp)
	{

	int nb = qp->dim->nb;

	CVT_VEC2STRVEC(nb, lb, qp->d, 0);
	VECSE(nb, 0.0, qp->m, 0);

	return;

	}



void DENSE_QCQP_SET_UB(REAL *ub, struct DENSE_QCQP *qp)
	{

	int nb = qp->dim->nb;
	int ng = qp->dim->ng;
	int nq = qp->dim->nq;

	CVT_VEC2STRVEC(nb, ub, qp->d, nb+ng+nq);
	VECSC(nb, -1.0, qp->d, nb+ng+nq);
	VECSE(nb, 0.0, qp->m, nb+ng+nq);

	return;

	}



void DENSE_QCQP_SET_C(REAL *C, struct DENSE_QCQP *qp)
	{

	int nv = qp->dim->nv;
	int ng = qp->dim->ng;

	CVT_TRAN_MAT2STRMAT(ng, nv, C, ng, qp->Ct, 0, 0);

	return;

	}



void DENSE_QCQP_SET_LG(REAL *lg, struct DENSE_QCQP *qp)
	{

	int nb = qp->dim->nb;
	int ng = qp->dim->ng;

	CVT_VEC2STRVEC(ng, lg, qp->d, nb);
	VECSE(ng, 0.0, qp->m, nb);

	return;

	}



void DENSE_QCQP_SET_UG(REAL *ug, struct DENSE_QCQP *qp)
	{

	int nb = qp->dim->nb;
	int ng = qp->dim->ng;
	int nq = qp->dim->nq;

	CVT_VEC2STRVEC(ng, ug, qp->d, 2*nb+ng+nq);
	VECSC(ng, -1.0, qp->d, 2*nb+ng+nq);
	VECSE(ng, 0.0, qp->m, 2*nb+ng+nq);

	return;

	}



void DENSE_QCQP_SET_HQ(REAL *HQ, struct DENSE_QCQP *qp)
	{

	int ii;

	int nv = qp->dim->nv;
	int nq = qp->dim->nq;

	for(ii=0; ii<nq; ii++)
		CVT_MAT2STRMAT(nv, nv, HQ+ii*nv*nv, nv, qp->Hq+ii, 0, 0);

	return;

	}



void DENSE_QCQP_SET_GQ(REAL *gq, struct DENSE_QCQP *qp)
	{

	int ii;

	int nv = qp->dim->nv;
	int nq = qp->dim->nq;

	for(ii=0; ii<nq; ii++)
		CVT_VEC2STRVEC(nv, gq+ii*nv, qp->gq+ii, 0);

	return;

	}



void DENSE_QCQP_SET_UQ(REAL *uq, struct DENSE_QCQP *qp)
	{

	int nb = qp->dim->nb;
	int ng = qp->dim->ng;
	int nq = qp->dim->nq;

	CVT_VEC2STRVEC(nq, uq, qp->d, 2*nb+2*ng+nq);
	VECSC(nq, -1.0, qp->d, 2*nb+2*ng+nq);
	VECSE(nq, 0.0, qp->m, 2*nb+2*ng+nq);

	// TODO put to a larger value, and check that it doesn't make convergence worse
	// TODO compute unconstr sol and set it just below
	REAL inf = 1e3;
	VECSE(nq, -inf, qp->d, nb+ng);
	VECSE(nq, 0.0, qp->m, nb+ng);

	return;

	}



void DENSE_QCQP_SET_IDXS(int *idxs, struct DENSE_QCQP *qp)
	{

	int ii;
	int ns = qp->dim->ns;

	for(ii=0; ii<ns; ii++) qp->idxs[ii] = idxs[ii];

	return;

	}



void DENSE_QCQP_SET_ZZL(REAL *Zl, struct DENSE_QCQP *qp)
	{

	int ns = qp->dim->ns;

	CVT_VEC2STRVEC(ns, Zl, qp->Z, 0);

	return;

	}



void DENSE_QCQP_SET_ZZU(REAL *Zu, struct DENSE_QCQP *qp)
	{

	int ns = qp->dim->ns;

	CVT_VEC2STRVEC(ns, Zu, qp->Z, ns);

	return;

	}




void DENSE_QCQP_SET_ZL(REAL *zl, struct DENSE_QCQP *qp)
	{

	int ns = qp->dim->ns;
	int nv = qp->dim->nv;

	CVT_VEC2STRVEC(ns, zl, qp->gz, nv);

	return;

	}



void DENSE_QCQP_SET_ZU(REAL *zu, struct DENSE_QCQP *qp)
	{

	int ns = qp->dim->ns;
	int nv = qp->dim->nv;

	CVT_VEC2STRVEC(ns, zu, qp->gz, nv+ns);

	return;

	}


void DENSE_QCQP_SET_LS(REAL *ls, struct DENSE_QCQP *qp)
	{

	int ns = qp->dim->ns;
	int nb = qp->dim->nb;
	int ng = qp->dim->ng;
	int nq = qp->dim->nq;

	CVT_VEC2STRVEC(ns, ls, qp->d, 2*nb+2*ng+2*nq);
	VECSE(ns, 0.0, qp->m, 2*nb+2*ng+2*nq);

	return;

	}



void DENSE_QCQP_SET_US(REAL *us, struct DENSE_QCQP *qp)
	{

	int ns = qp->dim->ns;
	int nb = qp->dim->nb;
	int ng = qp->dim->ng;
	int nq = qp->dim->nq;

	CVT_VEC2STRVEC(ns, us, qp->d, 2*nb+2*ng+2*nq+ns);
	VECSE(ns, 0.0, qp->m, 2*nb+2*ng+2*nq+ns);

	return;

	}



void DENSE_QCQP_GET_H(struct DENSE_QCQP *qp, REAL *H)
	{

	int nv = qp->dim->nv;

	CVT_STRMAT2MAT(nv, nv, qp->Hv, 0, 0, H, nv);

	return;

	}



void DENSE_QCQP_GET_G(struct DENSE_QCQP *qp, REAL *g)
	{

	int nv = qp->dim->nv;

	CVT_STRVEC2VEC(nv, qp->gz, 0, g);

	return;

	}



void DENSE_QCQP_GET_A(struct DENSE_QCQP *qp, REAL *A)
	{

	int nv = qp->dim->nv;
	int ne = qp->dim->ne;

	CVT_STRMAT2MAT(ne, nv, qp->A, 0, 0, A, ne);

	return;

	}



void DENSE_QCQP_GET_B(struct DENSE_QCQP *qp, REAL *b)
	{

	int ne = qp->dim->ne;

	CVT_STRVEC2VEC(ne, qp->b, 0, b);

	return;

	}



void DENSE_QCQP_GET_IDXB(struct DENSE_QCQP *qp, int *idxb)
	{

	int ii;
	int nb = qp->dim->nb;

	for(ii=0; ii<nb; ii++) idxb[ii] = qp->idxb[ii];

	return;

	}



void DENSE_QCQP_GET_LB(struct DENSE_QCQP *qp, REAL *lb)
	{

	int nb = qp->dim->nb;

	CVT_STRVEC2VEC(nb, qp->d, 0, lb);

	return;

	}



void DENSE_QCQP_GET_UB(struct DENSE_QCQP *qp, REAL *ub)
	{

	int ii;
	int nb = qp->dim->nb;
	int ng = qp->dim->ng;
	int nq = qp->dim->nq;

	CVT_STRVEC2VEC(nb, qp->d, nb+ng+nq, ub);
	for(ii=0; ii<nb; ii++) ub[ii] = - ub[ii];

	return;

	}



void DENSE_QCQP_GET_C(struct DENSE_QCQP *qp, REAL *C)
	{

	int nv = qp->dim->nv;
	int ng = qp->dim->ng;

	CVT_TRAN_STRMAT2MAT(nv, ng, qp->Ct, 0, 0, C, ng);

	return;

	}



void DENSE_QCQP_GET_LG(struct DENSE_QCQP *qp, REAL *lg)
	{

	int nb = qp->dim->nb;
	int ng = qp->dim->ng;

	CVT_STRVEC2VEC(ng, qp->d, nb, lg);

	return;

	}



void DENSE_QCQP_GET_UG(struct DENSE_QCQP *qp, REAL *ug)
	{

	int ii;
	int nb = qp->dim->nb;
	int ng = qp->dim->ng;
	int nq = qp->dim->nq;

	CVT_STRVEC2VEC(ng, qp->d, 2*nb+ng+nq, ug);
	for(ii=0; ii<ng; ii++) ug[ii] = - ug[ii];

	return;

	}



void DENSE_QCQP_GET_IDXS(struct DENSE_QCQP *qp, int *idxs)
	{

	int ii;
	int ns = qp->dim->ns;

	for(ii=0; ii<ns; ii++) idxs[ii] = qp->idxs[ii];

	return;

	}



void DENSE_QCQP_GET_ZZL(struct DENSE_QCQP *qp, REAL *Zl)
	{

	int ns = qp->dim->ns;

	CVT_STRVEC2VEC(ns, qp->Z, 0, Zl);

	return;

	}



void DENSE_QCQP_GET_ZZU(struct DENSE_QCQP *qp, REAL *Zu)
	{

	int ns = qp->dim->ns;

	CVT_STRVEC2VEC(ns, qp->Z, ns, Zu);

	return;

	}




void DENSE_QCQP_GET_ZL(struct DENSE_QCQP *qp, REAL *zl)
	{

	int ns = qp->dim->ns;
	int nv = qp->dim->nv;

	CVT_STRVEC2VEC(ns, qp->gz, nv, zl);

	return;

	}



void DENSE_QCQP_GET_ZU(struct DENSE_QCQP *qp, REAL *zu)
	{

	int ns = qp->dim->ns;
	int nv = qp->dim->nv;

	CVT_STRVEC2VEC(ns, qp->gz, nv+ns, zu);

	return;

	}


void DENSE_QCQP_GET_LS(struct DENSE_QCQP *qp, REAL *ls)
	{

	int ns = qp->dim->ns;
	int nb = qp->dim->nb;
	int ng = qp->dim->ng;
	int nq = qp->dim->nq;

	CVT_STRVEC2VEC(ns, qp->d, 2*nb+2*ng+2*nq, ls);

	return;

	}



void DENSE_QCQP_GET_US(struct DENSE_QCQP *qp, REAL *us)
	{

	int ns = qp->dim->ns;
	int nb = qp->dim->nb;
	int ng = qp->dim->ng;
	int nq = qp->dim->nq;

	CVT_STRVEC2VEC(ns, qp->d, 2*nb+2*ng+2*nq+ns, us);

	return;

	}




