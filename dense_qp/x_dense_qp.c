/**************************************************************************************************
*                                                                                                 *
* This file is part of HPIPM.                                                                     *
*                                                                                                 *
* HPIPM -- High Performance Interior Point Method.                                                *
* Copyright (C) 2017 by Gianluca Frison.                                                          *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
*                                                                                                 *
* HPIPM is free software; you can redistribute it and/or                                          *
* modify it under the terms of the GNU Lesser General Public                                      *
* License as published by the Free Software Foundation; either                                    *
* version 2.1 of the License, or (at your option) any later version.                              *
*                                                                                                 *
* HPIPM is distributed in the hope that it will be useful,                                        *
* but WITHOUT ANY WARRANTY; without even the implied warranty of                                  *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                                            *
* See the GNU Lesser General Public License for more details.                                     *
*                                                                                                 *
* You should have received a copy of the GNU Lesser General Public                                *
* License along with HPIPM; if not, write to the Free Software                                    *
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA                  *
*                                                                                                 *
* Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             *
*                                                                                                 *
**************************************************************************************************/



int MEMSIZE_DENSE_QP(struct DENSE_QP_DIM *dim)
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



void CREATE_DENSE_QP(struct DENSE_QP_DIM *dim, struct DENSE_QP *qp, void *mem)
	{

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

	qp->memsize = MEMSIZE_DENSE_QP(dim);


#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + qp->memsize)
		{
		printf("\nCreate_ocp_qp: outsize memory bounds!\n\n");
		exit(1);
		}
#endif


	return;

	}



void CVT_COLMAJ_TO_DENSE_QP(REAL *H, REAL *g, REAL *A, REAL *b, int *idxb, REAL *d_lb, REAL *d_ub, REAL *C, REAL *d_lg, REAL *d_ug, REAL *Zl, REAL *Zu, REAL *zl, REAL *zu, int *idxs, REAL *d_ls, REAL *d_us, struct DENSE_QP *qp)
	{

	int ii;

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



void CVT_DENSE_QP_TO_COLMAJ(struct DENSE_QP *qp, REAL *H, REAL *g, REAL *A, REAL *b, int *idxb, REAL *d_lb, REAL *d_ub, REAL *C, REAL *d_lg, REAL *d_ug, REAL *Zl, REAL *Zu, REAL *zl, REAL *zu, int *idxs, REAL *d_ls, REAL *d_us)
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
		for(ii=0; ii<ns; ii++) qp->idxs[ii] = idxs[ii];
		CVT_STRVEC2VEC(ns, qp->Z, 0, Zl);
		CVT_STRVEC2VEC(ns, qp->Z, ns, Zu);
		CVT_STRVEC2VEC(ns, qp->gz, nv, zl);
		CVT_STRVEC2VEC(ns, qp->gz, nv+ns, zu);
		CVT_STRVEC2VEC(ns, qp->d, 2*nb+2*ng, d_ls);
		CVT_STRVEC2VEC(ns, qp->d, 2*nb+2*ng+ns, d_us);
		}

	return;

	}



void CVT_ROWMAJ_TO_DENSE_QP(REAL *H, REAL *g, REAL *A, REAL *b, int *idxb, REAL *d_lb, REAL *d_ub, REAL *C, REAL *d_lg, REAL *d_ug, REAL *Zl, REAL *Zu, REAL *zl, REAL *zu, int *idxs, REAL *d_ls, REAL *d_us, struct DENSE_QP *qp)
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



void CVT_DENSE_QP_TO_ROWMAJ(struct DENSE_QP *qp, REAL *H, REAL *g, REAL *A, REAL *b, int *idxb, REAL *d_lb, REAL *d_ub, REAL *C, REAL *d_lg, REAL *d_ug, REAL *Zl, REAL *Zu, REAL *zl, REAL *zu, int *idxs, REAL *d_ls, REAL *d_us)
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

