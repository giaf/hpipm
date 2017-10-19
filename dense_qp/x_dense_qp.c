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



int MEMSIZE_DENSE_QP(int nv, int ne, int nb, int ng, int ns)
	{

	int size = 0;

	size += 5*sizeof(struct STRVEC); // g b d Z z
	size += 3*sizeof(struct STRMAT); // Hg A Ct

	size += 1*SIZE_STRVEC(nv); // g
	size += 1*SIZE_STRVEC(ne); // b
	size += 1*SIZE_STRVEC(2*nb+2*ng); // d
	size += 2*SIZE_STRVEC(2*ns); // Z z
	size += 1*nb*sizeof(int); // idxb
	size += 1*ns*sizeof(int); // idxb

	size += 1*SIZE_STRMAT(nv+1, nv); // Hg
	size += 1*SIZE_STRMAT(ne, nv); // A
	size += 1*SIZE_STRMAT(nv, ng); // Ct

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size
	
	return size;

	}



void CREATE_DENSE_QP(int nv, int ne, int nb, int ng, int ns, struct DENSE_QP *qp, void *mem)
	{


	// problem size
	qp->nv = nv;
	qp->ne = ne;
	qp->nb = nb;
	qp->ng = ng;
	qp->ns = ns;


	// matrix struct stuff
	struct STRMAT *sm_ptr = (struct STRMAT *) mem;

	qp->Hg = sm_ptr;
	sm_ptr += 1;

	qp->A = sm_ptr;
	sm_ptr += 1;

	qp->Ct = sm_ptr;
	sm_ptr += 1;


	// vector struct stuff
	struct STRVEC *sv_ptr = (struct STRVEC *) sm_ptr;

	qp->g = sv_ptr;
	sv_ptr += 1;

	qp->b = sv_ptr;
	sv_ptr += 1;

	qp->d = sv_ptr;
	sv_ptr += 1;

	qp->Z = sv_ptr;
	sv_ptr += 1;

	qp->z = sv_ptr;
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

	CREATE_STRMAT(nv+1, nv, qp->Hg, c_ptr);
	c_ptr += qp->Hg->memory_size;

	CREATE_STRMAT(ne, nv, qp->A, c_ptr);
	c_ptr += qp->A->memory_size;

	CREATE_STRMAT(nv, ng, qp->Ct, c_ptr);
	c_ptr += qp->Ct->memory_size;

	CREATE_STRVEC(nv, qp->g, c_ptr);
	c_ptr += qp->g->memory_size;

	CREATE_STRVEC(ne, qp->b, c_ptr);
	c_ptr += qp->b->memory_size;

	CREATE_STRVEC(2*nb+2*ng, qp->d, c_ptr);
	c_ptr += qp->d->memory_size;

	CREATE_STRVEC(2*ns, qp->Z, c_ptr);
	c_ptr += qp->Z->memory_size;

	CREATE_STRVEC(2*ns, qp->z, c_ptr);
	c_ptr += qp->z->memory_size;


	qp->memsize = MEMSIZE_DENSE_QP(nv, ne, nb, ng, ns);


#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + qp->memsize)
		{
		printf("\nCreate_ocp_qp: outsize memory bounds!\n\n");
		exit(1);
		}
#endif


	return;

	}



void CVT_COLMAJ_TO_DENSE_QP(REAL *H, REAL *g, REAL *A, REAL *b, int *idxb, REAL *d_lb, REAL *d_ub, REAL *C, REAL *d_lg, REAL *d_ug, REAL *Zl, REAL *Zu, REAL *zl, REAL *zu, int *idxs, struct DENSE_QP *qp)
	{

	int ii;

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;
	int ns = qp->ns;

	CVT_MAT2STRMAT(nv, nv, H, nv, qp->Hg, 0, 0);
	CVT_TRAN_MAT2STRMAT(nv, 1, g, nv, qp->Hg, nv, 0);
	CVT_VEC2STRVEC(nv, g, qp->g, 0);
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
		}
	if(ng>0)
		{
		CVT_TRAN_MAT2STRMAT(ng, nv, C, ng, qp->Ct, 0, 0);
		CVT_VEC2STRVEC(ng, d_lg, qp->d, nb);
		CVT_VEC2STRVEC(ng, d_ug, qp->d, 2*nb+ng);
		}
	if(ns>0)
		{
		for(ii=0; ii<ns; ii++) qp->idxs[ii] = idxs[ii];
		CVT_VEC2STRVEC(ns, Zl, qp->Z, 0);
		CVT_VEC2STRVEC(ns, Zu, qp->Z, ns);
		CVT_VEC2STRVEC(ns, zl, qp->z, 0);
		CVT_VEC2STRVEC(ns, zu, qp->z, ns);
		}

	return;

	}



void CVT_DENSE_QP_TO_COLMAJ(struct DENSE_QP *qp, REAL *H, REAL *g, REAL *A, REAL *b, int *idxb, REAL *d_lb, REAL *d_ub, REAL *C, REAL *d_lg, REAL *d_ug, REAL *Zl, REAL *Zu, REAL *zl, REAL *zu, int *idxs)
	{

	int ii;

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;
	int ns = qp->ns;

	CVT_STRMAT2MAT(nv, nv, qp->Hg, 0, 0, H, nv);
	CVT_STRVEC2VEC(nv, qp->g, 0, g);
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
		}
	if(ng>0)
		{
		CVT_TRAN_STRMAT2MAT(nv, ng, qp->Ct, 0, 0, C, ng);
		CVT_STRVEC2VEC(ng, qp->d, nb, d_lg);
		CVT_STRVEC2VEC(ng, qp->d, 2*nb+ng, d_ug);
		}
	if(ns>0)
		{
		for(ii=0; ii<ns; ii++) qp->idxs[ii] = idxs[ii];
		CVT_STRVEC2VEC(ns, qp->Z, 0, Zl);
		CVT_STRVEC2VEC(ns, qp->Z, ns, Zu);
		CVT_STRVEC2VEC(ns, qp->z, 0, zl);
		CVT_STRVEC2VEC(ns, qp->z, ns, zu);
		}

	return;

	}



void CVT_ROWMAJ_TO_DENSE_QP(REAL *H, REAL *g, REAL *A, REAL *b, int *idxb, REAL *d_lb, REAL *d_ub, REAL *C, REAL *d_lg, REAL *d_ug, REAL *Zl, REAL *Zu, REAL *zl, REAL *zu, int *idxs, struct DENSE_QP *qp)
	{

	int ii;

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;
	int ns = qp->ns;

	CVT_TRAN_MAT2STRMAT(nv, nv, H, nv, qp->Hg, 0, 0);
	CVT_TRAN_MAT2STRMAT(nv, 1, g, nv, qp->Hg, nv, 0);
	CVT_VEC2STRVEC(nv, g, qp->g, 0);
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
		}
	if(ng>0)
		{
		CVT_MAT2STRMAT(nv, ng, C, nv, qp->Ct, 0, 0);
		CVT_VEC2STRVEC(ng, d_lg, qp->d, nb);
		CVT_VEC2STRVEC(ng, d_ug, qp->d, 2*nb+ng);
		}
	if(ns>0)
		{
		for(ii=0; ii<ns; ii++) qp->idxs[ii] = idxs[ii];
		CVT_VEC2STRVEC(ns, Zl, qp->Z, 0);
		CVT_VEC2STRVEC(ns, Zu, qp->Z, ns);
		CVT_VEC2STRVEC(ns, zl, qp->z, 0);
		CVT_VEC2STRVEC(ns, zu, qp->z, ns);
		}

	return;

	}



void CVT_DENSE_QP_TO_ROWMAJ(struct DENSE_QP *qp, REAL *H, REAL *g, REAL *A, REAL *b, int *idxb, REAL *d_lb, REAL *d_ub, REAL *C, REAL *d_lg, REAL *d_ug, REAL *Zl, REAL *Zu, REAL *zl, REAL *zu, int *idxs)
	{

	int ii;

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;
	int ns = qp->ns;

	CVT_TRAN_STRMAT2MAT(nv, nv, qp->Hg, 0, 0, H, nv);
	CVT_STRVEC2VEC(nv, qp->g, 0, g);
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
		}
	if(ng>0)
		{
		CVT_STRMAT2MAT(nv, ng, qp->Ct, 0, 0, C, nv);
		CVT_STRVEC2VEC(ng, qp->d, nb, d_lg);
		CVT_STRVEC2VEC(ng, qp->d, 2*nb+ng, d_ug);
		}
	if(ns>0)
		{
		for(ii=0; ii<ns; ii++) qp->idxs[ii] = idxs[ii];
		CVT_STRVEC2VEC(ns, qp->Z, 0, Zl);
		CVT_STRVEC2VEC(ns, qp->Z, ns, Zu);
		CVT_STRVEC2VEC(ns, qp->z, 0, zl);
		CVT_STRVEC2VEC(ns, qp->z, ns, zu);
		}

	return;

	}



void CVT_LIBSTR_TO_DENSE_QP(struct STRMAT *H, struct STRMAT *A, struct STRMAT *C, struct STRVEC *g, struct STRVEC *b, struct STRVEC *d_lb, struct STRVEC *d_ub, struct STRVEC *d_lg, struct STRVEC *d_ug, int *idxb, struct STRVEC *Zl, struct STRVEC *Zu, struct STRVEC *zl, struct STRVEC *zu, int *idxs, struct DENSE_QP *qp)
	{

	int ii;

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;
	int ns = qp->ns;

	GECP_LIBSTR(nv, nv, H, 0, 0, qp->Hg, 0, 0);
	ROWIN_LIBSTR(nv, 1.0, g, 0, qp->Hg, nv, 0);
	VECCP_LIBSTR(nv, g, 0, qp->g, 0);
	if(ne>0)
		{
		GECP_LIBSTR(ne, nv, A, 0, 0, qp->A, 0, 0);
		VECCP_LIBSTR(ne, b, 0, qp->b, 0);
		}
	if(nb>0)
		{
		for(ii=0; ii<nb; ii++) qp->idxb[ii] = idxb[ii];
		VECCP_LIBSTR(nb, d_lb, 0, qp->d, 0);
		VECCP_LIBSTR(nb, d_ub, 0, qp->d, nb+ng);
		}
	if(ng>0)
		{
		GETR_LIBSTR(ng, nv, C, 0, 0, qp->Ct, 0, 0);
		VECCP_LIBSTR(ng, d_lg, 0, qp->d, nb);
		VECCP_LIBSTR(ng, d_ug, 0, qp->d, 2*nb+ng);
		}
	if(ns>0)
		{
		for(ii=0; ii<ns; ii++) qp->idxs[ii] = idxs[ii];
		VECCP_LIBSTR(ns, Zl, 0, qp->Z, 0);
		VECCP_LIBSTR(ns, Zu, 0, qp->Z, ns);
		VECCP_LIBSTR(ns, zl, 0, qp->z, 0);
		VECCP_LIBSTR(ns, zu, 0, qp->z, ns);
		}

	return;

	}



void CVT_DENSE_QP_TO_LIBSTR(struct DENSE_QP *qp, struct STRMAT *H, struct STRMAT *A, struct STRMAT *C, struct STRVEC *g, struct STRVEC *b, struct STRVEC *d_lb, struct STRVEC *d_ub, struct STRVEC *d_lg, struct STRVEC *d_ug, int *idxb, struct STRVEC *Zl, struct STRVEC *Zu, struct STRVEC *zl, struct STRVEC *zu, int *idxs)
	{

	int ii;

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;
	int ns = qp->ns;

	GECP_LIBSTR(nv, nv, qp->Hg, 0, 0, H, 0, 0);
	VECCP_LIBSTR(nv, qp->g, 0, g, 0);
	if(ne>0)
		{
		GECP_LIBSTR(ne, nv, qp->A, 0, 0, A, 0, 0);
		VECCP_LIBSTR(ne, qp->b, 0, b, 0);
		}
	if(nb>0)
		{
		for(ii=0; ii<nb; ii++) idxb[ii] = qp->idxb[ii];
		VECCP_LIBSTR(nb, qp->d, 0, d_lb, 0);
		VECCP_LIBSTR(nb, qp->d, nb+ng, d_ub, 0);
		}
	if(ng>0)
		{
		GETR_LIBSTR(nv, ng, qp->Ct, 0, 0, C, 0, 0);
		VECCP_LIBSTR(ng, qp->d, nb, d_lg, 0);
		VECCP_LIBSTR(ng, qp->d, 2*nb+ng, d_ug, 0);
		}
	if(ns>0)
		{
		for(ii=0; ii<ns; ii++) qp->idxs[ii] = idxs[ii];
		VECCP_LIBSTR(ns, qp->Z, 0, Zl, 0);
		VECCP_LIBSTR(ns, qp->Z, ns, Zu, 0);
		VECCP_LIBSTR(ns, qp->z, 0, zl, 0);
		VECCP_LIBSTR(ns, qp->z, ns, zu, 0);
		}

	return;

	}


