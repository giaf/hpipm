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



int MEMSIZE_DENSE_QP(int nv, int ne, int nb, int ng)
	{

	int size = 0;

	size += 3*sizeof(struct STRVEC); // g b d
	size += 3*sizeof(struct STRMAT); // Hg A Ct

	size += 1*SIZE_STRVEC(nv); // g
	size += 1*SIZE_STRVEC(ne); // b
	size += 1*SIZE_STRVEC(2*nb+2*ng); // d
	size += 1*nb*sizeof(int); // idxb

	size += 1*SIZE_STRMAT(nv+1, nv); // Hg
	size += 1*SIZE_STRMAT(ne, nv); // A
	size += 1*SIZE_STRMAT(nv, ng); // Ct

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size
	
	return size;

	}



void CREATE_DENSE_QP(int nv, int ne, int nb, int ng, struct DENSE_QP *qp, void *memory)
	{

	qp->memsize = MEMSIZE_DENSE_QP(nv, ne, nb, ng);


	// problem size
	qp->nv = nv;
	qp->ne = ne;
	qp->nb = nb;
	qp->ng = ng;


	// matrix struct stuff
	struct STRMAT *sm_ptr = (struct STRMAT *) memory;

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


	// int stuff
	int *i_ptr;
	i_ptr = (int *) sv_ptr;

	// idxb
	qp->idxb = i_ptr;
	i_ptr += nb;


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

	return;

	}



void CVT_COLMAJ_TO_DENSE_QP(REAL *H, REAL *g, REAL *A, REAL *b, int *idxb, REAL *d_lb, REAL *d_ub, REAL *C, REAL *d_lg, REAL *d_ug, struct DENSE_QP *qp)
	{

	int ii;

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;

	CVT_MAT2STRMAT(nv, nv, H, nv, qp->Hg, 0, 0);
	CVT_TRAN_MAT2STRMAT(nv, 1, g, nv, qp->Hg, nv, 0);
	CVT_MAT2STRMAT(ne, nv, A, ne, qp->A, 0, 0);
	CVT_TRAN_MAT2STRMAT(ng, nv, C, ng, qp->Ct, 0, 0);
	CVT_VEC2STRVEC(nv, g, qp->g, 0);
	CVT_VEC2STRVEC(ne, b, qp->b, 0);
	CVT_VEC2STRVEC(nb, d_lb, qp->d, 0);
	CVT_VEC2STRVEC(nb, d_ub, qp->d, nb+ng);
	CVT_VEC2STRVEC(ng, d_lg, qp->d, nb);
	CVT_VEC2STRVEC(ng, d_ug, qp->d, 2*nb+ng);
	for(ii=0; ii<nb; ii++) qp->idxb[ii] = idxb[ii];

	return;

	}



void CVT_DENSE_QP_TO_COLMAJ(struct DENSE_QP *qp, REAL *H, REAL *g, REAL *A, REAL *b, int *idxb, REAL *d_lb, REAL *d_ub, REAL *C, REAL *d_lg, REAL *d_ug)
	{

	int ii;

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;

	CVT_STRMAT2MAT(nv, nv, qp->Hg, 0, 0, H, nv);
	CVT_STRMAT2MAT(ne, nv, qp->A, 0, 0, A, ne);
	CVT_TRAN_STRMAT2MAT(nv, ng, qp->Ct, 0, 0, C, ng);
	CVT_STRVEC2VEC(nv, qp->g, 0, g);
	CVT_STRVEC2VEC(ne, qp->b, 0, b);
	CVT_STRVEC2VEC(nb, qp->d, 0, d_lb);
	CVT_STRVEC2VEC(nb, qp->d, nb+ng, d_ub);
	CVT_STRVEC2VEC(ng, qp->d, nb, d_lg);
	CVT_STRVEC2VEC(ng, qp->d, 2*nb+ng, d_ug);
	for(ii=0; ii<nb; ii++) idxb[ii] = qp->idxb[ii];

	return;

	}



void CVT_ROWMAJ_TO_DENSE_QP(REAL *H, REAL *g, REAL *A, REAL *b, int *idxb, REAL *d_lb, REAL *d_ub, REAL *C, REAL *d_lg, REAL *d_ug, struct DENSE_QP *qp)
	{

	int ii;

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;

	CVT_TRAN_MAT2STRMAT(nv, nv, H, nv, qp->Hg, 0, 0);
	CVT_TRAN_MAT2STRMAT(nv, 1, g, nv, qp->Hg, nv, 0);
	CVT_TRAN_MAT2STRMAT(nv, ne, A, nv, qp->A, 0, 0);
	CVT_MAT2STRMAT(nv, ng, C, nv, qp->Ct, 0, 0);
	CVT_VEC2STRVEC(nv, g, qp->g, 0);
	CVT_VEC2STRVEC(ne, b, qp->b, 0);
	CVT_VEC2STRVEC(nb, d_lb, qp->d, 0);
	CVT_VEC2STRVEC(nb, d_ub, qp->d, nb+ng);
	CVT_VEC2STRVEC(ng, d_lg, qp->d, nb);
	CVT_VEC2STRVEC(ng, d_ug, qp->d, 2*nb+ng);
	for(ii=0; ii<nb; ii++) qp->idxb[ii] = idxb[ii];

	return;

	}



void CVT_DENSE_QP_TO_ROWMAJ(struct DENSE_QP *qp, REAL *H, REAL *g, REAL *A, REAL *b, int *idxb, REAL *d_lb, REAL *d_ub, REAL *C, REAL *d_lg, REAL *d_ug)
	{

	int ii;

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;

	CVT_TRAN_STRMAT2MAT(nv, nv, qp->Hg, 0, 0, H, nv);
	CVT_TRAN_STRMAT2MAT(ne, nv, qp->A, 0, 0, A, nv);
	CVT_STRMAT2MAT(nv, ng, qp->Ct, 0, 0, C, nv);
	CVT_STRVEC2VEC(nv, qp->g, 0, g);
	CVT_STRVEC2VEC(ne, qp->b, 0, b);
	CVT_STRVEC2VEC(nb, qp->d, 0, d_lb);
	CVT_STRVEC2VEC(nb, qp->d, nb+ng, d_ub);
	CVT_STRVEC2VEC(ng, qp->d, nb, d_lg);
	CVT_STRVEC2VEC(ng, qp->d, 2*nb+ng, d_ug);
	for(ii=0; ii<nb; ii++) idxb[ii] = qp->idxb[ii];

	return;

	}



void CVT_LIBSTR_TO_DENSE_QP(struct STRMAT *H, struct STRMAT *A, struct STRMAT *C, struct STRVEC *g, struct STRVEC *b, struct STRVEC *d_lb, struct STRVEC *d_ub, struct STRVEC *d_lg, struct STRVEC *d_ug, int *idxb, struct DENSE_QP *qp)
	{

	int ii;

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;

	GECP_LIBSTR(nv, nv, H, 0, 0, qp->Hg, 0, 0);
	ROWIN_LIBSTR(nv, 1.0, g, 0, qp->Hg, nv, 0);
	GECP_LIBSTR(ne, nv, A, 0, 0, qp->A, 0, 0);
	GETR_LIBSTR(ng, nv, C, 0, 0, qp->Ct, 0, 0);
	VECCP_LIBSTR(nv, g, 0, qp->g, 0);
	VECCP_LIBSTR(ne, b, 0, qp->b, 0);
	VECCP_LIBSTR(nb, d_lb, 0, qp->d, 0);
	VECCP_LIBSTR(nb, d_ub, 0, qp->d, nb+ng);
	VECCP_LIBSTR(ng, d_lg, 0, qp->d, nb);
	VECCP_LIBSTR(ng, d_ug, 0, qp->d, 2*nb+ng);
	for(ii=0; ii<nb; ii++) qp->idxb[ii] = idxb[ii];

	return;

	}



void CVT_DENSE_QP_TO_LIBSTR(struct DENSE_QP *qp, struct STRMAT *H, struct STRMAT *A, struct STRMAT *C, struct STRVEC *g, struct STRVEC *b, struct STRVEC *d_lb, struct STRVEC *d_ub, struct STRVEC *d_lg, struct STRVEC *d_ug, int *idxb)
	{

	int ii;

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;

	GECP_LIBSTR(nv, nv, qp->Hg, 0, 0, H, 0, 0);
	GECP_LIBSTR(ne, nv, qp->A, 0, 0, A, 0, 0);
	GETR_LIBSTR(nv, ng, qp->Ct, 0, 0, C, 0, 0);
	VECCP_LIBSTR(nv, qp->g, 0, g, 0);
	VECCP_LIBSTR(ne, qp->b, 0, b, 0);
	VECCP_LIBSTR(nb, qp->d, 0, d_lb, 0);
	VECCP_LIBSTR(nb, qp->d, nb+ng, d_ub, 0);
	VECCP_LIBSTR(ng, qp->d, nb, d_lg, 0);
	VECCP_LIBSTR(ng, qp->d, 2*nb+ng, d_ug, 0);
	for(ii=0; ii<nb; ii++) idxb[ii] = qp->idxb[ii];

	return;

	}


