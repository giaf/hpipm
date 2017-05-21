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



int SIZE_DENSE_QP_DIM(struct DENSE_QP_DIM *qp_dim)
	{

	int size = 0;

	return size;
	}



int SIZE_DENSE_QP_VEC(struct DENSE_QP_DIM *qp_dim)
	{

	// extract problem dimension
	int nv = qp_dim->nv;
	int ne = qp_dim->ne;
	int nb = qp_dim->nb;
	int nc = qp_dim->nc;

	int size = 0;

	size += 6*sizeof(struct STRVEC);

	size += 1*SIZE_STRVEC(nv); // g
	size += 1*SIZE_STRVEC(ne); // be
	size += 2*SIZE_STRVEC(nb); // lb ub
	size += 2*SIZE_STRVEC(nc); // lc uc
	size += 1*nb*sizeof(int); // idxb

	size = (size+63)/64*64; // make multiple of typical cache line size

	return size;
	}



int SIZE_DENSE_QP_MAT(struct DENSE_QP_DIM *qp_dim)
	{

	// extract problem dimension
	int nv = qp_dim->nv;
	int ne = qp_dim->ne;
	int nb = qp_dim->nb;
	int nc = qp_dim->nc;

	int size = 0;

	size += 3*sizeof(struct STRMAT);

	size += 1*SIZE_STRMAT(nv, nv); // H
	size += 1*SIZE_STRMAT(ne, nv); // A
	size += 1*SIZE_STRMAT(nv, nc); // Ct

	size = (size+63)/64*64; // make multiple of typical cache line size

	return size;
	}



int SIZE_DENSE_QP(struct DENSE_QP_DIM *qp_dim)
	{

	int nv = qp_dim->nv;
	int ne = qp_dim->ne;
	int nb = qp_dim->nb;
	int nc = qp_dim->nc;

	int ii;

	int size = 0;

	size += sizeof(struct DENSE_QP_DIM);
	size += SIZE_DENSE_QP_DIM(qp_dim);
	size += sizeof(struct DENSE_QP_VEC);
	size += SIZE_DENSE_QP_VEC(qp_dim);
	size += sizeof(struct DENSE_QP_MAT);
	size += SIZE_DENSE_QP_MAT(qp_dim);
	size += 1*sizeof(void (*)(void));
	
	return size;

	}



void CREATE_DENSE_QP_DIM(struct DENSE_QP_DIM *qp_dim, struct DENSE_QP_DIM *qp_dim_out, void *memory)
	{

	// memory size
	qp_dim_out->mem_size = SIZE_DENSE_QP_DIM(qp_dim);

	return;

	}



void CREATE_DENSE_QP_VEC(struct DENSE_QP_DIM *qp_dim, struct DENSE_QP_VEC *qp_vec, void *memory)
	{

	int nv = qp_dim->nv;
	int ne = qp_dim->ne;
	int nb = qp_dim->nb;
	int nc = qp_dim->nc;

	// vector struct stuff
	struct STRVEC *sv_ptr = (struct STRVEC *) memory;

	// g
	qp_vec->g = sv_ptr;
	sv_ptr += 1;

	// be
	qp_vec->be = sv_ptr;
	sv_ptr += 1;

	// lb
	qp_vec->lb = sv_ptr;
	sv_ptr += 1;

	// ub
	qp_vec->ub = sv_ptr;
	sv_ptr += 1;

	// lc
	qp_vec->lc = sv_ptr;
	sv_ptr += 1;

	// uc
	qp_vec->uc = sv_ptr;
	sv_ptr += 1;


	// int stuff
	int *i_ptr;
	i_ptr = (int *) sv_ptr;

	// idxb
	qp_vec->idxb = i_ptr;
	i_ptr += nb;


	// align to typical cache line size
	long long l_ptr = (long long) i_ptr;
	l_ptr = (l_ptr+63)/64*64;


	// double stuff
	void *v_ptr;
	v_ptr = (void *) l_ptr;

	// q
	CREATE_STRVEC(nv, qp_vec->g, v_ptr);
	v_ptr += qp_vec->g->memory_size;

	// b
	CREATE_STRVEC(ne, qp_vec->be, v_ptr);
	v_ptr += qp_vec->be->memory_size;

	// lb
	CREATE_STRVEC(nb, qp_vec->lb, v_ptr);
	v_ptr += qp_vec->lb->memory_size;

	// ub
	CREATE_STRVEC(nb, qp_vec->ub, v_ptr);
	v_ptr += qp_vec->ub->memory_size;

	// lg
	CREATE_STRVEC(nc, qp_vec->lc, v_ptr);
	v_ptr += qp_vec->lc->memory_size;

	// ug
	CREATE_STRVEC(nc, qp_vec->uc, v_ptr);
	v_ptr += qp_vec->uc->memory_size;

	// memory size
	qp_vec->mem_size = SIZE_DENSE_QP_VEC(qp_dim);

	return;

	}



void CREATE_DENSE_QP_MAT(struct DENSE_QP_DIM *qp_dim, struct DENSE_QP_MAT *qp_mat, void *memory)
	{

	int nv = qp_dim->nv;
	int ne = qp_dim->ne;
	int nb = qp_dim->nb;
	int nc = qp_dim->nc;

	// vector struct stuff
	struct STRMAT *sm_ptr = (struct STRMAT *) memory;

	// H
	qp_mat->H = sm_ptr;
	sm_ptr += 1;

	// A
	qp_mat->A = sm_ptr;
	sm_ptr += 1;

	// Ct
	qp_mat->Ct = sm_ptr;
	sm_ptr += 1;



	// align to typical cache line size
	long long l_ptr = (long long) sm_ptr;
	l_ptr = (l_ptr+63)/64*64;


	// double stuff
	void *v_ptr;
	v_ptr = (void *) l_ptr;

	// H
	CREATE_STRMAT(nv, nv, qp_mat->H, v_ptr);
	v_ptr += qp_mat->H->memory_size;

	// A
	CREATE_STRMAT(ne, nv, qp_mat->A, v_ptr);
	v_ptr += qp_mat->A->memory_size;

	// H
	CREATE_STRMAT(nv, nc, qp_mat->Ct, v_ptr);
	v_ptr += qp_mat->Ct->memory_size;

	// memory size
	qp_mat->mem_size = SIZE_DENSE_QP_MAT(qp_dim);

	return;

	}



void CREATE_DENSE_QP(struct DENSE_QP_DIM *qp_dim, struct DENSE_QP *qp, void *memory)
	{

	int nv = qp_dim->nv;
	int ne = qp_dim->ne;
	int nb = qp_dim->nb;
	int nc = qp_dim->nc;

	qp->dim = memory;
	memory += sizeof(struct DENSE_QP_DIM);

	qp->vec = memory;
	memory += sizeof(struct DENSE_QP_VEC);

	qp->mat = memory;
	memory += sizeof(struct DENSE_QP_MAT);

	CREATE_DENSE_QP_DIM(qp_dim, qp->dim, memory);
	memory += qp->dim->mem_size;

	CREATE_DENSE_QP_VEC(qp_dim, qp->vec, memory);
	memory += qp->vec->mem_size;

	CREATE_DENSE_QP_MAT(qp_dim, qp->mat, memory);
	memory += qp->mat->mem_size;

	return;
	}



void INIT_DENSE_QP_DIM(int nv, int ne, int nb, int nc, struct DENSE_QP_DIM *qp_dim_out)
	{

	qp_dim_out->nv = nv;
	qp_dim_out->ne = ne;
	qp_dim_out->nb = nb;
	qp_dim_out->nc = nc;

	return;

	}



void INIT_DENSE_QP_VEC(struct DENSE_QP_DIM *qp_dim, struct STRVEC *g, struct STRVEC *be, struct STRVEC *lb, struct STRVEC *ub, struct STRVEC *lc, struct STRVEC *uc, struct DENSE_QP_VEC *qp_vec)
	{

	int nv = qp_dim->nv;
	int ne = qp_dim->ne;
	int nb = qp_dim->nb;
	int nc = qp_dim->nc;

	dveccp_libstr(nv, g, 0, qp_vec->g, 0);
	dveccp_libstr(ne, be, 0, qp_vec->be, 0);
	dveccp_libstr(nb, lb, 0, qp_vec->lb, 0);
	dveccp_libstr(nb, ub, 0, qp_vec->ub, 0);
	dveccp_libstr(nc, lc, 0, qp_vec->lc, 0);
	dveccp_libstr(nc, uc, 0, qp_vec->uc, 0);

	return;

	}



void INIT_DENSE_QP_MAT(struct DENSE_QP_DIM *qp_dim, struct STRMAT *H, struct STRMAT *A, struct STRMAT *Ct, struct DENSE_QP_MAT *qp_mat)
	{

	int nv = qp_dim->nv;
	int ne = qp_dim->ne;
	int nb = qp_dim->nb;
	int nc = qp_dim->nc;

	dgecp_libstr(nv, nv, H, 0, 0, qp_mat->H, 0, 0);
	dgecp_libstr(ne, nv, A, 0, 0, qp_mat->A, 0, 0);
	dgecp_libstr(nv, nc, Ct, 0, 0, qp_mat->Ct, 0, 0);

	return;

	}



void INIT_DENSE_QP(struct DENSE_QP_DIM *qp_dim, struct STRVEC *g, struct STRVEC *be, struct STRVEC *lb, struct STRVEC *ub, struct STRVEC *lc, struct STRVEC *uc, struct STRMAT *H, struct STRMAT *A, struct STRMAT *Ct, struct DENSE_QP *qp)
	{

	INIT_DENSE_QP_DIM(qp_dim->nv, qp_dim->ne, qp_dim->nb, qp_dim->nc, qp->dim);
	INIT_DENSE_QP_VEC(qp->dim, g, be, lb, ub, lc, uc, qp->vec);
	INIT_DENSE_QP_MAT(qp->dim, H, A, Ct, qp->mat);

	return;

	}



void CAST_DENSE_QP_DIM(int nv, int ne, int nb, int nc, struct DENSE_QP_DIM *qp_dim_out)
	{

	qp_dim_out->nv = nv;
	qp_dim_out->ne = ne;
	qp_dim_out->nb = nb;
	qp_dim_out->nc = nc;

	return;

	}



#if 0
void CREATE_DENSE_QP(struct DENSE_QP_DIM *qp_dim, struct DENSE_QP *qp_out, void *memory)
	{

	int ii;


	// int stuff
	int *i_ptr;
	i_ptr = (int *) memory;

	// idxb
	str_out->idxb = i_ptr;
	i_ptr += nb;


	// align to typical cache line size
	long long l_ptr = (long long) i_ptr;
	l_ptr = (l_ptr+63)/64*64;


	// double stuff
	void *v_ptr;
	v_ptr = (void *) l_ptr;

	// Q
	CREATE_STRMAT(nv, nv, &(str_out->sQ), v_ptr);
	v_ptr += str_out->sQ.memory_size;

	// A
	CREATE_STRMAT(ne, nv, &(str_out->sA), v_ptr);
	v_ptr += str_out->sA.memory_size;

	// Ct
	CREATE_STRMAT(nv, ng, &(str_out->sCt), v_ptr);
	v_ptr += str_out->sCt.memory_size;

	// q
	CREATE_STRVEC(nv, &(str_out->sq), v_ptr);
	v_ptr += str_out->sq.memory_size;

	// b
	CREATE_STRVEC(ne, &(str_out->sb), v_ptr);
	v_ptr += str_out->sb.memory_size;

	// lb
	CREATE_STRVEC(nb, &(str_out->slb), v_ptr);
	v_ptr += str_out->slb.memory_size;

	// ub
	CREATE_STRVEC(nb, &(str_out->sub), v_ptr);
	v_ptr += str_out->sub.memory_size;

	// lg
	CREATE_STRVEC(ng, &(str_out->slg), v_ptr);
	v_ptr += str_out->slg.memory_size;

	// ug
	CREATE_STRVEC(ng, &(str_out->sug), v_ptr);
	v_ptr += str_out->sug.memory_size;

	return;

	}
#endif



#if 0
void COPY_DENSE_QP(struct DENSE_QP *str_in, struct DENSE_QP *str_out)
	{

	int ii;

#if defined(RUNTIME_CHECKS)
	if(str_out->nv != str_in->nv)
		{
		printf("\nError : d_copy_dense_qp : str_out->nv != str_out->nv : %d != %d\n\n", str_out->nv, str_in->nv);
		exit(1);
		}
	if(str_out->ne != str_in->ne)
		{
		printf("\nError : d_copy_dense_qp : str_out->ne != str_out->ne : %d != %d\n\n", str_out->ne, str_in->ne);
		exit(1);
		}
	if(str_out->nb != str_in->nb)
		{
		printf("\nError : d_copy_dense_qp : str_out->nb != str_out->nb : %d != %d\n\n", str_out->nb, str_in->nb);
		exit(1);
		}
	if(str_out->ng != str_in->ng)
		{
		printf("\nError : d_copy_dense_qp : str_out->ng != str_out->ng : %d != %d\n\n", str_out->ng, str_in->ng);
		exit(1);
		}
#endif

	for(ii=0; ii<str_in->nb; ii++) str_out->idxb[ii] = str_in->idxb[ii];
	GECP_LIBSTR(str_in->nv, str_in->nv, &(str_out->sQ), 0, 0, &(str_in->sQ), 0, 0);
	VECCP_LIBSTR(str_in->nv, &(str_out->sq), 0, &(str_in->sq), 0);
	GECP_LIBSTR(str_in->ne, str_in->nv, &(str_out->sA), 0, 0, &(str_in->sA), 0, 0);
	VECCP_LIBSTR(str_in->ne, &(str_out->sb), 0, &(str_in->sb), 0);
	GECP_LIBSTR(str_in->nv, str_in->ng, &(str_out->sCt), 0, 0, &(str_in->sCt), 0, 0);
	VECCP_LIBSTR(str_in->nb, &(str_out->slb), 0, &(str_in->slb), 0);
	VECCP_LIBSTR(str_in->nb, &(str_out->sub), 0, &(str_in->sub), 0);
	VECCP_LIBSTR(str_in->ng, &(str_out->slg), 0, &(str_in->slg), 0);
	VECCP_LIBSTR(str_in->ng, &(str_out->sug), 0, &(str_in->sug), 0);

	return;

	}
#endif


