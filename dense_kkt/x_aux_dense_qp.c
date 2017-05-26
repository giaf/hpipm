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



int SIZE_DENSE_QP(int nv, int ne, int nb, int ng)
	{

	int size = 0;

	size += 6*sizeof(struct STRVEC); // g be lb ub lg ug
	size += 3*sizeof(struct STRMAT); // H A Ct

	size += 1*SIZE_STRVEC(nv); // g
	size += 1*SIZE_STRVEC(ne); // be
	size += 2*SIZE_STRVEC(nb); // lb ub
	size += 2*SIZE_STRVEC(ng); // lg ug
	size += 1*nb*sizeof(int); // idxb

	size += 1*SIZE_STRMAT(nv, nv); // H
	size += 1*SIZE_STRMAT(ne, nv); // A
	size += 1*SIZE_STRMAT(nv, ng); // Ct

	size = (size+63)/64*64; // make multiple of typical cache line size
	
	return size;

	}



void CREATE_DENSE_QP(int nv, int ne, int nb, int ng, struct DENSE_QP *qp, void *memory)
	{

	// problem size
	qp->nv = nv;
	qp->ne = ne;
	qp->nb = nb;
	qp->ng = ng;


	// matrix struct stuff
	struct STRMAT *sm_ptr = (struct STRMAT *) memory;

	// H
	qp->H = sm_ptr;
	sm_ptr += 1;

	// A
	qp->A = sm_ptr;
	sm_ptr += 1;

	// Ct
	qp->Ct = sm_ptr;
	sm_ptr += 1;


	// vector struct stuff
	struct STRVEC *sv_ptr = (struct STRVEC *) sm_ptr;

	// g
	qp->g = sv_ptr;
	sv_ptr += 1;

	// be
	qp->be = sv_ptr;
	sv_ptr += 1;

	// lb
	qp->lb = sv_ptr;
	sv_ptr += 1;

	// ub
	qp->ub = sv_ptr;
	sv_ptr += 1;

	// lg
	qp->lg = sv_ptr;
	sv_ptr += 1;

	// ug
	qp->ug = sv_ptr;
	sv_ptr += 1;


	// int stuff
	int *i_ptr;
	i_ptr = (int *) sv_ptr;

	// idxb
	qp->idxb = i_ptr;
	i_ptr += nb;


	// align to typical cache line size
	long long l_ptr = (long long) i_ptr;
	l_ptr = (l_ptr+63)/64*64;


	// double stuff
	void *v_ptr;
	v_ptr = (void *) l_ptr;

	// H
	CREATE_STRMAT(nv, nv, qp->H, v_ptr);
	v_ptr += qp->H->memory_size;

	// A
	CREATE_STRMAT(ne, nv, qp->A, v_ptr);
	v_ptr += qp->A->memory_size;

	// H
	CREATE_STRMAT(nv, ng, qp->Ct, v_ptr);
	v_ptr += qp->Ct->memory_size;

	// q
	CREATE_STRVEC(nv, qp->g, v_ptr);
	v_ptr += qp->g->memory_size;

	// b
	CREATE_STRVEC(ne, qp->be, v_ptr);
	v_ptr += qp->be->memory_size;

	// lb
	CREATE_STRVEC(nb, qp->lb, v_ptr);
	v_ptr += qp->lb->memory_size;

	// ub
	CREATE_STRVEC(nb, qp->ub, v_ptr);
	v_ptr += qp->ub->memory_size;

	// lg
	CREATE_STRVEC(ng, qp->lg, v_ptr);
	v_ptr += qp->lg->memory_size;

	// ug
	CREATE_STRVEC(ng, qp->ug, v_ptr);
	v_ptr += qp->ug->memory_size;

	// memory size
	qp->mem_size = SIZE_DENSE_QP(nv, ne, nb, ng);

	return;

	}



void INIT_DENSE_QP(struct STRMAT *H, struct STRMAT *A, struct STRMAT *C, struct STRVEC *g, struct STRVEC *be, struct STRVEC *lb, struct STRVEC *ub, struct STRVEC *lg, struct STRVEC *ug, int *idxb, struct DENSE_QP *qp)
	{

	int ii;

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;

	GECP_LIBSTR(nv, nv, H, 0, 0, qp->H, 0, 0);
	GECP_LIBSTR(ne, nv, A, 0, 0, qp->A, 0, 0);
	GECP_LIBSTR(ng, nv, C, 0, 0, qp->Ct, 0, 0);
	VECCP_LIBSTR(nv, g, 0, qp->g, 0);
	VECCP_LIBSTR(ne, be, 0, qp->be, 0);
	VECCP_LIBSTR(nb, lb, 0, qp->lb, 0);
	VECCP_LIBSTR(nb, ub, 0, qp->ub, 0);
	VECCP_LIBSTR(ng, lg, 0, qp->lg, 0);
	VECCP_LIBSTR(ng, ug, 0, qp->ug, 0);
	for(ii=0; ii<nb; ii++) qp->idxb[ii] = idxb[ii];

	return;

	}



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


