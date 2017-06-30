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



int MEMSIZE_DENSE_QP_SOL(int nv, int ne, int nb, int ng)
	{

	int size = 0;

	size += (10)*sizeof(struct STRVEC); // v pi lam_lb lam_ub lam_lg lam_ug t_lb t_ub t_lg t_ug

	size += 1*SIZE_STRVEC(nv); // ux
	size += 1*SIZE_STRVEC(ne); // pi
	size += 8*SIZE_STRVEC(2*nb+2*ng); // lam_lb lam_ub lam_lg lam_ug t_lb t_ub t_lg t_ug

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 64; // align to typical cache line size
	
	return size;

	}



void CREATE_DENSE_QP_SOL(int nv, int ne, int nb, int ng, struct DENSE_QP_SOL *qp_sol, void *memory)
	{

	qp_sol->memsize = MEMSIZE_DENSE_QP_SOL(nv, ne, nb, ng);


	// vector struct stuff
	struct STRVEC *sv_ptr = (struct STRVEC *) memory;

	qp_sol->v = sv_ptr;
	sv_ptr += 1;
	qp_sol->pi = sv_ptr;
	sv_ptr += 1;
	qp_sol->lam_lb = sv_ptr;
	sv_ptr += 1;
	qp_sol->lam_ub = sv_ptr;
	sv_ptr += 1;
	qp_sol->lam_lg = sv_ptr;
	sv_ptr += 1;
	qp_sol->lam_ug = sv_ptr;
	sv_ptr += 1;
	qp_sol->t_lb = sv_ptr;
	sv_ptr += 1;
	qp_sol->t_ub = sv_ptr;
	sv_ptr += 1;
	qp_sol->t_lg = sv_ptr;
	sv_ptr += 1;
	qp_sol->t_ug = sv_ptr;
	sv_ptr += 1;


	// align to typical cache line size
	long long l_ptr = (long long) sv_ptr;
	l_ptr = (l_ptr+63)/64*64;


	// double stuff
	void *v_ptr;
	v_ptr = (void *) l_ptr;

	void *tmp_ptr;

	// v
	CREATE_STRVEC(nv, qp_sol->v, v_ptr);
	v_ptr += qp_sol->v->memory_size;
	// pi
	CREATE_STRVEC(ne, qp_sol->pi, v_ptr);
	v_ptr += qp_sol->pi->memory_size;
	// lam
	tmp_ptr = v_ptr;
	v_ptr += SIZE_STRVEC(2*nb+2*ng);
	// lam_lb
	CREATE_STRVEC(nb, qp_sol->lam_lb, tmp_ptr);
	tmp_ptr += (nb)*sizeof(REAL);
	// lam_lg
	CREATE_STRVEC(ng, qp_sol->lam_lg, tmp_ptr);
	tmp_ptr += (ng)*sizeof(REAL);
	// lam_ub
	CREATE_STRVEC(nb, qp_sol->lam_ub, tmp_ptr);
	tmp_ptr += (nb)*sizeof(REAL);
	// lam_ug
	CREATE_STRVEC(ng, qp_sol->lam_ug, tmp_ptr);
	tmp_ptr += (ng)*sizeof(REAL);
	// t_lb
	CREATE_STRVEC(nb, qp_sol->t_lb, tmp_ptr);
	tmp_ptr += (nb)*sizeof(REAL);
	// t_lg
	CREATE_STRVEC(ng, qp_sol->t_lg, tmp_ptr);
	tmp_ptr += (ng)*sizeof(REAL);
	// t_ub
	CREATE_STRVEC(nb, qp_sol->t_ub, tmp_ptr);
	tmp_ptr += (nb)*sizeof(REAL);
	// t_ug
	CREATE_STRVEC(ng, qp_sol->t_ug, tmp_ptr);
	tmp_ptr += (ng)*sizeof(REAL);

	return;

	}



void CVT_DENSE_QP_SOL_TO_COLMAJ(struct DENSE_QP *qp, struct DENSE_QP_SOL *qp_sol, REAL *v, REAL *pi, REAL *lam_lb, REAL *lam_ub, REAL *lam_lg, REAL *lam_ug)
	{

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;

	CVT_STRVEC2VEC(nv, qp_sol->v, 0, v);
	CVT_STRVEC2VEC(ne, qp_sol->pi, 0, pi);
	CVT_STRVEC2VEC(nb, qp_sol->lam_lb, 0, lam_lb);
	CVT_STRVEC2VEC(nb, qp_sol->lam_ub, 0, lam_ub);
	CVT_STRVEC2VEC(ng, qp_sol->lam_lg, 0, lam_lg);
	CVT_STRVEC2VEC(ng, qp_sol->lam_ug, 0, lam_ug);

	return;

	}



void CVT_DENSE_QP_SOL_TO_ROWMAJ(struct DENSE_QP *qp, struct DENSE_QP_SOL *qp_sol, REAL *v, REAL *pi, REAL *lam_lb, REAL *lam_ub, REAL *lam_lg, REAL *lam_ug)
	{

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;

	CVT_STRVEC2VEC(nv, qp_sol->v, 0, v);
	CVT_STRVEC2VEC(ne, qp_sol->pi, 0, pi);
	CVT_STRVEC2VEC(nb, qp_sol->lam_lb, 0, lam_lb);
	CVT_STRVEC2VEC(nb, qp_sol->lam_ub, 0, lam_ub);
	CVT_STRVEC2VEC(ng, qp_sol->lam_lg, 0, lam_lg);
	CVT_STRVEC2VEC(ng, qp_sol->lam_ug, 0, lam_ug);

	return;

	}



void CVT_DENSE_QP_SOL_TO_LIBSTR(struct DENSE_QP *qp, struct DENSE_QP_SOL *qp_sol, struct STRVEC *v, struct STRVEC *pi, struct STRVEC *lam_lb, struct STRVEC *lam_ub, struct STRVEC *lam_lg, struct STRVEC *lam_ug)
	{

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;

	VECCP_LIBSTR(nv, qp_sol->v, 0, v, 0);
	VECCP_LIBSTR(ne, qp_sol->pi, 0, pi, 0);
	VECCP_LIBSTR(nb, qp_sol->lam_lb, 0, lam_lb, 0);
	VECCP_LIBSTR(nb, qp_sol->lam_ub, 0, lam_ub, 0);
	VECCP_LIBSTR(ng, qp_sol->lam_lg, 0, lam_lg, 0);
	VECCP_LIBSTR(ng, qp_sol->lam_ug, 0, lam_ug, 0);

	return;

	}
