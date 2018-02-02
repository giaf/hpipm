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



int MEMSIZE_DENSE_QP_SOL(struct DENSE_QP_DIM *dim)
	{

	int nv = dim->nv;
	int ne = dim->ne;
	int nb = dim->nb;
	int ng = dim->ng;
	int ns = dim->ns;

	int size = 0;

	size += 4*sizeof(struct STRVEC); // v pi lam t

	size += 1*SIZE_STRVEC(nv+2*ns); // ux
	size += 1*SIZE_STRVEC(ne); // pi
	size += 2*SIZE_STRVEC(2*nb+2*ng+2*ns); // lam t

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 64; // align to typical cache line size
	
	return size;

	}



void CREATE_DENSE_QP_SOL(struct DENSE_QP_DIM *dim, struct DENSE_QP_SOL *qp_sol, void *mem)
	{

	int nv = dim->nv;
	int ne = dim->ne;
	int nb = dim->nb;
	int ng = dim->ng;
	int ns = dim->ns;


	// vector struct stuff
	struct STRVEC *sv_ptr = (struct STRVEC *) mem;

	qp_sol->v = sv_ptr;
	sv_ptr += 1;
	qp_sol->pi = sv_ptr;
	sv_ptr += 1;
	qp_sol->lam = sv_ptr;
	sv_ptr += 1;
	qp_sol->t = sv_ptr;
	sv_ptr += 1;


	// align to typical cache line size
	long long l_ptr = (long long) sv_ptr;
	l_ptr = (l_ptr+63)/64*64;


	// double stuff
	char *c_ptr;
	c_ptr = (char *) l_ptr;

	char *tmp_ptr;

	// v
	CREATE_STRVEC(nv+2*ns, qp_sol->v, c_ptr);
	c_ptr += qp_sol->v->memsize;
	// pi
	CREATE_STRVEC(ne, qp_sol->pi, c_ptr);
	c_ptr += qp_sol->pi->memsize;
	// lam
	CREATE_STRVEC(2*nb+2*ng+2*ns, qp_sol->lam, c_ptr);
	c_ptr += qp_sol->lam->memsize;
	// t
	CREATE_STRVEC(2*nb+2*ng+2*ns, qp_sol->t, c_ptr);
	c_ptr += qp_sol->t->memsize;


	qp_sol->dim = dim;

	qp_sol->memsize = MEMSIZE_DENSE_QP_SOL(dim);


#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + qp_sol->memsize)
		{
		printf("\nCreate_ocp_qp_sol: outsize memory bounds!\n\n");
		exit(1);
		}
#endif


	return;

	}



void CVT_DENSE_QP_SOL_TO_COLMAJ(struct DENSE_QP_SOL *qp_sol, REAL *v, REAL *ls, REAL *us, REAL *pi, REAL *lam_lb, REAL *lam_ub, REAL *lam_lg, REAL *lam_ug, REAL *lam_ls, REAL *lam_us)
	{

	int nv = qp_sol->dim->nv;
	int ne = qp_sol->dim->ne;
	int nb = qp_sol->dim->nb;
	int ng = qp_sol->dim->ng;
	int ns = qp_sol->dim->ns;

	CVT_STRVEC2VEC(nv, qp_sol->v, 0, v);
	if(ne>0)
		{
		CVT_STRVEC2VEC(ne, qp_sol->pi, 0, pi);
		}
	if(nb>0)
		{
		CVT_STRVEC2VEC(nb, qp_sol->lam, 0, lam_lb);
		CVT_STRVEC2VEC(nb, qp_sol->lam, nb+ng, lam_ub);
		}
	if(ng>0)
		{
		CVT_STRVEC2VEC(ng, qp_sol->lam, nb, lam_lg);
		CVT_STRVEC2VEC(ng, qp_sol->lam, 2*nb+ng, lam_ug);
		}
	if(ns>0)
		{
		CVT_STRVEC2VEC(ns, qp_sol->v, nv, ls);
		CVT_STRVEC2VEC(ns, qp_sol->v, nv+ns, us);
		CVT_STRVEC2VEC(ns, qp_sol->lam, 2*nb+2*ng, lam_ls);
		CVT_STRVEC2VEC(ns, qp_sol->lam, 2*nb+2*ng+ns, lam_us);
		}

	return;

	}



void CVT_DENSE_QP_SOL_TO_ROWMAJ(struct DENSE_QP_SOL *qp_sol, REAL *v, REAL *ls, REAL *us, REAL *pi, REAL *lam_lb, REAL *lam_ub, REAL *lam_lg, REAL *lam_ug, REAL *lam_ls, REAL *lam_us)
	{

	int nv = qp_sol->dim->nv;
	int ne = qp_sol->dim->ne;
	int nb = qp_sol->dim->nb;
	int ng = qp_sol->dim->ng;
	int ns = qp_sol->dim->ns;

	CVT_STRVEC2VEC(nv, qp_sol->v, 0, v);
	if(ne>0)
		{
		CVT_STRVEC2VEC(ne, qp_sol->pi, 0, pi);
		}
	if(nb>0)
		{
		CVT_STRVEC2VEC(nb, qp_sol->lam, 0, lam_lb);
		CVT_STRVEC2VEC(nb, qp_sol->lam, nb+ng, lam_ub);
		}
	if(ng>0)
		{
		CVT_STRVEC2VEC(ng, qp_sol->lam, nb, lam_lg);
		CVT_STRVEC2VEC(ng, qp_sol->lam, 2*nb+ng, lam_ug);
		}
	if(ns>0)
		{
		CVT_STRVEC2VEC(ns, qp_sol->v, nv, ls);
		CVT_STRVEC2VEC(ns, qp_sol->v, nv+ns, us);
		CVT_STRVEC2VEC(ns, qp_sol->lam, 2*nb+2*ng, lam_ls);
		CVT_STRVEC2VEC(ns, qp_sol->lam, 2*nb+2*ng+ns, lam_us);
		}

	return;

	}



void CVT_DENSE_QP_SOL_TO_LIBSTR(struct DENSE_QP_SOL *qp_sol, struct STRVEC *v, struct STRVEC *ls, struct STRVEC *us, struct STRVEC *pi, struct STRVEC *lam_lb, struct STRVEC *lam_ub, struct STRVEC *lam_lg, struct STRVEC *lam_ug, struct STRVEC *lam_ls, struct STRVEC *lam_us)
	{

	int nv = qp_sol->dim->nv;
	int ne = qp_sol->dim->ne;
	int nb = qp_sol->dim->nb;
	int ng = qp_sol->dim->ng;
	int ns = qp_sol->dim->ns;

	VECCP_LIBSTR(nv, qp_sol->v, 0, v, 0);
	if(ne>0)
		{
		VECCP_LIBSTR(ne, qp_sol->pi, 0, pi, 0);
		}
	if(nb>0)
		{
		VECCP_LIBSTR(nb, qp_sol->lam, 0, lam_lb, 0);
		VECCP_LIBSTR(nb, qp_sol->lam, nb+ng, lam_ub, 0);
		}
	if(ng>0)
		{
		VECCP_LIBSTR(ng, qp_sol->lam, nb, lam_lg, 0);
		VECCP_LIBSTR(ng, qp_sol->lam, 2*nb+ng, lam_ug, 0);
		}
	if(ns>0)
		{
		VECCP_LIBSTR(ns, qp_sol->v, nv, ls, 0);
		VECCP_LIBSTR(ns, qp_sol->v, nv+ns, us, 0);
		VECCP_LIBSTR(ns, qp_sol->lam, 2*nb+2*ng, lam_ls, 0);
		VECCP_LIBSTR(ns, qp_sol->lam, 2*nb+2*ng+ns, lam_us, 0);
		}

	return;

	}
