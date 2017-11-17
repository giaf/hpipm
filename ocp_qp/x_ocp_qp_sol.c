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



int MEMSIZE_OCP_QP_SOL(struct OCP_QP_DIM *dim)
	{

	// extract dim
	int N = dim->N;
	int *nx = dim->nx;
	int *nu = dim->nu;
	int *nb = dim->nb;
	int *ng = dim->ng;
	int *ns = dim->ns;

	// loop index
	int ii;

	int nvt = 0;
	int net = 0;
	int nct = 0;
	for(ii=0; ii<N; ii++)
		{
		nvt += nu[ii]+nx[ii]+2*ns[ii];
		net += nx[ii+1];
		nct += nb[ii]+ng[ii]+ns[ii];
		}
	nvt += nu[ii]+nx[ii]+2*ns[ii];
	nct += nb[ii]+ng[ii]+ns[ii];

	int size = 0;

	size += 3*(N+1)*sizeof(struct STRVEC); // ux lam t
	size += 1*N*sizeof(struct STRVEC); // pi

	size += 1*SIZE_STRVEC(nvt); // ux
	size += 1*SIZE_STRVEC(net); // pi
	size += 2*SIZE_STRVEC(2*nct); // lam t

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 64; // align to typical cache line size

	return size;

	}



void CREATE_OCP_QP_SOL(struct OCP_QP_DIM *dim, struct OCP_QP_SOL *qp_sol, void *mem)
	{

	// extract dim
	int N = dim->N;
	int *nx = dim->nx;
	int *nu = dim->nu;
	int *nb = dim->nb;
	int *ng = dim->ng;
	int *ns = dim->ns;

	// loop index
	int ii;

	int nvt = 0;
	int net = 0;
	int nct = 0;
	for(ii=0; ii<N; ii++)
		{
		nvt += nu[ii]+nx[ii]+2*ns[ii];
		net += nx[ii+1];
		nct += nb[ii]+ng[ii]+ns[ii];
		}
	nvt += nu[ii]+nx[ii]+2*ns[ii];
	nct += nb[ii]+ng[ii]+ns[ii];


	// vector struct stuff
	struct STRVEC *sv_ptr = (struct STRVEC *) mem;

	qp_sol->ux = sv_ptr;
	sv_ptr += N+1;
	qp_sol->pi = sv_ptr;
	sv_ptr += N;
	qp_sol->lam = sv_ptr;
	sv_ptr += N+1;
	qp_sol->t = sv_ptr;
	sv_ptr += N+1;


	// align to typical cache line size
	long long l_ptr = (long long) sv_ptr;
	l_ptr = (l_ptr+63)/64*64;


	// double stuff
	char *c_ptr;
	c_ptr = (char *) l_ptr;

	char *tmp_ptr;

	// ux
	tmp_ptr = c_ptr;
	c_ptr += SIZE_STRVEC(nvt);
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(nu[ii]+nx[ii]+2*ns[ii], qp_sol->ux+ii, tmp_ptr);
		tmp_ptr += nu[ii]*sizeof(REAL); // u
		tmp_ptr += nx[ii]*sizeof(REAL); // x
		tmp_ptr += ns[ii]*sizeof(REAL); // s_ls
		tmp_ptr += ns[ii]*sizeof(REAL); // s_us
		}
	// pi
	tmp_ptr = c_ptr;
	c_ptr += SIZE_STRVEC(net);
	for(ii=0; ii<N; ii++)
		{
		CREATE_STRVEC(nx[ii+1], qp_sol->pi+ii, tmp_ptr);
		tmp_ptr += (nx[ii+1])*sizeof(REAL); // pi
		}
	// lam
	tmp_ptr = c_ptr;
	c_ptr += SIZE_STRVEC(2*nct);
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(2*nb[ii]+2*ng[ii]+2*ns[ii], qp_sol->lam+ii, tmp_ptr);
		tmp_ptr += nb[ii]*sizeof(REAL); // lb
		tmp_ptr += ng[ii]*sizeof(REAL); // lg
		tmp_ptr += nb[ii]*sizeof(REAL); // ub
		tmp_ptr += ng[ii]*sizeof(REAL); // ug
		tmp_ptr += ns[ii]*sizeof(REAL); // ls
		tmp_ptr += ns[ii]*sizeof(REAL); // us
		}
	// t
	tmp_ptr = c_ptr;
	c_ptr += SIZE_STRVEC(2*nct);
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(2*nb[ii]+2*ng[ii]+2*ns[ii], qp_sol->t+ii, tmp_ptr);
		tmp_ptr += nb[ii]*sizeof(REAL); // lb
		tmp_ptr += ng[ii]*sizeof(REAL); // lg
		tmp_ptr += nb[ii]*sizeof(REAL); // ub
		tmp_ptr += ng[ii]*sizeof(REAL); // ug
		tmp_ptr += ns[ii]*sizeof(REAL); // ls
		tmp_ptr += ns[ii]*sizeof(REAL); // us
		}

	qp_sol->dim = dim;

	qp_sol->memsize = MEMSIZE_OCP_QP_SOL(dim);


#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + qp_sol->memsize)
		{
		printf("\nCreate_ocp_qp_sol: outsize memory bounds!\n\n");
		exit(1);
		}
#endif


	return;

	}



void CVT_OCP_QP_SOL_TO_COLMAJ(struct OCP_QP_SOL *qp_sol, REAL **u, REAL **x, REAL **ls, REAL **us, REAL **pi, REAL **lam_lb, REAL **lam_ub, REAL **lam_lg, REAL **lam_ug, REAL **lam_ls, REAL **lam_us)
	{

	int N = qp_sol->dim->N;
	int *nx = qp_sol->dim->nx;
	int *nu = qp_sol->dim->nu;
	int *nb = qp_sol->dim->nb;
	int *ng = qp_sol->dim->ng;
	int *ns = qp_sol->dim->ns;

	int ii;

	for(ii=0; ii<N; ii++)
		{
		CVT_STRVEC2VEC(nx[ii+1], qp_sol->pi+ii, 0, pi[ii]);
		}

	for(ii=0; ii<=N; ii++)
		{
		CVT_STRVEC2VEC(nu[ii], qp_sol->ux+ii, 0, u[ii]);
		CVT_STRVEC2VEC(nx[ii], qp_sol->ux+ii, nu[ii], x[ii]);
		if(nb[ii]>0)
			{
			CVT_STRVEC2VEC(nb[ii], qp_sol->lam+ii, 0, lam_lb[ii]);
			CVT_STRVEC2VEC(nb[ii], qp_sol->lam+ii, nb[ii]+ng[ii], lam_ub[ii]);
			}
		if(ng[ii]>0)
			{
			CVT_STRVEC2VEC(ng[ii], qp_sol->lam+ii, nb[ii], lam_lg[ii]);
			CVT_STRVEC2VEC(ng[ii], qp_sol->lam+ii, 2*nb[ii]+ng[ii], lam_ug[ii]);
			}
		if(ns[ii]>0)
			{
			CVT_STRVEC2VEC(ns[ii], qp_sol->ux+ii, nu[ii]+nx[ii], ls[ii]);
			CVT_STRVEC2VEC(ns[ii], qp_sol->ux+ii, nu[ii]+nx[ii]+ns[ii], us[ii]);
			CVT_STRVEC2VEC(ns[ii], qp_sol->lam+ii, 2*nb[ii]+2*ng[ii], lam_ls[ii]);
			CVT_STRVEC2VEC(ns[ii], qp_sol->lam+ii, 2*nb[ii]+2*ng[ii]+ns[ii], lam_us[ii]);
			}
		}

	return;

	}



void CVT_COLMAJ_TO_OCP_QP_SOL(REAL **u, REAL **x, REAL **ls, REAL **us, REAL **pi, REAL **lam_lb, REAL **lam_ub, REAL **lam_lg, REAL **lam_ug, REAL **lam_ls, REAL **lam_us, struct OCP_QP_SOL *qp_sol)
	{

	int N = qp_sol->dim->N;
	int *nx = qp_sol->dim->nx;
	int *nu = qp_sol->dim->nu;
	int *nb = qp_sol->dim->nb;
	int *ng = qp_sol->dim->ng;
	int *ns = qp_sol->dim->ns;

	int ii;

	for(ii=0; ii<N; ii++)
		{
		CVT_VEC2STRVEC(nx[ii+1], pi[ii], qp_sol->pi+ii, 0);
		}

	for(ii=0; ii<=N; ii++)
		{
		CVT_VEC2STRVEC(nu[ii], u[ii], qp_sol->ux+ii, 0);
		CVT_VEC2STRVEC(nx[ii], x[ii], qp_sol->ux+ii, nu[ii]);
		if(nb[ii]>0)
			{
			CVT_VEC2STRVEC(nb[ii], lam_lb[ii], qp_sol->lam+ii, 0);
			CVT_VEC2STRVEC(nb[ii], lam_ub[ii], qp_sol->lam+ii, nb[ii]+ng[ii]);
			}
		if(ng[ii]>0)
			{
			CVT_VEC2STRVEC(ng[ii], lam_lg[ii], qp_sol->lam+ii, nb[ii]);
			CVT_VEC2STRVEC(ng[ii], lam_ug[ii], qp_sol->lam+ii, 2*nb[ii]+ng[ii]);
			}
		if(ns[ii]>0)
			{
			CVT_VEC2STRVEC(ns[ii], ls[ii], qp_sol->ux+ii, nu[ii]+nx[ii]);
			CVT_VEC2STRVEC(ns[ii], us[ii], qp_sol->ux+ii, nu[ii]+nx[ii]+ns[ii]);
			CVT_VEC2STRVEC(ns[ii], lam_ls[ii], qp_sol->lam+ii, 2*nb[ii]+2*ng[ii]);
			CVT_VEC2STRVEC(ns[ii], lam_us[ii], qp_sol->lam+ii, 2*nb[ii]+2*ng[ii]+ns[ii]);
			}
		}

	return;

	}



void CVT_OCP_QP_SOL_TO_ROWMAJ(struct OCP_QP_SOL *qp_sol, REAL **u, REAL **x, REAL **ls, REAL **us, REAL **pi, REAL **lam_lb, REAL **lam_ub, REAL **lam_lg, REAL **lam_ug, REAL **lam_ls, REAL **lam_us)
	{

	int N = qp_sol->dim->N;
	int *nx = qp_sol->dim->nx;
	int *nu = qp_sol->dim->nu;
	int *nb = qp_sol->dim->nb;
	int *ng = qp_sol->dim->ng;
	int *ns = qp_sol->dim->ns;

	int ii;

	for(ii=0; ii<N; ii++)
		{
		CVT_STRVEC2VEC(nx[ii+1], qp_sol->pi+ii, 0, pi[ii]);
		}

	for(ii=0; ii<=N; ii++)
		{
		CVT_STRVEC2VEC(nu[ii], qp_sol->ux+ii, 0, u[ii]);
		CVT_STRVEC2VEC(nx[ii], qp_sol->ux+ii, nu[ii], x[ii]);
		if(nb[ii]>0)
			{
			CVT_STRVEC2VEC(nb[ii], qp_sol->lam+ii, 0, lam_lb[ii]);
			CVT_STRVEC2VEC(nb[ii], qp_sol->lam+ii, nb[ii]+ng[ii], lam_ub[ii]);
			}
		if(ng[ii]>0)
			{
			CVT_STRVEC2VEC(ng[ii], qp_sol->lam+ii, nb[ii], lam_lg[ii]);
			CVT_STRVEC2VEC(ng[ii], qp_sol->lam+ii, 2*nb[ii]+ng[ii], lam_ug[ii]);
			}
		if(ns[ii]>0)
			{
			CVT_STRVEC2VEC(ns[ii], qp_sol->ux+ii, nu[ii]+nx[ii], ls[ii]);
			CVT_STRVEC2VEC(ns[ii], qp_sol->ux+ii, nu[ii]+nx[ii]+ns[ii], us[ii]);
			CVT_STRVEC2VEC(ns[ii], qp_sol->lam+ii, 2*nb[ii]+2*ng[ii], lam_ls[ii]);
			CVT_STRVEC2VEC(ns[ii], qp_sol->lam+ii, 2*nb[ii]+2*ng[ii]+ns[ii], lam_us[ii]);
			}
		}

	return;

	}



void CVT_OCP_QP_SOL_TO_LIBSTR(struct OCP_QP_SOL *qp_sol, struct STRVEC *u, struct STRVEC *ls, struct STRVEC *us, struct STRVEC *x, struct STRVEC *pi, struct STRVEC *lam_lb, struct STRVEC *lam_ub, struct STRVEC *lam_lg, struct STRVEC *lam_ug, struct STRVEC *lam_ls, struct STRVEC *lam_us)
	{

	int N = qp_sol->dim->N;
	int *nx = qp_sol->dim->nx;
	int *nu = qp_sol->dim->nu;
	int *nb = qp_sol->dim->nb;
	int *ng = qp_sol->dim->ng;
	int *ns = qp_sol->dim->ns;

	int ii;

	for(ii=0; ii<N; ii++)
		{
		VECCP_LIBSTR(nx[ii+1], qp_sol->pi+ii, 0, pi+ii, 0);
		}

	for(ii=0; ii<=N; ii++)
		{
		VECCP_LIBSTR(nu[ii], qp_sol->ux+ii, 0, u+ii, 0);
		VECCP_LIBSTR(nx[ii], qp_sol->ux+ii, nu[ii], x+ii, 0);
		if(nb[ii]>0)
			{
			VECCP_LIBSTR(nb[ii], qp_sol->lam+ii, 0, lam_lb+ii, 0);
			VECCP_LIBSTR(nb[ii], qp_sol->lam+ii, nb[ii]+ng[ii], lam_ub+ii, 0);
			}
		if(ng[ii]>0)
			{
			VECCP_LIBSTR(ng[ii], qp_sol->lam+ii, nb[ii], lam_lg+ii, 0);
			VECCP_LIBSTR(ng[ii], qp_sol->lam+ii, 2*nb[ii]+ng[ii], lam_ug+ii, 0);
			}
		if(ns[ii]>0)
			{
			VECCP_LIBSTR(ns[ii], qp_sol->ux+ii, nu[ii]+nx[ii], ls+ii, 0);
			VECCP_LIBSTR(ns[ii], qp_sol->ux+ii, nu[ii]+nx[ii]+ns[ii], us+ii, 0);
			VECCP_LIBSTR(ns[ii], qp_sol->lam+ii, 2*nb[ii]+2*ng[ii], lam_ls+ii, 0);
			VECCP_LIBSTR(ns[ii], qp_sol->lam+ii, 2*nb[ii]+2*ng[ii]+ns[ii], lam_us+ii, 0);
			}
		}

	return;

	}
