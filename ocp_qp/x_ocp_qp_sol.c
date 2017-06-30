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



int MEMSIZE_OCP_QP_SOL(int N, int *nx, int *nu, int *nb, int *ng)
	{

	int ii;

	int nvt = 0;
	int net = 0;
	int nbt = 0;
	int ngt = 0;
	for(ii=0; ii<N; ii++)
		{
		nvt += nu[ii]+nx[ii];
		net += nx[ii+1];
		nbt += nb[ii];
		ngt += ng[ii];
		}
	nvt += nu[ii]+nx[ii];
	nbt += nb[ii];
	ngt += ng[ii];

	int size = 0;

	size += (0+10*(N+1))*sizeof(struct STRVEC); // ux pi lam_lb lam_ub lam_lg lam_ug t_lb t_ub t_lg t_ug

	size += 1*SIZE_STRVEC(nvt); // ux
	size += 1*SIZE_STRVEC(net); // pi
	size += 8*SIZE_STRVEC(2*nbt+2*ngt); // lam_lb lam_ub lam_lg lam_ug t_lb t_ub t_lg t_ug

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 64; // align to typical cache line size
	
	return size;

	}



void CREATE_OCP_QP_SOL(int N, int *nx, int *nu, int *nb, int *ng, struct OCP_QP_SOL *qp_sol, void *memory)
	{

	int ii;


	// memsize
	qp_sol->memsize = MEMSIZE_OCP_QP_SOL(N, nx, nu, nb, ng);


	int nvt = 0;
	int net = 0;
	int nbt = 0;
	int ngt = 0;
	for(ii=0; ii<N; ii++)
		{
		nvt += nu[ii]+nx[ii];
		net += nx[ii+1];
		nbt += nb[ii];
		ngt += ng[ii];
		}
	nvt += nu[ii]+nx[ii];
	nbt += nb[ii];
	ngt += ng[ii];


	// vector struct stuff
	struct STRVEC *sv_ptr = (struct STRVEC *) memory;

	qp_sol->ux = sv_ptr;
	sv_ptr += N+1;
	qp_sol->pi = sv_ptr;
	sv_ptr += N+1;
	qp_sol->lam_lb = sv_ptr;
	sv_ptr += N+1;
	qp_sol->lam_ub = sv_ptr;
	sv_ptr += N+1;
	qp_sol->lam_lg = sv_ptr;
	sv_ptr += N+1;
	qp_sol->lam_ug = sv_ptr;
	sv_ptr += N+1;
	qp_sol->t_lb = sv_ptr;
	sv_ptr += N+1;
	qp_sol->t_ub = sv_ptr;
	sv_ptr += N+1;
	qp_sol->t_lg = sv_ptr;
	sv_ptr += N+1;
	qp_sol->t_ug = sv_ptr;
	sv_ptr += N+1;


	// align to typical cache line size
	long long l_ptr = (long long) sv_ptr;
	l_ptr = (l_ptr+63)/64*64;


	// double stuff
	void *v_ptr;
	v_ptr = (void *) l_ptr;

	void *tmp_ptr;

	// ux
	tmp_ptr = v_ptr;
	v_ptr += SIZE_STRVEC(nvt);
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(nu[ii]+nx[ii], qp_sol->ux+ii, tmp_ptr);
		tmp_ptr += (nu[ii]+nx[ii])*sizeof(REAL);
		}
	// pi
	tmp_ptr = v_ptr;
	v_ptr += SIZE_STRVEC(net);
	for(ii=0; ii<N; ii++)
		{
		CREATE_STRVEC(nx[ii+1], qp_sol->pi+ii, tmp_ptr);
		tmp_ptr += (nx[ii+1])*sizeof(REAL);
		}
	// lam
	tmp_ptr = v_ptr;
	v_ptr += SIZE_STRVEC(2*nbt+2*ngt);
	// lam_lb
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(nb[ii], qp_sol->lam_lb+ii, tmp_ptr);
		tmp_ptr += (nb[ii])*sizeof(REAL);
		}
	// lam_lg
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(ng[ii], qp_sol->lam_lg+ii, tmp_ptr);
		tmp_ptr += (ng[ii])*sizeof(REAL);
		}
	// lam_ub
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(nb[ii], qp_sol->lam_ub+ii, tmp_ptr);
		tmp_ptr += (nb[ii])*sizeof(REAL);
		}
	// lam_ug
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(ng[ii], qp_sol->lam_ug+ii, tmp_ptr);
		tmp_ptr += (ng[ii])*sizeof(REAL);
		}
	// t
	tmp_ptr = v_ptr;
	v_ptr += SIZE_STRVEC(2*nbt+2*ngt);
	// t_lb
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(nb[ii], qp_sol->t_lb+ii, tmp_ptr);
		tmp_ptr += (nb[ii])*sizeof(REAL);
		}
	// t_lg
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(ng[ii], qp_sol->t_lg+ii, tmp_ptr);
		tmp_ptr += (ng[ii])*sizeof(REAL);
		}
	// t_ub
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(nb[ii], qp_sol->t_ub+ii, tmp_ptr);
		tmp_ptr += (nb[ii])*sizeof(REAL);
		}
	// t_ug
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(ng[ii], qp_sol->t_ug+ii, tmp_ptr);
		tmp_ptr += (ng[ii])*sizeof(REAL);
		}
	
	return;

	}



void CVT_OCP_QP_SOL_TO_COLMAJ(struct OCP_QP *qp, struct OCP_QP_SOL *qp_sol, REAL **u, REAL **x, REAL **pi, REAL **lam_lb, REAL **lam_ub, REAL **lam_lg, REAL **lam_ug)
	{

	int N = qp->N;
	int *nx = qp->nx;
	int *nu = qp->nu;
	int *nb = qp->nb;
	int *ng = qp->ng;
	
	int ii;

	for(ii=0; ii<N; ii++)
		{
		CVT_STRVEC2VEC(nu[ii], qp_sol->ux+ii, 0, u[ii]);
		CVT_STRVEC2VEC(nx[ii], qp_sol->ux+ii, nu[ii], x[ii]);
		CVT_STRVEC2VEC(nx[ii+1], qp_sol->pi+ii, 0, pi[ii]);
		CVT_STRVEC2VEC(nb[ii], qp_sol->lam_lb+ii, 0, lam_lb[ii]);
		CVT_STRVEC2VEC(nb[ii], qp_sol->lam_ub+ii, 0, lam_ub[ii]);
		CVT_STRVEC2VEC(ng[ii], qp_sol->lam_lg+ii, 0, lam_lg[ii]);
		CVT_STRVEC2VEC(ng[ii], qp_sol->lam_ug+ii, 0, lam_ug[ii]);
		}
	CVT_STRVEC2VEC(nu[ii], qp_sol->ux+ii, 0, u[ii]);
	CVT_STRVEC2VEC(nx[ii], qp_sol->ux+ii, nu[ii], x[ii]);
	CVT_STRVEC2VEC(nb[ii], qp_sol->lam_lb+ii, 0, lam_lb[ii]);
	CVT_STRVEC2VEC(nb[ii], qp_sol->lam_ub+ii, 0, lam_ub[ii]);
	CVT_STRVEC2VEC(ng[ii], qp_sol->lam_lg+ii, 0, lam_lg[ii]);
	CVT_STRVEC2VEC(ng[ii], qp_sol->lam_ug+ii, 0, lam_ug[ii]);

	return;

	}



void CVT_OCP_QP_SOL_TO_ROWMAJ(struct OCP_QP *qp, struct OCP_QP_SOL *qp_sol, REAL **u, REAL **x, REAL **pi, REAL **lam_lb, REAL **lam_ub, REAL **lam_lg, REAL **lam_ug)
	{

	int N = qp->N;
	int *nx = qp->nx;
	int *nu = qp->nu;
	int *nb = qp->nb;
	int *ng = qp->ng;
	
	int ii;

	for(ii=0; ii<N; ii++)
		{
		CVT_STRVEC2VEC(nu[ii], qp_sol->ux+ii, 0, u[ii]);
		CVT_STRVEC2VEC(nx[ii], qp_sol->ux+ii, nu[ii], x[ii]);
		CVT_STRVEC2VEC(nx[ii+1], qp_sol->pi+ii, 0, pi[ii]);
		CVT_STRVEC2VEC(nb[ii], qp_sol->lam_lb+ii, 0, lam_lb[ii]);
		CVT_STRVEC2VEC(nb[ii], qp_sol->lam_ub+ii, 0, lam_ub[ii]);
		CVT_STRVEC2VEC(ng[ii], qp_sol->lam_lg+ii, 0, lam_lg[ii]);
		CVT_STRVEC2VEC(ng[ii], qp_sol->lam_ug+ii, 0, lam_ug[ii]);
		}
	CVT_STRVEC2VEC(nu[ii], qp_sol->ux+ii, 0, u[ii]);
	CVT_STRVEC2VEC(nx[ii], qp_sol->ux+ii, nu[ii], x[ii]);
	CVT_STRVEC2VEC(nb[ii], qp_sol->lam_lb+ii, 0, lam_lb[ii]);
	CVT_STRVEC2VEC(nb[ii], qp_sol->lam_ub+ii, 0, lam_ub[ii]);
	CVT_STRVEC2VEC(ng[ii], qp_sol->lam_lg+ii, 0, lam_lg[ii]);
	CVT_STRVEC2VEC(ng[ii], qp_sol->lam_ug+ii, 0, lam_ug[ii]);

	return;

	}



void CVT_OCP_QP_SOL_TO_LIBSTR(struct OCP_QP *qp, struct OCP_QP_SOL *qp_sol, struct STRVEC *u, struct STRVEC *x, struct STRVEC *pi, struct STRVEC *lam_lb, struct STRVEC *lam_ub, struct STRVEC *lam_lg, struct STRVEC *lam_ug)
	{

	int N = qp->N;
	int *nx = qp->nx;
	int *nu = qp->nu;
	int *nb = qp->nb;
	int *ng = qp->ng;
	
	int ii;

	for(ii=0; ii<N; ii++)
		{
		VECCP_LIBSTR(nu[ii], qp_sol->ux+ii, 0, u+ii, 0);
		VECCP_LIBSTR(nx[ii], qp_sol->ux+ii, nu[ii], x+ii, 0);
		VECCP_LIBSTR(nx[ii+1], qp_sol->pi+ii, 0, pi+ii, 0);
		VECCP_LIBSTR(nb[ii], qp_sol->lam_lb+ii, 0, lam_lb+ii, 0);
		VECCP_LIBSTR(nb[ii], qp_sol->lam_ub+ii, 0, lam_ub+ii, 0);
		VECCP_LIBSTR(ng[ii], qp_sol->lam_lg+ii, 0, lam_lg+ii, 0);
		VECCP_LIBSTR(ng[ii], qp_sol->lam_ug+ii, 0, lam_ug+ii, 0);
		}
	VECCP_LIBSTR(nu[ii], qp_sol->ux+ii, 0, u+ii, 0);
	VECCP_LIBSTR(nx[ii], qp_sol->ux+ii, nu[ii], x+ii, 0);
	VECCP_LIBSTR(nb[ii], qp_sol->lam_lb+ii, 0, lam_lb+ii, 0);
	VECCP_LIBSTR(nb[ii], qp_sol->lam_ub+ii, 0, lam_ub+ii, 0);
	VECCP_LIBSTR(ng[ii], qp_sol->lam_lg+ii, 0, lam_lg+ii, 0);
	VECCP_LIBSTR(ng[ii], qp_sol->lam_ug+ii, 0, lam_ug+ii, 0);

	return;

	}
