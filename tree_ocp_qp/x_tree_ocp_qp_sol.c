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

int MEMSIZE_TREE_OCP_QP_SOL(struct TREE_OCP_QP_DIM *dim)
	{

	// extract dim
	struct tree *ttree = dim->ttree;
	int Nn = dim->Nn;
	int *nx = dim->nx;
	int *nu = dim->nu;
	int *nb = dim->nb;
	int *ng = dim->ng;
	int *ns = dim->ns;

	int ii, idx, idxdad;

	int nvt = 0;
	int net = 0;
	int nct = 0;
	for(ii=0; ii<Nn; ii++)
		{
		nvt += nu[ii]+nx[ii]+2*ns[ii];
		nct += 2*nb[ii]+2*ng[ii]+2*ns[ii];
		}
	for(ii=0; ii<Nn-1; ii++)
		{
		idx = ii+1;
		net += nx[idx];
		}

	int size = 0;

	size += 3*Nn*sizeof(struct STRVEC); // ux lam t
	size += 1*(Nn-1)*sizeof(struct STRVEC); // pi

	size += 1*SIZE_STRVEC(nvt); // ux
	size += 1*SIZE_STRVEC(net); // pi
	size += 2*SIZE_STRVEC(nct); // lam t

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 64; // align to typical cache line size
	
	return size;

	}



void CREATE_TREE_OCP_QP_SOL(struct TREE_OCP_QP_DIM *dim, struct TREE_OCP_QP_SOL *qp_sol, void *mem)
	{

	// extract dim
	struct tree *ttree = dim->ttree;
	int Nn = dim->Nn;
	int *nx = dim->nx;
	int *nu = dim->nu;
	int *nb = dim->nb;
	int *ng = dim->ng;
	int *ns = dim->ns;

	int ii, idx, idxdad;


	int nvt = 0;
	int net = 0;
	int nct = 0;
	for(ii=0; ii<Nn-1; ii++)
		{
		nvt += nu[ii]+nx[ii]+2*ns[ii];
		net += nx[ii+1];
		nct += 2*nb[ii]+2*ng[ii]+2*ns[ii];
		}
	nvt += nu[ii]+nx[ii]+2*ns[ii];
	nct += 2*nb[ii]+2*ng[ii]+2*ns[ii];


	// vector struct stuff
	struct STRVEC *sv_ptr = (struct STRVEC *) mem;

	qp_sol->ux = sv_ptr;
	sv_ptr += Nn;
	qp_sol->pi = sv_ptr;
	sv_ptr += Nn-1;
	qp_sol->lam = sv_ptr;
	sv_ptr += Nn;
	qp_sol->t = sv_ptr;
	sv_ptr += Nn;


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
	for(ii=0; ii<Nn; ii++)
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
	for(ii=0; ii<Nn-1; ii++)
		{
		idx = ii+1;
		CREATE_STRVEC(nx[idx], qp_sol->pi+ii, tmp_ptr);
		tmp_ptr += (nx[idx])*sizeof(REAL); // pi
		}
	// lam
	tmp_ptr = c_ptr;
	c_ptr += SIZE_STRVEC(nct);
	for(ii=0; ii<Nn; ii++)
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
	c_ptr += SIZE_STRVEC(nct);
	for(ii=0; ii<Nn; ii++)
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

	qp_sol->memsize = MEMSIZE_TREE_OCP_QP_SOL(dim);


#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + qp_sol->memsize)
		{
		printf("\nCreate_tree_ocp_qp_sol: outsize memory bounds!\n\n");
		exit(1);
		}
#endif


	return;

	}



void CVT_TREE_OCP_QP_SOL_TO_COLMAJ(struct TREE_OCP_QP *qp, struct TREE_OCP_QP_SOL *qp_sol, REAL **u, REAL **x, REAL **ls, REAL **us, REAL **pi, REAL **lam_lb, REAL **lam_ub, REAL **lam_lg, REAL **lam_ug, REAL **lam_ls, REAL **lam_us)
	{

	int Nn = qp->dim->Nn;
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;
	int *nb = qp->dim->nb;
	int *ng = qp->dim->ng;
	int *ns = qp->dim->ns;

	int ii;

	for(ii=0; ii<Nn-1; ii++)
		{
		CVT_STRVEC2VEC(nx[ii+1], qp_sol->pi+ii, 0, pi[ii]);
		}

	for(ii=0; ii<Nn; ii++)
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



void CVT_TREE_OCP_QP_SOL_TO_ROWMAJ(struct TREE_OCP_QP *qp, struct TREE_OCP_QP_SOL *qp_sol, REAL **u, REAL **x, REAL **ls, REAL **us, REAL **pi, REAL **lam_lb, REAL **lam_ub, REAL **lam_lg, REAL **lam_ug, REAL **lam_ls, REAL **lam_us)
	{

	int Nn = qp->dim->Nn;
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;
	int *nb = qp->dim->nb;
	int *ng = qp->dim->ng;
	int *ns = qp->dim->ns;

	int ii;

	for(ii=0; ii<Nn-1; ii++)
		{
		CVT_STRVEC2VEC(nx[ii+1], qp_sol->pi+ii, 0, pi[ii]);
		}

	for(ii=0; ii<Nn-1; ii++)
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


