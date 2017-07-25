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



#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_d_aux.h>

#include "../include/hpipm_tree.h"
#include "../include/hpipm_scenario_tree.h"
#include "../include/hpipm_d_tree_ocp_qp.h"

#define CREATE_STRMAT d_create_strmat
#define CREATE_STRVEC d_create_strvec
#define CVT_MAT2STRMAT d_cvt_mat2strmat
#define CVT_TRAN_MAT2STRMAT d_cvt_tran_mat2strmat
#define CVT_VEC2STRVEC d_cvt_vec2strvec
#define REAL double
#define STRMAT d_strmat
#define STRVEC d_strvec
#define TREE_OCP_QP d_tree_ocp_qp

#define MEMSIZE_TREE_OCP_QP d_memsize_tree_ocp_qp
#define CREATE_TREE_OCP_QP d_create_tree_ocp_qp
#define CVT_COLMAJ_TO_TREE_OCP_QP d_cvt_colmaj_to_tree_ocp_qp



int MEMSIZE_TREE_OCP_QP(struct tree *ttree, int *nx, int *nu, int *nb, int *ng)
	{

	int ii, idx, idxdad;

	int Nn = ttree->Nn;

	int size = 0;

	size += 4*Nn*sizeof(int); // nx nu nb ng
	size += 1*Nn*sizeof(int *); // idxb
	size += (1*(Nn-1)+2*Nn)*sizeof(struct d_strmat); // BAbt ; RSQrq DCt
	size += (1*(Nn-1)+5*Nn)*sizeof(struct d_strvec); // b ; rq d_lb d_ub d_lg d_ug

	for(ii=0; ii<Nn-1; ii++)
		{
		idx = ii+1;
		idxdad = (ttree->root+idx)->dad;
		size += d_size_strmat(nu[idxdad]+nx[idxdad]+1, nx[idx]); // BAbt
		size += d_size_strvec(nx[idx]); // b
		}
	
	for(ii=0; ii<Nn; ii++)
		{
		size += nb[ii]*sizeof(int); // idxb
		size += d_size_strmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii]); // RSQrq
		size += d_size_strvec(nu[ii]+nx[ii]); // rq
		size += d_size_strmat(nu[ii]+nx[ii], ng[ii]); // DCt
		size += 2*d_size_strvec(nb[ii]); // d_lb d_ub
		size += 2*d_size_strvec(ng[ii]); // d_lg d_ug
		}
	
	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 64; // align to typical cache line size

	return size;

	}



void CREATE_TREE_OCP_QP(struct tree *ttree, int *nx, int *nu, int *nb, int *ng, struct d_tree_ocp_qp *qp, void *memory)
	{

	int ii, idx, idxdad;


	// memsize
	qp->memsize = d_memsize_tree_ocp_qp(ttree, nx, nu, nb, ng);


	// tree
	qp->ttree = ttree;

	// number of nodes
	int Nn = ttree->Nn;
	qp->Nn = Nn;


	// int pointer stuff
	int **ip_ptr;
	ip_ptr = (int **) memory;

	// idxb
	qp->idxb = ip_ptr;
	ip_ptr += Nn;


	// matrix struct stuff
	struct STRMAT *sm_ptr = (struct STRMAT *) ip_ptr;

	// BAbt
	qp->BAbt = sm_ptr;
	sm_ptr += Nn-1;

	// RSQrq
	qp->RSQrq = sm_ptr;
	sm_ptr += Nn;

	// DCt
	qp->DCt = sm_ptr;
	sm_ptr += Nn;


	// vector struct stuff
	struct STRVEC *sv_ptr = (struct STRVEC *) sm_ptr;

	// b
	qp->b = sv_ptr;
	sv_ptr += Nn-1;

	// rq
	qp->rq = sv_ptr;
	sv_ptr += Nn;

	// d_lb
	qp->d_lb = sv_ptr;
	sv_ptr += Nn;

	// d_ub
	qp->d_ub = sv_ptr;
	sv_ptr += Nn;

	// d_lg
	qp->d_lg = sv_ptr;
	sv_ptr += Nn;

	// d_ug
	qp->d_ug = sv_ptr;
	sv_ptr += Nn;


	// integer stuff
	int *i_ptr;
	i_ptr = (int *) sv_ptr;

	// nx
	qp->nx = i_ptr;
	for(ii=0; ii<Nn; ii++)
		{
		i_ptr[ii] = nx[ii];
		}
	i_ptr += Nn;
	
	// nu
	qp->nu = i_ptr;
	for(ii=0; ii<Nn; ii++)
		{
		i_ptr[ii] = nu[ii];
		}
	i_ptr += Nn;
	
	// nb
	qp->nb = i_ptr;
	for(ii=0; ii<Nn; ii++)
		{
		i_ptr[ii] = nb[ii];
		}
	i_ptr += Nn;

	// ng
	qp->ng = i_ptr;
	for(ii=0; ii<Nn; ii++)
		{
		i_ptr[ii] = ng[ii];
		}
	i_ptr += Nn;
	
	// idxb
	for(ii=0; ii<Nn; ii++)
		{
		(qp->idxb)[ii] = i_ptr;
		i_ptr += nb[ii];
		}


	// align to typical cache line size
	long long l_ptr = (long long) i_ptr;
	l_ptr = (l_ptr+63)/64*64;


	// double stuff
	char *c_ptr;
	c_ptr = (char *) l_ptr;

	// BAbt
	for(ii=0; ii<Nn-1; ii++)
		{
		idx = ii+1;
		idxdad = (ttree->root+idx)->dad;
		CREATE_STRMAT(nu[idxdad]+nx[idxdad]+1, nx[idx], qp->BAbt+ii, c_ptr);
		c_ptr += (qp->BAbt+ii)->memory_size;
		}

	// RSQrq
	for(ii=0; ii<Nn; ii++)
		{
		CREATE_STRMAT(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], qp->RSQrq+ii, c_ptr);
		c_ptr += (qp->RSQrq+ii)->memory_size;
		}

	// DCt
	for(ii=0; ii<Nn; ii++)
		{
		CREATE_STRMAT(nu[ii]+nx[ii], ng[ii], qp->DCt+ii, c_ptr);
		c_ptr += (qp->DCt+ii)->memory_size;
		}

	// b
	for(ii=0; ii<Nn-1; ii++)
		{
		idx = ii+1;
		CREATE_STRVEC(nx[idx], qp->b+ii, c_ptr);
		c_ptr += (qp->b+ii)->memory_size;
		}

	// rq
	for(ii=0; ii<Nn; ii++)
		{
		CREATE_STRVEC(nu[ii]+nx[ii], qp->rq+ii, c_ptr);
		c_ptr += (qp->rq+ii)->memory_size;
		}

	// d_lb
	for(ii=0; ii<Nn; ii++)
		{
		CREATE_STRVEC(nb[ii], qp->d_lb+ii, c_ptr);
		c_ptr += (qp->d_lb+ii)->memory_size;
		}

	// d_ub
	for(ii=0; ii<Nn; ii++)
		{
		CREATE_STRVEC(nb[ii], qp->d_ub+ii, c_ptr);
		c_ptr += (qp->d_ub+ii)->memory_size;
		}

	// d_lg
	for(ii=0; ii<Nn; ii++)
		{
		CREATE_STRVEC(ng[ii], qp->d_lg+ii, c_ptr);
		c_ptr += (qp->d_lg+ii)->memory_size;
		}

	// d_ug
	for(ii=0; ii<Nn; ii++)
		{
		CREATE_STRVEC(ng[ii], qp->d_ug+ii, c_ptr);
		c_ptr += (qp->d_ug+ii)->memory_size;
		}

	return;

	}



void CVT_COLMAJ_TO_TREE_OCP_QP(REAL **A, REAL **B, REAL **b, REAL **Q, REAL **S, REAL **R, REAL **q, REAL **r, int **idxb, REAL **d_lb, REAL **d_ub, REAL **C, REAL **D, REAL **d_lg, REAL **d_ug, struct TREE_OCP_QP *qp)
	{

	int Nn = qp->Nn;
	int *nx = qp->nx;
	int *nu = qp->nu;
	int *nb = qp->nb;
	int *ng = qp->ng;
	struct tree *ttree = qp->ttree;

	int ii, jj, idx, idxdad;

	for(ii=0; ii<Nn-1; ii++)
		{
		idx = ii+1;
		idxdad = (ttree->root+idx)->dad;
		CVT_TRAN_MAT2STRMAT(nx[idx], nu[idxdad], B[ii], nx[idx], qp->BAbt+ii, 0, 0);
		CVT_TRAN_MAT2STRMAT(nx[idx], nx[idxdad], A[ii], nx[idx], qp->BAbt+ii, nu[idxdad], 0);
		CVT_TRAN_MAT2STRMAT(nx[idx], 1, b[ii], nx[idx], qp->BAbt+ii, nu[idxdad]+nx[idxdad], 0);
		CVT_VEC2STRVEC(nx[idx], b[ii], qp->b+ii, 0);
		}
	
	for(ii=0; ii<Nn; ii++)
		{
		CVT_MAT2STRMAT(nu[ii], nu[ii], R[ii], nu[ii], qp->RSQrq+ii, 0, 0);
		CVT_TRAN_MAT2STRMAT(nu[ii], nx[ii], S[ii], nu[ii], qp->RSQrq+ii, nu[ii], 0);
		CVT_MAT2STRMAT(nx[ii], nx[ii], Q[ii], nx[ii], qp->RSQrq+ii, nu[ii], nu[ii]);
		CVT_TRAN_MAT2STRMAT(nu[ii], 1, r[ii], nu[ii], qp->RSQrq+ii, nu[ii]+nx[ii], 0);
		CVT_TRAN_MAT2STRMAT(nx[ii], 1, q[ii], nx[ii], qp->RSQrq+ii, nu[ii]+nx[ii], nu[ii]);
		CVT_VEC2STRVEC(nu[ii], r[ii], qp->rq+ii, 0);
		CVT_VEC2STRVEC(nx[ii], q[ii], qp->rq+ii, nu[ii]);
		}
	
	for(ii=0; ii<Nn; ii++)
		{
		for(jj=0; jj<nb[ii]; jj++)
			qp->idxb[ii][jj] = idxb[ii][jj];
		CVT_VEC2STRVEC(nb[ii], d_lb[ii], qp->d_lb+ii, 0);
		CVT_VEC2STRVEC(nb[ii], d_ub[ii], qp->d_ub+ii, 0);
		}
	
	for(ii=0; ii<Nn; ii++)
		{
		CVT_TRAN_MAT2STRMAT(ng[ii], nu[ii], D[ii], ng[ii], qp->DCt+ii, 0, 0);
		CVT_TRAN_MAT2STRMAT(ng[ii], nx[ii], C[ii], ng[ii], qp->DCt+ii, nu[ii], 0);
		CVT_VEC2STRVEC(ng[ii], d_lg[ii], qp->d_lg+ii, 0);
		CVT_VEC2STRVEC(ng[ii], d_ug[ii], qp->d_ug+ii, 0);
		}

	return;

	}


