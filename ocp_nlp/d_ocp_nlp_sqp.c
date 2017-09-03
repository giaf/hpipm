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



#if defined(RUNTIME_CHECKS)
#include <stdlib.h>
#include <stdio.h>
#endif

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_d_aux.h>
#include <blasfeo_d_blas.h>

#include "../include/hpipm_d_rk_int.h"
#include "../include/hpipm_d_erk_int.h"
#include "../include/hpipm_d_ocp_qp.h"
#include "../include/hpipm_d_ocp_qp_sol.h"
#include "../include/hpipm_d_ocp_qp_ipm.h"
#include "../include/hpipm_d_ocp_qp_kkt.h"
#include "../include/hpipm_d_ocp_qp_sim.h"
#include "../include/hpipm_d_ocp_nlp.h"
#include "../include/hpipm_d_ocp_nlp_sol.h"
#include "../include/hpipm_d_ocp_nlp_sqp.h"



#define CREATE_ERK_INT d_create_erk_int
#define CREATE_OCP_QP d_create_ocp_qp
#define CREATE_OCP_QP_IPM d_create_ocp_qp_ipm
#define CREATE_OCP_QP_SOL d_create_ocp_qp_sol
#define CREATE_STRVEC d_create_strvec
#define ERK_ARG d_erk_args
#define ERK_WORKSPACE d_erk_workspace
#define MEMSIZE_ERK_INT d_memsize_erk_int
#define MEMSIZE_OCP_NLP_SQP d_memsize_ocp_nlp_sqp
#define MEMSIZE_OCP_QP d_memsize_ocp_qp
#define MEMSIZE_OCP_QP_IPM d_memsize_ocp_qp_ipm
#define MEMSIZE_OCP_QP_SOL d_memsize_ocp_qp_sol
#define OCP_NLP d_ocp_nlp
#define OCP_NLP_SQP_ARG d_ocp_nlp_sqp_arg
#define OCP_NLP_SOL d_ocp_nlp_sol
#define OCP_NLP_SQP_WORKSPACE d_ocp_nlp_sqp_workspace
#define OCP_QP d_ocp_qp
#define OCP_QP_IPM_ARG d_ocp_qp_ipm_arg
#define OCP_QP_IPM_WORKSPACE d_ocp_qp_ipm_workspace
#define OCP_QP_SOL d_ocp_qp_sol
#define REAL double
#define SIZE_STRVEC d_size_strvec
#define STRVEC d_strvec

#define MEMSIZE_OCP_NLP_SQP d_memsize_ocp_nlp_sqp
#define CREATE_OCP_NLP_SQP d_create_ocp_nlp_sqp
#define SOLVE_OCP_NLP_SQP d_solve_ocp_nlp_sqp



// TODO eliminate x0 in QP !!!
int MEMSIZE_OCP_NLP_SQP(struct OCP_NLP *nlp, struct OCP_NLP_SQP_ARG *arg)
	{

	int ii;

	int N = nlp->N;
	int *nx = nlp->nx;
	int *nu = nlp->nu;
	int *nb = nlp->nb;
	int *ng = nlp->ng;
	int *ns = nlp->ns;

	int nuxM = 0;
	int nbgM = 0;
	for(ii=0; ii<=N; ii++)
		{
		nuxM = nu[ii]+nx[ii]>nuxM ? nu[ii]+nx[ii] : nuxM;
		nbgM = nb[ii]+ng[ii]>nuxM ? nb[ii]+ng[ii] : nbgM;
		}

	int *i_ptr;

	int size = 0;

	size += 1*sizeof(struct OCP_QP);
	size += 2*sizeof(struct OCP_QP_SOL);
	size += 1*sizeof(struct OCP_QP_IPM_WORKSPACE);
	size += N*sizeof(struct ERK_WORKSPACE);
	size += 2*sizeof(struct STRVEC); // tmp_nuxM tmp_nbgM
	size += 1*(N+1)*sizeof(struct STRVEC); // rq
	size += N*sizeof(REAL *); // fs

	size += MEMSIZE_OCP_QP(N, nx, nu, nb, ng, ns);
	size += MEMSIZE_OCP_QP_SOL(N, nx, nu, nb, ng, ns);
	size += SIZE_STRVEC(nuxM); // tmp_nuxM
	size += SIZE_STRVEC(nbgM); // tmp_nbgM

	for(ii=0; ii<=N; ii++)
		{
		size += SIZE_STRVEC(nu[ii]+nx[ii]); // rq
		}

	struct OCP_QP qp;
	qp.N = nlp->N;
	qp.nx = nlp->nx;
	qp.nu = nlp->nu;
	qp.nb = nlp->nb;
	qp.ng = nlp->ng;
	qp.ns = nlp->ns;
	qp.idxb = nlp->idxb;
	qp.idxs = nlp->idxs;

	size += MEMSIZE_OCP_QP_IPM(&qp, arg->ipm_arg);

	for(ii=0; ii<N; ii++)
		{
		size += MEMSIZE_ERK_INT(arg->rk_data, nx[ii], nx[ii]+nu[ii], nu[ii]);
		size += N*nx[ii]*(nu[ii]+nx[ii])*sizeof(REAL); // fs
		}

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

	}



// TODO eliminate x0 in QP !!!
void CREATE_OCP_NLP_SQP(struct OCP_NLP *nlp, struct OCP_NLP_SQP_ARG *arg, struct OCP_NLP_SQP_WORKSPACE *ws, void *mem)
	{

	int ii, jj;

	int N = nlp->N;
	int *nx = nlp->nx;
	int *nu = nlp->nu;
	int *nb = nlp->nb;
	int *ng = nlp->ng;
	int *ns = nlp->ns;

	int nuxM = 0;
	int nbgM = 0;
	for(ii=0; ii<=N; ii++)
		{
		nuxM = nu[ii]+nx[ii]>nuxM ? nu[ii]+nx[ii] : nuxM;
		nbgM = nb[ii]+ng[ii]>nuxM ? nb[ii]+ng[ii] : nbgM;
		}

	// ocp qp
	struct OCP_QP *qp_ptr = mem;
	//
	ws->qp = qp_ptr;
	qp_ptr += 1;

	// ocp qp sol
	struct OCP_QP_SOL *qp_sol_ptr = (struct OCP_QP_SOL *) qp_ptr;
	//
	ws->qp_sol = qp_sol_ptr;
	qp_sol_ptr += 2;

	// ocp qp ipm ws
	struct OCP_QP_IPM_WORKSPACE *ipm_ws_ptr = (struct OCP_QP_IPM_WORKSPACE *) qp_sol_ptr;
	//
	ws->ipm_workspace = ipm_ws_ptr;
	ipm_ws_ptr += 1;

	// erk ws
	struct ERK_WORKSPACE *erk_ws_ptr = (struct ERK_WORKSPACE *) ipm_ws_ptr;
	//
	ws->erk_workspace = erk_ws_ptr;
	erk_ws_ptr += N;

	// void stuf
	char *c_ptr = (char *) erk_ws_ptr;
	//
	CREATE_OCP_QP(N, nx, nu, nb, ng, ns, ws->qp, c_ptr);
	c_ptr += ws->qp->memsize;
	//
	CREATE_OCP_QP_SOL(N, nx, nu, nb, ng, ns, ws->qp_sol, c_ptr);
	c_ptr += ws->qp_sol->memsize;
	//
	CREATE_OCP_QP_IPM(ws->qp, arg->ipm_arg, ws->ipm_workspace, c_ptr);
	c_ptr += ws->ipm_workspace->memsize;
	//
	for(ii=0; ii<N; ii++)
		{
		CREATE_ERK_INT(arg->rk_data, nx[ii], nx[ii]+nu[ii], nu[ii], ws->erk_workspace+ii, c_ptr);
		c_ptr += (ws->erk_workspace+ii)->memsize;
		}
	
	//
	REAL **dp_ptr = (REAL **) c_ptr;
	//
	ws->fs = dp_ptr;
	dp_ptr += N;

	//
	REAL *d_ptr = (REAL *) dp_ptr;
	//
	for(ii=0; ii<N; ii++)
		{
		ws->fs[ii] = d_ptr;
		d_ptr += nx[ii]*(nu[ii]+nx[ii]);
		for(jj=0; jj<nx[ii]*(nu[ii]+nx[ii]); jj++)
			ws->fs[ii][jj] = 0.0;
		for(jj=0; jj<nx[ii]; jj++)
			ws->fs[ii][nu[ii]*nx[ii]+jj*(nx[ii]+1)] = 1.0;
		}
	
	//
	struct STRVEC *sv_ptr = (struct STRVEC *) d_ptr;
	//
	ws->rq = sv_ptr;
	sv_ptr += N+1;
	//
	ws->tmp_nuxM = sv_ptr;
	sv_ptr += 1;
	//
	ws->tmp_nbgM = sv_ptr;
	sv_ptr += 1;


	// align to typicl cache line size
	size_t s_ptr = (size_t) sv_ptr;
	s_ptr = (s_ptr+63)/64*64;


	// void stuf
	c_ptr = (char *) s_ptr;
	//
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(nu[ii]+nx[ii], ws->rq+ii, c_ptr);
		c_ptr += (ws->rq+ii)->memory_size;
		}
	//
	CREATE_STRVEC(nuxM, ws->tmp_nuxM+0, c_ptr);
	c_ptr += (ws->tmp_nuxM+0)->memory_size;
	//
	CREATE_STRVEC(nbgM, ws->tmp_nbgM+0, c_ptr);
	c_ptr += (ws->tmp_nbgM+0)->memory_size;


	ws->memsize = MEMSIZE_OCP_NLP_SQP(nlp, arg);


#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + ws->memsize)
		{
		printf("\nCreate_ocp_nlp_sqp: outsize memory bounds!\n\n");
		exit(1);
		}
#endif


	return;

	}



int SOLVE_OCP_NLP_SQP(struct OCP_NLP *nlp, struct OCP_NLP_SOL *nlp_sol, struct OCP_NLP_SQP_ARG *arg, struct OCP_NLP_SQP_WORKSPACE *ws)
	{

	struct OCP_QP *qp = ws->qp;
	struct OCP_QP_SOL *qp_sol = ws->qp_sol;
	struct OCP_QP_SOL *tmp_qp_sol = ws->qp_sol+1;
	struct OCP_QP_IPM_WORKSPACE *ipm_ws = ws->ipm_workspace;
	struct ERK_WORKSPACE *erk_ws = ws->erk_workspace;
	struct STRVEC *tmp_nuxM = ws->tmp_nuxM;
	struct STRVEC *tmp_nbgM = ws->tmp_nbgM;

	struct OCP_QP_IPM_ARG *ipm_arg = arg->ipm_arg;
	struct ERK_ARG *erk_arg = arg->erk_arg;

	int ss, nn, ii;

	int N = nlp->N;
	int *nx = nlp->nx;
	int *nu = nlp->nu;
	int *nb = nlp->nb;
	int *ng = nlp->ng;
	int *ns = nlp->ns;

	double *x, *xn, *u;

	// compute rq
	for(nn=0; nn<=N; nn++)
		{
		dsymv_l_libstr(nu[nn]+nx[nn], nu[nn]+nx[nn], 1.0, nlp->RSQ+nn, 0, 0, nlp->ux_ref, 0, 0.0, ws->rq+nn, 0, ws->rq+nn, 0);
		}

	// copy nlp into qp
	for(nn=0; nn<=N; nn++)
		{
		dgecp_libstr(nu[nn]+nx[nn], nu[nn]+nx[nn], nlp->RSQ+nn, 0, 0, qp->RSQrq+nn, 0, 0);
		dgecp_libstr(nu[nn]+nx[nn], ng[nn], nlp->DCt+nn, 0, 0, qp->DCt+nn, 0, 0);
		dveccp_libstr(nu[nn]+nx[nn], ws->rq+nn, 0, qp->rq+nn, 0); // XXX
		drowin_libstr(nu[nn]+nx[nn], 1.0, qp->rq+nn, 0, qp->RSQrq+nn, nu[nn]+nx[nn], 0); // XXX
		dveccp_libstr(2*nb[nn]+2*ng[nn]+2*ns[nn], nlp->d+nn, 0, qp->d+nn, 0); // XXX
		for(ii=0; ii<nb[nn]; ii++) qp->idxb[nn][ii] = nlp->idxb[nn][ii];
		for(ii=0; ii<ns[nn]; ii++) qp->idxs[nn][ii] = nlp->idxs[nn][ii];
		}

	// initialize solution (to zero atm)
	for(nn=0; nn<=N; nn++)
		dvecse_libstr(nu[nn]+nx[nn], 0.0, nlp_sol->ux+nn, 0);

	int sqp_steps = 2;
	for(ss=0; ss<sqp_steps; ss++)	
		{

		// simulation & sensitivity propagation
		for(nn=0; nn<N; nn++)
			{
			x  = (nlp_sol->ux+nn)->pa+nu[nn];
			xn = (nlp_sol->ux+nn+1)->pa+nu[nn+1];
			u  = (nlp_sol->ux+nn)->pa;
			d_init_erk_int(x, ws->fs[nn], u, (nlp->model+nn)->expl_vde, (nlp->model+nn)->arg, erk_ws+nn);
			d_erk_int(erk_arg+nn, erk_ws+nn);
			d_cvt_erk_int_to_ocp_qp(nn, erk_ws+nn, xn, qp);
			}

	
		// setup qp
		for(nn=0; nn<=N; nn++)
			dveccp_libstr(nu[nn]+nx[nn]+2*ns[nn], ws->rq+nn, 0, qp->rq+nn, 0);
		for(nn=0; nn<N; nn++)
			dvecse_libstr(nx[nn+1], 0.0, qp->b+nn, 0);
		for(nn=0; nn<=N; nn++)
			dveccp_libstr(2*nb[nn]+2*ng[nn]+2*ns[nn], nlp->d+nn, 0, qp->d+nn, 0);

		// copy nlp_sol into qp_sol
		tmp_qp_sol->ux = nlp_sol->ux;
		tmp_qp_sol->pi = nlp_sol->pi;
		tmp_qp_sol->lam = nlp_sol->lam;
		tmp_qp_sol->t = nlp_sol->t;

		// compute residuals // XXX use adjoing sensitivities to avoid A'*pi ???
		d_compute_res_ocp_qp(qp, tmp_qp_sol, ipm_ws);

		// copy residuals into qp rhs
		for(nn=0; nn<=N; nn++)
			{
			dveccp_libstr(nu[nn]+nx[nn]+2*ns[nn], ipm_ws->res_g+nn, 0, qp->rq+nn, 0);
			drowin_libstr(nu[nn]+nx[nn], 1.0, qp->rq+nn, 0, qp->RSQrq+nn, nu[nn]+nx[nn], 0);
			}
		for(nn=0; nn<N; nn++)
			{
			dveccp_libstr(nx[nn+1], ipm_ws->res_b+nn, 0, qp->b+nn, 0);
			drowin_libstr(nx[nn+1], 1.0, qp->b+nn, 0, qp->BAbt+nn, nu[nn]+nx[nn], 0);
			}
		for(nn=0; nn<=N; nn++)
			{
			dveccp_libstr(2*nb[nn]+2*ng[nn], ipm_ws->res_d+nn, 0, qp->d+nn, 0);
			dvecsc_libstr(nb[nn]+ng[nn], -1.0, qp->d+nn, nb[nn]+ng[nn]);
			}


#if 0
		for(nn=0; nn<=N; nn++)
			d_print_strmat(nu[nn]+nx[nn]+1, nu[nn]+nx[nn], qp->RSQrq+nn, 0, 0);
		for(nn=0; nn<N; nn++)
			d_print_strmat(nu[nn]+nx[nn]+1, nx[nn+1], qp->BAbt+nn, 0, 0);
		for(nn=0; nn<=N; nn++)
			d_print_tran_strvec(2*nb[nn]+2*ng[nn]+2*ns[nn], qp->d+nn, 0);
//		exit(1);
#endif

		// solve qp
		d_solve_ocp_qp_ipm(qp, qp_sol, ipm_arg, ipm_ws);

#if 0
		for(nn=0; nn<=N; nn++)
			d_print_tran_strvec(nu[nn]+nx[nn]+2*ns[nn], qp_sol->ux+nn, 0);
		for(nn=0; nn<N; nn++)
			d_print_tran_strvec(nx[nn+1], qp_sol->pi+nn, 0);
		for(nn=0; nn<=N; nn++)
			d_print_tran_strvec(2*nb[nn]+2*ng[nn]+2*ns[nn], qp_sol->lam+nn, 0);
		for(nn=0; nn<=N; nn++)
			d_print_tran_strvec(2*nb[nn]+2*ng[nn]+2*ns[nn], qp_sol->t+nn, 0);
		d_print_e_tran_mat(5, ipm_ws->iter, ipm_ws->stat, 5);
#endif

		// update variables (full step)
		for(nn=0; nn<=N; nn++)
			daxpy_libstr(nu[nn]+nx[nn]+2*ns[ii], 1.0, qp_sol->ux+nn, 0, nlp_sol->ux+nn, 0, nlp_sol->ux+nn, 0);
		for(nn=0; nn<N; nn++)
			daxpy_libstr(nx[nn+1], 1.0, qp_sol->pi+nn, 0, nlp_sol->pi+nn, 0, nlp_sol->pi+nn, 0);
		for(nn=0; nn<=N; nn++)
			daxpy_libstr(2*nb[nn]+2*ng[nn]+2*ns[ii], 1.0, qp_sol->lam+nn, 0, nlp_sol->lam+nn, 0, nlp_sol->lam+nn, 0);
		for(nn=0; nn<=N; nn++)
			daxpy_libstr(2*nb[nn]+2*ng[nn]+2*ns[ii], 1.0, qp_sol->t+nn, 0, nlp_sol->t+nn, 0, nlp_sol->t+nn, 0);

		}
	
	return 0;

	}
