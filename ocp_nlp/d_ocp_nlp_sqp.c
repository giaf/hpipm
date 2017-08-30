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
#include <blasfeo_d_blas.h>

#include "../include/hpipm_d_rk_int.h"
#include "../include/hpipm_d_erk_int.h"
#include "../include/hpipm_d_ocp_qp.h"
#include "../include/hpipm_d_ocp_qp_sol.h"
#include "../include/hpipm_d_ocp_qp_ipm.h"
#include "../include/hpipm_d_ocp_nlp.h"
#include "../include/hpipm_d_ocp_nlp_sol.h"
#include "../include/hpipm_d_ocp_nlp_sqp.h"



#define CREATE_ERK_INT d_create_erk_int
#define CREATE_OCP_QP d_create_ocp_qp
#define CREATE_OCP_QP_IPM d_create_ocp_qp_ipm
#define CREATE_OCP_QP_SOL d_create_ocp_qp_sol
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
#define OCP_QP_IPM_WORKSPACE d_ocp_qp_ipm_workspace
#define OCP_QP_SOL d_ocp_qp_sol

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

	int *i_ptr;

	int size = 0;

	size += 1*sizeof(struct OCP_QP);
	size += 1*sizeof(struct OCP_QP_SOL);
	size += 1*sizeof(struct OCP_QP_IPM_WORKSPACE);
	size += N*sizeof(struct ERK_WORKSPACE);
	size += N*sizeof(struct ERK_ARG);
	size += N*sizeof(double *);

	size += MEMSIZE_OCP_QP(N, nx, nu, nb, ng, ns);
	size += MEMSIZE_OCP_QP_SOL(N, nx, nu, nb, ng, ns);

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
		size += N*nx[ii]*(nu[ii]+nx[ii])*sizeof(double);
		}

	size = (size+63)/64*64; // make multiple of typical cache line size

	return size;

	}



// TODO eliminate x0 in QP !!!
void CREATE_OCP_NLP_SQP(struct OCP_NLP *nlp, struct OCP_NLP_SQP_ARG *arg, struct OCP_NLP_SQP_WORKSPACE *ws, void *mem)
	{

	ws->memsize = MEMSIZE_OCP_NLP_SQP(nlp, arg);

	int ii, jj;

	int N = nlp->N;
	int *nx = nlp->nx;
	int *nu = nlp->nu;
	int *nb = nlp->nb;
	int *ng = nlp->ng;
	int *ns = nlp->ns;

	// ocp qp
	struct OCP_QP *qp_ptr = mem;
	//
	ws->qp = qp_ptr;
	qp_ptr += 1;

	// ocp qp sol
	struct OCP_QP_SOL *qp_sol_ptr = (struct OCP_QP_SOL *) qp_ptr;
	//
	ws->qp_sol = qp_sol_ptr;
	qp_sol_ptr += 1;

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

	// erk ws
	struct ERK_ARG *erk_arg_ptr = (struct ERK_ARG *) erk_ws_ptr;
	//
	ws->erk_arg = erk_arg_ptr;
	erk_arg_ptr += N;

	// void stuf
	char *c_ptr = (char *) erk_arg_ptr;
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
	double **dp_ptr = (double **) c_ptr;
	//
	ws->fs = dp_ptr;
	dp_ptr += N;

	//
	double *d_ptr = (double *) dp_ptr;
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
	for(ii=0; ii<N; ii++)
		ws->erk_arg[ii] = arg->erk_arg[ii];

	return;

	}



int SOLVE_OCP_NLP_SQP(struct OCP_NLP *nlp, struct OCP_NLP_SOL *nlp_sol, struct OCP_NLP_SQP_WORKSPACE *ws)
	{

	struct OCP_QP *qp = ws->qp;
	struct OCP_QP_SOL *qp_sol = ws->qp_sol;
	struct OCP_QP_IPM_WORKSPACE *ipm_ws = ws->ipm_workspace;
	struct ERK_WORKSPACE *erk_ws = ws->erk_workspace;
	struct ERK_ARG *erk_arg = ws->erk_arg;

	int ss, nn, ii;

	int N = nlp->N;
	int *nx = nlp->nx;
	int *nu = nlp->nu;

	// initialize solution to zero
	for(nn=0; nn<=N; nn++)
		for(ii=0; ii<nu[nn]+nx[nn]; ii++)
			(nlp_sol->ux+nn)->pa[ii] = 0.0;

//	double *fs1; d_zeros(&fs1, nx_*nf1, 1); // TODO
//	for(ii=0; ii<nx_; ii++) fs1[nu_*nx_+ii*(nx_+1)] = 1.0; // TODO

	int sqp_steps = 1;
	for(ss=0; ss<sqp_steps; ss++)	
		{

		// initial stage
		// XXX x0 in the QP is zero since x0 in the nlp is initialized to x0 !!!
//		nn = 0;
		// XXX it does not need the sensitivities wrt x here
//		d_init_erk_int(x0, fs0, u[nn], &d_linear_vde0, &ls, &erk_workspace0);
//		d_erk_int(&erk_arg, &erk_workspace0);
//		d_cvt_erk_int_to_ocp_qp(nn, &erk_workspace0, x[nn+1], &qp);

		// other stages
		for(nn=0; nn<N; nn++)
			{
//			d_print_mat(nx[nn], nu[nn]+nx[nn], ws->fs[nn], nx[nn]);
//			d_init_erk_int((nlp_sol->ux+nn)->pa[0]+nu[0], ws->fs[nn], (nlp_sol->ux+nn)->pa[0], nlp->expl_vde[nn], &ls, erk_ws+nn);
//			d_erk_int(erk_arg+nn, erk_ws+nn);
//			d_cvt_erk_int_to_ocp_qp(nn, erk_ws+nn, x[nn+1], qp);
			}

//		for(nn=0; nn<=N; nn++)
//			d_print_strmat(nu[nn]+nx[nn]+1, nu[nn]+nx[nn], qp.RSQrq+nn, 0, 0);
//		for(nn=0; nn<N; nn++)
//			d_print_strmat(nu[nn]+nx[nn]+1, nx[nn+1], qp.BAbt+nn, 0, 0);
//		for(nn=0; nn<=N; nn++)
//			d_print_tran_strvec(nb[nn], qp.d+nn, 0);

		d_solve_ocp_qp_ipm2(qp, qp_sol, ipm_ws);

//		d_print_e_tran_mat(5, workspace.iter, workspace.stat, 5);

//		d_cvt_ocp_qp_sol_to_colmaj(&qp, &qp_sol, du, dx, dls, dus, dpi, dlam_lb, dlam_ub, dlam_lg, dlam_ug, dlam_ls, dlam_us);


		for(nn=0; nn<=N; nn++)
			daxpy_libstr(nu[nn]+nx[nn], 1.0, qp_sol->ux+nn, 0, nlp_sol->ux+nn, 0, nlp_sol->ux+nn, 0);

		}
	
//	free(fs1); // TODO

	return 0;

	}
