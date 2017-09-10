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

#include "../include/hpipm_d_core_qp_ipm.h"
#include "../include/hpipm_d_core_qp_ipm_aux.h"
#include "../include/hpipm_d_rk_int.h"
#include "../include/hpipm_d_erk_int.h"
#include "../include/hpipm_d_ocp_qp.h"
#include "../include/hpipm_d_ocp_qp_sol.h"
#include "../include/hpipm_d_ocp_qp_ipm.h"
#include "../include/hpipm_d_ocp_qp_kkt.h"
#include "../include/hpipm_d_ocp_nlp.h"
#include "../include/hpipm_d_ocp_nlp_sol.h"
#include "../include/hpipm_d_ocp_nlp_ipm.h"
#include "../include/hpipm_d_ocp_qp_sim.h"



#define COMPUTE_ALPHA_QP d_compute_alpha_qp
#define COMPUTE_CENTERING_CORRECTION_QP d_compute_centering_correction_qp
#define COMPUTE_MU_AFF_QP d_compute_mu_aff_qp
#define COMPUTE_RES_OCP_QP d_compute_res_ocp_qp
#define CORE_QP_IPM_WORKSPACE d_core_qp_ipm_workspace
#define CREATE_ERK_INT d_create_erk_int
#define CREATE_OCP_QP d_create_ocp_qp
#define CREATE_OCP_QP_IPM d_create_ocp_qp_ipm
#define CREATE_OCP_QP_SOL d_create_ocp_qp_sol
#define CREATE_STRVEC d_create_strvec
#define ERK_ARG d_erk_args
#define ERK_WORKSPACE d_erk_workspace
#define FACT_SOLVE_KKT_STEP_OCP_QP d_fact_solve_kkt_step_ocp_qp
//#define FACT_SOLVE_KKT_UNCONSTR_OCP_QP d_fact_solve_kkt_unconstr_ocp_qp
#define INIT_VAR_OCP_QP d_init_var_ocp_qp
#define MEMSIZE_ERK_INT d_memsize_erk_int
#define MEMSIZE_OCP_NLP_IPM d_memsize_ocp_nlp_ipm
#define MEMSIZE_OCP_QP d_memsize_ocp_qp
#define MEMSIZE_OCP_QP_IPM d_memsize_ocp_qp_ipm
#define MEMSIZE_OCP_QP_SOL d_memsize_ocp_qp_sol
#define OCP_NLP d_ocp_nlp
#define OCP_NLP_IPM_ARG d_ocp_nlp_ipm_arg
#define OCP_NLP_SOL d_ocp_nlp_sol
#define OCP_NLP_IPM_WORKSPACE d_ocp_nlp_ipm_workspace
#define OCP_QP d_ocp_qp
#define OCP_QP_IPM_ARG d_ocp_qp_ipm_arg
#define OCP_QP_IPM_WORKSPACE d_ocp_qp_ipm_workspace
#define OCP_QP_SOL d_ocp_qp_sol
#define REAL double
#define SIZE_STRVEC d_size_strvec
#define SOLVE_KKT_STEP_OCP_QP d_solve_kkt_step_ocp_qp
#define STRVEC d_strvec
#define UPDATE_VAR_QP d_update_var_qp

#define MEMSIZE_OCP_NLP_IPM d_memsize_ocp_nlp_ipm
#define CREATE_OCP_NLP_IPM d_create_ocp_nlp_ipm
#define SOLVE_OCP_NLP_IPM d_solve_ocp_nlp_ipm



// TODO eliminate x0 in QP !!!
int MEMSIZE_OCP_NLP_IPM(struct OCP_NLP *nlp, struct OCP_NLP_IPM_ARG *arg)
	{

	int ii;

	int N = nlp->N;
	int *nx = nlp->nx;
	int *nu = nlp->nu;
	int *nb = nlp->nb;
	int *ng = nlp->ng;
	int *ns = nlp->ns;
	int ne0 = nlp->ne0;

	int nuxM = 0;
	int nbgM = 0;
	for(ii=0; ii<=N; ii++)
		{
		nuxM = nu[ii]+nx[ii]>nuxM ? nu[ii]+nx[ii] : nuxM;
		nbgM = nb[ii]+ng[ii]>nuxM ? nb[ii]+ng[ii] : nbgM;
		}
	
//	int tmp_nx0 = 0;

	int *i_ptr;

	int size = 0;

	size += 1*sizeof(struct OCP_QP);
	size += 2*sizeof(struct OCP_QP_SOL);
	size += 1*sizeof(struct OCP_QP_IPM_WORKSPACE);
	size += N*sizeof(struct ERK_WORKSPACE);
	size += 2*sizeof(struct STRVEC); // tmp_nuxM tmp_nbgM
	size += N*sizeof(REAL *); // fs

	size += SIZE_STRVEC(nuxM); // tmp_nuxM
	size += SIZE_STRVEC(nbgM); // tmp_nbgM

	// remove initial state form optimization variables
//	if(ne0>0) // TODO what if 0 < ne0 < nx[0] ???
//		{
//		nx0_tmp = nx[0];
//		nx[0] = 0;
//		}

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

	// restore initial state into optimization variables
//	if(ne0>0) // TODO what if 0 < ne0 < nx[0] ???
//		{
//		nx[0] = nx0_tmp;
//		}

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
void CREATE_OCP_NLP_IPM(struct OCP_NLP *nlp, struct OCP_NLP_IPM_ARG *arg, struct OCP_NLP_IPM_WORKSPACE *ws, void *mem)
	{

	int ii, jj;

	int N = nlp->N;
	int *nx = nlp->nx;
	int *nu = nlp->nu;
	int *nb = nlp->nb;
	int *ng = nlp->ng;
	int *ns = nlp->ns;
	int ne0 = nlp->ne0;

	int nuxM = 0;
	int nbgM = 0;
	for(ii=0; ii<=N; ii++)
		{
		nuxM = nu[ii]+nx[ii]>nuxM ? nu[ii]+nx[ii] : nuxM;
		nbgM = nb[ii]+ng[ii]>nuxM ? nb[ii]+ng[ii] : nbgM;
		}

	int tmp_nx0 = 0;

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

	// remove initial state form optimization variables
//	if(ne0>0) // TODO what if 0 < ne0 < nx[0] ???
//		{
//		nx0_tmp = nx[0];
//		nx[0] = 0;
//		}

	//
	CREATE_OCP_QP(N, nx, nu, nb, ng, ns, ws->qp, c_ptr);
	c_ptr += ws->qp->memsize;
	//
	CREATE_OCP_QP_SOL(N, nx, nu, nb, ng, ns, ws->qp_sol, c_ptr);
	c_ptr += ws->qp_sol->memsize;
	//
	CREATE_OCP_QP_IPM(ws->qp, arg->ipm_arg, ws->ipm_workspace, c_ptr);
	c_ptr += ws->ipm_workspace->memsize;

	// restore initial state into optimization variables
//	if(ne0>0) // TODO what if 0 < ne0 < nx[0] ???
//		{
//		nx[0] = nx0_tmp;
//		}

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
	CREATE_STRVEC(nuxM, ws->tmp_nuxM+0, c_ptr);
	c_ptr += (ws->tmp_nuxM+0)->memory_size;
	//
	CREATE_STRVEC(nbgM, ws->tmp_nbgM+0, c_ptr);
	c_ptr += (ws->tmp_nbgM+0)->memory_size;


	ws->memsize = MEMSIZE_OCP_NLP_IPM(nlp, arg);


#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + ws->memsize)
		{
		printf("\nCreate_ocp_nlp_sqp: outsize memory bounds!\n\n");
		exit(1);
		}
#endif


	return;

	}



int SOLVE_OCP_NLP_IPM(struct OCP_NLP *nlp, struct OCP_NLP_SOL *nlp_sol, struct OCP_NLP_IPM_ARG *arg, struct OCP_NLP_IPM_WORKSPACE *ws)
	{

	struct OCP_QP *qp = ws->qp;
	struct OCP_QP_SOL *qp_sol = ws->qp_sol;
	struct OCP_QP_SOL *tmp_qp_sol = ws->qp_sol+1;
	struct OCP_QP_IPM_WORKSPACE *ipm_ws = ws->ipm_workspace;
	struct ERK_WORKSPACE *erk_ws = ws->erk_workspace;
	struct STRVEC *tmp_nuxM = ws->tmp_nuxM;
	struct STRVEC *tmp_nbgM = ws->tmp_nbgM;

	struct CORE_QP_IPM_WORKSPACE *cws = ipm_ws->core_workspace;

	struct OCP_QP_IPM_ARG *ipm_arg = arg->ipm_arg;
	struct ERK_ARG *erk_arg = arg->erk_arg;

	// alias qp vectors into qp_sol
	cws->v = qp_sol->ux->pa;
	cws->pi = qp_sol->pi->pa;
	cws->lam = qp_sol->lam->pa;
	cws->t = qp_sol->t->pa;

	ipm_ws->mu0 = ipm_arg->mu0;

	int ss, nn, ii;

	int kk;

	// qp size
	int N = qp->N;
	int *nx = qp->nx;
	int *nu = qp->nu;
	int *nb = qp->nb;
	int *ng = qp->ng;
	int *ns = qp->ns;

	int ne0 = nlp->ne0;

	double *x, *xn, *u;

	double tmp;

	struct STRVEC str_res_g;
	struct STRVEC str_res_b;
	struct STRVEC str_res_d;
	struct STRVEC str_res_m;
	str_res_g.m = cws->nv;
	str_res_b.m = cws->ne;
	str_res_d.m = cws->nc;
	str_res_m.m = cws->nc;
	str_res_g.pa = cws->res_g;
	str_res_b.pa = cws->res_b;
	str_res_d.pa = cws->res_d;
	str_res_m.pa = cws->res_m;

	double nlp_res[4];
	double qp_res[4];


	// initialize nlp sol
	for(nn=0; nn<=N; nn++)
		dvecse_libstr(nlp->nu[nn]+nlp->nx[nn], 0.0, nlp_sol->ux+nn, 0);
	for(nn=0; nn<N; nn++)
		dvecse_libstr(nlp->nx[nn+1], 0.0, nlp_sol->pi+nn, 0);
	for(nn=0; nn<=N; nn++)
		dvecse_libstr(2*nlp->nb[nn]+2*nlp->ng[nn], 0.0, nlp_sol->lam+nn, 0);
	for(nn=0; nn<=N; nn++)
		dvecse_libstr(2*nlp->nb[nn]+2*nlp->ng[nn], 0.0, nlp_sol->t+nn, 0);
	// fix initial state to e0
	dveccp_libstr(ne0, nlp->e0, 0, nlp_sol->ux+0, nlp->nu[0]);


	// copy nlp into qp
	nn = 0;
	if(ne0>0)
		{
		// TODO what if 0 < ne0 < nx[0] ???
		// remove initial state form QP variables
		nx[0] = 0;
		// cost function
		dgecp_libstr(nu[0], nu[0], nlp->RSQ+0, 0, 0, qp->RSQrq+0, 0, 0);
		dgemv_t_libstr(nlp->nx[0], nlp->nu[0], 1.0, nlp->RSQ+0, nlp->nu[0], 0, nlp->e0, 0, 1.0, nlp->rq+0, 0, qp->rq+0, 0);
		// box constraints
		// XXX assume that there are not box constraints on x0 !!!!!
		dveccp_libstr(nb[0], nlp->d+0, 0, qp->d+0, 0);
		dveccp_libstr(nb[0], nlp->d+0, nlp->nb[0]+nlp->ng[0], qp->d+0, nb[0]+ng[0]);
		for(ii=0; ii<nb[0]; ii++) qp->idxb[0][ii] = nlp->idxb[0][ii];
		// general constraints
		dgecp_libstr(nu[0], ng[0], nlp->DCt+0, 0, 0, qp->DCt+0, 0, 0);
		dgemv_t_libstr(nlp->nx[0], nlp->ng[0], 1.0, nlp->DCt+0, nlp->nu[0], 0, nlp->e0, 0, 0.0, tmp_nbgM, 0, tmp_nbgM, 0);
		daxpy_libstr(ng[0], -1.0, tmp_nbgM, 0, nlp->d+0, nlp->nb[0], qp->d+0, nb[0]);
		daxpy_libstr(ng[0], -1.0, tmp_nbgM, 0, nlp->d+0, 2*(nlp->nb[0])+nlp->ng[0], qp->d+0, 2*nb[0]+ng[0]);
		// softed constraints
		dveccp_libstr(2*ns[0], nlp->d+0, 2*(nlp->nb[0]+nlp->ng[0]), qp->d+0, 2*nb[0]+2*ng[0]);
		for(ii=0; ii<ns[0]; ii++) qp->idxs[0][ii] = nlp->idxs[0][ii];
		nn++;
		}
	for(; nn<=N; nn++)
		{
		dgecp_libstr(nu[nn]+nx[nn], nu[nn]+nx[nn], nlp->RSQ+nn, 0, 0, qp->RSQrq+nn, 0, 0);
		dgecp_libstr(nu[nn]+nx[nn], ng[nn], nlp->DCt+nn, 0, 0, qp->DCt+nn, 0, 0);
		dveccp_libstr(nu[nn]+nx[nn], nlp->rq+nn, 0, qp->rq+nn, 0); // XXX
		drowin_libstr(nu[nn]+nx[nn], 1.0, qp->rq+nn, 0, qp->RSQrq+nn, nu[nn]+nx[nn], 0); // XXX
		dveccp_libstr(2*nb[nn]+2*ng[nn]+2*ns[nn], nlp->d+nn, 0, qp->d+nn, 0); // XXX
		for(ii=0; ii<nb[nn]; ii++) qp->idxb[nn][ii] = nlp->idxb[nn][ii];
		for(ii=0; ii<ns[nn]; ii++) qp->idxs[nn][ii] = nlp->idxs[nn][ii];
		}


	// initialize solution
	INIT_VAR_OCP_QP(qp, qp_sol, ipm_ws);


	// nlp loop
	for(ss=0; ss<arg->nlp_iter_max; ss++)	
		{

		// simulation & sensitivity propagation
		for(nn=0; nn<N; nn++)
			{
			x  = (nlp_sol->ux+nn)->pa+nu[nn];
			xn = (nlp_sol->ux+nn+1)->pa+nu[nn+1];
			u  = (nlp_sol->ux+nn)->pa;
			d_init_erk_int(x, ws->fs[nn], u, (nlp->model+nn)->expl_vde, (nlp->model+nn)->arg, erk_ws+nn);
			d_erk_int(erk_arg+nn, erk_ws+nn);
			d_cvt_erk_int_to_ocp_qp(nn, erk_ws+nn, xn, qp, nlp_sol);
			}

//for(ii=0; ii<N; ii++)
//	d_print_e_strmat(nlp->nu[ii]+nlp->nx[ii]+1, nlp->nx[ii+1], qp->BAbt+ii, 0, 0);
	

		// eliminate x0 from optimization variables
		dgemv_t_libstr(nlp->nx[0], nlp->nx[1], 1.0, qp->BAbt+0, nu[0], 0, nlp->e0, 0, 1.0, qp->b+0, 0, qp->b+0, 0);
		// TODO r0 d_lg0 d_ug0

#if 0
printf("\n%d %d\n", nu[0], nx[0]);
d_print_tran_strvec(nu[0]+nx[0], qp->rq+0, 0);
d_print_tran_strvec(nx[1], qp->b+0, 0);
exit(1);
#endif

		// compute residuals
		COMPUTE_RES_OCP_QP(qp, qp_sol, ipm_ws);
		cws->mu = ipm_ws->res_mu;

		// compute infinity norm of residuals
		dvecnrm_inf_libstr(cws->nv, &str_res_g, 0, &nlp_res[0]);
		dvecnrm_inf_libstr(cws->ne, &str_res_b, 0, &nlp_res[1]);
		dvecnrm_inf_libstr(cws->nc, &str_res_d, 0, &nlp_res[2]);
		dvecnrm_inf_libstr(cws->nc, &str_res_m, 0, &nlp_res[3]);

#if 0
printf("\nresiduals\n");
d_print_e_mat(1, cws->nv, cws->res_g, 1);
d_print_e_mat(1, cws->ne, cws->res_b, 1);
d_print_e_mat(1, cws->nc, cws->res_d, 1);
d_print_e_mat(1, cws->nc, cws->res_m, 1);
#endif

#if 0
printf("\n\niter %d nlp inf norm res %e %e %e %e\n", ss, nlp_res[0], nlp_res[1], nlp_res[2], nlp_res[3]);
#endif

		// exit condition on residuals
		if(!(nlp_res[0]>arg->nlp_res_g_max | nlp_res[1]>arg->nlp_res_b_max | nlp_res[2]>arg->nlp_res_d_max | nlp_res[3]>arg->nlp_res_m_max))
			{
			ws->iter = ss;
			ws->nlp_res_g = nlp_res[0];
			ws->nlp_res_b = nlp_res[1];
			ws->nlp_res_d = nlp_res[2];
			ws->nlp_res_m = nlp_res[3];
			return 0;
			}

		qp_res[0] = 1.0;
		qp_res[1] = 1.0;
		qp_res[2] = 1.0;
		qp_res[3] = 1.0;

		// XXX
//		cws->mu = 1.0;

		// qp loop
//		for(kk=0; kk<ipm_arg->iter_max & cws->mu>ipm_arg->mu_max; kk++)
		for(kk=0; kk<ipm_arg->iter_max & (qp_res[0]>ipm_arg->res_g_max | qp_res[1]>ipm_arg->res_b_max | qp_res[2]>ipm_arg->res_d_max | qp_res[3]>ipm_arg->res_m_max); kk++)
			{

			// fact and solve kkt
			FACT_SOLVE_KKT_STEP_OCP_QP(qp, ipm_ws);

			// alpha
			COMPUTE_ALPHA_QP(cws);
			if(kk<ipm_ws->stat_max)
				ipm_ws->stat[5*kk+0] = cws->alpha;

			// Mehrotra's corrector
			if(ipm_arg->pred_corr==1)
				{
				// mu_aff
				COMPUTE_MU_AFF_QP(cws);
				if(kk<ipm_ws->stat_max)
					ipm_ws->stat[5*kk+1] = cws->mu_aff;

				tmp = cws->mu_aff/cws->mu;
				cws->sigma = tmp*tmp*tmp;
				if(kk<ipm_ws->stat_max)
					ipm_ws->stat[5*kk+2] = cws->sigma;

				COMPUTE_CENTERING_CORRECTION_QP(cws);

				// fact and solve kkt
				SOLVE_KKT_STEP_OCP_QP(qp, ipm_ws);

				// alpha
				COMPUTE_ALPHA_QP(cws);
				if(kk<ipm_ws->stat_max)
					ipm_ws->stat[5*kk+3] = cws->alpha;
				}

#if 0
printf("\nstep\n");
d_print_e_mat(1, cws->nv, cws->dv, 1);
d_print_e_mat(1, cws->ne, cws->dpi, 1);
d_print_e_mat(1, cws->nc, cws->dlam, 1);
d_print_e_mat(1, cws->nc, cws->dt, 1);
#endif
			//
			UPDATE_VAR_QP(cws);

			// compute residuals
			COMPUTE_RES_OCP_QP(qp, qp_sol, ipm_ws);
			cws->mu = ipm_ws->res_mu;
			if(kk<ipm_ws->stat_max)
				ipm_ws->stat[5*kk+4] = ipm_ws->res_mu;

			// compute infinity norm of residuals
			dvecnrm_inf_libstr(cws->nv, &str_res_g, 0, &qp_res[0]);
			dvecnrm_inf_libstr(cws->ne, &str_res_b, 0, &qp_res[1]);
			dvecnrm_inf_libstr(cws->nc, &str_res_d, 0, &qp_res[2]);
			dvecnrm_inf_libstr(cws->nc, &str_res_m, 0, &qp_res[3]);

#if 0
printf("\nqp inf norm res %e %e %e %e\n", qp_res[0], qp_res[1], qp_res[2], qp_res[3]);
#endif

			}
		
		ipm_ws->iter = kk;

#if 0
d_print_e_tran_mat(5, kk, ipm_ws->stat, 5);
#endif

		// update NLP variables
		for(nn=0; nn<=N; nn++)
			dveccp_libstr(nu[nn]+nx[nn]+2*ns[ii], qp_sol->ux+nn, 0, nlp_sol->ux+nn, 0);
		for(nn=0; nn<N; nn++)
			dveccp_libstr(nx[nn+1], qp_sol->pi+nn, 0, nlp_sol->pi+nn, 0);
		for(nn=0; nn<=N; nn++)
			dveccp_libstr(2*nb[nn]+2*ng[nn]+2*ns[ii], qp_sol->lam+nn, 0, nlp_sol->lam+nn, 0);
		for(nn=0; nn<=N; nn++)
			dveccp_libstr(2*nb[nn]+2*ng[nn]+2*ns[ii], qp_sol->t+nn, 0, nlp_sol->t+nn, 0);

		}
	
	// maximum iteration number reached
	ws->iter = ss;
	ws->nlp_res_g = nlp_res[0];
	ws->nlp_res_b = nlp_res[1];
	ws->nlp_res_d = nlp_res[2];
	ws->nlp_res_m = nlp_res[3];

	return 0;

	}

