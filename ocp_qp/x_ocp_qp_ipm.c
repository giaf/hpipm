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



int MEMSIZE_OCP_QP_IPM(struct OCP_QP *qp, struct OCP_QP_IPM_ARG *arg)
	{

	// loop index
	int ii;

	// extract ocp qp size
	int N = qp->N;
	int *nx = qp->nx;
	int *nu = qp->nu;
	int *nb = qp->nb;
	int *ng = qp->ng;
	int *ns = qp->ns;

	// compute core qp size and max size
	int nvt = 0;
	int net = 0;
	int nct = 0;
	int nxM = 0;
	int nuM = 0;
	int nbM = 0;
	int ngM = 0;
	int nsM = 0;
	for(ii=0; ii<N; ii++)
		{
		nvt += nx[ii]+nu[ii]+2*ns[ii];
		net += nx[ii+1];
		nct += nb[ii]+ng[ii]+ns[ii];
		nxM = nx[ii]>nxM ? nx[ii] : nxM;
		nuM = nu[ii]>nuM ? nu[ii] : nuM;
		nbM = nb[ii]>nbM ? nb[ii] : nbM;
		ngM = ng[ii]>ngM ? ng[ii] : ngM;
		nsM = ns[ii]>nsM ? ns[ii] : nsM;
		}
	nvt += nx[ii]+nu[ii]+2*ns[ii];
	nct += nb[ii]+ng[ii]+ns[ii];
	nxM = nx[ii]>nxM ? nx[ii] : nxM;
	nuM = nu[ii]>nuM ? nu[ii] : nuM;
	nbM = nb[ii]>nbM ? nb[ii] : nbM;
	ngM = ng[ii]>ngM ? ng[ii] : ngM;
	nsM = ns[ii]>nsM ? ns[ii] : nsM;

	int size = 0;

	size += 8*(N+1)*sizeof(struct STRVEC); // dux dt res_g res_d res_m Gamma gamma Zs_inv
	size += 3*N*sizeof(struct STRVEC); // dpi res_b Pb
	size += 6*sizeof(struct STRVEC); // tmp_nxM 4*tmp_nbgM tmp_nsM

	size += 1*(N+1)*sizeof(struct STRMAT); // L
	size += 2*sizeof(struct STRMAT); // AL

	size += 1*SIZE_STRVEC(nxM); // tmp_nxM
	size += 4*SIZE_STRVEC(nbM+ngM); // tmp_nbgM
	size += 1*SIZE_STRVEC(nsM); // tmp_nsM
	for(ii=0; ii<N; ii++) size += 1*SIZE_STRVEC(nx[ii+1]); // Pb
	for(ii=0; ii<=N; ii++) size += 1*SIZE_STRVEC(2*ns[ii]); // Zs_inv
	for(ii=0; ii<=N; ii++) size += 1*SIZE_STRMAT(nu[ii]+nx[ii]+1, nu[ii]+nx[ii]); // L
	size += 2*SIZE_STRMAT(nuM+nxM+1, nxM+ngM); // AL

	size += 1*sizeof(struct CORE_QP_IPM_WORKSPACE);
	size += 1*MEMSIZE_CORE_QP_IPM(nvt, net, nct, arg->iter_max);

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

	}



void CREATE_OCP_QP_IPM(struct OCP_QP *qp, struct OCP_QP_IPM_ARG *arg, struct OCP_QP_IPM_WORKSPACE *workspace, void *mem)
	{

	// loop index
	int ii;

	// extract ocp qp size
	int N = qp->N;
	int *nx = qp->nx;
	int *nu = qp->nu;
	int *nb = qp->nb;
	int *ng = qp->ng;
	int *ns = qp->ns;


	workspace->memsize = MEMSIZE_OCP_QP_IPM(qp, arg);


	// compute core qp size and max size
	int nvt = 0;
	int net = 0;
	int nct = 0;
	int nxM = 0;
	int nuM = 0;
	int nbM = 0;
	int ngM = 0;
	int nsM = 0;
	for(ii=0; ii<N; ii++)
		{
		nvt += nx[ii]+nu[ii]+2*ns[ii];
		net += nx[ii+1];
		nct += nb[ii]+ng[ii]+ns[ii];
		nxM = nx[ii]>nxM ? nx[ii] : nxM;
		nuM = nu[ii]>nuM ? nu[ii] : nuM;
		nbM = nb[ii]>nbM ? nb[ii] : nbM;
		ngM = ng[ii]>ngM ? ng[ii] : ngM;
		nsM = ns[ii]>nsM ? ns[ii] : nsM;
		}
	nvt += nx[ii]+nu[ii]+2*ns[ii];
	nct += nb[ii]+ng[ii]+ns[ii];
	nxM = nx[ii]>nxM ? nx[ii] : nxM;
	nuM = nu[ii]>nuM ? nu[ii] : nuM;
	nbM = nb[ii]>nbM ? nb[ii] : nbM;
	ngM = ng[ii]>ngM ? ng[ii] : ngM;
	nsM = ns[ii]>nsM ? ns[ii] : nsM;


	// core struct
	struct CORE_QP_IPM_WORKSPACE *sr_ptr = mem;

	// core workspace
	workspace->core_workspace = sr_ptr;
	sr_ptr += 1;
	struct CORE_QP_IPM_WORKSPACE *cws = workspace->core_workspace;


	// matrix struct
	struct STRMAT *sm_ptr = (struct STRMAT *) sr_ptr;

	workspace->L = sm_ptr;
	sm_ptr += N+1;
	workspace->AL = sm_ptr;
	sm_ptr += 2;


	// vector struct
	struct STRVEC *sv_ptr = (struct STRVEC *) sm_ptr;

	workspace->dux = sv_ptr;
	sv_ptr += N+1;
	workspace->dpi = sv_ptr;
	sv_ptr += N;
	workspace->dt = sv_ptr;
	sv_ptr += N+1;
	workspace->res_g = sv_ptr;
	sv_ptr += N+1;
	workspace->res_b = sv_ptr;
	sv_ptr += N;
	workspace->res_d = sv_ptr;
	sv_ptr += N+1;
	workspace->res_m = sv_ptr;
	sv_ptr += N+1;
	workspace->Gamma = sv_ptr;
	sv_ptr += N+1;
	workspace->gamma = sv_ptr;
	sv_ptr += N+1;
	workspace->Pb = sv_ptr;
	sv_ptr += N;
	workspace->Zs_inv = sv_ptr;
	sv_ptr += N+1;
	workspace->tmp_nxM = sv_ptr;
	sv_ptr += 1;
	workspace->tmp_nbgM = sv_ptr;
	sv_ptr += 4;
	workspace->tmp_nsM = sv_ptr;
	sv_ptr += 1;


	// align to typicl cache line size
	size_t s_ptr = (size_t) sv_ptr;
	s_ptr = (s_ptr+63)/64*64;


	// void stuf
	char *c_ptr = (char *) s_ptr;

	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRMAT(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], workspace->L+ii, c_ptr);
		c_ptr += (workspace->L+ii)->memory_size;
		}

	CREATE_STRMAT(nuM+nxM+1, nxM+ngM, workspace->AL+0, c_ptr);
	c_ptr += (workspace->AL+0)->memory_size;

	CREATE_STRMAT(nuM+nxM+1, nxM+ngM, workspace->AL+1, c_ptr);
	c_ptr += (workspace->AL+1)->memory_size;

	for(ii=0; ii<N; ii++)
		{
		CREATE_STRVEC(nx[ii+1], workspace->Pb+ii, c_ptr);
		c_ptr += (workspace->Pb+ii)->memory_size;
		}

	for(ii=0; ii<N+1; ii++)
		{
		CREATE_STRVEC(2*ns[ii], workspace->Zs_inv+ii, c_ptr);
		c_ptr += (workspace->Zs_inv+ii)->memory_size;
		}

	CREATE_STRVEC(nxM, workspace->tmp_nxM, c_ptr);
	c_ptr += workspace->tmp_nxM->memory_size;

	CREATE_STRVEC(nbM+ngM, workspace->tmp_nbgM+0, c_ptr);
	c_ptr += (workspace->tmp_nbgM+0)->memory_size;

	CREATE_STRVEC(nbM+ngM, workspace->tmp_nbgM+1, c_ptr);
	c_ptr += (workspace->tmp_nbgM+1)->memory_size;

	CREATE_STRVEC(nbM+ngM, workspace->tmp_nbgM+2, c_ptr);
	c_ptr += (workspace->tmp_nbgM+2)->memory_size;

	CREATE_STRVEC(nbM+ngM, workspace->tmp_nbgM+3, c_ptr);
	c_ptr += (workspace->tmp_nbgM+3)->memory_size;

	CREATE_STRVEC(nsM, workspace->tmp_nsM+0, c_ptr);
	c_ptr += (workspace->tmp_nsM+0)->memory_size;



	cws->nv = nvt;
	cws->ne = net;
	cws->nc = nct;
	cws->iter_max = arg->iter_max;
	CREATE_CORE_QP_IPM(cws, c_ptr);
	c_ptr += workspace->core_workspace->memsize;

	cws->alpha_min = arg->alpha_min;
	cws->mu_max = arg->mu_max;
	cws->mu0 = arg->mu0;
	cws->nt_inv = 1.0/(2*nct); // TODO avoid computation if nt=0 XXX


	// alias members of workspace and core_workspace
	//
	c_ptr = (char *) cws->dv;
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(nu[ii]+nx[ii], workspace->dux+ii, c_ptr);
		c_ptr += (nu[ii]+nx[ii])*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		}
	//
	c_ptr = (char *) cws->dpi;
	for(ii=0; ii<N; ii++)
		{
		CREATE_STRVEC(nx[ii+1], workspace->dpi+ii, c_ptr);
		c_ptr += (nx[ii+1])*sizeof(REAL);
		}
	//
	c_ptr = (char *) cws->dt;
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(nb[ii], workspace->dt+ii, c_ptr);
		c_ptr += nb[ii]*sizeof(REAL);
		c_ptr += ng[ii]*sizeof(REAL);
		c_ptr += nb[ii]*sizeof(REAL);
		c_ptr += ng[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		}
	//
	c_ptr = (char *) cws->res_g;
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(nu[ii]+nx[ii], workspace->res_g+ii, c_ptr);
		c_ptr += (nu[ii]+nx[ii])*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		}
	//
	c_ptr = (char *) cws->res_b;
	for(ii=0; ii<N; ii++)
		{
		CREATE_STRVEC(nx[ii+1], workspace->res_b+ii, c_ptr);
		c_ptr += (nx[ii+1])*sizeof(REAL);
		}
	//
	c_ptr = (char *) cws->res_d;
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(nb[ii], workspace->res_d+ii, c_ptr);
		c_ptr += nb[ii]*sizeof(REAL);
		c_ptr += ng[ii]*sizeof(REAL);
		c_ptr += nb[ii]*sizeof(REAL);
		c_ptr += ng[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		}
	//
	c_ptr = (char *) cws->res_m;
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(nb[ii], workspace->res_m+ii, c_ptr);
		c_ptr += nb[ii]*sizeof(REAL);
		c_ptr += ng[ii]*sizeof(REAL);
		c_ptr += nb[ii]*sizeof(REAL);
		c_ptr += ng[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		}
	//
	c_ptr = (char *) cws->Gamma;
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(nb[ii], workspace->Gamma+ii, c_ptr);
		c_ptr += nb[ii]*sizeof(REAL);
		c_ptr += ng[ii]*sizeof(REAL);
		c_ptr += nb[ii]*sizeof(REAL);
		c_ptr += ng[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		}
	//
	c_ptr = (char *) cws->gamma;
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(nb[ii], workspace->gamma+ii, c_ptr);
		c_ptr += nb[ii]*sizeof(REAL);
		c_ptr += ng[ii]*sizeof(REAL);
		c_ptr += nb[ii]*sizeof(REAL);
		c_ptr += ng[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		}
	//
	workspace->stat = cws->stat;
	//
	workspace->iter = cws->iter_max;

	return;

	}



int SOLVE_OCP_QP_IPM(struct OCP_QP *qp, struct OCP_QP_SOL *qp_sol, struct OCP_QP_IPM_WORKSPACE *ws)
	{

	struct CORE_QP_IPM_WORKSPACE *cws = ws->core_workspace;

	// alias qp vectors into qp_sol
	cws->v = qp_sol->ux->pa;
	cws->pi = qp_sol->pi->pa;
	cws->lam = qp_sol->lam->pa;
	cws->t = qp_sol->t->pa;

//	printf("\n%e\n", cws->nt_inv);
//	exit(1);
//cws->nt_inv = 1.063830e-2;
//cws->nt_inv = 5.263158e-3;

	int kk = 0;

	if(cws->nc==0)
		{
		FACT_SOLVE_KKT_UNCONSTR_OCP_QP(qp, qp_sol, ws);
		COMPUTE_RES_OCP_QP(qp, qp_sol, ws);
		cws->mu = ws->res_mu;
		ws->iter = 0;
		return 0;
		}

	// init solver
	INIT_VAR_OCP_QP(qp, qp_sol, ws);

#if 0
int ii;
double *ptr;
printf("\nuxs\n");
ptr = cws->v;
for(ii=0; ii<=qp->N; ii++)
	{
	d_print_e_mat(1, qp->nu[ii]+qp->nx[ii]+2*qp->ns[ii], ptr, 1);
	ptr += qp->nu[ii]+qp->nx[ii]+2*qp->ns[ii];
	}
printf("\npi\n");
ptr = cws->pi;
for(ii=0; ii<qp->N; ii++)
	{
	d_print_e_mat(1, qp->nx[ii+1], ptr, 1);
	ptr += qp->nx[ii+1];
	}
printf("\nlam\n");
ptr = cws->lam;
for(ii=0; ii<=qp->N; ii++)
	{
	d_print_e_mat(1, 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii], ptr, 1);
	ptr += 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii];
	}
printf("\nt\n");
ptr = cws->t;
for(ii=0; ii<=qp->N; ii++)
	{
	d_print_e_mat(1, 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii], ptr, 1);
	ptr += 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii];
	}
exit(1);
#endif

	// compute residuals
	COMPUTE_RES_OCP_QP(qp, qp_sol, ws);
	cws->mu = ws->res_mu;

	for(kk=0; kk<cws->iter_max & cws->mu>cws->mu_max; kk++)
		{

#if 0
int ii;
double *ptr;
printf("\nres_g\n");
ptr = cws->res_g;
for(ii=0; ii<=qp->N; ii++)
	{
	d_print_e_mat(1, qp->nu[ii]+qp->nx[ii]+2*qp->ns[ii], ptr, 1);
	ptr += qp->nu[ii]+qp->nx[ii]+2*qp->ns[ii];
	}
printf("\nres_b\n");
ptr = cws->res_b;
for(ii=0; ii<qp->N; ii++)
	{
	d_print_e_mat(1, qp->nx[ii+1], ptr, 1);
	ptr += qp->nx[ii+1];
	}
printf("\nres_d\n");
ptr = cws->res_d;
for(ii=0; ii<=qp->N; ii++)
	{
	d_print_e_mat(1, 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii], ptr, 1);
	ptr += 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii];
	}
printf("\nres_m\n");
ptr = cws->res_m;
for(ii=0; ii<=qp->N; ii++)
	{
	d_print_e_mat(1, 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii], ptr, 1);
	ptr += 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii];
	}
exit(1);
#endif

		// fact and solve kkt
		FACT_SOLVE_KKT_STEP_OCP_QP(qp, ws);

#if 0
int ii;
double *ptr;
printf("\nduxs\n");
ptr = cws->dv;
for(ii=0; ii<=qp->N; ii++)
	{
	d_print_e_mat(1, qp->nu[ii]+qp->nx[ii]+2*qp->ns[ii], ptr, 1);
	ptr += qp->nu[ii]+qp->nx[ii]+2*qp->ns[ii];
	}
printf("\ndpi\n");
ptr = cws->dpi;
for(ii=0; ii<qp->N; ii++)
	{
	d_print_e_mat(1, qp->nx[ii+1], ptr, 1);
	ptr += qp->nx[ii+1];
	}
printf("\ndlam\n");
ptr = cws->dlam;
for(ii=0; ii<=qp->N; ii++)
	{
	d_print_e_mat(1, 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii], ptr, 1);
	ptr += 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii];
	}
printf("\ndt\n");
ptr = cws->dt;
for(ii=0; ii<=qp->N; ii++)
	{
	d_print_e_mat(1, 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii], ptr, 1);
	ptr += 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii];
	}
exit(1);
#endif

		// alpha
		COMPUTE_ALPHA_QP(cws);
		cws->stat[5*kk+0] = cws->alpha;

//printf("\nalpha = %e\n", cws->alpha);
//exit(1);
//cws->alpha = 4.768102e-1;
		//
		UPDATE_VAR_QP(cws);

#if 0
int ii;
double *ptr;
printf("\nuxs\n");
ptr = cws->v;
for(ii=0; ii<=qp->N; ii++)
	{
	d_print_e_mat(1, qp->nu[ii]+qp->nx[ii]+2*qp->ns[ii], ptr, 1);
	ptr += qp->nu[ii]+qp->nx[ii]+2*qp->ns[ii];
	}
printf("\npi\n");
ptr = cws->pi;
for(ii=0; ii<qp->N; ii++)
	{
	d_print_e_mat(1, qp->nx[ii+1], ptr, 1);
	ptr += qp->nx[ii+1];
	}
printf("\nlam\n");
ptr = cws->lam;
for(ii=0; ii<=qp->N; ii++)
	{
	d_print_e_mat(1, 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii], ptr, 1);
	ptr += 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii];
	}
printf("\nt\n");
ptr = cws->t;
for(ii=0; ii<=qp->N; ii++)
	{
	d_print_e_mat(1, 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii], ptr, 1);
	ptr += 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii];
	}
exit(1);
#endif

		// compute residuals
		COMPUTE_RES_OCP_QP(qp, qp_sol, ws);
		cws->mu = ws->res_mu;
		cws->stat[5*kk+1] = ws->res_mu;

#if 0
int ii;
double *ptr;
printf("\nres_g\n");
ptr = cws->res_g;
for(ii=0; ii<=qp->N; ii++)
	{
	d_print_e_mat(1, qp->nu[ii]+qp->nx[ii]+2*qp->ns[ii], ptr, 1);
	ptr += qp->nu[ii]+qp->nx[ii]+2*qp->ns[ii];
	}
printf("\nres_b\n");
ptr = cws->res_b;
for(ii=0; ii<qp->N; ii++)
	{
	d_print_e_mat(1, qp->nx[ii+1], ptr, 1);
	ptr += qp->nx[ii+1];
	}
printf("\nres_d\n");
ptr = cws->res_d;
for(ii=0; ii<=qp->N; ii++)
	{
	d_print_e_mat(1, 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii], ptr, 1);
	ptr += 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii];
	}
printf("\nres_m\n");
ptr = cws->res_m;
for(ii=0; ii<=qp->N; ii++)
	{
	d_print_e_mat(1, 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii], ptr, 1);
	ptr += 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii];
	}
exit(1);
#endif

//	ws->iter = kk+1;
//	return;
		}
	
	ws->iter = kk;
	
	// max iteration number reached
	if(kk==cws->iter_max)
		return 1;

	// normal return
	return 0;

	}



int SOLVE_OCP_QP_IPM2(struct OCP_QP *qp, struct OCP_QP_SOL *qp_sol, struct OCP_QP_IPM_WORKSPACE *ws)
	{

	struct CORE_QP_IPM_WORKSPACE *cws = ws->core_workspace;

	// alias qp vectors into qp_sol
	cws->v = qp_sol->ux->pa;
	cws->pi = qp_sol->pi->pa;
	cws->lam = qp_sol->lam->pa;
	cws->t = qp_sol->t->pa;

	REAL tmp;

	if(cws->nc==0)
		{
		FACT_SOLVE_KKT_UNCONSTR_OCP_QP(qp, qp_sol, ws);
		COMPUTE_RES_OCP_QP(qp, qp_sol, ws);
		cws->mu = ws->res_mu;
		ws->iter = 0;
		return 0;
		}

	// init solver
	INIT_VAR_OCP_QP(qp, qp_sol, ws);

#if 0
int ii;
double *ptr;
printf("\nuxs\n");
ptr = cws->v;
for(ii=0; ii<=qp->N; ii++)
	{
	d_print_e_mat(1, qp->nu[ii]+qp->nx[ii]+2*qp->ns[ii], ptr, 1);
	ptr += qp->nu[ii]+qp->nx[ii]+2*qp->ns[ii];
	}
printf("\npi\n");
ptr = cws->pi;
for(ii=0; ii<qp->N; ii++)
	{
	d_print_e_mat(1, qp->nx[ii+1], ptr, 1);
	ptr += qp->nx[ii+1];
	}
printf("\nlam\n");
ptr = cws->lam;
for(ii=0; ii<=qp->N; ii++)
	{
	d_print_e_mat(1, 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii], ptr, 1);
	ptr += 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii];
	}
printf("\nt\n");
ptr = cws->t;
for(ii=0; ii<=qp->N; ii++)
	{
	d_print_e_mat(1, 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii], ptr, 1);
	ptr += 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii];
	}
exit(1);
#endif

	// compute residuals
	COMPUTE_RES_OCP_QP(qp, qp_sol, ws);
	cws->mu = ws->res_mu;

#if 0
	int ii;
	for(ii=0; ii<1; ii++)
		{
		cws->sigma = 1.0;
		cws->stat[5*kk+2] = cws->sigma;
		COMPUTE_CENTERING_CORRECTION_QP(cws);
		FACT_SOLVE_KKT_STEP_OCP_QP(qp, ws);
		COMPUTE_ALPHA_QP(cws);
		cws->stat[5*kk+3] = cws->alpha;
		UPDATE_VAR_QP(cws);
		COMPUTE_RES_OCP_QP(qp, qp_sol, ws);
		cws->mu = ws->res_mu;
		cws->stat[5*kk+4] = ws->res_mu;
		kk++;
		}
//	ws->iter = kk;
//		return;
#endif

	int kk = 0;
	for(; kk<cws->iter_max & cws->mu>cws->mu_max; kk++)
		{

#if 0
int ii;
double *ptr;
printf("\nres_g\n");
ptr = cws->res_g;
for(ii=0; ii<=qp->N; ii++)
	{
	d_print_e_mat(1, qp->nu[ii]+qp->nx[ii]+2*qp->ns[ii], ptr, 1);
	ptr += qp->nu[ii]+qp->nx[ii]+2*qp->ns[ii];
	}
printf("\nres_b\n");
ptr = cws->res_b;
for(ii=0; ii<qp->N; ii++)
	{
	d_print_e_mat(1, qp->nx[ii+1], ptr, 1);
	ptr += qp->nx[ii+1];
	}
printf("\nres_d\n");
ptr = cws->res_d;
for(ii=0; ii<=qp->N; ii++)
	{
	d_print_e_mat(1, 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii], ptr, 1);
	ptr += 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii];
	}
printf("\nres_m\n");
ptr = cws->res_m;
for(ii=0; ii<=qp->N; ii++)
	{
	d_print_e_mat(1, 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii], ptr, 1);
	ptr += 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii];
	}
exit(1);
#endif

		// fact and solve kkt
		FACT_SOLVE_KKT_STEP_OCP_QP(qp, ws);

#if 0
int ii;
double *ptr;
printf("\nduxs\n");
ptr = cws->dv;
for(ii=0; ii<=qp->N; ii++)
	{
	d_print_e_mat(1, qp->nu[ii]+qp->nx[ii]+2*qp->ns[ii], ptr, 1);
	ptr += qp->nu[ii]+qp->nx[ii]+2*qp->ns[ii];
	}
printf("\ndpi\n");
ptr = cws->dpi;
for(ii=0; ii<qp->N; ii++)
	{
	d_print_e_mat(1, qp->nx[ii+1], ptr, 1);
	ptr += qp->nx[ii+1];
	}
printf("\ndlam\n");
ptr = cws->dlam;
for(ii=0; ii<=qp->N; ii++)
	{
	d_print_e_mat(1, 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii], ptr, 1);
	ptr += 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii];
	}
printf("\ndt\n");
ptr = cws->dt;
for(ii=0; ii<=qp->N; ii++)
	{
	d_print_e_mat(1, 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii], ptr, 1);
	ptr += 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii];
	}
exit(1);
#endif

		// alpha
		COMPUTE_ALPHA_QP(cws);
		cws->stat[5*kk+0] = cws->alpha;

		// mu_aff
		COMPUTE_MU_AFF_QP(cws);
		cws->stat[5*kk+1] = cws->mu_aff;

		tmp = cws->mu_aff/cws->mu;
		cws->sigma = tmp*tmp*tmp;
		cws->stat[5*kk+2] = cws->sigma;

		COMPUTE_CENTERING_CORRECTION_QP(cws);

		// fact and solve kkt
		SOLVE_KKT_STEP_OCP_QP(qp, ws);

#if 0
int ii;
double *ptr;
printf("\nduxs\n");
ptr = cws->dv;
for(ii=0; ii<=qp->N; ii++)
	{
	d_print_e_mat(1, qp->nu[ii]+qp->nx[ii]+2*qp->ns[ii], ptr, 1);
	ptr += qp->nu[ii]+qp->nx[ii]+2*qp->ns[ii];
	}
printf("\ndpi\n");
ptr = cws->dpi;
for(ii=0; ii<qp->N; ii++)
	{
	d_print_e_mat(1, qp->nx[ii+1], ptr, 1);
	ptr += qp->nx[ii+1];
	}
printf("\ndlam\n");
ptr = cws->dlam;
for(ii=0; ii<=qp->N; ii++)
	{
	d_print_e_mat(1, 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii], ptr, 1);
	ptr += 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii];
	}
printf("\ndt\n");
ptr = cws->dt;
for(ii=0; ii<=qp->N; ii++)
	{
	d_print_e_mat(1, 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii], ptr, 1);
	ptr += 2*qp->nb[ii]+2*qp->ng[ii]+2*qp->ns[ii];
	}
exit(1);
#endif

#if 0
int ii;
for(ii=0; ii<=qp->N; ii++)
	d_print_tran_strvec(qp->nu[ii]+qp->nx[ii], ws->dux+ii, 0);
for(ii=0; ii<qp->N; ii++)
	d_print_tran_strvec(qp->nx[ii+1], ws->dpi+ii, 0);
exit(1);
#endif
		// alpha
		COMPUTE_ALPHA_QP(cws);
		cws->stat[5*kk+3] = cws->alpha;

		//
		UPDATE_VAR_QP(cws);

		// compute residuals
		COMPUTE_RES_OCP_QP(qp, qp_sol, ws);
		cws->mu = ws->res_mu;
		cws->stat[5*kk+4] = ws->res_mu;

		}
	
	ws->iter = kk;
	
	// max iteration number reached
	if(kk==cws->iter_max)
		return 1;

	// normal return
	return 0;

	}



void CVT_OCP_QP_RES_TO_COLMAJ(struct OCP_QP *qp, struct OCP_QP_IPM_WORKSPACE *ws, REAL **res_r, REAL **res_q, REAL **res_ls, REAL **res_us, REAL **res_b, REAL **res_d_lb, REAL **res_d_ub, REAL **res_d_lg, REAL **res_d_ug, REAL **res_d_ls, REAL **res_d_us, REAL **res_m_lb, REAL **res_m_ub, REAL **res_m_lg, REAL **res_m_ug, REAL **res_m_ls, REAL **res_m_us)
	{

	int N = qp->N;
	int *nx = qp->nx;
	int *nu = qp->nu;
	int *nb = qp->nb;
	int *ng = qp->ng;
	int *ns = qp->ns;
	
	int ii;

	for(ii=0; ii<N; ii++)
		{
		CVT_STRVEC2VEC(nu[ii], ws->res_g+ii, 0, res_r[ii]);
		CVT_STRVEC2VEC(nx[ii], ws->res_g+ii, nu[ii], res_q[ii]);
		CVT_STRVEC2VEC(ns[ii], ws->res_g+ii, nu[ii]+nx[ii], res_ls[ii]);
		CVT_STRVEC2VEC(ns[ii], ws->res_g+ii, nu[ii]+nx[ii]+ns[ii], res_us[ii]);

		CVT_STRVEC2VEC(nx[ii+1], ws->res_b+ii, 0, res_b[ii]);
		CVT_STRVEC2VEC(nb[ii], ws->res_d+ii, 0, res_d_lb[ii]);
		CVT_STRVEC2VEC(ng[ii], ws->res_d+ii, nb[ii], res_d_lg[ii]);
		CVT_STRVEC2VEC(nb[ii], ws->res_d+ii, nb[ii]+ng[ii], res_d_ub[ii]);
		CVT_STRVEC2VEC(ng[ii], ws->res_d+ii, 2*nb[ii]+ng[ii], res_d_ug[ii]);
		CVT_STRVEC2VEC(ns[ii], ws->res_d+ii, 2*nb[ii]+2*ng[ii], res_d_ls[ii]);
		CVT_STRVEC2VEC(ns[ii], ws->res_d+ii, 2*nb[ii]+2*ng[ii]+ns[ii], res_d_us[ii]);
		CVT_STRVEC2VEC(nb[ii], ws->res_m+ii, 0, res_m_lb[ii]);
		CVT_STRVEC2VEC(ng[ii], ws->res_m+ii, nb[ii], res_m_lg[ii]);
		CVT_STRVEC2VEC(nb[ii], ws->res_m+ii, nb[ii]+ng[ii], res_m_ub[ii]);
		CVT_STRVEC2VEC(ng[ii], ws->res_m+ii, 2*nb[ii]+ng[ii], res_m_ug[ii]);
		CVT_STRVEC2VEC(ns[ii], ws->res_m+ii, 2*nb[ii]+2*ng[ii], res_m_ls[ii]);
		CVT_STRVEC2VEC(ns[ii], ws->res_m+ii, 2*nb[ii]+2*ng[ii]+ns[ii], res_m_us[ii]);
		}
	CVT_STRVEC2VEC(nu[ii], ws->res_g+ii, 0, res_r[ii]);
	CVT_STRVEC2VEC(nx[ii], ws->res_g+ii, nu[ii], res_q[ii]);
	CVT_STRVEC2VEC(ns[ii], ws->res_g+ii, nu[ii]+nx[ii], res_ls[ii]);
	CVT_STRVEC2VEC(ns[ii], ws->res_g+ii, nu[ii]+nx[ii]+ns[ii], res_us[ii]);

	CVT_STRVEC2VEC(nb[ii], ws->res_d+ii, 0, res_d_lb[ii]);
	CVT_STRVEC2VEC(ng[ii], ws->res_d+ii, nb[ii], res_d_lg[ii]);
	CVT_STRVEC2VEC(nb[ii], ws->res_d+ii, nb[ii]+ng[ii], res_d_ub[ii]);
	CVT_STRVEC2VEC(ng[ii], ws->res_d+ii, 2*nb[ii]+ng[ii], res_d_ug[ii]);
	CVT_STRVEC2VEC(ns[ii], ws->res_d+ii, 2*nb[ii]+2*ng[ii], res_d_ls[ii]);
	CVT_STRVEC2VEC(ns[ii], ws->res_d+ii, 2*nb[ii]+2*ng[ii]+ns[ii], res_d_us[ii]);
	CVT_STRVEC2VEC(nb[ii], ws->res_m+ii, 0, res_m_lb[ii]);
	CVT_STRVEC2VEC(ng[ii], ws->res_m+ii, nb[ii], res_m_lg[ii]);
	CVT_STRVEC2VEC(nb[ii], ws->res_m+ii, nb[ii]+ng[ii], res_m_ub[ii]);
	CVT_STRVEC2VEC(ng[ii], ws->res_m+ii, 2*nb[ii]+ng[ii], res_m_ug[ii]);
	CVT_STRVEC2VEC(ns[ii], ws->res_m+ii, 2*nb[ii]+2*ng[ii], res_m_ls[ii]);
	CVT_STRVEC2VEC(ns[ii], ws->res_m+ii, 2*nb[ii]+2*ng[ii]+ns[ii], res_m_us[ii]);

	return;

	}



void CVT_OCP_QP_RES_TO_ROWMAJ(struct OCP_QP *qp, struct OCP_QP_IPM_WORKSPACE *ws, REAL **res_r, REAL **res_q, REAL **res_ls, REAL **res_us, REAL **res_b, REAL **res_d_lb, REAL **res_d_ub, REAL **res_d_lg, REAL **res_d_ug, REAL **res_d_ls, REAL **res_d_us, REAL **res_m_lb, REAL **res_m_ub, REAL **res_m_lg, REAL **res_m_ug, REAL **res_m_ls, REAL **res_m_us)
	{

	int N = qp->N;
	int *nx = qp->nx;
	int *nu = qp->nu;
	int *nb = qp->nb;
	int *ng = qp->ng;
	int *ns = qp->ns;
	
	int ii;

	for(ii=0; ii<N; ii++)
		{
		CVT_STRVEC2VEC(nu[ii], ws->res_g+ii, 0, res_r[ii]);
		CVT_STRVEC2VEC(nx[ii], ws->res_g+ii, nu[ii], res_q[ii]);
		CVT_STRVEC2VEC(ns[ii], ws->res_g+ii, nu[ii]+nx[ii], res_ls[ii]);
		CVT_STRVEC2VEC(ns[ii], ws->res_g+ii, nu[ii]+nx[ii]+ns[ii], res_us[ii]);

		CVT_STRVEC2VEC(nx[ii+1], ws->res_b+ii, 0, res_b[ii]);
		CVT_STRVEC2VEC(nb[ii], ws->res_d+ii, 0, res_d_lb[ii]);
		CVT_STRVEC2VEC(ng[ii], ws->res_d+ii, nb[ii], res_d_lg[ii]);
		CVT_STRVEC2VEC(nb[ii], ws->res_d+ii, nb[ii]+ng[ii], res_d_ub[ii]);
		CVT_STRVEC2VEC(ng[ii], ws->res_d+ii, 2*nb[ii]+ng[ii], res_d_ug[ii]);
		CVT_STRVEC2VEC(ns[ii], ws->res_d+ii, 2*nb[ii]+2*ng[ii], res_d_ls[ii]);
		CVT_STRVEC2VEC(ns[ii], ws->res_d+ii, 2*nb[ii]+2*ng[ii]+ns[ii], res_d_us[ii]);
		CVT_STRVEC2VEC(nb[ii], ws->res_m+ii, 0, res_m_lb[ii]);
		CVT_STRVEC2VEC(ng[ii], ws->res_m+ii, nb[ii], res_m_lg[ii]);
		CVT_STRVEC2VEC(nb[ii], ws->res_m+ii, nb[ii]+ng[ii], res_m_ub[ii]);
		CVT_STRVEC2VEC(ng[ii], ws->res_m+ii, 2*nb[ii]+ng[ii], res_m_ug[ii]);
		CVT_STRVEC2VEC(ns[ii], ws->res_m+ii, 2*nb[ii]+2*ng[ii], res_m_ls[ii]);
		CVT_STRVEC2VEC(ns[ii], ws->res_m+ii, 2*nb[ii]+2*ng[ii]+ns[ii], res_m_us[ii]);
		}
	CVT_STRVEC2VEC(nu[ii], ws->res_g+ii, 0, res_r[ii]);
	CVT_STRVEC2VEC(nx[ii], ws->res_g+ii, nu[ii], res_q[ii]);
	CVT_STRVEC2VEC(ns[ii], ws->res_g+ii, nu[ii]+nx[ii], res_ls[ii]);
	CVT_STRVEC2VEC(ns[ii], ws->res_g+ii, nu[ii]+nx[ii]+ns[ii], res_us[ii]);

	CVT_STRVEC2VEC(nb[ii], ws->res_d+ii, 0, res_d_lb[ii]);
	CVT_STRVEC2VEC(ng[ii], ws->res_d+ii, nb[ii], res_d_lg[ii]);
	CVT_STRVEC2VEC(nb[ii], ws->res_d+ii, nb[ii]+ng[ii], res_d_ub[ii]);
	CVT_STRVEC2VEC(ng[ii], ws->res_d+ii, 2*nb[ii]+ng[ii], res_d_ug[ii]);
	CVT_STRVEC2VEC(ns[ii], ws->res_d+ii, 2*nb[ii]+2*ng[ii], res_d_ls[ii]);
	CVT_STRVEC2VEC(ns[ii], ws->res_d+ii, 2*nb[ii]+2*ng[ii]+ns[ii], res_d_us[ii]);
	CVT_STRVEC2VEC(nb[ii], ws->res_m+ii, 0, res_m_lb[ii]);
	CVT_STRVEC2VEC(ng[ii], ws->res_m+ii, nb[ii], res_m_lg[ii]);
	CVT_STRVEC2VEC(nb[ii], ws->res_m+ii, nb[ii]+ng[ii], res_m_ub[ii]);
	CVT_STRVEC2VEC(ng[ii], ws->res_m+ii, 2*nb[ii]+ng[ii], res_m_ug[ii]);
	CVT_STRVEC2VEC(ns[ii], ws->res_m+ii, 2*nb[ii]+2*ng[ii], res_m_ls[ii]);
	CVT_STRVEC2VEC(ns[ii], ws->res_m+ii, 2*nb[ii]+2*ng[ii]+ns[ii], res_m_us[ii]);

	return;

	}




