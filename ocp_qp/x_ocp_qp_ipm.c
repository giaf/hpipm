/**************************************************************************************************
*                                                                                                 *
* This file is part of HPIPM.                                                                     *
*                                                                                                 *
* HPIPM -- High-Performance Interior Point Method.                                                *
* Copyright (C) 2019 by Gianluca Frison.                                                          *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
*                                                                                                 *
* The 2-Clause BSD License                                                                        *
*                                                                                                 *
* Redistribution and use in source and binary forms, with or without                              *
* modification, are permitted provided that the following conditions are met:                     *
*                                                                                                 *
* 1. Redistributions of source code must retain the above copyright notice, this                  *
*    list of conditions and the following disclaimer.                                             *
* 2. Redistributions in binary form must reproduce the above copyright notice,                    *
*    this list of conditions and the following disclaimer in the documentation                    *
*    and/or other materials provided with the distribution.                                       *
*                                                                                                 *
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND                 *
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED                   *
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE                          *
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR                 *
* ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES                  *
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;                    *
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND                     *
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT                      *
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS                   *
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                                    *
*                                                                                                 *
* Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             *
*                                                                                                 *
**************************************************************************************************/



int OCP_QP_IPM_ARG_STRSIZE()
	{
	return sizeof(struct OCP_QP_IPM_ARG);
	}



int OCP_QP_IPM_ARG_MEMSIZE(struct OCP_QP_DIM *dim)
	{

	return 0;

	}



void OCP_QP_IPM_ARG_CREATE(struct OCP_QP_DIM *dim, struct OCP_QP_IPM_ARG *arg, void *mem)
	{

#if 0
	// loop index
	int ii;

	// zero memory (to avoid corrupted memory like e.g. NaN)
	int memsize = OCP_QP_IPM_ARG_MEMSIZE(dim);
	int memsize_m8 = memsize/8; // sizeof(double) is 8
//	int memsize_r8 = memsize - 8*memsize_m8;
	double *double_ptr = mem;
	// XXX exploit that it is multiple of 64 bytes !!!!!
	for(ii=0; ii<memsize_m8-7; ii+=8)
		{
		double_ptr[ii+0] = 0.0;
		double_ptr[ii+1] = 0.0;
		double_ptr[ii+2] = 0.0;
		double_ptr[ii+3] = 0.0;
		double_ptr[ii+4] = 0.0;
		double_ptr[ii+5] = 0.0;
		double_ptr[ii+6] = 0.0;
		double_ptr[ii+7] = 0.0;
		}
//	for(; ii<memsize_m8; ii++)
//		{
//		double_ptr[ii] = 0.0;
//		}
//	char *char_ptr = (char *) (&double_ptr[ii]);
//	for(ii=0; ii<memsize_r8; ii++)
//		{
//		char_ptr[ii] = 0;
//		}
#endif

	arg->memsize = 0;

	return;

	}



void OCP_QP_IPM_ARG_SET_DEFAULT(enum HPIPM_MODE mode, struct OCP_QP_IPM_ARG *arg)
	{

	if(mode==SPEED_ABS)
		{
		arg->mu0 = 1e1;
		arg->alpha_min = 1e-12;
		arg->res_g_max = 1e0; // not used
		arg->res_b_max = 1e0; // not used
		arg->res_d_max = 1e0; // not used
		arg->res_m_max = 1e-8;
		arg->iter_max = 15;
		arg->stat_max = 15;
		arg->pred_corr = 1;
		arg->cond_pred_corr = 0; // not used
		arg->itref_pred_max = 0; // not used
		arg->itref_corr_max = 0; // not used
		arg->reg_prim = 1e-15;
		arg->square_root_alg = 1;
		arg->lq_fact = 0; // not used
		arg->lam_min = 1e-30;
		arg->t_min = 1e-30;
		arg->warm_start = 0;
		arg->abs_form = 1;
		arg->comp_dual_sol = 0;
		arg->comp_res_exit = 0;
		arg->comp_res_pred = 0;
		arg->mode = mode;
		}
	else if(mode==SPEED)
		{
		arg->mu0 = 1e1;
		arg->alpha_min = 1e-12;
		arg->res_g_max = 1e-6;
		arg->res_b_max = 1e-8;
		arg->res_d_max = 1e-8;
		arg->res_m_max = 1e-8;
		arg->iter_max = 15;
		arg->stat_max = 15;
		arg->pred_corr = 1;
		arg->cond_pred_corr = 1;
		arg->itref_pred_max = 0;
		arg->itref_corr_max = 0;
		arg->reg_prim = 1e-15;
		arg->square_root_alg = 1;
		arg->lq_fact = 0;
		arg->lam_min = 1e-30;
		arg->t_min = 1e-30;
		arg->warm_start = 0;
		arg->abs_form = 0;
		arg->comp_dual_sol = 1;
		arg->comp_res_exit = 1;
		arg->comp_res_pred = 1;
		arg->mode = mode;
		}
	else if(mode==BALANCE)
		{
		arg->mu0 = 1e1;
		arg->alpha_min = 1e-12;
		arg->res_g_max = 1e-6;
		arg->res_b_max = 1e-8;
		arg->res_d_max = 1e-8;
		arg->res_m_max = 1e-8;
		arg->iter_max = 30;
		arg->stat_max = 30;
		arg->pred_corr = 1;
		arg->cond_pred_corr = 1;
		arg->itref_pred_max = 0;
		arg->itref_corr_max = 2;
		arg->reg_prim = 1e-15;
		arg->square_root_alg = 1;
		arg->lq_fact = 1;
		arg->lam_min = 1e-30;
		arg->t_min = 1e-30;
		arg->warm_start = 0;
		arg->abs_form = 0;
		arg->comp_dual_sol = 1;
		arg->comp_res_exit = 1;
		arg->comp_res_pred = 1;
		arg->mode = mode;
		}
	else if(mode==ROBUST)
		{
		arg->mu0 = 1e2;
		arg->alpha_min = 1e-12;
		arg->res_g_max = 1e-6;
		arg->res_b_max = 1e-8;
		arg->res_d_max = 1e-8;
		arg->res_m_max = 1e-8;
		arg->iter_max = 100;
		arg->stat_max = 100;
		arg->pred_corr = 1;
		arg->cond_pred_corr = 1;
		arg->itref_pred_max = 0;
		arg->itref_corr_max = 4;
		arg->reg_prim = 1e-15;
		arg->square_root_alg = 1;
		arg->lq_fact = 2;
		arg->lam_min = 1e-30;
		arg->t_min = 1e-30;
		arg->warm_start = 0;
		arg->abs_form = 0;
		arg->comp_dual_sol = 1;
		arg->comp_res_exit = 1;
		arg->comp_res_pred = 1;
		arg->mode = mode;
		}
	else
		{
		printf("\nerror: OCP_QP_IPM_ARG_SET_DEFAULT: wrong set default mode\n");
		exit(1);
		}

	return;

	}



void OCP_QP_IPM_ARG_SET(char *field, void *value, struct OCP_QP_IPM_ARG *arg)
	{
	if(hpipm_strcmp(field, "iter_max")) 
		{
		OCP_QP_IPM_ARG_SET_ITER_MAX(value, arg);
		}
	else if(hpipm_strcmp(field, "alpha_min")) 
		{
		OCP_QP_IPM_ARG_SET_ALPHA_MIN(value, arg);
		}
	else if(hpipm_strcmp(field, "mu0")) 
		{
		OCP_QP_IPM_ARG_SET_MU0(value, arg);
		}
	else if(hpipm_strcmp(field, "tol_stat")) 
		{
		OCP_QP_IPM_ARG_SET_TOL_STAT(value, arg);
		}
	else if(hpipm_strcmp(field, "tol_eq")) 
		{
		OCP_QP_IPM_ARG_SET_TOL_EQ(value, arg);
		}
	else if(hpipm_strcmp(field, "tol_ineq")) 
		{
		OCP_QP_IPM_ARG_SET_TOL_INEQ(value, arg);
		}
	else if(hpipm_strcmp(field, "tol_comp")) 
		{
		OCP_QP_IPM_ARG_SET_TOL_COMP(value, arg);
		}
	else if(hpipm_strcmp(field, "reg_prim")) 
		{
		OCP_QP_IPM_ARG_SET_REG_PRIM(value, arg);
		}
	else if(hpipm_strcmp(field, "warm_start")) 
		{
		OCP_QP_IPM_ARG_SET_WARM_START(value, arg);
		}
	else if(hpipm_strcmp(field, "pred_corr")) 
		{
		OCP_QP_IPM_ARG_SET_PRED_CORR(value, arg);
		}
	else if(hpipm_strcmp(field, "ric_alg")) 
		{
		OCP_QP_IPM_ARG_SET_RIC_ALG(value, arg);
		}
	else if(hpipm_strcmp(field, "comp_res_pred")) 
		{
		OCP_QP_IPM_ARG_SET_COMP_RES_PRED(value, arg);
		}
	else
		{
		printf("error: OCP_QP_IPM_ARG_SET: wrong field %s\n", field);
		exit(1);	
		}
	return;
	}



void OCP_QP_IPM_ARG_SET_ITER_MAX(int *iter_max, struct OCP_QP_IPM_ARG *arg)
	{
	arg->iter_max = *iter_max;
	return;
	}



void OCP_QP_IPM_ARG_SET_ALPHA_MIN(REAL *alpha_min, struct OCP_QP_IPM_ARG *arg)
	{
	arg->alpha_min = *alpha_min;
	return;
	}



void OCP_QP_IPM_ARG_SET_MU0(REAL *mu0, struct OCP_QP_IPM_ARG *arg)
	{
	arg->mu0 = *mu0;
	return;
	}



void OCP_QP_IPM_ARG_SET_TOL_STAT(REAL *tol_stat, struct OCP_QP_IPM_ARG *arg)
	{
	arg->res_g_max = *tol_stat;
	return;
	}



void OCP_QP_IPM_ARG_SET_TOL_EQ(REAL *tol_eq, struct OCP_QP_IPM_ARG *arg)
	{
	arg->res_b_max = *tol_eq;
	return;
	}



void OCP_QP_IPM_ARG_SET_TOL_INEQ(REAL *tol_ineq, struct OCP_QP_IPM_ARG *arg)
	{
	arg->res_d_max = *tol_ineq;
	return;
	}



void OCP_QP_IPM_ARG_SET_TOL_COMP(REAL *tol_comp, struct OCP_QP_IPM_ARG *arg)
	{
	arg->res_m_max = *tol_comp;
	return;
	}



void OCP_QP_IPM_ARG_SET_REG_PRIM(REAL *reg, struct OCP_QP_IPM_ARG *arg)
	{
	arg->reg_prim = *reg;
	return;
	}



void OCP_QP_IPM_ARG_SET_WARM_START(int *warm_start, struct OCP_QP_IPM_ARG *arg)
	{
	arg->warm_start = *warm_start;
	return;
	}



void OCP_QP_IPM_ARG_SET_PRED_CORR(int *pred_corr, struct OCP_QP_IPM_ARG *arg)
	{
	arg->pred_corr = *pred_corr;
	return;
	}



void OCP_QP_IPM_ARG_SET_RIC_ALG(int *ric_alg, struct OCP_QP_IPM_ARG *arg)
	{
	arg->square_root_alg = *ric_alg;
	return;
	}



void OCP_QP_IPM_ARG_SET_COMP_RES_PRED(int *comp_res_pred, struct OCP_QP_IPM_ARG *arg)
	{
	arg->comp_res_pred = *comp_res_pred;
	return;
	}



int OCP_QP_IPM_WS_STRSIZE()
	{
	return sizeof(struct OCP_QP_IPM_WS);
	}



int OCP_QP_IPM_WS_MEMSIZE(struct OCP_QP_DIM *dim, struct OCP_QP_IPM_ARG *arg)
	{

	// stat_max is at least as big as iter_max
	if(arg->iter_max > arg->stat_max)
		arg->stat_max = arg->iter_max;

	// loop index
	int ii;

	// extract ocp qp size
	int N = dim->N;
	int *nx = dim->nx;
	int *nu = dim->nu;
	int *nb = dim->nb;
	int *ng = dim->ng;
	int *ns = dim->ns;

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
		nct += 2*nb[ii]+2*ng[ii]+2*ns[ii];
		nxM = nx[ii]>nxM ? nx[ii] : nxM;
		nuM = nu[ii]>nuM ? nu[ii] : nuM;
		nbM = nb[ii]>nbM ? nb[ii] : nbM;
		ngM = ng[ii]>ngM ? ng[ii] : ngM;
		nsM = ns[ii]>nsM ? ns[ii] : nsM;
		}
	nvt += nx[ii]+nu[ii]+2*ns[ii];
	nct += 2*nb[ii]+2*ng[ii]+2*ns[ii];
	nxM = nx[ii]>nxM ? nx[ii] : nxM;
	nuM = nu[ii]>nuM ? nu[ii] : nuM;
	nbM = nb[ii]>nbM ? nb[ii] : nbM;
	ngM = ng[ii]>ngM ? ng[ii] : ngM;
	nsM = ns[ii]>nsM ? ns[ii] : nsM;

	int size = 0;

	size += 1*sizeof(struct CORE_QP_IPM_WORKSPACE);
	size += 1*MEMSIZE_CORE_QP_IPM(nvt, net, nct);

	size += 1*sizeof(struct OCP_QP_RES_WS); // res_workspace

	size += 2*sizeof(struct OCP_QP); // qp_step qp_itref

	size += 2*sizeof(struct OCP_QP_SOL); // sol_step sol_itref
	size += 1*OCP_QP_SOL_MEMSIZE(dim); // sol_itref

	size += 2*sizeof(struct OCP_QP_RES); // res res_itref
	size += 1*OCP_QP_RES_MEMSIZE(dim); // res_itref

	size += 9*(N+1)*sizeof(struct STRVEC); // res_g res_d res_m Gamma gamma Zs_inv sol_step(v,lam,t)
	size += 3*N*sizeof(struct STRVEC); // res_b Pb sol_step(pi) 
	size += 10*sizeof(struct STRVEC); // tmp_nxM (4+2)*tmp_nbgM (1+1)*tmp_nsM tmp_m

	size += 1*(N+1)*sizeof(struct STRMAT); // L
	if(!arg->square_root_alg)
		{
		size += 1*(N+1)*sizeof(struct STRMAT); // P
		size += 1*sizeof(struct STRMAT); // Ls
		}
	if(arg->lq_fact>0)
		{
		size += 1*(N+1)*sizeof(struct STRMAT); // Lh
		}
	size += 2*sizeof(struct STRMAT); // AL
	if(arg->lq_fact>0)
		{
		size += 1*sizeof(struct STRMAT); // lq0
		}

	size += 1*SIZE_STRVEC(nxM); // tmp_nxM
	size += 4*SIZE_STRVEC(nbM+ngM); // tmp_nbgM
	size += 1*SIZE_STRVEC(nsM); // tmp_nsM
	for(ii=0; ii<N; ii++) size += 1*SIZE_STRVEC(nx[ii+1]); // Pb
	for(ii=0; ii<=N; ii++) size += 1*SIZE_STRVEC(2*ns[ii]); // Zs_inv

	for(ii=0; ii<=N; ii++) size += 1*SIZE_STRMAT(nu[ii]+nx[ii]+1, nu[ii]+nx[ii]); // L
	if(!arg->square_root_alg)
		{
		for(ii=0; ii<=N; ii++) size += 1*SIZE_STRMAT(nx[ii]+1, nx[ii]); // P
		size += 1*SIZE_STRMAT(nxM+1, nuM); // Ls
		}
	if(arg->lq_fact>0)
		{
		for(ii=0; ii<=N; ii++) size += 1*SIZE_STRMAT(nu[ii]+nx[ii]+1, nu[ii]+nx[ii]); // Lh
		}
	size += 2*SIZE_STRMAT(nuM+nxM+1, nxM+ngM); // AL
	if(arg->lq_fact>0)
		{
		size += 1*SIZE_STRMAT(nuM+nxM, 2*nuM+3*nxM+ngM); // lq0
		}
	size += 1*SIZE_STRVEC(nct); // tmp_m

	if(arg->lq_fact>0)
		{
		size += 1*GELQF_WORKSIZE(nuM+nxM, 2*nuM+3*nxM+ngM); // lq_work0
		}

	size += 9*(1+arg->stat_max)*sizeof(REAL); // stat

	size += (N+1)*sizeof(int); // use_hess_fact

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

	}



void OCP_QP_IPM_WS_CREATE(struct OCP_QP_DIM *dim, struct OCP_QP_IPM_ARG *arg, struct OCP_QP_IPM_WS *workspace, void *mem)
	{

	// loop index
	int ii;

	// zero memory (to avoid corrupted memory like e.g. NaN)
	int memsize = OCP_QP_IPM_WS_MEMSIZE(dim, arg);
	int memsize_m8 = memsize/8; // sizeof(double) is 8
//	int memsize_r8 = memsize - 8*memsize_m8;
	double *double_ptr = mem;
	// XXX exploit that it is multiple of 64 bytes !!!!!
	for(ii=0; ii<memsize_m8-7; ii+=8)
		{
		double_ptr[ii+0] = 0.0;
		double_ptr[ii+1] = 0.0;
		double_ptr[ii+2] = 0.0;
		double_ptr[ii+3] = 0.0;
		double_ptr[ii+4] = 0.0;
		double_ptr[ii+5] = 0.0;
		double_ptr[ii+6] = 0.0;
		double_ptr[ii+7] = 0.0;
		}
//	for(; ii<memsize_m8; ii++)
//		{
//		double_ptr[ii] = 0.0;
//		}
//	char *char_ptr = (char *) (&double_ptr[ii]);
//	for(ii=0; ii<memsize_r8; ii++)
//		{
//		char_ptr[ii] = 0;
//		}

	// extract ocp qp size
	int N = dim->N;
	int *nx = dim->nx;
	int *nu = dim->nu;
	int *nb = dim->nb;
	int *ng = dim->ng;
	int *ns = dim->ns;


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
		nct += 2*nb[ii]+2*ng[ii]+2*ns[ii];
		nxM = nx[ii]>nxM ? nx[ii] : nxM;
		nuM = nu[ii]>nuM ? nu[ii] : nuM;
		nbM = nb[ii]>nbM ? nb[ii] : nbM;
		ngM = ng[ii]>ngM ? ng[ii] : ngM;
		nsM = ns[ii]>nsM ? ns[ii] : nsM;
		}
	nvt += nx[ii]+nu[ii]+2*ns[ii];
	nct += 2*nb[ii]+2*ng[ii]+2*ns[ii];
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


	// res struct
	struct OCP_QP_RES *res_ptr = (struct OCP_QP_RES *) sr_ptr;
	workspace->res = res_ptr;
	res_ptr += 1;
	workspace->res_itref = res_ptr;
	res_ptr += 1;


	// res workspace struct
	struct OCP_QP_RES_WS *res_ws_ptr = (struct OCP_QP_RES_WS *) res_ptr;
	workspace->res_workspace = res_ws_ptr;
	res_ws_ptr += 1;


	// qp sol struct
	struct OCP_QP_SOL *qp_sol_ptr = (struct OCP_QP_SOL *) res_ws_ptr;

	workspace->sol_step = qp_sol_ptr;
	qp_sol_ptr += 1;
	workspace->sol_itref = qp_sol_ptr;
	qp_sol_ptr += 1;


	// qp struct
	struct OCP_QP *qp_ptr = (struct OCP_QP *) qp_sol_ptr;

	workspace->qp_step = qp_ptr;
	qp_ptr += 1;
	workspace->qp_itref = qp_ptr;
	qp_ptr += 1;


	// matrix struct
	struct STRMAT *sm_ptr = (struct STRMAT *) qp_ptr;

	workspace->L = sm_ptr;
	sm_ptr += N+1;
	if(!arg->square_root_alg)
		{
		workspace->P = sm_ptr;
		sm_ptr += N+1;
		workspace->Ls = sm_ptr;
		sm_ptr += 1;
		}
	if(arg->lq_fact>0)
		{
		workspace->Lh = sm_ptr;
		sm_ptr += N+1;
		}
	workspace->AL = sm_ptr;
	sm_ptr += 2;
	if(arg->lq_fact>0)
		{
		workspace->lq0 = sm_ptr;
		sm_ptr += 1;
		}


	// vector struct
	struct STRVEC *sv_ptr = (struct STRVEC *) sm_ptr;

	workspace->sol_step->ux = sv_ptr;
	sv_ptr += N+1;
	workspace->sol_step->pi = sv_ptr;
	sv_ptr += N;
	workspace->sol_step->lam = sv_ptr;
	sv_ptr += N+1;
	workspace->sol_step->t = sv_ptr;
	sv_ptr += N+1;
	workspace->res->res_g = sv_ptr;
	sv_ptr += N+1;
	workspace->res->res_b = sv_ptr;
	sv_ptr += N;
	workspace->res->res_d = sv_ptr;
	sv_ptr += N+1;
	workspace->res->res_m = sv_ptr;
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
	workspace->res_workspace->tmp_nbgM = sv_ptr;
	sv_ptr += 2;
	workspace->tmp_nsM = sv_ptr;
	sv_ptr += 1;
	workspace->res_workspace->tmp_nsM = sv_ptr;
	sv_ptr += 1;
	workspace->tmp_m = sv_ptr;
	sv_ptr += 1;


	// double/float stuff
	REAL *d_ptr = (REAL *) sv_ptr;

	workspace->stat = d_ptr;
	d_ptr += 9*(1+arg->stat_max);

	// int stuff
	int *i_ptr = (int *) d_ptr;

	workspace->use_hess_fact = i_ptr;
	i_ptr += N+1;


	// align to typicl cache line size
	size_t s_ptr = (size_t) i_ptr;
	s_ptr = (s_ptr+63)/64*64;


	// void stuf
	char *c_ptr = (char *) s_ptr;

	OCP_QP_SOL_CREATE(dim, workspace->sol_itref, c_ptr);
	c_ptr += workspace->sol_itref->memsize;

	OCP_QP_RES_CREATE(dim, workspace->res_itref, c_ptr);
	c_ptr += workspace->res_itref->memsize;

	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRMAT(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], workspace->L+ii, c_ptr);
		c_ptr += (workspace->L+ii)->memsize;
		}
	if(!arg->square_root_alg)
		{
		for(ii=0; ii<=N; ii++)
			{
			CREATE_STRMAT(nx[ii]+1, nx[ii], workspace->P+ii, c_ptr);
			c_ptr += (workspace->P+ii)->memsize;
			}
		CREATE_STRMAT(nxM+1, nuM, workspace->Ls, c_ptr);
		c_ptr += (workspace->Ls)->memsize;
		}

	if(arg->lq_fact>0)
		{
		for(ii=0; ii<=N; ii++)
			{
			CREATE_STRMAT(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], workspace->Lh+ii, c_ptr);
			c_ptr += (workspace->Lh+ii)->memsize;
			}
		}

	CREATE_STRMAT(nuM+nxM+1, nxM+ngM, workspace->AL+0, c_ptr);
	c_ptr += (workspace->AL+0)->memsize;

	CREATE_STRMAT(nuM+nxM+1, nxM+ngM, workspace->AL+1, c_ptr);
	c_ptr += (workspace->AL+1)->memsize;

	if(arg->lq_fact>0)
		{
		CREATE_STRMAT(nuM+nxM, 2*nuM+3*nxM+ngM, workspace->lq0, c_ptr);
		c_ptr += (workspace->lq0)->memsize;
		}

	for(ii=0; ii<N; ii++)
		{
		CREATE_STRVEC(nx[ii+1], workspace->Pb+ii, c_ptr);
		c_ptr += (workspace->Pb+ii)->memsize;
		}

	for(ii=0; ii<N+1; ii++)
		{
		CREATE_STRVEC(2*ns[ii], workspace->Zs_inv+ii, c_ptr);
		c_ptr += (workspace->Zs_inv+ii)->memsize;
		}

	CREATE_STRVEC(nxM, workspace->tmp_nxM, c_ptr);
	c_ptr += workspace->tmp_nxM->memsize;

	CREATE_STRVEC(nbM+ngM, workspace->tmp_nbgM+0, c_ptr);
	CREATE_STRVEC(nbM+ngM, workspace->res_workspace->tmp_nbgM+0, c_ptr);
	c_ptr += (workspace->tmp_nbgM+0)->memsize;

	CREATE_STRVEC(nbM+ngM, workspace->tmp_nbgM+1, c_ptr);
	CREATE_STRVEC(nbM+ngM, workspace->res_workspace->tmp_nbgM+1, c_ptr);
	c_ptr += (workspace->tmp_nbgM+1)->memsize;

	CREATE_STRVEC(nbM+ngM, workspace->tmp_nbgM+2, c_ptr);
	c_ptr += (workspace->tmp_nbgM+2)->memsize;

	CREATE_STRVEC(nbM+ngM, workspace->tmp_nbgM+3, c_ptr);
	c_ptr += (workspace->tmp_nbgM+3)->memsize;

	CREATE_STRVEC(nsM, workspace->tmp_nsM+0, c_ptr);
	CREATE_STRVEC(nsM, workspace->res_workspace->tmp_nsM+0, c_ptr);
	c_ptr += (workspace->tmp_nsM+0)->memsize;

	CREATE_STRVEC(nct, workspace->tmp_m, c_ptr);
	c_ptr += SIZE_STRVEC(nct);

	CREATE_CORE_QP_IPM(nvt, net, nct, cws, c_ptr);
	c_ptr += workspace->core_workspace->memsize;

	if(arg->lq_fact>0)
		{
		workspace->lq_work0 = c_ptr;
		c_ptr += GELQF_WORKSIZE(nuM+nxM, 2*nuM+3*nxM+ngM);
		}


	// alias members of workspace and core_workspace
	//
	c_ptr = (char *) cws->dv;
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(nu[ii]+nx[ii]+2*ns[ii], workspace->sol_step->ux+ii, c_ptr);
		c_ptr += (nu[ii]+nx[ii])*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		}
	//
	c_ptr = (char *) cws->dpi;
	for(ii=0; ii<N; ii++)
		{
		CREATE_STRVEC(nx[ii+1], workspace->sol_step->pi+ii, c_ptr);
		c_ptr += (nx[ii+1])*sizeof(REAL);
		}
	//
	c_ptr = (char *) cws->dlam;
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(2*nb[ii]+2*ng[ii]+2*ns[ii], workspace->sol_step->lam+ii, c_ptr);
		c_ptr += nb[ii]*sizeof(REAL);
		c_ptr += ng[ii]*sizeof(REAL);
		c_ptr += nb[ii]*sizeof(REAL);
		c_ptr += ng[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		}
	//
	c_ptr = (char *) cws->dt;
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(2*nb[ii]+2*ng[ii]+2*ns[ii], workspace->sol_step->t+ii, c_ptr);
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
		CREATE_STRVEC(nu[ii]+nx[ii]+2*ns[ii], workspace->res->res_g+ii, c_ptr);
		c_ptr += (nu[ii]+nx[ii])*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		}
	//
	c_ptr = (char *) cws->res_b;
	for(ii=0; ii<N; ii++)
		{
		CREATE_STRVEC(nx[ii+1], workspace->res->res_b+ii, c_ptr);
		c_ptr += (nx[ii+1])*sizeof(REAL);
		}
	//
	c_ptr = (char *) cws->res_d;
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(2*nb[ii]+2*ng[ii]+2*ns[ii], workspace->res->res_d+ii, c_ptr);
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
		CREATE_STRVEC(2*nb[ii]+2*ng[ii]+2*ns[ii], workspace->res->res_m+ii, c_ptr);
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
		CREATE_STRVEC(2*nb[ii]+2*ng[ii]+2*ns[ii], workspace->Gamma+ii, c_ptr);
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
		CREATE_STRVEC(2*nb[ii]+2*ng[ii]+2*ns[ii], workspace->gamma+ii, c_ptr);
		c_ptr += nb[ii]*sizeof(REAL);
		c_ptr += ng[ii]*sizeof(REAL);
		c_ptr += nb[ii]*sizeof(REAL);
		c_ptr += ng[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		}


	workspace->res->dim = dim;

	workspace->stat_max = arg->stat_max;

	workspace->stat_m = 9;

	for(ii=0; ii<=N; ii++)
		workspace->use_hess_fact[ii] = 0;
	
	workspace->use_Pb = 0;

	workspace->memsize = memsize; //OCP_QP_IPM_WS_MEMSIZE(dim, arg);


#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + workspace->memsize)
		{
		printf("\nCreate_ocp_qp_ipm: outside memory bounds!\n\n");
		exit(1);
		}
#endif


	return;

	}


void OCP_QP_IPM_GET(char *field, struct OCP_QP_IPM_WS *ws, void *value)
	{
	if(hpipm_strcmp(field, "status"))
		{ 
		OCP_QP_IPM_GET_STATUS(ws, value);
		}
	else if(hpipm_strcmp(field, "iter"))
		{ 
		OCP_QP_IPM_GET_ITER(ws, value);
		}
	else if(hpipm_strcmp(field, "max_res_stat"))
		{ 
		OCP_QP_IPM_GET_MAX_RES_STAT(ws, value);
		}
	else if(hpipm_strcmp(field, "max_res_eq"))
		{ 
		OCP_QP_IPM_GET_MAX_RES_EQ(ws, value);
		}
	else if(hpipm_strcmp(field, "max_res_ineq"))
		{ 
		OCP_QP_IPM_GET_MAX_RES_INEQ(ws, value);
		}
	else if(hpipm_strcmp(field, "max_res_comp"))
		{ 
		OCP_QP_IPM_GET_MAX_RES_COMP(ws, value);
		}
	else if(hpipm_strcmp(field, "stat"))
		{ 
		OCP_QP_IPM_GET_STAT(ws, value);
		}
	else if(hpipm_strcmp(field, "stat_m"))
		{ 
		OCP_QP_IPM_GET_STAT_M(ws, value);
		}
	else 
		{
		printf("error: OCP_QP_IPM_GET: wrong field %s\n", field);
		exit(1);
		}
	return;
	}



void OCP_QP_IPM_GET_STATUS(struct OCP_QP_IPM_WS *ws, int *status)
	{
	*status = ws->status;
	return;
	}



void OCP_QP_IPM_GET_ITER(struct OCP_QP_IPM_WS *ws, int *iter)
	{
	*iter = ws->iter;
	return;
	}



void OCP_QP_IPM_GET_MAX_RES_STAT(struct OCP_QP_IPM_WS *ws, REAL *res_stat)
	{
	*res_stat = ws->qp_res[0];
	return;
	}



void OCP_QP_IPM_GET_MAX_RES_EQ(struct OCP_QP_IPM_WS *ws, REAL *res_eq)
	{
	*res_eq = ws->qp_res[1];
	return;
	}



void OCP_QP_IPM_GET_MAX_RES_INEQ(struct OCP_QP_IPM_WS *ws, REAL *res_ineq)
	{
	*res_ineq = ws->qp_res[2];
	return;
	}



void OCP_QP_IPM_GET_MAX_RES_COMP(struct OCP_QP_IPM_WS *ws, REAL *res_comp)
	{
	*res_comp = ws->qp_res[3];
	return;
	}



void OCP_QP_IPM_GET_STAT(struct OCP_QP_IPM_WS *ws, REAL **stat)
	{
	*stat = ws->stat;
	}



void OCP_QP_IPM_GET_STAT_M(struct OCP_QP_IPM_WS *ws, int *stat_m)
	{
	*stat_m = ws->stat_m;
	}



void OCP_QP_IPM_SOLVE(struct OCP_QP *qp, struct OCP_QP_SOL *qp_sol, struct OCP_QP_IPM_ARG *arg, struct OCP_QP_IPM_WS *ws)
	{

#if 0
	OCP_QP_DIM_PRINT(qp->dim);
	OCP_QP_PRINT(qp->dim, qp);
#endif

	struct CORE_QP_IPM_WORKSPACE *cws = ws->core_workspace;

	// arg to core workspace
	cws->lam_min = arg->lam_min;
	cws->t_min = arg->t_min;

	// alias qp vectors into qp_sol
	cws->v = qp_sol->ux->pa;
	cws->pi = qp_sol->pi->pa;
	cws->lam = qp_sol->lam->pa;
	cws->t = qp_sol->t->pa;

	// alias members of qp_step
	ws->qp_step->dim = qp->dim;
	ws->qp_step->RSQrq = qp->RSQrq;
	ws->qp_step->BAbt = qp->BAbt;
	ws->qp_step->DCt = qp->DCt;
	ws->qp_step->Z = qp->Z;
	ws->qp_step->idxb = qp->idxb;
	ws->qp_step->idxs = qp->idxs;
	ws->qp_step->rqz = ws->res->res_g;
	ws->qp_step->b = ws->res->res_b;
	ws->qp_step->d = ws->res->res_d;
	ws->qp_step->m = ws->res->res_m;

	// alias members of qp_itref
	ws->qp_itref->dim = qp->dim;
	ws->qp_itref->RSQrq = qp->RSQrq;
	ws->qp_itref->BAbt = qp->BAbt;
	ws->qp_itref->DCt = qp->DCt;
	ws->qp_itref->Z = qp->Z;
	ws->qp_itref->idxb = qp->idxb;
	ws->qp_itref->idxs = qp->idxs;
	ws->qp_itref->rqz = ws->res_itref->res_g;
	ws->qp_itref->b = ws->res_itref->res_b;
	ws->qp_itref->d = ws->res_itref->res_d;
	ws->qp_itref->m = ws->res_itref->res_m;

	// blasfeo alias for residuals
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

	REAL *qp_res = ws->qp_res;
	qp_res[0] = 0;
	qp_res[1] = 0;
	qp_res[2] = 0;
	qp_res[3] = 0;

	// no constraints
	if(cws->nc==0)
		{
		FACT_SOLVE_KKT_UNCONSTR_OCP_QP(qp, qp_sol, arg, ws);
		OCP_QP_RES_COMPUTE(qp, qp_sol, ws->res, ws->res_workspace);
		// compute infinity norm of residuals
		VECNRM_INF(cws->nv, &str_res_g, 0, &qp_res[0]);
		VECNRM_INF(cws->ne, &str_res_b, 0, &qp_res[1]);
		VECNRM_INF(cws->nc, &str_res_d, 0, &qp_res[2]);
		VECNRM_INF(cws->nc, &str_res_m, 0, &qp_res[3]);
		ws->stat[5] = qp_res[0];
		ws->stat[6] = qp_res[1];
		ws->stat[7] = qp_res[2];
		ws->stat[8] = qp_res[3];
		cws->mu = ws->res->res_mu;
		ws->iter = 0;
		ws->status = 0;
		return;
		}

	int N = qp->dim->N;
	int *nx = qp->dim->nx;
	int *nu = qp->dim->nu;
	int *nb = qp->dim->nb;
	int *ng = qp->dim->ng;
	int *ns = qp->dim->ns;

	int kk, ii, itref0=0, itref1=0, iter_ref_step;
	REAL tmp;
	REAL mu_aff0, mu;

	// init solver
	INIT_VAR_OCP_QP(qp, qp_sol, arg, ws);

	cws->alpha = 1.0;



	// absolute IPM formulation

	if(arg->abs_form)
		{

		// alias members of qp_step
		ws->qp_step->dim = qp->dim;
		ws->qp_step->RSQrq = qp->RSQrq;
		ws->qp_step->BAbt = qp->BAbt;
		ws->qp_step->DCt = qp->DCt;
		ws->qp_step->Z = qp->Z;
		ws->qp_step->idxb = qp->idxb;
		ws->qp_step->idxs = qp->idxs;
		ws->qp_step->rqz = qp->rqz;
		ws->qp_step->b = qp->b;
		ws->qp_step->d = qp->d;
		ws->qp_step->m = ws->tmp_m;

		// alias core workspace
		cws->res_m = ws->qp_step->m->pa;
		cws->res_m_bkp = ws->qp_step->m->pa;

		mu = VECMULDOT(cws->nc, qp_sol->lam, 0, qp_sol->t, 0, ws->tmp_m, 0);
		mu /= cws->nc;
		cws->mu = mu;

		// IPM loop (absolute formulation)
		for(kk=0; \
				kk<arg->iter_max & \
				cws->alpha>arg->alpha_min & \
				mu>arg->res_m_max; kk++)
			{

			VECSC(cws->nc, -1.0, ws->tmp_m, 0);

			// fact solve
			FACT_SOLVE_KKT_STEP_OCP_QP(ws->qp_step, ws->sol_step, arg, ws);
//blasfeo_print_tran_dvec(cws->nv, ws->sol_step->ux, 0);

			// compute step
			AXPY(cws->nv, -1.0, qp_sol->ux, 0, ws->sol_step->ux, 0, ws->sol_step->ux, 0);
			AXPY(cws->ne, -1.0, qp_sol->pi, 0, ws->sol_step->pi, 0, ws->sol_step->pi, 0);
			AXPY(cws->nc, -1.0, qp_sol->lam, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
			AXPY(cws->nc, -1.0, qp_sol->t, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);

			// alpha
			COMPUTE_ALPHA_QP(cws);
			if(kk<ws->stat_max)
				ws->stat[ws->stat_m*(kk+1)+0] = cws->alpha;

			// Mehrotra's predictor-corrector
			if(arg->pred_corr==1)
				{
				// mu_aff
				COMPUTE_MU_AFF_QP(cws);
				if(kk<ws->stat_max)
					ws->stat[ws->stat_m*(kk+1)+1] = cws->mu_aff;

				tmp = cws->mu_aff/cws->mu;
				cws->sigma = tmp*tmp*tmp;
				if(kk<ws->stat_max)
					ws->stat[ws->stat_m*(kk+1)+2] = cws->sigma;

				COMPUTE_CENTERING_CORRECTION_QP(cws);

				// fact and solve kkt
				ws->use_Pb = 1;
				SOLVE_KKT_STEP_OCP_QP(ws->qp_step, ws->sol_step, arg, ws);

				// compute step
				AXPY(cws->nv, -1.0, qp_sol->ux, 0, ws->sol_step->ux, 0, ws->sol_step->ux, 0);
				AXPY(cws->ne, -1.0, qp_sol->pi, 0, ws->sol_step->pi, 0, ws->sol_step->pi, 0);
				AXPY(cws->nc, -1.0, qp_sol->lam, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
				AXPY(cws->nc, -1.0, qp_sol->t, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);

				// alpha
				COMPUTE_ALPHA_QP(cws);
				if(kk<ws->stat_max)
					ws->stat[ws->stat_m*(kk+1)+3] = cws->alpha;

				}

			//
			UPDATE_VAR_QP(cws);

			// compute mu
			mu = VECMULDOT(cws->nc, qp_sol->lam, 0, qp_sol->t, 0, ws->tmp_m, 0);
			mu /= cws->nc;
			cws->mu = mu;
			if(kk<ws->stat_max)
				ws->stat[ws->stat_m*(kk+1)+4] = mu;

	//		exit(1);

			}

		if(arg->comp_res_exit & arg->comp_dual_sol)
			{
			// compute residuals
			OCP_QP_RES_COMPUTE(qp, qp_sol, ws->res, ws->res_workspace);

			// compute infinity norm of residuals
			VECNRM_INF(cws->nv, &str_res_g, 0, &qp_res[0]);
			VECNRM_INF(cws->ne, &str_res_b, 0, &qp_res[1]);
			VECNRM_INF(cws->nc, &str_res_d, 0, &qp_res[2]);
			VECNRM_INF(cws->nc, &str_res_m, 0, &qp_res[3]);
			}

		ws->iter = kk;

		if(kk == arg->iter_max)
			{
			// max iteration number reached
			ws->status = 1;
			}
		else if(cws->alpha <= arg->alpha_min)
			{
			// min step lenght
			ws->status = 2;
			}
#ifdef USE_C99_MATH
		else if(isnan(cws->mu))
			{
			// NaN in the solution
			ws->status = 3;
			}
#else
		else if(cws->mu != cws->mu)
			{
			// NaN in the solution
			ws->status = 3;
			}
#endif
		else
			{
			// normal return
			ws->status = 0;
			}

		return;

		}



	// compute residuals
	OCP_QP_RES_COMPUTE(qp, qp_sol, ws->res, ws->res_workspace);
	BACKUP_RES_M(cws);
	cws->mu = ws->res->res_mu;

	// compute infinity norm of residuals
	VECNRM_INF(cws->nv, &str_res_g, 0, &qp_res[0]);
	VECNRM_INF(cws->ne, &str_res_b, 0, &qp_res[1]);
	VECNRM_INF(cws->nc, &str_res_d, 0, &qp_res[2]);
	VECNRM_INF(cws->nc, &str_res_m, 0, &qp_res[3]);

	ws->stat[ws->stat_m*(0)+5] = qp_res[0];
	ws->stat[ws->stat_m*(0)+6] = qp_res[1];
	ws->stat[ws->stat_m*(0)+7] = qp_res[2];
	ws->stat[ws->stat_m*(0)+8] = qp_res[3];

//printf("\niter %d\t%e\t%e\t%e\t%e\n", -1, qp_res[0], qp_res[1], qp_res[2], qp_res[3]);

	REAL itref_qp_norm[4] = {0,0,0,0};
	REAL itref_qp_norm0[4] = {0,0,0,0};
	int ndp0, ndp1;

	int force_lq = 0;



	// relative IPM formulation

	// IPM loop
	for(kk=0; kk<arg->iter_max & \
			cws->alpha>arg->alpha_min & \
			(qp_res[0]>arg->res_g_max | \
			qp_res[1]>arg->res_b_max | \
			qp_res[2]>arg->res_d_max | \
			qp_res[3]>arg->res_m_max); kk++)
		{

		// fact and solve kkt
		if(arg->lq_fact==0)
			{

			// syrk+cholesky
			FACT_SOLVE_KKT_STEP_OCP_QP(ws->qp_step, ws->sol_step, arg, ws);

			}
		else if(arg->lq_fact==1 & force_lq==0)
			{

			// syrk+chol, switch to lq when needed
			FACT_SOLVE_KKT_STEP_OCP_QP(ws->qp_step, ws->sol_step, arg, ws);

			// compute res of linear system
			OCP_QP_RES_COMPUTE_LIN(ws->qp_step, qp_sol, ws->sol_step, ws->res_itref, ws->res_workspace);
			VECNRM_INF(cws->nv, ws->res_itref->res_g, 0, &itref_qp_norm[0]);
			VECNRM_INF(cws->ne, ws->res_itref->res_b, 0, &itref_qp_norm[1]);
			VECNRM_INF(cws->nc, ws->res_itref->res_d, 0, &itref_qp_norm[2]);
			VECNRM_INF(cws->nc, ws->res_itref->res_m, 0, &itref_qp_norm[3]);

//printf("\n%e\t%e\t%e\t%e\n", itref_qp_norm[0], itref_qp_norm[1], itref_qp_norm[2], itref_qp_norm[3]);

			// inaccurate factorization: switch to lq
			if(
#ifdef USE_C99_MATH
				( itref_qp_norm[0]==0.0 & isnan(BLASFEO_DVECEL(ws->res_itref->res_g+0, 0)) ) |
#else
				( itref_qp_norm[0]==0.0 & BLASFEO_DVECEL(ws->res_itref->res_g+0, 0)!=BLASFEO_DVECEL(ws->res_itref->res_g+0, 0) ) |
#endif
				itref_qp_norm[0]>1e-5 |
				itref_qp_norm[1]>1e-5 |
				itref_qp_norm[2]>1e-5 |
				itref_qp_norm[3]>1e-5 )
				{

#if 0
blasfeo_print_tran_dvec(cws->nv, ws->sol_step->v, 0);
blasfeo_print_tran_dvec(cws->ne, ws->sol_step->pi, 0);
blasfeo_print_tran_dvec(cws->nc, ws->sol_step->lam, 0);
blasfeo_print_tran_dvec(cws->nc, ws->sol_step->t, 0);
#endif

				// refactorize using lq
				FACT_LQ_SOLVE_KKT_STEP_OCP_QP(ws->qp_step, ws->sol_step, arg, ws);

				// switch to lq
				force_lq = 1;

#if 0
blasfeo_print_tran_dvec(cws->nv, ws->sol_step->v, 0);
blasfeo_print_tran_dvec(cws->ne, ws->sol_step->pi, 0);
blasfeo_print_tran_dvec(cws->nc, ws->sol_step->lam, 0);
blasfeo_print_tran_dvec(cws->nc, ws->sol_step->t, 0);
#endif

				}

			}
		else // arg->lq_fact==2
			{

			FACT_LQ_SOLVE_KKT_STEP_OCP_QP(ws->qp_step, ws->sol_step, arg, ws);

			}

		// iterative refinement on prediction step
		for(itref0=0; itref0<arg->itref_pred_max; itref0++)
			{

			OCP_QP_RES_COMPUTE_LIN(ws->qp_step, qp_sol, ws->sol_step, ws->res_itref, ws->res_workspace);

			VECNRM_INF(cws->nv, ws->res_itref->res_g, 0, &itref_qp_norm[0]);
			VECNRM_INF(cws->ne, ws->res_itref->res_b, 0, &itref_qp_norm[1]);
			VECNRM_INF(cws->nc, ws->res_itref->res_d, 0, &itref_qp_norm[2]);
			VECNRM_INF(cws->nc, ws->res_itref->res_m, 0, &itref_qp_norm[3]);

			if(itref0==0)
				{
				itref_qp_norm0[0] = itref_qp_norm[0];
				itref_qp_norm0[1] = itref_qp_norm[1];
				itref_qp_norm0[2] = itref_qp_norm[2];
				itref_qp_norm0[3] = itref_qp_norm[3];
				}

			if( \
					(itref_qp_norm[0]<1e0*arg->res_g_max | itref_qp_norm[0]<1e-3*qp_res[0]) & \
					(itref_qp_norm[1]<1e0*arg->res_b_max | itref_qp_norm[1]<1e-3*qp_res[1]) & \
					(itref_qp_norm[2]<1e0*arg->res_d_max | itref_qp_norm[2]<1e-3*qp_res[2]) & \
					(itref_qp_norm[3]<1e0*arg->res_m_max | itref_qp_norm[3]<1e-3*qp_res[3]) )
//					(itref_qp_norm[0]<=arg->res_g_max) & \
					(itref_qp_norm[1]<=arg->res_b_max) & \
					(itref_qp_norm[2]<=arg->res_d_max) & \
					(itref_qp_norm[3]<=arg->res_m_max) )
				{
				break;
				}

			ws->use_Pb = 0;
			SOLVE_KKT_STEP_OCP_QP(ws->qp_itref, ws->sol_itref, arg, ws);

			for(ii=0; ii<=N; ii++)
				AXPY(nu[ii]+nx[ii]+2*ns[ii], 1.0, ws->sol_itref->ux+ii, 0, ws->sol_step->ux+ii, 0, ws->sol_step->ux+ii, 0);
			for(ii=0; ii<N; ii++)
				AXPY(nx[ii+1], 1.0, ws->sol_itref->pi+ii, 0, ws->sol_step->pi+ii, 0, ws->sol_step->pi+ii, 0);
			for(ii=0; ii<=N; ii++)
				AXPY(2*nb[ii]+2*ng[ii]+2*ns[ii], 1.0, ws->sol_itref->lam+ii, 0, ws->sol_step->lam+ii, 0, ws->sol_step->lam+ii, 0);
			for(ii=0; ii<=N; ii++)
				AXPY(2*nb[ii]+2*ng[ii]+2*ns[ii], 1.0, ws->sol_itref->t+ii, 0, ws->sol_step->t+ii, 0, ws->sol_step->t+ii, 0);

			}

#if 0
int N = qp->dim->N;
int *nx = qp->dim->nx;
int *nu = qp->dim->nu;
int *nb = qp->dim->nb;
int *ng = qp->dim->ng;
int *ns = qp->dim->ns;

int ii;

//exit(1);

//for(ii=0; ii<=N; ii++)
//	blasfeo_print_dmat(nu[ii]+nx[ii], nu[ii]+nx[ii], ws->L+ii, 0, 0);
//exit(1);

printf("\nux\n");
for(ii=0; ii<=N; ii++)
	blasfeo_print_tran_dvec(nu[ii]+nx[ii]+2*ns[ii], ws->dux+ii, 0);
printf("\npi\n");
for(ii=0; ii<N; ii++)
	blasfeo_print_tran_dvec(nx[ii+1], ws->dpi+ii, 0);
//printf("\nlam\n");
//for(ii=0; ii<=N; ii++)
//	blasfeo_print_tran_dvec(2*nb[ii]+2*ng[ii]+2*ns[ii], ws->dlam+ii, 0);
printf("\nt\n");
for(ii=0; ii<=N; ii++)
	blasfeo_print_tran_dvec(2*nb[ii]+2*ng[ii]+2*ns[ii], ws->dt+ii, 0);

SOLVE_KKT_STEP_OCP_QP(qp, ws);

printf("\nux\n");
for(ii=0; ii<=N; ii++)
	blasfeo_print_tran_dvec(nu[ii]+nx[ii]+2*ns[ii], ws->dux+ii, 0);
printf("\npi\n");
for(ii=0; ii<N; ii++)
	blasfeo_print_tran_dvec(nx[ii+1], ws->dpi+ii, 0);
//printf("\nlam\n");
//for(ii=0; ii<=N; ii++)
//	blasfeo_print_tran_dvec(2*nb[ii]+2*ng[ii]+2*ns[ii], ws->dlam+ii, 0);
printf("\nt\n");
for(ii=0; ii<=N; ii++)
	blasfeo_print_tran_dvec(2*nb[ii]+2*ng[ii]+2*ns[ii], ws->dt+ii, 0);

exit(1);
#endif

		// alpha
		COMPUTE_ALPHA_QP(cws);
		if(kk<ws->stat_max)
			ws->stat[ws->stat_m*(kk+1)+0] = cws->alpha;

		// Mehrotra's predictor-corrector
		if(arg->pred_corr==1)
			{
			// mu_aff
			COMPUTE_MU_AFF_QP(cws);
			if(kk<ws->stat_max)
				ws->stat[ws->stat_m*(kk+1)+1] = cws->mu_aff;

			tmp = cws->mu_aff/cws->mu;
			cws->sigma = tmp*tmp*tmp;
			if(kk<ws->stat_max)
				ws->stat[ws->stat_m*(kk+1)+2] = cws->sigma;

			COMPUTE_CENTERING_CORRECTION_QP(cws);

			// fact and solve kkt
			ws->use_Pb = 1;
			SOLVE_KKT_STEP_OCP_QP(ws->qp_step, ws->sol_step, arg, ws);

			// alpha
			COMPUTE_ALPHA_QP(cws);
			if(kk<ws->stat_max)
				ws->stat[ws->stat_m*(kk+1)+3] = cws->alpha;

			// conditional Mehrotra's predictor-corrector
			if(arg->cond_pred_corr==1)
				{

				// save mu_aff (from prediction step)
				mu_aff0 = cws->mu_aff;

				// compute mu for predictor-corrector-centering
				COMPUTE_MU_AFF_QP(cws);

//				if(cws->mu_aff > 2.0*cws->mu)
				if(cws->mu_aff > 2.0*mu_aff0)
					{

					// centering direction
					COMPUTE_CENTERING_QP(cws);

					// solve kkt
					ws->use_Pb = 1;
					SOLVE_KKT_STEP_OCP_QP(ws->qp_step, ws->sol_step, arg, ws);

					// alpha
					COMPUTE_ALPHA_QP(cws);
					if(kk<ws->stat_max)
						ws->stat[ws->stat_m*(kk+1)+3] = cws->alpha;

					}

				}

			iter_ref_step = 0;
			for(itref1=0; itref1<arg->itref_corr_max; itref1++)
				{

				OCP_QP_RES_COMPUTE_LIN(ws->qp_step, qp_sol, ws->sol_step, ws->res_itref, ws->res_workspace);

//for(ii=0; ii<=N; ii++)
//	blasfeo_dvecse(nu[ii]+nx[ii], 0.0, ws->res_itref->res_g+ii, 0);
//for(ii=0; ii<N; ii++)
//	blasfeo_dvecse(nx[ii+1], 0.0, ws->res_itref->res_b+ii, 0);
//for(ii=0; ii<=N; ii++)
//	blasfeo_dvecse(2*nb[ii]+2*ng[ii], 0.0, ws->res_itref->res_d+ii, 0);
//for(ii=0; ii<=N; ii++)
//	blasfeo_dvecse(2*nb[ii]+2*ng[ii], 0.0, ws->res_itref->res_m+ii, 0);
				VECNRM_INF(cws->nv, ws->res_itref->res_g, 0, &itref_qp_norm[0]);
				VECNRM_INF(cws->ne, ws->res_itref->res_b, 0, &itref_qp_norm[1]);
				VECNRM_INF(cws->nc, ws->res_itref->res_d, 0, &itref_qp_norm[2]);
				VECNRM_INF(cws->nc, ws->res_itref->res_m, 0, &itref_qp_norm[3]);

				if(itref1==0)
					{
					itref_qp_norm0[0] = itref_qp_norm[0];
					itref_qp_norm0[1] = itref_qp_norm[1];
					itref_qp_norm0[2] = itref_qp_norm[2];
					itref_qp_norm0[3] = itref_qp_norm[3];
					}

//printf("\nitref1 %d\t%e\t%e\t%e\t%e\n", itref1, itref_qp_norm[0], itref_qp_norm[1], itref_qp_norm[2], itref_qp_norm[3]);

				if( \
						(itref_qp_norm[0]<1e0*arg->res_g_max | itref_qp_norm[0]<1e-3*qp_res[0]) & \
						(itref_qp_norm[1]<1e0*arg->res_b_max | itref_qp_norm[1]<1e-3*qp_res[1]) & \
						(itref_qp_norm[2]<1e0*arg->res_d_max | itref_qp_norm[2]<1e-3*qp_res[2]) & \
						(itref_qp_norm[3]<1e0*arg->res_m_max | itref_qp_norm[3]<1e-3*qp_res[3]) )
//						(itref_qp_norm[0]<=arg->res_g_max) & \
						(itref_qp_norm[1]<=arg->res_b_max) & \
						(itref_qp_norm[2]<=arg->res_d_max) & \
						(itref_qp_norm[3]<=arg->res_m_max) )
					{
					break;
					}

//printf("\nres_g\n");
//ii = 0;
//blasfeo_print_exp_tran_dvec(nu[ii]+nx[ii], ws->res_itref->res_g+ii, 0);
//for(ii=0; ii<=N; ii++)
//	blasfeo_print_exp_tran_dvec(nu[ii]+nx[ii], ws->res_itref->res_g+ii, 0);
//printf("\nres_b\n");
//for(ii=0; ii<N; ii++)
//	blasfeo_print_exp_tran_dvec(nx[ii+1], ws->res_itref->res_b+ii, 0);
//printf("\nres_d\n");
//for(ii=0; ii<=N; ii++)
//	blasfeo_print_exp_tran_dvec(2*nb[ii]+2*ng[ii], ws->res_itref->res_d+ii, 0);
//printf("\nres_m\n");
//for(ii=0; ii<=N; ii++)
//	blasfeo_print_exp_tran_dvec(2*nb[ii]+2*ng[ii], ws->res_itref->res_m+ii, 0);

				ws->use_Pb = 0;
				SOLVE_KKT_STEP_OCP_QP(ws->qp_itref, ws->sol_itref, arg, ws);
//				FACT_SOLVE_LQ_KKT_STEP_OCP_QP(ws->qp_itref, ws->sol_itref, arg, ws);
				iter_ref_step = 1;

//printf("\nux_corr\n");
//ii = 0;
//blasfeo_print_exp_tran_dvec(nu[ii]+nx[ii], ws->sol_itref->ux+ii, 0);
//for(ii=0; ii<=N; ii++)
//	blasfeo_print_exp_tran_dvec(nu[ii]+nx[ii], ws->sol_itref->ux+ii, 0);
//printf("\npi_corr\n");
//for(ii=0; ii<N; ii++)
//	blasfeo_print_exp_tran_dvec(nx[ii+1], ws->sol_itref->pi+ii, 0);
//printf("\nlam_corr\n");
//for(ii=0; ii<=N; ii++)
//	blasfeo_print_exp_tran_dvec(2*nb[ii]+2*ng[ii], ws->sol_itref->lam+ii, 0);
//printf("\nt_corr\n");
//for(ii=0; ii<=N; ii++)
//	blasfeo_print_exp_tran_dvec(2*nb[ii]+2*ng[ii], ws->sol_itref->t+ii, 0);

				for(ii=0; ii<=N; ii++)
					AXPY(nu[ii]+nx[ii]+2*ns[ii], 1.0, ws->sol_itref->ux+ii, 0, ws->sol_step->ux+ii, 0, ws->sol_step->ux+ii, 0);
				for(ii=0; ii<N; ii++)
					AXPY(nx[ii+1], 1.0, ws->sol_itref->pi+ii, 0, ws->sol_step->pi+ii, 0, ws->sol_step->pi+ii, 0);
				for(ii=0; ii<=N; ii++)
					AXPY(2*nb[ii]+2*ng[ii]+2*ns[ii], 1.0, ws->sol_itref->lam+ii, 0, ws->sol_step->lam+ii, 0, ws->sol_step->lam+ii, 0);
				for(ii=0; ii<=N; ii++)
					AXPY(2*nb[ii]+2*ng[ii]+2*ns[ii], 1.0, ws->sol_itref->t+ii, 0, ws->sol_step->t+ii, 0, ws->sol_step->t+ii, 0);

				}

			if(iter_ref_step)
				{
				// alpha
				COMPUTE_ALPHA_QP(cws);
				if(kk<ws->stat_max)
					ws->stat[ws->stat_m*(kk+1)+3] = cws->alpha;
				}

			}

		//
		UPDATE_VAR_QP(cws);

		// compute residuals
		OCP_QP_RES_COMPUTE(qp, qp_sol, ws->res, ws->res_workspace);
		BACKUP_RES_M(cws);
		cws->mu = ws->res->res_mu;
		if(kk<ws->stat_max)
			ws->stat[ws->stat_m*(kk+1)+4] = ws->res->res_mu;

		// compute infinity norm of residuals
		VECNRM_INF(cws->nv, &str_res_g, 0, &qp_res[0]);
		VECNRM_INF(cws->ne, &str_res_b, 0, &qp_res[1]);
		VECNRM_INF(cws->nc, &str_res_d, 0, &qp_res[2]);
		VECNRM_INF(cws->nc, &str_res_m, 0, &qp_res[3]);

		ws->stat[ws->stat_m*(kk+1)+5] = qp_res[0];
		ws->stat[ws->stat_m*(kk+1)+6] = qp_res[1];
		ws->stat[ws->stat_m*(kk+1)+7] = qp_res[2];
		ws->stat[ws->stat_m*(kk+1)+8] = qp_res[3];

//printf("\niter %d\t%e\t%e\t%e\t%e\n", kk, qp_res[0], qp_res[1], qp_res[2], qp_res[3]);

		}

	ws->iter = kk;

#if 0
	printf("\nux\n");
	for(ii=0; ii<=N; ii++)
		blasfeo_print_tran_dvec(nu[ii]+nx[ii], qp_sol->ux+ii, 0);
#endif

	if(kk == arg->iter_max)
		{
		// max iteration number reached
		ws->status = 1;
		}
	else if(cws->alpha <= arg->alpha_min)
		{
		// min step lenght
		ws->status = 2;
		}
#ifdef USE_C99_MATH
	else if(isnan(cws->mu))
		{
		// NaN in the solution
		ws->status = 3;
		}
#else
	else if(cws->mu != cws->mu)
		{
		// NaN in the solution
		ws->status = 3;
		}
#endif
	else
		{
		// normal return
		ws->status = 0;
		}

	return;

	}



void OCP_QP_IPM_PREDICT(struct OCP_QP *qp, struct OCP_QP_SOL *qp_sol, struct OCP_QP_IPM_ARG *arg, struct OCP_QP_IPM_WS *ws)
	{

#if 0
	OCP_QP_DIM_PRINT(qp->dim);
	OCP_QP_PRINT(qp->dim, qp);
#endif

	int ii;

	struct CORE_QP_IPM_WORKSPACE *cws = ws->core_workspace;

	// arg to core workspace
	cws->lam_min = arg->lam_min;
	cws->t_min = arg->t_min;

	// alias qp vectors into qp_sol
	cws->v = qp_sol->ux->pa;
	cws->pi = qp_sol->pi->pa;
	cws->lam = qp_sol->lam->pa;
	cws->t = qp_sol->t->pa;

	// alias members of qp_step
	ws->qp_step->dim = qp->dim;
	ws->qp_step->RSQrq = qp->RSQrq;
	ws->qp_step->BAbt = qp->BAbt;
	ws->qp_step->DCt = qp->DCt;
	ws->qp_step->Z = qp->Z;
	ws->qp_step->idxb = qp->idxb;
	ws->qp_step->idxs = qp->idxs;
	ws->qp_step->rqz = ws->res->res_g;
	ws->qp_step->b = ws->res->res_b;
	ws->qp_step->d = ws->res->res_d;
	ws->qp_step->m = ws->res->res_m;

	// TODO ?

	// blasfeo alias for residuals
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

	REAL *qp_res = ws->qp_res;
	qp_res[0] = 0;
	qp_res[1] = 0;
	qp_res[2] = 0;
	qp_res[3] = 0;

#if 0
// compute residuals
OCP_QP_RES_COMPUTE(qp, qp_sol, ws->res, ws->res_workspace);

// compute infinity norm of residuals
VECNRM_INF(cws->nv, &str_res_g, 0, &qp_res[0]);
VECNRM_INF(cws->ne, &str_res_b, 0, &qp_res[1]);
VECNRM_INF(cws->nc, &str_res_d, 0, &qp_res[2]);
VECNRM_INF(cws->nc, &str_res_m, 0, &qp_res[3]);

printf("\npredict\t%e\t%e\t%e\t%e\n", qp_res[0], qp_res[1], qp_res[2], qp_res[3]);
#endif

	// load sol from bkp
	for(ii=0; ii<cws->nv; ii++)
		cws->v[ii] = cws->v_bkp[ii];
	for(ii=0; ii<cws->ne; ii++)
		cws->pi[ii] = cws->pi_bkp[ii];
	for(ii=0; ii<cws->nc; ii++)
		cws->lam[ii] = cws->lam_bkp[ii];
	for(ii=0; ii<cws->nc; ii++)
		cws->t[ii] = cws->t_bkp[ii];

	// TODO absolute formulation !!!!!

	// TODO robust formulation !!!!!

	// compute residuals
	OCP_QP_RES_COMPUTE(qp, qp_sol, ws->res, ws->res_workspace);

	// compute infinity norm of residuals
	VECNRM_INF(cws->nv, &str_res_g, 0, &qp_res[0]);
	VECNRM_INF(cws->ne, &str_res_b, 0, &qp_res[1]);
	VECNRM_INF(cws->nc, &str_res_d, 0, &qp_res[2]);
	VECNRM_INF(cws->nc, &str_res_m, 0, &qp_res[3]);

//printf("\npredict\t%e\t%e\t%e\t%e\n", qp_res[0], qp_res[1], qp_res[2], qp_res[3]);

	// solve kkt
	ws->use_Pb = 0;
	SOLVE_KKT_STEP_OCP_QP(ws->qp_step, ws->sol_step, arg, ws);

	// alpha TODO fix alpha=1 !!!!!
//	COMPUTE_ALPHA_QP(cws);
	cws->alpha = 1.0;

	//
	UPDATE_VAR_QP(cws);

	if(arg->comp_res_pred)
		{
		// compute residuals in exit
		OCP_QP_RES_COMPUTE(qp, qp_sol, ws->res, ws->res_workspace);

		// compute infinity norm of residuals
		VECNRM_INF(cws->nv, &str_res_g, 0, &qp_res[0]);
		VECNRM_INF(cws->ne, &str_res_b, 0, &qp_res[1]);
		VECNRM_INF(cws->nc, &str_res_d, 0, &qp_res[2]);
		VECNRM_INF(cws->nc, &str_res_m, 0, &qp_res[3]);
		}

//printf("\npredict\t%e\t%e\t%e\t%e\n", qp_res[0], qp_res[1], qp_res[2], qp_res[3]);

	// TODO

	// do not change status

	return;

	}



void OCP_QP_IPM_SENS(struct OCP_QP *qp, struct OCP_QP_SOL *qp_sol, struct OCP_QP_IPM_ARG *arg, struct OCP_QP_IPM_WS *ws)
	{

#if 0
	OCP_QP_DIM_PRINT(qp->dim);
	OCP_QP_PRINT(qp->dim, qp);
#endif

	int ii;

	struct CORE_QP_IPM_WORKSPACE *cws = ws->core_workspace;

	// arg to core workspace
	cws->lam_min = arg->lam_min;
	cws->t_min = arg->t_min;

	// alias qp vectors into qp_sol
	cws->v = qp_sol->ux->pa;
	cws->pi = qp_sol->pi->pa;
	cws->lam = qp_sol->lam->pa;
	cws->t = qp_sol->t->pa;

#if 0
	// alias members of qp_step
	ws->qp_step->dim = qp->dim;
	ws->qp_step->RSQrq = qp->RSQrq;
	ws->qp_step->BAbt = qp->BAbt;
	ws->qp_step->DCt = qp->DCt;
	ws->qp_step->Z = qp->Z;
	ws->qp_step->idxb = qp->idxb;
	ws->qp_step->idxs = qp->idxs;
	ws->qp_step->rqz = res->res_g;
	ws->qp_step->b = res->res_b;
	ws->qp_step->d = res->res_d;
	ws->qp_step->m = res->res_m;
#endif

	// TODO ?

#if 0
	// blasfeo alias for residuals
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

	REAL *qp_res = ws->qp_res;
	qp_res[0] = 0;
	qp_res[1] = 0;
	qp_res[2] = 0;
	qp_res[3] = 0;
#endif

#if 0
// compute residuals
OCP_QP_RES_COMPUTE(qp, qp_sol, ws->res, ws->res_workspace);

// compute infinity norm of residuals
VECNRM_INF(cws->nv, &str_res_g, 0, &qp_res[0]);
VECNRM_INF(cws->ne, &str_res_b, 0, &qp_res[1]);
VECNRM_INF(cws->nc, &str_res_d, 0, &qp_res[2]);
VECNRM_INF(cws->nc, &str_res_m, 0, &qp_res[3]);

printf("\npredict\t%e\t%e\t%e\t%e\n", qp_res[0], qp_res[1], qp_res[2], qp_res[3]);
#endif

	// load sol from bkp
	for(ii=0; ii<cws->nv; ii++)
		cws->v[ii] = cws->v_bkp[ii];
	for(ii=0; ii<cws->ne; ii++)
		cws->pi[ii] = cws->pi_bkp[ii];
	for(ii=0; ii<cws->nc; ii++)
		cws->lam[ii] = cws->lam_bkp[ii];
	for(ii=0; ii<cws->nc; ii++)
		cws->t[ii] = cws->t_bkp[ii];

//printf("\npredict\t%e\t%e\t%e\t%e\n", qp_res[0], qp_res[1], qp_res[2], qp_res[3]);

	// solve kkt
	ws->use_Pb = 0;
//	SOLVE_KKT_STEP_OCP_QP(ws->qp_step, ws->sol_step, arg, ws);
	SOLVE_KKT_STEP_OCP_QP(qp, qp_sol, arg, ws);

#if 0
	// alpha TODO fix alpha=1 !!!!!
//	COMPUTE_ALPHA_QP(cws);
	cws->alpha = 1.0;

	//
	UPDATE_VAR_QP(cws);

	if(arg->comp_res_pred)
		{
		// compute residuals in exit
		OCP_QP_RES_COMPUTE(qp, qp_sol, ws->res, ws->res_workspace);

		// compute infinity norm of residuals
		VECNRM_INF(cws->nv, &str_res_g, 0, &qp_res[0]);
		VECNRM_INF(cws->ne, &str_res_b, 0, &qp_res[1]);
		VECNRM_INF(cws->nc, &str_res_d, 0, &qp_res[2]);
		VECNRM_INF(cws->nc, &str_res_m, 0, &qp_res[3]);
		}
#endif

//printf("\npredict\t%e\t%e\t%e\t%e\n", qp_res[0], qp_res[1], qp_res[2], qp_res[3]);

	// TODO

	// do not change status

	return;

	}
