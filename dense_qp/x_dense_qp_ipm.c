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



int MEMSIZE_DENSE_QP_IPM_ARG(struct DENSE_QP_DIM *dim)
	{

	return 0;

	}



void CREATE_DENSE_QP_IPM_ARG(struct DENSE_QP_DIM *dim, struct DENSE_QP_IPM_ARG *arg, void *mem)
	{

	arg->memsize = 0;

	return;

	}



void SET_DEFAULT_DENSE_QP_IPM_ARG(struct DENSE_QP_IPM_ARG *arg)
	{

	arg->mu0 = 1e1;
	arg->alpha_min = 1e-12;
	arg->res_g_max = 1e-6;
	arg->res_b_max = 1e-8;
	arg->res_d_max = 1e-8;
	arg->res_m_max = 1e-8;
	arg->iter_max = 20;
	arg->stat_max = 20;
	arg->pred_corr = 1;
	arg->cond_pred_corr = 1;
	arg->scale = 1;
	arg->itref_pred_max = 0;
	arg->itref_corr_max = 4;
	arg->reg_prim = 1e-15;
	arg->reg_dual = 1e-15;
	arg->warm_start = 0;

	// TODO if(performance_mode) {} else /* reliability_mode */ {}
	return;

	}



int MEMSIZE_DENSE_QP_IPM(struct DENSE_QP_DIM *dim, struct DENSE_QP_IPM_ARG *arg)
	{

	int nv = dim->nv;
	int ne = dim->ne;
	int nb = dim->nb;
	int ng = dim->ng;
	int ns = dim->ns;

	int size = 0;

	size += 1*sizeof(struct CORE_QP_IPM_WORKSPACE);
	size += 1*MEMSIZE_CORE_QP_IPM(nv+2*ns, ne, 2*nb+2*ng+2*ns);

	size += 1*sizeof(struct DENSE_QP_RES); // res
	size += 1*sizeof(struct DENSE_QP_RES_WORKSPACE); // res_workspace

	size += 2*sizeof(struct DENSE_QP_SOL); // sol_step(v,pi,lam,t) sol_itref
	size += 1*MEMSIZE_DENSE_QP_SOL(dim); // sol_itref

	size += 2*sizeof(struct DENSE_QP); // qp_step qp_itref

	size += 1*sizeof(struct DENSE_QP_RES); // res_itref
	size += 1*MEMSIZE_DENSE_QP_RES(dim); // res_itref

	size += 22*sizeof(struct STRVEC); // sol_step(v,pi,lam,t) res_g res_b res_d res_m lv (4+2)*tmp_nbg (1+1)*tmp_ns Gamma gamma Zs_inv sv se
	size += 7*sizeof(struct STRMAT); // 2*Lv AL Le Ctx lq0 lq1

	size += 4*SIZE_STRVEC(nb+ng); // 4*tmp_nbg
	size += 1*SIZE_STRVEC(ns); // tmp_ns
	size += 2*SIZE_STRVEC(nv); // lv sv
	size += 1*SIZE_STRVEC(ne); // se
	size += 1*SIZE_STRVEC(2*ns); // Zs_inv
	size += 2*SIZE_STRMAT(nv+1, nv); // Lv
	size += 1*SIZE_STRMAT(ne, nv); // AL
	size += 1*SIZE_STRMAT(ne, ne); // Le
	size += 1*SIZE_STRMAT(nv+1, ng); // Ctx
	size += 1*SIZE_STRMAT(ne, ne+nv); // lq0
	size += 1*SIZE_STRMAT(nv, nv+nv+ng); // lq1

	size += nv*sizeof(int); // ipiv_v
	size += ne*sizeof(int); // ipiv_e

	size += 1*GELQF_WORKSIZE(ne, nv); // lq_work0
	size += 1*GELQF_WORKSIZE(nv, nv+nv+ng); // lq_work1

	size += 5*arg->stat_max*sizeof(REAL);

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

	}



void CREATE_DENSE_QP_IPM(struct DENSE_QP_DIM *dim, struct DENSE_QP_IPM_ARG *arg, struct DENSE_QP_IPM_WORKSPACE *workspace, void *mem)
	{

	int nv = dim->nv;
	int ne = dim->ne;
	int nb = dim->nb;
	int ng = dim->ng;
	int ns = dim->ns;


	// core struct
	struct CORE_QP_IPM_WORKSPACE *sr_ptr = mem;

	// core workspace
	workspace->core_workspace = sr_ptr;
	sr_ptr += 1;
	struct CORE_QP_IPM_WORKSPACE *cws = workspace->core_workspace;


	// res struct
	struct DENSE_QP_RES *res_ptr = (struct DENSE_QP_RES *) sr_ptr;

	workspace->res = res_ptr;
	res_ptr += 1;
	workspace->res_itref = res_ptr;
	res_ptr += 1;


	// res workspace struct
	struct DENSE_QP_RES_WORKSPACE *res_ws_ptr = (struct DENSE_QP_RES_WORKSPACE *) res_ptr;

	workspace->res_workspace = res_ws_ptr;
	res_ws_ptr += 1;


	// qp sol struct
	struct DENSE_QP_SOL *qp_sol_ptr = (struct DENSE_QP_SOL *) res_ws_ptr;

	workspace->sol_step = qp_sol_ptr;
	qp_sol_ptr += 1;
	workspace->sol_itref = qp_sol_ptr;
	qp_sol_ptr += 1;


	// qp struct
	struct DENSE_QP *qp_ptr = (struct DENSE_QP *) qp_sol_ptr;

	workspace->qp_step = qp_ptr;
	qp_ptr += 1;
	workspace->qp_itref = qp_ptr;
	qp_ptr += 1;


	// matrix struct
	struct STRMAT *sm_ptr = (struct STRMAT *) qp_ptr;

	workspace->Lv = sm_ptr;
	sm_ptr += 2;
	workspace->AL = sm_ptr;
	sm_ptr += 1;
	workspace->Le = sm_ptr;
	sm_ptr += 1;
	workspace->Ctx = sm_ptr;
	sm_ptr += 1;
	workspace->lq0 = sm_ptr;
	sm_ptr += 1;
	workspace->lq1 = sm_ptr;
	sm_ptr += 1;


	// vector struct
	struct STRVEC *sv_ptr = (struct STRVEC *) sm_ptr;

	workspace->sol_step->v = sv_ptr;
	sv_ptr += 1;
	workspace->sol_step->pi = sv_ptr;
	sv_ptr += 1;
	workspace->sol_step->lam = sv_ptr;
	sv_ptr += 1;
	workspace->sol_step->t = sv_ptr;
	sv_ptr += 1;
	workspace->res->res_g = sv_ptr;
	sv_ptr += 1;
	workspace->res->res_b = sv_ptr;
	sv_ptr += 1;
	workspace->res->res_d = sv_ptr;
	sv_ptr += 1;
	workspace->res->res_m = sv_ptr;
	sv_ptr += 1;
	workspace->Gamma = sv_ptr;
	sv_ptr += 1;
	workspace->gamma = sv_ptr;
	sv_ptr += 1;
	workspace->Zs_inv = sv_ptr;
	sv_ptr += 1;
	workspace->lv = sv_ptr;
	sv_ptr += 1;
	workspace->sv = sv_ptr;
	sv_ptr += 1;
	workspace->se = sv_ptr;
	sv_ptr += 1;
	workspace->tmp_nbg = sv_ptr;
	sv_ptr += 4;
	workspace->res_workspace->tmp_nbg = sv_ptr;
	sv_ptr += 2;
	workspace->tmp_ns = sv_ptr;
	sv_ptr += 1;
	workspace->res_workspace->tmp_ns = sv_ptr;
	sv_ptr += 1;


	// double/float stuff
	REAL *d_ptr = (REAL *) sv_ptr;
	
	workspace->stat = d_ptr;
	d_ptr += 5*arg->stat_max;

	workspace->stat_max = arg->stat_max;
	workspace->warm_start = arg->warm_start;


	// int suff
	int *i_ptr = (int *) d_ptr;

	workspace->ipiv_v = i_ptr;
	i_ptr += nv;

	workspace->ipiv_e = i_ptr;
	i_ptr += ne;


	// align to typicl cache line size
	size_t s_ptr = (size_t) i_ptr;
	s_ptr = (s_ptr+63)/64*64;


	// void stuf
	char *c_ptr = (char *) s_ptr;

	CREATE_DENSE_QP_SOL(dim, workspace->sol_itref, c_ptr);
	c_ptr += workspace->sol_itref->memsize;

	CREATE_DENSE_QP_RES(dim, workspace->res_itref, c_ptr);
	c_ptr += workspace->res_itref->memsize;

	CREATE_STRMAT(nv+1, nv, workspace->Lv, c_ptr);
	c_ptr += workspace->Lv->memsize;

	CREATE_STRMAT(nv+1, nv, workspace->Lv+1, c_ptr);
	c_ptr += workspace->Lv[1].memsize;

	CREATE_STRMAT(ne, nv, workspace->AL, c_ptr);
	c_ptr += workspace->AL->memsize;

	CREATE_STRMAT(ne, ne, workspace->Le, c_ptr);
	c_ptr += workspace->Le->memsize;

	CREATE_STRMAT(nv+1, ng, workspace->Ctx, c_ptr);
	c_ptr += workspace->Ctx->memsize;

	CREATE_STRMAT(ne, ne+nv, workspace->lq0, c_ptr);
	c_ptr += workspace->lq0->memsize;

	CREATE_STRMAT(nv, nv+nv+ng, workspace->lq1, c_ptr);
	c_ptr += workspace->lq1->memsize;

	CREATE_STRVEC(nv, workspace->lv, c_ptr);
	c_ptr += workspace->lv->memsize;

	CREATE_STRVEC(nv, workspace->sv, c_ptr);
	c_ptr += workspace->sv->memsize;

	CREATE_STRVEC(ne, workspace->se, c_ptr);
	c_ptr += workspace->se->memsize;

	CREATE_STRVEC(2*ns, workspace->Zs_inv, c_ptr);
	c_ptr += workspace->Zs_inv->memsize;

	CREATE_STRVEC(nb+ng, workspace->tmp_nbg+0, c_ptr);
	CREATE_STRVEC(nb+ng, workspace->res_workspace->tmp_nbg+0, c_ptr);
	c_ptr += (workspace->tmp_nbg+0)->memsize;

	CREATE_STRVEC(nb+ng, workspace->tmp_nbg+1, c_ptr);
	CREATE_STRVEC(nb+ng, workspace->res_workspace->tmp_nbg+1, c_ptr);
	c_ptr += (workspace->tmp_nbg+1)->memsize;

	CREATE_STRVEC(nb+ng, workspace->tmp_nbg+2, c_ptr);
	c_ptr += (workspace->tmp_nbg+2)->memsize;

	CREATE_STRVEC(nb+ng, workspace->tmp_nbg+3, c_ptr);
	c_ptr += (workspace->tmp_nbg+3)->memsize;

	CREATE_STRVEC(ns, workspace->tmp_ns+0, c_ptr);
	CREATE_STRVEC(ns, workspace->res_workspace->tmp_ns+0, c_ptr);
	c_ptr += (workspace->tmp_ns+0)->memsize;

	CREATE_CORE_QP_IPM(nv+2*ns, ne, 2*nb+2*ng+2*ns, cws, c_ptr);
	c_ptr += workspace->core_workspace->memsize;

	workspace->lq_work0 = c_ptr;
	c_ptr += GELQF_WORKSIZE(ne, nv);

	workspace->lq_work1 = c_ptr;
	c_ptr += GELQF_WORKSIZE(nv, nv+nv+ng);


	// alias members of workspace and core_workspace
	//
	CREATE_STRVEC(nv+2*ns, workspace->sol_step->v, cws->dv);
	//
	CREATE_STRVEC(ne, workspace->sol_step->pi, cws->dpi);
	//
	CREATE_STRVEC(2*nb+2*ng+2*ns, workspace->sol_step->lam, cws->dlam);
	//
	CREATE_STRVEC(2*nb+2*ng+2*ns, workspace->sol_step->t, cws->dt);
	//
	CREATE_STRVEC(nv+2*ns, workspace->res->res_g, cws->res_g);
	//
	CREATE_STRVEC(ne, workspace->res->res_b, cws->res_b);
	//
	CREATE_STRVEC(2*nb+2*ng+2*ns, workspace->res->res_d, cws->res_d);
	//
	CREATE_STRVEC(2*nb+2*ng+2*ns, workspace->res->res_m, cws->res_m);
	//
	CREATE_STRVEC(2*nb+2*ng+2*ns, workspace->Gamma, cws->Gamma);
	//
	CREATE_STRVEC(2*nb+2*ng+2*ns, workspace->gamma, cws->gamma);

	//
	workspace->sol_step->dim = dim;

	//
	workspace->memsize = MEMSIZE_DENSE_QP_IPM(dim, arg);


#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + workspace->memsize)
		{
		printf("\nCreate_dense_qp_ipm: outsize memory bounds!\n\n");
		exit(1);
		}
#endif


	return;

	}



int SOLVE_DENSE_QP_IPM(struct DENSE_QP *qp, struct DENSE_QP_SOL *qp_sol, struct DENSE_QP_IPM_ARG *arg, struct DENSE_QP_IPM_WORKSPACE *ws)
	{

	struct CORE_QP_IPM_WORKSPACE *cws = ws->core_workspace;

	// dims
	int nv = qp->dim->nv;
	int ne = qp->dim->ne;
	int nb = qp->dim->nb;
	int ng = qp->dim->ng;
	int ns = qp->dim->ns;

	// alias qp vectors into qp_sol
	cws->v = qp_sol->v->pa;
	cws->pi = qp_sol->pi->pa;
	cws->lam = qp_sol->lam->pa;
	cws->t = qp_sol->t->pa;

	// alias members of qp_step
	ws->qp_step->dim = qp->dim;
	ws->qp_step->Hv = qp->Hv;
	ws->qp_step->A = qp->A;
	ws->qp_step->Ct = qp->Ct;
	ws->qp_step->Z = qp->Z;
	ws->qp_step->idxb = qp->idxb;
	ws->qp_step->idxs = qp->idxs;
	ws->qp_step->gz = ws->res->res_g;
	ws->qp_step->b = ws->res->res_b;
	ws->qp_step->d = ws->res->res_d;
	ws->qp_step->m = ws->res->res_m;

	// alias members of qp_itref
	ws->qp_itref->dim = qp->dim;
	ws->qp_itref->Hv = qp->Hv;
	ws->qp_itref->A = qp->A;
	ws->qp_itref->Ct = qp->Ct;
	ws->qp_itref->Z = qp->Z;
	ws->qp_itref->idxb = qp->idxb;
	ws->qp_itref->idxs = qp->idxs;
	ws->qp_itref->gz = ws->res_itref->res_g;
	ws->qp_itref->b = ws->res_itref->res_b;
	ws->qp_itref->d = ws->res_itref->res_d;
	ws->qp_itref->m = ws->res_itref->res_m;

	// no constraints
	if(cws->nc==0)
		{
		FACT_SOLVE_KKT_UNCONSTR_DENSE_QP(qp, qp_sol, ws);
		COMPUTE_RES_DENSE_QP(qp, qp_sol, ws->res, ws->res_workspace);
		cws->mu = ws->res->res_mu;
		ws->iter = 0;
		return 0;
		}
	
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

	ws->mu0 = arg->mu0;

	int kk, ii, itref0=0, itref1=0;
	REAL tmp;
	REAL mu_aff0;
	int iter_ref_step;

	// init solver
	INIT_VAR_DENSE_QP(qp, qp_sol, ws);

	// compute residuals
	COMPUTE_RES_DENSE_QP(qp, qp_sol, ws->res, ws->res_workspace);
	BACKUP_RES_M(cws);
	cws->mu = ws->res->res_mu;

	cws->alpha = 1.0;

	// compute infinity norm of residuals
	VECNRM_INF_LIBSTR(cws->nv, &str_res_g, 0, &qp_res[0]);
	VECNRM_INF_LIBSTR(cws->ne, &str_res_b, 0, &qp_res[1]);
	VECNRM_INF_LIBSTR(cws->nc, &str_res_d, 0, &qp_res[2]);
	VECNRM_INF_LIBSTR(cws->nc, &str_res_m, 0, &qp_res[3]);

//	REAL sigma_min = 1e9;
//	sigma_min = arg->res_g_max<sigma_min ? arg->res_g_max : sigma_min;
//	sigma_min = arg->res_b_max<sigma_min ? arg->res_b_max : sigma_min;
//	sigma_min = arg->res_d_max<sigma_min ? arg->res_d_max : sigma_min;
//	sigma_min = arg->res_m_max<sigma_min ? arg->res_m_max : sigma_min;
//	sigma_min *= 0.1;

	REAL itref_qp_norm[4] = {0,0,0,0};
	REAL itref_qp_norm0[4] = {0,0,0,0};
	int ndp0, ndp1;

	for(kk=0; \
			kk<arg->iter_max & \
			cws->alpha>arg->alpha_min & \
			(qp_res[0]>arg->res_g_max | \
			qp_res[1]>arg->res_b_max | \
			qp_res[2]>arg->res_d_max | \
			qp_res[3]>arg->res_m_max); kk++)
		{

		ws->scale = arg->scale;

		// fact and solve kkt
//		FACT_SOLVE_KKT_STEP_DENSE_QP(ws->qp_step, ws->sol_step, arg, ws);
		FACT_SOLVE_LQ_KKT_STEP_DENSE_QP(ws->qp_step, ws->sol_step, arg, ws);
//		FACT_SOLVE_LU_KKT_STEP_DENSE_QP(ws->qp_step, ws->sol_step, arg, ws);

		for(itref0=0; itref0<arg->itref_pred_max; itref0++)
			{

			COMPUTE_LIN_RES_DENSE_QP(ws->qp_step, qp_sol, ws->sol_step, ws->res_itref, ws->res_workspace);

			VECNRM_INF_LIBSTR(cws->nv, ws->res_itref->res_g, 0, &itref_qp_norm[0]);
			VECNRM_INF_LIBSTR(cws->ne, ws->res_itref->res_b, 0, &itref_qp_norm[1]);
			VECNRM_INF_LIBSTR(cws->nc, ws->res_itref->res_d, 0, &itref_qp_norm[2]);
			VECNRM_INF_LIBSTR(cws->nc, ws->res_itref->res_m, 0, &itref_qp_norm[3]);

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

			SOLVE_KKT_STEP_DENSE_QP(ws->qp_itref, ws->sol_itref, arg, ws);

			AXPY_LIBSTR(nv+2*ns, 1.0, ws->sol_itref->v, 0, ws->sol_step->v, 0, ws->sol_step->v, 0);
			AXPY_LIBSTR(ne, 1.0, ws->sol_itref->pi, 0, ws->sol_step->pi, 0, ws->sol_step->pi, 0);
			AXPY_LIBSTR(2*nb+2*ng+2*ns, 1.0, ws->sol_itref->lam, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
			AXPY_LIBSTR(2*nb+2*ng+2*ns, 1.0, ws->sol_itref->t, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);

			}

#if 0
		COMPUTE_LIN_RES_DENSE_QP(ws->qp_step, qp_sol, ws->sol_step, ws->res_itref, ws->res_workspace);
		VECNRM_INF_LIBSTR(cws->nv, ws->res_itref->res_g, 0, &itref_qp_norm0[0]);
		VECNRM_INF_LIBSTR(cws->ne, ws->res_itref->res_b, 0, &itref_qp_norm0[1]);
		VECNRM_INF_LIBSTR(cws->nc, ws->res_itref->res_d, 0, &itref_qp_norm0[2]);
		VECNRM_INF_LIBSTR(cws->nc, ws->res_itref->res_m, 0, &itref_qp_norm0[3]);
//		printf("\nkk = %d\n", kk);
//		blasfeo_print_exp_tran_dvec(qp->dim->nv, ws->res_itref->res_g, 0);
//		blasfeo_print_exp_tran_dvec(qp->dim->ne, ws->res_itref->res_b, 0);
//		blasfeo_print_exp_tran_dvec(2*qp->dim->nb+2*qp->dim->ng, ws->res_itref->res_d, 0);
#endif

#if 0
		COMPUTE_RES_DENSE_QP(ws->qp_step, ws->sol_step, ws->res_itref, ws->res_workspace);
		VECNRM_INF_LIBSTR(cws->nv, ws->res_itref->res_g, 0, &itref_qp_norm0[0]);
		VECNRM_INF_LIBSTR(cws->ne, ws->res_itref->res_b, 0, &itref_qp_norm0[1]);
		VECNRM_INF_LIBSTR(cws->nc, ws->res_itref->res_d, 0, &itref_qp_norm0[2]);
		VECNRM_INF_LIBSTR(cws->nc, ws->res_itref->res_m, 0, &itref_qp_norm0[3]);
//		printf("\nkk = %d\n", kk);
//		blasfeo_print_exp_tran_dvec(qp->dim->nv, ws->res_itref->res_g, 0);
//		blasfeo_print_exp_tran_dvec(qp->dim->ne, ws->res_itref->res_b, 0);
//		blasfeo_print_exp_tran_dvec(2*qp->dim->nb+2*qp->dim->ng, ws->res_itref->res_d, 0);
#endif

#if 1
		ndp0 = 0;
		for(ii=0; ii<qp->dim->nv; ii++)
			{
			if(ws->Lv->dA[ii]<=0)
				{
				ndp0 = ii;
				break;
				}
			}
		ndp1 = 0;
		for(ii=0; ii<qp->dim->ne; ii++)
			{
			if(ws->Le->dA[ii]<=0)
				{
				ndp1 = ii;
				break;
				}
			}
#endif

		// alpha
		COMPUTE_ALPHA_QP(cws);
		if(kk<ws->stat_max)
			ws->stat[5*kk+0] = cws->alpha;

		// Mehrotra's predictor-corrector
		if(arg->pred_corr==1)
			{
			// mu_aff
			COMPUTE_MU_AFF_QP(cws);
			if(kk<ws->stat_max)
				ws->stat[5*kk+1] = cws->mu_aff;

			// compute centering parameter
			tmp = cws->mu_aff/cws->mu;
			cws->sigma = tmp*tmp*tmp;
//			cws->sigma = sigma_min>cws->sigma ? sigma_min : cws->sigma;
			if(kk<ws->stat_max)
				ws->stat[5*kk+2] = cws->sigma;

			COMPUTE_CENTERING_CORRECTION_QP(cws);

			// solve kkt
			SOLVE_KKT_STEP_DENSE_QP(ws->qp_step, ws->sol_step, arg, ws);

			// alpha
			COMPUTE_ALPHA_QP(cws);
			if(kk<ws->stat_max)
				ws->stat[5*kk+3] = cws->alpha;

			// conditional Mehrotra's predictor-corrector
			if(arg->cond_pred_corr==1)
				{

				// save mu_aff (from prediction sol_step)
				mu_aff0 = cws->mu_aff;

				// compute mu for predictor-corrector-centering
				COMPUTE_MU_AFF_QP(cws);

//				if(cws->mu_aff > 2.0*cws->mu)
				if(cws->mu_aff > 2.0*mu_aff0)
					{

					// centering direction
					COMPUTE_CENTERING_QP(cws);

					// solve kkt
					SOLVE_KKT_STEP_DENSE_QP(ws->qp_step, ws->sol_step, arg, ws);

					// alpha
					COMPUTE_ALPHA_QP(cws);
					if(kk<ws->stat_max)
						ws->stat[5*kk+3] = cws->alpha;

					}

				}

			iter_ref_step = 0;
			for(itref1=0; itref1<arg->itref_corr_max; itref1++)
				{

				COMPUTE_LIN_RES_DENSE_QP(ws->qp_step, qp_sol, ws->sol_step, ws->res_itref, ws->res_workspace);

				VECNRM_INF_LIBSTR(cws->nv, ws->res_itref->res_g, 0, &itref_qp_norm[0]);
				VECNRM_INF_LIBSTR(cws->ne, ws->res_itref->res_b, 0, &itref_qp_norm[1]);
				VECNRM_INF_LIBSTR(cws->nc, ws->res_itref->res_d, 0, &itref_qp_norm[2]);
				VECNRM_INF_LIBSTR(cws->nc, ws->res_itref->res_m, 0, &itref_qp_norm[3]);

				if(itref1==0)
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
//						(itref_qp_norm[0]<=arg->res_g_max) & \
						(itref_qp_norm[1]<=arg->res_b_max) & \
						(itref_qp_norm[2]<=arg->res_d_max) & \
						(itref_qp_norm[3]<=arg->res_m_max) )
					{
					break;
					}

				SOLVE_KKT_STEP_DENSE_QP(ws->qp_itref, ws->sol_itref, arg, ws);
				iter_ref_step = 1;

				AXPY_LIBSTR(nv+2*ns, 1.0, ws->sol_itref->v, 0, ws->sol_step->v, 0, ws->sol_step->v, 0);
				AXPY_LIBSTR(ne, 1.0, ws->sol_itref->pi, 0, ws->sol_step->pi, 0, ws->sol_step->pi, 0);
				AXPY_LIBSTR(2*nb+2*ng+2*ns, 1.0, ws->sol_itref->lam, 0, ws->sol_step->lam, 0, ws->sol_step->lam, 0);
				AXPY_LIBSTR(2*nb+2*ng+2*ns, 1.0, ws->sol_itref->t, 0, ws->sol_step->t, 0, ws->sol_step->t, 0);

				}

			if(iter_ref_step)
				{
				// alpha
				COMPUTE_ALPHA_QP(cws);
				if(kk<ws->stat_max)
					ws->stat[5*kk+3] = cws->alpha;
				}

			}

		//
		UPDATE_VAR_QP(cws);

		// compute residuals
		COMPUTE_RES_DENSE_QP(qp, qp_sol, ws->res, ws->res_workspace);
		BACKUP_RES_M(cws);
		cws->mu = ws->res->res_mu;
		if(kk<ws->stat_max)
			ws->stat[5*kk+4] = ws->res->res_mu;

		// compute infinity norm of residuals
		VECNRM_INF_LIBSTR(cws->nv, &str_res_g, 0, &qp_res[0]);
		VECNRM_INF_LIBSTR(cws->ne, &str_res_b, 0, &qp_res[1]);
		VECNRM_INF_LIBSTR(cws->nc, &str_res_d, 0, &qp_res[2]);
		VECNRM_INF_LIBSTR(cws->nc, &str_res_m, 0, &qp_res[3]);

#if 0
printf("%e %e %e\n", cws->alpha, cws->alpha_prim, cws->alpha_dual);
#endif

#if 0
printf("%e %e %e %e %e\n", ws->stat[5*kk+0], ws->stat[5*kk+1], ws->stat[5*kk+2], ws->stat[5*kk+3], ws->stat[5*kk+4]);
#endif

#if 0
printf("%e %e %e %e\n", qp_res[0], qp_res[1], qp_res[2], qp_res[3]);
#endif

#if 0
printf("iter %3d   alpha %1.3e %1.3e   sigma %1.3e   ndp %3d %3d   itref %d %d   res_kkt %1.3e %1.3e %1.3e %1.3e   res_kkt %1.3e %1.3e %1.3e %1.3e   res_qp %1.3e %1.3e %1.3e %1.3e\n", kk, cws->alpha_prim, cws->alpha_dual, cws->sigma, ndp0, ndp1, itref0, itref1, itref_qp_norm0[0], itref_qp_norm0[1], itref_qp_norm0[2], itref_qp_norm0[3], itref_qp_norm[0], itref_qp_norm[1], itref_qp_norm[2], itref_qp_norm[3], qp_res[0], qp_res[1], qp_res[2], qp_res[3]);
#endif

		}

	ws->iter = kk;

	// max iteration number reached
	if(kk == arg->iter_max)
		return 1;

	// min step lenght
	if(cws->alpha <= arg->alpha_min)
		return 2;
	
	// NaN in the solution
#ifdef USE_C99_MATH
	if(isnan(cws->mu))
		return 3;
#else
	if(cws->mu != cws->mu)
		return 3;
#endif

	// normal return
	return 0;

	}


