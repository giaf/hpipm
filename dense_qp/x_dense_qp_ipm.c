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

	arg->mu0 = 100;
	arg->alpha_min = 1e-12;
	arg->res_g_max = 1e-6;
	arg->res_b_max = 1e-8;
	arg->res_d_max = 1e-8;
	arg->res_m_max = 1e-8;
	arg->iter_max = 20;
	arg->stat_max = 20;
	arg->pred_corr = 1;
	arg->warm_start = 0;

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

	size += 1*sizeof(struct DENSE_QP_SOL); // step(v,pi,lam,t)

	size += 1*sizeof(struct DENSE_QP); // itref_qp

	size += 1*sizeof(struct DENSE_QP_RES); // itref_res
	size += 1*MEMSIZE_DENSE_QP_RES(dim); // itref_res

	size += 20*sizeof(struct STRVEC); // step(v,pi,lam,t) res_g res_b res_d res_m lv (4+2)*tmp_nbg (1+1)*tmp_ns Gamma gamma Zs_inv
	size += 4*sizeof(struct STRMAT); // Lv AL Le Ctx

	size += 4*SIZE_STRVEC(nb+ng); // 4*tmp_nbg
	size += 1*SIZE_STRVEC(ns); // tmp_ns
	size += 1*SIZE_STRVEC(nv); // lv
	size += 1*SIZE_STRVEC(2*ns); // Zs_inv
	size += 1*SIZE_STRMAT(nv+1, nv); // Lv
	size += 1*SIZE_STRMAT(ne, nv); // AL
	size += 1*SIZE_STRMAT(ne, ne); // Le
	size += 1*SIZE_STRMAT(nv+1, ng); // Ctx

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
	workspace->itref_res = res_ptr;
	res_ptr += 1;


	// res workspace struct
	struct DENSE_QP_RES_WORKSPACE *res_ws_ptr = (struct DENSE_QP_RES_WORKSPACE *) res_ptr;
	workspace->res_workspace = res_ws_ptr;
	res_ws_ptr += 1;


	// qp sol struct
	struct DENSE_QP_SOL *qp_sol_ptr = (struct DENSE_QP_SOL *) res_ws_ptr;
	workspace->step = qp_sol_ptr;
	qp_sol_ptr += 1;


	// qp struct
	struct DENSE_QP *qp_ptr = (struct DENSE_QP *) qp_sol_ptr;
	workspace->itref_qp = qp_ptr;
	qp_ptr += 1;


	// matrix struct
	struct STRMAT *sm_ptr = (struct STRMAT *) qp_ptr;

	workspace->Lv = sm_ptr;
	sm_ptr += 1;
	workspace->AL = sm_ptr;
	sm_ptr += 1;
	workspace->Le = sm_ptr;
	sm_ptr += 1;
	workspace->Ctx = sm_ptr;
	sm_ptr += 1;


	// vector struct
	struct STRVEC *sv_ptr = (struct STRVEC *) sm_ptr;

	workspace->step->v = sv_ptr;
	sv_ptr += 1;
	workspace->step->pi = sv_ptr;
	sv_ptr += 1;
	workspace->step->lam = sv_ptr;
	sv_ptr += 1;
	workspace->step->t = sv_ptr;
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


	// align to typicl cache line size
	size_t s_ptr = (size_t) d_ptr;
	s_ptr = (s_ptr+63)/64*64;


	// void stuf
	char *c_ptr = (char *) s_ptr;

	CREATE_DENSE_QP_RES(dim, workspace->itref_res, c_ptr);
	c_ptr += workspace->itref_res->memsize;

	CREATE_STRMAT(nv+1, nv, workspace->Lv, c_ptr);
	c_ptr += workspace->Lv->memory_size;

	CREATE_STRMAT(ne, nv, workspace->AL, c_ptr);
	c_ptr += workspace->AL->memory_size;

	CREATE_STRMAT(ne, ne, workspace->Le, c_ptr);
	c_ptr += workspace->Le->memory_size;

	CREATE_STRMAT(nv+1, ng, workspace->Ctx, c_ptr);
	c_ptr += workspace->Ctx->memory_size;

	CREATE_STRVEC(nv, workspace->lv, c_ptr);
	c_ptr += workspace->lv->memory_size;

	CREATE_STRVEC(2*ns, workspace->Zs_inv, c_ptr);
	c_ptr += workspace->Zs_inv->memory_size;

	CREATE_STRVEC(nb+ng, workspace->tmp_nbg+0, c_ptr);
	CREATE_STRVEC(nb+ng, workspace->res_workspace->tmp_nbg+0, c_ptr);
	c_ptr += (workspace->tmp_nbg+0)->memory_size;

	CREATE_STRVEC(nb+ng, workspace->tmp_nbg+1, c_ptr);
	CREATE_STRVEC(nb+ng, workspace->res_workspace->tmp_nbg+1, c_ptr);
	c_ptr += (workspace->tmp_nbg+1)->memory_size;

	CREATE_STRVEC(nb+ng, workspace->tmp_nbg+2, c_ptr);
	c_ptr += (workspace->tmp_nbg+2)->memory_size;

	CREATE_STRVEC(nb+ng, workspace->tmp_nbg+3, c_ptr);
	c_ptr += (workspace->tmp_nbg+3)->memory_size;

	CREATE_STRVEC(ns, workspace->tmp_ns+0, c_ptr);
	CREATE_STRVEC(ns, workspace->res_workspace->tmp_ns+0, c_ptr);
	c_ptr += (workspace->tmp_ns+0)->memory_size;

	CREATE_CORE_QP_IPM(nv+2*ns, ne, 2*nb+2*ng+2*ns, cws, c_ptr);
	c_ptr += workspace->core_workspace->memsize;


	// alias members of workspace and core_workspace
	//
	CREATE_STRVEC(nv+2*ns, workspace->step->v, cws->dv);
	//
	CREATE_STRVEC(ne, workspace->step->pi, cws->dpi);
	//
	CREATE_STRVEC(2*nb+2*ng+2*ns, workspace->step->lam, cws->dlam);
	//
	CREATE_STRVEC(2*nb+2*ng+2*ns, workspace->step->t, cws->dt);
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
	workspace->step->dim = dim;

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

	// alias qp vectors into qp_sol
	cws->v = qp_sol->v->pa;
	cws->pi = qp_sol->pi->pa;
	cws->lam = qp_sol->lam->pa;
	cws->t = qp_sol->t->pa;

	// alias members of itref_qp
	ws->itref_qp->dim = qp->dim;
	ws->itref_qp->Hv = qp->Hv;
	ws->itref_qp->A = qp->A;
	ws->itref_qp->Ct = qp->Ct;
	ws->itref_qp->Z = qp->Z; // TODO ???
	ws->itref_qp->z = qp->z; // TODO ???
	ws->itref_qp->idxb = qp->idxb;
	ws->itref_qp->idxs = qp->idxs;
	ws->itref_qp->g = ws->res->res_g;
	ws->itref_qp->b = ws->res->res_b;
	ws->itref_qp->d = ws->res->res_d;
	// res_m ???

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

	int kk = 0;
	REAL tmp;

	// init solver
	INIT_VAR_DENSE_QP(qp, qp_sol, ws);

	// compute residuals
	COMPUTE_RES_DENSE_QP(qp, qp_sol, ws->res, ws->res_workspace);
	cws->mu = ws->res->res_mu;

	cws->alpha = 1.0;

	// compute infinity norm of residuals
	VECNRM_INF_LIBSTR(cws->nv, &str_res_g, 0, &qp_res[0]);
	VECNRM_INF_LIBSTR(cws->ne, &str_res_b, 0, &qp_res[1]);
	VECNRM_INF_LIBSTR(cws->nc, &str_res_d, 0, &qp_res[2]);
	VECNRM_INF_LIBSTR(cws->nc, &str_res_m, 0, &qp_res[3]);

	for(kk=0; kk<arg->iter_max & cws->alpha>arg->alpha_min & (qp_res[0]>arg->res_g_max | qp_res[1]>arg->res_b_max | qp_res[2]>arg->res_d_max | qp_res[3]>arg->res_m_max); kk++)
		{

		// fact and solve kkt
		FACT_SOLVE_KKT_STEP_DENSE_QP(qp, ws);

#if 0
		COMPUTE_RES_DENSE_QP(ws->itref_qp, ws->step, ws->itref_res, ws->res_workspace);
		printf("\nkk = %d\n", kk);
		d_print_e_tran_strvec(qp->dim->nv, ws->itref_res->res_g, 0);
		d_print_e_tran_strvec(qp->dim->ne, ws->itref_res->res_b, 0);
		d_print_e_tran_strvec(2*qp->dim->nb+2*qp->dim->ng, ws->itref_res->res_d, 0);
#endif

		// alpha
		COMPUTE_ALPHA_QP(cws);
		if(kk<ws->stat_max)
			ws->stat[5*kk+0] = cws->alpha;

		if(arg->pred_corr==1)
			{
			// mu_aff
			COMPUTE_MU_AFF_QP(cws);
			if(kk<ws->stat_max)
				ws->stat[5*kk+1] = cws->mu_aff;

			tmp = cws->mu_aff/cws->mu;
			cws->sigma = tmp*tmp*tmp;
			if(kk<ws->stat_max)
				ws->stat[5*kk+2] = cws->sigma;

			COMPUTE_CENTERING_CORRECTION_QP(cws);
//			COMPUTE_CENTERING_QP(cws);

			// solve kkt
			SOLVE_KKT_STEP_DENSE_QP(qp, ws);

#if 0
			COMPUTE_RES_DENSE_QP(ws->itref_qp, ws->step, ws->itref_res, ws->res_workspace);
			double itref_qp_norm[4] = {0,0,0,0};
			VECNRM_INF_LIBSTR(cws->nv, ws->itref_res->res_g, 0, &itref_qp_norm[0]);
			VECNRM_INF_LIBSTR(cws->ne, ws->itref_res->res_b, 0, &itref_qp_norm[1]);
			VECNRM_INF_LIBSTR(cws->nc, ws->itref_res->res_d, 0, &itref_qp_norm[2]);
			VECNRM_INF_LIBSTR(cws->nc, ws->itref_res->res_m, 0, &itref_qp_norm[3]);
			printf("\nkk = %d %e %e %e %e\n", kk, itref_qp_norm[0], itref_qp_norm[1], itref_qp_norm[2], itref_qp_norm[3]);
//			d_print_e_tran_strvec(qp->dim->nv, ws->itref_res->res_g, 0);
//			d_print_e_tran_strvec(qp->dim->ne, ws->itref_res->res_b, 0);
//			d_print_e_tran_strvec(2*qp->dim->nb+2*qp->dim->ng, ws->itref_res->res_d, 0);
#endif

			// alpha
			COMPUTE_ALPHA_QP(cws);
			if(kk<ws->stat_max)
				ws->stat[5*kk+3] = cws->alpha;

			// conditional predictor-corrector
#if 1

			// save mu_aff (from prediction step)
			REAL mu_aff0 = cws->mu_aff;

			// compute mu for predictor-corrector-centering
			COMPUTE_MU_AFF_QP(cws);

//			if(cws->mu_aff > cws->mu)
			if(cws->mu_aff > 2.0*mu_aff0)
				{

				// centering direction
				COMPUTE_CENTERING_QP(cws);

				// solve kkt
				SOLVE_KKT_STEP_DENSE_QP(qp, ws);

				// alpha
				COMPUTE_ALPHA_QP(cws);
				if(kk<ws->stat_max)
					ws->stat[5*kk+3] = cws->alpha;

				}

#endif // conditional predictor-corrector

			}

		//
		UPDATE_VAR_QP(cws);

		// compute residuals
		COMPUTE_RES_DENSE_QP(qp, qp_sol, ws->res, ws->res_workspace);
		cws->mu = ws->res->res_mu;
		if(kk<ws->stat_max)
			ws->stat[5*kk+4] = ws->res->res_mu;

		// compute infinity norm of residuals
		VECNRM_INF_LIBSTR(cws->nv, &str_res_g, 0, &qp_res[0]);
		VECNRM_INF_LIBSTR(cws->ne, &str_res_b, 0, &qp_res[1]);
		VECNRM_INF_LIBSTR(cws->nc, &str_res_d, 0, &qp_res[2]);
		VECNRM_INF_LIBSTR(cws->nc, &str_res_m, 0, &qp_res[3]);

#if 0
printf("%e %e %e %e %e\n", ws->stat[5*kk+0], ws->stat[5*kk+1], ws->stat[5*kk+2], ws->stat[5*kk+3], ws->stat[5*kk+4]);
#endif

#if 0
printf("%e %e %e %e\n", qp_res[0], qp_res[1], qp_res[2], qp_res[3]);
#endif

		}

	ws->iter = kk;

	// max iteration number reached
	if(kk==arg->iter_max)
		return 1;

	if(cws->alpha <= arg->alpha_min)
		return 2;

	return 0;

	}


