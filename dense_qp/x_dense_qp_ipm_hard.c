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



int MEMSIZE_IPM_HARD_DENSE_QP(struct DENSE_QP *qp, struct IPM_HARD_DENSE_QP_ARG *arg)
	{

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;

	int size = 0;

	size += 22*sizeof(struct STRVEC); // dv dpi dlam dt dt_lb dt_ub dt_lg dt_ug res_g res_b res_d res_d_lb res_d_ub res_d_lg res_d_ug res_m Qx qx lv tmp_nb tmp_ng0 tmp_ng1
	size += 4*sizeof(struct STRMAT); // Lv AL Le Ctx

	size += 1*SIZE_STRVEC(nb); // tmp_nb
	size += 2*SIZE_STRVEC(ng); // tmp_ng0 tmp_ng1
	size += 1*SIZE_STRVEC(nv); // lv
	size += 1*SIZE_STRMAT(nv+1, nv); // Lv
	size += 1*SIZE_STRMAT(ne, nv); // AL
	size += 1*SIZE_STRMAT(ne, ne); // Le
	size += 1*SIZE_STRMAT(nv+1, ng); // Ctx

	size += 1*sizeof(struct IPM_HARD_CORE_QP_WORKSPACE);
	size += 1*MEMSIZE_IPM_HARD_CORE_QP(qp->nv, qp->ne, qp->nb, qp->ng, arg->iter_max);

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

	}



void CREATE_IPM_HARD_DENSE_QP(struct DENSE_QP *qp, struct IPM_HARD_DENSE_QP_ARG *arg, struct IPM_HARD_DENSE_QP_WORKSPACE *workspace, void *mem)
	{

	workspace->memsize = MEMSIZE_IPM_HARD_DENSE_QP(qp, arg);


	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;


	// core struct
	struct IPM_HARD_CORE_QP_WORKSPACE *sr_ptr = mem;

	// core workspace
	workspace->core_workspace = sr_ptr;
	sr_ptr += 1;
	struct IPM_HARD_CORE_QP_WORKSPACE *rwork = workspace->core_workspace;


	// matrix struct
	struct STRMAT *sm_ptr = (struct STRMAT *) sr_ptr;

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

	workspace->dv = sv_ptr;
	sv_ptr += 1;
	workspace->dpi = sv_ptr;
	sv_ptr += 1;
	workspace->dlam = sv_ptr;
	sv_ptr += 1;
	workspace->dt = sv_ptr;
	sv_ptr += 1;
	workspace->dt_lb = sv_ptr;
	sv_ptr += 1;
	workspace->dt_ub = sv_ptr;
	sv_ptr += 1;
	workspace->dt_lg = sv_ptr;
	sv_ptr += 1;
	workspace->dt_ug = sv_ptr;
	sv_ptr += 1;
	workspace->res_g = sv_ptr;
	sv_ptr += 1;
	workspace->res_b = sv_ptr;
	sv_ptr += 1;
	workspace->res_d = sv_ptr;
	sv_ptr += 1;
	workspace->res_d_lb = sv_ptr;
	sv_ptr += 1;
	workspace->res_d_ub = sv_ptr;
	sv_ptr += 1;
	workspace->res_d_lg = sv_ptr;
	sv_ptr += 1;
	workspace->res_d_ug = sv_ptr;
	sv_ptr += 1;
	workspace->res_m = sv_ptr;
	sv_ptr += 1;
	workspace->Qx = sv_ptr;
	sv_ptr += 1;
	workspace->qx = sv_ptr;
	sv_ptr += 1;
	workspace->lv = sv_ptr;
	sv_ptr += 1;
	workspace->tmp_nb = sv_ptr;
	sv_ptr += 1;
	workspace->tmp_ng0 = sv_ptr;
	sv_ptr += 1;
	workspace->tmp_ng1 = sv_ptr;
	sv_ptr += 1;


	// align to typicl cache line size
	size_t s_ptr = (size_t) sv_ptr;
	s_ptr = (s_ptr+63)/64*64;


	// void stuf
	void *v_ptr = (void *) s_ptr;

	CREATE_STRMAT(nv+1, nv, workspace->Lv, v_ptr);
	v_ptr += workspace->Lv->memory_size;

	CREATE_STRMAT(ne, nv, workspace->AL, v_ptr);
	v_ptr += workspace->AL->memory_size;

	CREATE_STRMAT(ne, ne, workspace->Le, v_ptr);
	v_ptr += workspace->Le->memory_size;

	CREATE_STRMAT(nv+1, ng, workspace->Ctx, v_ptr);
	v_ptr += workspace->Ctx->memory_size;

	CREATE_STRVEC(nv, workspace->lv, v_ptr);
	v_ptr += workspace->lv->memory_size;

	CREATE_STRVEC(nb, workspace->tmp_nb, v_ptr);
	v_ptr += workspace->tmp_nb->memory_size;

	CREATE_STRVEC(ng, workspace->tmp_ng0, v_ptr);
	v_ptr += workspace->tmp_ng0->memory_size;

	CREATE_STRVEC(ng, workspace->tmp_ng1, v_ptr);
	v_ptr += workspace->tmp_ng1->memory_size;

	rwork->nv = nv;
	rwork->ne = ne;
	rwork->nb = nb;
	rwork->ng = ng;
	rwork->iter_max = arg->iter_max;
	CREATE_IPM_HARD_CORE_QP(rwork, v_ptr);
	v_ptr += workspace->core_workspace->memsize;

	rwork->alpha_min = arg->alpha_min;
	rwork->mu_max = arg->mu_max;
	rwork->mu0 = arg->mu0;
	rwork->nt_inv = 1.0/(2*nb+2*ng);


	// alias members of workspace and core_workspace
	CREATE_STRVEC(nv, workspace->dv, rwork->dv);
	CREATE_STRVEC(ne, workspace->dpi, rwork->dpi);
	CREATE_STRVEC(2*nb+2*ng, workspace->dlam, rwork->dlam);
	CREATE_STRVEC(2*nb+2*ng, workspace->dt, rwork->dt);
	CREATE_STRVEC(nb, workspace->dt_lb, rwork->dt_lb);
	CREATE_STRVEC(nb, workspace->dt_ub, rwork->dt_ub);
	CREATE_STRVEC(ng, workspace->dt_lg, rwork->dt_lg);
	CREATE_STRVEC(ng, workspace->dt_ug, rwork->dt_ug);
	CREATE_STRVEC(nv, workspace->res_g, rwork->res_g);
	CREATE_STRVEC(ne, workspace->res_b, rwork->res_b);
	CREATE_STRVEC(2*nb+2*ng, workspace->res_d, rwork->res_d);
	CREATE_STRVEC(nb, workspace->res_d_lb, rwork->res_d_lb);
	CREATE_STRVEC(nb, workspace->res_d_ub, rwork->res_d_ub);
	CREATE_STRVEC(ng, workspace->res_d_lg, rwork->res_d_lg);
	CREATE_STRVEC(ng, workspace->res_d_ug, rwork->res_d_ug);
	CREATE_STRVEC(2*nb+2*ng, workspace->res_m, rwork->res_m);
	CREATE_STRVEC(nb+ng, workspace->Qx, rwork->Qx);
	CREATE_STRVEC(nb+ng, workspace->qx, rwork->qx);
	workspace->stat = rwork->stat;

	return;

	}



void SOLVE_IPM_HARD_DENSE_QP(struct DENSE_QP *qp, struct DENSE_QP_SOL *qp_sol, struct IPM_HARD_DENSE_QP_WORKSPACE *ws)
	{

	struct IPM_HARD_CORE_QP_WORKSPACE *cws = ws->core_workspace;

	// alias qp vectors into qp
	cws->d = qp->d->pa; // TODO REMOVE
	cws->d_lb = qp->d_lb->pa;
	cws->d_ub = qp->d_ub->pa;
	cws->d_lg = qp->d_lg->pa;
	cws->d_ug = qp->d_ug->pa;

	// alias qp vectors into qp_sol
	cws->v = qp_sol->v->pa;
	cws->pi = qp_sol->pi->pa;
	cws->lam = qp_sol->lam_lb->pa;
	cws->lam_lb = qp_sol->lam_lb->pa;
	cws->lam_ub = qp_sol->lam_ub->pa;
	cws->lam_lg = qp_sol->lam_lg->pa;
	cws->lam_ug = qp_sol->lam_ug->pa;
	cws->t = qp_sol->t_lb->pa;
	cws->t_lb = qp_sol->t_lb->pa;
	cws->t_ub = qp_sol->t_ub->pa;
	cws->t_lg = qp_sol->t_lg->pa;
	cws->t_ug = qp_sol->t_ug->pa;

	if(cws->nb+cws->ng==0)
		{
		FACT_SOLVE_KKT_UNCONSTR_DENSE_QP(qp, qp_sol, ws);
		COMPUTE_RES_HARD_DENSE_QP(qp, qp_sol, ws);
		cws->mu = ws->res_mu;
		ws->iter = 0;
		return;
		}

	// init solver
	INIT_VAR_HARD_DENSE_QP(qp, qp_sol, ws);

	// compute residuals
	COMPUTE_RES_HARD_DENSE_QP(qp, qp_sol, ws);
	cws->mu = ws->res_mu;

	int kk;
	for(kk=0; kk<cws->iter_max & cws->mu>cws->mu_max; kk++)
		{

		// fact and solve kkt
		FACT_SOLVE_KKT_STEP_HARD_DENSE_QP(qp, ws);

		// alpha
		COMPUTE_ALPHA_HARD_QP(cws);
		cws->stat[5*kk+0] = cws->alpha;

		//
		UPDATE_VAR_HARD_QP(cws);

		// compute residuals
		COMPUTE_RES_HARD_DENSE_QP(qp, qp_sol, ws);
		cws->mu = ws->res_mu;
		cws->stat[5*kk+1] = ws->res_mu;

		}

	ws->iter = kk;

	return;

	}



void SOLVE_IPM2_HARD_DENSE_QP(struct DENSE_QP *qp, struct DENSE_QP_SOL *qp_sol, struct IPM_HARD_DENSE_QP_WORKSPACE *ws)
	{

	struct IPM_HARD_CORE_QP_WORKSPACE *cws = ws->core_workspace;

	// alias qp vectors into qp
	cws->d = qp->d->pa; // TODO REMOVE
	cws->d_lb = qp->d_lb->pa;
	cws->d_ub = qp->d_ub->pa;
	cws->d_lg = qp->d_lg->pa;
	cws->d_ug = qp->d_ug->pa;

	// alias qp vectors into qp_sol
	cws->v = qp_sol->v->pa;
	cws->pi = qp_sol->pi->pa;
	cws->lam = qp_sol->lam_lb->pa;
	cws->lam_lb = qp_sol->lam_lb->pa;
	cws->lam_ub = qp_sol->lam_ub->pa;
	cws->lam_lg = qp_sol->lam_lg->pa;
	cws->lam_ug = qp_sol->lam_ug->pa;
	cws->t = qp_sol->t_lb->pa;
	cws->t_lb = qp_sol->t_lb->pa;
	cws->t_ub = qp_sol->t_ub->pa;
	cws->t_lg = qp_sol->t_lg->pa;
	cws->t_ug = qp_sol->t_ug->pa;

	REAL tmp;

	if(cws->nb+cws->ng==0)
		{
		FACT_SOLVE_KKT_UNCONSTR_DENSE_QP(qp, qp_sol, ws);
		COMPUTE_RES_HARD_DENSE_QP(qp, qp_sol, ws);
		cws->mu = ws->res_mu;
		ws->iter = 0;
		return;
		}

	// init solver
	INIT_VAR_HARD_DENSE_QP(qp, qp_sol, ws);

	// compute residuals
	COMPUTE_RES_HARD_DENSE_QP(qp, qp_sol, ws);
	cws->mu = ws->res_mu;

	int kk;
	for(kk=0; kk<cws->iter_max & cws->mu>cws->mu_max; kk++)
		{

		// fact and solve kkt
		FACT_SOLVE_KKT_STEP_HARD_DENSE_QP(qp, ws);

		// alpha
		COMPUTE_ALPHA_HARD_QP(cws);
		cws->stat[5*kk+0] = cws->alpha;

		// mu_aff
		COMPUTE_MU_AFF_HARD_QP(cws);
		cws->stat[5*kk+1] = cws->mu_aff;

		tmp = cws->mu_aff/cws->mu;
		cws->sigma = tmp*tmp*tmp;
		cws->stat[5*kk+2] = cws->sigma;

		COMPUTE_CENTERING_CORRECTION_HARD_QP(cws);

		// fact and solve kkt
		SOLVE_KKT_STEP_HARD_DENSE_QP(qp, ws);

		// alpha
		COMPUTE_ALPHA_HARD_QP(cws);
		cws->stat[5*kk+3] = cws->alpha;

		//
		UPDATE_VAR_HARD_QP(cws);

		// compute residuals
		COMPUTE_RES_HARD_DENSE_QP(qp, qp_sol, ws);
		cws->mu = ws->res_mu;
		cws->stat[5*kk+4] = ws->res_mu;

		}

	ws->iter = kk;

	return;

	}

