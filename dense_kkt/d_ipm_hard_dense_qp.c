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



#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_d_aux.h>

#include "../include/hpipm_d_dense_qp.h"
#include "../include/hpipm_d_dense_qp_sol.h"
#include "../include/hpipm_d_ipm_hard_dense_qp.h"
#include "../include/hpipm_d_ipm_hard_core_qp.h"
#include "../include/hpipm_d_aux_ipm_hard.h"
#include "../include/hpipm_d_kkt_dense_qp.h"



int d_memsize_ipm_hard_dense_qp(struct d_dense_qp *qp, struct d_ipm_hard_dense_qp_arg *arg)
	{

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;

	int size = 0;

	size += 22*sizeof(struct d_strvec); // dv dpi dlam dt dt_lb dt_ub dt_lg dt_ug res_g res_b res_d res_d_lb res_d_ub res_d_lg res_d_ug res_m Qx qx lv tmp_nb tmp_ng0 tmp_ng1
	size += 4*sizeof(struct d_strmat); // Lv AL Le Ctx

	size += 1*d_size_strvec(nb); // tmp_nb
	size += 2*d_size_strvec(ng); // tmp_ng0 tmp_ng1
	size += 1*d_size_strvec(nv); // lv
	size += 1*d_size_strmat(nv, nv); // Lv
	size += 1*d_size_strmat(ne, nv); // AL
	size += 1*d_size_strmat(ne, ne); // Le
	size += 1*d_size_strmat(nv, ng); // Ctx

	size += 1*sizeof(struct d_ipm_hard_core_qp_workspace);
	size += 1*d_memsize_ipm_hard_core_qp(qp->nv, qp->ne, qp->nb, qp->ng, arg->iter_max);

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

	}



void d_create_ipm_hard_dense_qp(struct d_dense_qp *qp, struct d_ipm_hard_dense_qp_arg *arg, struct d_ipm_hard_dense_qp_workspace *workspace, void *mem)
	{

	int nv = qp->nv;
	int ne = qp->ne;
	int nb = qp->nb;
	int ng = qp->ng;


	// core struct
	struct d_ipm_hard_core_qp_workspace *sr_ptr = mem;

	// core workspace
	workspace->core_workspace = sr_ptr;
	sr_ptr += 1;
	struct d_ipm_hard_core_qp_workspace *rwork = workspace->core_workspace;


	// matrix struct
	struct d_strmat *sm_ptr = (struct d_strmat *) sr_ptr;

	workspace->Lv = sm_ptr;
	sm_ptr += 1;
	workspace->AL = sm_ptr;
	sm_ptr += 1;
	workspace->Le = sm_ptr;
	sm_ptr += 1;
	workspace->Ctx = sm_ptr;
	sm_ptr += 1;


	// vector struct
	struct d_strvec *sv_ptr = (struct d_strvec *) sm_ptr;

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

	d_create_strmat(nv, nv, workspace->Lv, v_ptr);
	v_ptr += workspace->Lv->memory_size;

	d_create_strmat(ne, nv, workspace->AL, v_ptr);
	v_ptr += workspace->AL->memory_size;

	d_create_strmat(ne, ne, workspace->Le, v_ptr);
	v_ptr += workspace->Le->memory_size;

	d_create_strmat(nv, ng, workspace->Ctx, v_ptr);
	v_ptr += workspace->Ctx->memory_size;

	d_create_strvec(nv, workspace->lv, v_ptr);
	v_ptr += workspace->lv->memory_size;

	d_create_strvec(nb, workspace->tmp_nb, v_ptr);
	v_ptr += workspace->tmp_nb->memory_size;

	d_create_strvec(ng, workspace->tmp_ng0, v_ptr);
	v_ptr += workspace->tmp_ng0->memory_size;

	d_create_strvec(ng, workspace->tmp_ng1, v_ptr);
	v_ptr += workspace->tmp_ng1->memory_size;

	rwork->nv = nv;
	rwork->ne = ne;
	rwork->nb = nb;
	rwork->ng = ng;
	rwork->iter_max = arg->iter_max;
	d_create_ipm_hard_core_qp(rwork, v_ptr);
	v_ptr += workspace->core_workspace->memsize;

	rwork->alpha_min = arg->alpha_min;
	rwork->mu_max = arg->mu_max;
	rwork->mu0 = arg->mu0;


	// alias members of workspace and core_workspace
	d_create_strvec(nv, workspace->dv, rwork->dv);
	d_create_strvec(ne, workspace->dpi, rwork->dpi);
	d_create_strvec(2*nb+2*ng, workspace->dlam, rwork->dlam);
	d_create_strvec(2*nb+2*ng, workspace->dt, rwork->dt);
	d_create_strvec(nb, workspace->dt_lb, rwork->dt_lb);
	d_create_strvec(nb, workspace->dt_ub, rwork->dt_ub);
	d_create_strvec(ng, workspace->dt_lg, rwork->dt_lg);
	d_create_strvec(ng, workspace->dt_ug, rwork->dt_ug);
	d_create_strvec(nv, workspace->res_g, rwork->res_g);
	d_create_strvec(ne, workspace->res_b, rwork->res_b);
	d_create_strvec(2*nb+2*ng, workspace->res_d, rwork->res_d);
	d_create_strvec(nb, workspace->res_d_lb, rwork->res_d_lb);
	d_create_strvec(nb, workspace->res_d_ub, rwork->res_d_ub);
	d_create_strvec(ng, workspace->res_d_lg, rwork->res_d_lg);
	d_create_strvec(ng, workspace->res_d_ug, rwork->res_d_ug);
	d_create_strvec(2*nb+2*ng, workspace->res_m, rwork->res_m);
	d_create_strvec(nb+ng, workspace->Qx, rwork->Qx);
	d_create_strvec(nb+ng, workspace->qx, rwork->qx);
	workspace->stat = rwork->stat;


	workspace->nt_inv = 1.0/(2*nb+2*ng);

	return;

	}



void d_solve_ipm_hard_dense_qp(struct d_dense_qp *qp, struct d_dense_qp_sol *qp_sol, struct d_ipm_hard_dense_qp_workspace *ws)
	{

	struct d_ipm_hard_core_qp_workspace *cws = ws->core_workspace;

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

	// init solver
	d_init_var_hard_dense_qp(qp, qp_sol, ws);

#if 0
	// XXX hard-code solution for debug
	ws->v->pa[0] = 0.25;
	ws->v->pa[1] = 0.75;
	ws->pi->pa[0] = 2.75;
	ws->lam_lb->pa[0] = 0.0;
	ws->lam_lb->pa[1] = 0.0;
	ws->lam_ub->pa[0] = 0.0;
	ws->lam_ub->pa[1] = 0.0;
	ws->t_lb->pa[0] = 1.25;
	ws->t_lb->pa[1] = 1.75;
	ws->t_ub->pa[0] = 1.75;
	ws->t_ub->pa[1] = 1.25;
#endif

	// compute residuals
	d_compute_res_hard_dense_qp(qp, qp_sol, ws);
	cws->mu = ws->res_mu;

#if 0
	printf("\nres_g\n");
	d_print_tran_strvec(nv, ws->res_g, 0);
	printf("\nres_b\n");
	d_print_tran_strvec(ne, ws->res_b, 0);
	printf("\nres_d\n");
	d_print_tran_strvec(2*nb+2*ng, ws->res_d, 0);
	printf("\nres_m\n");
	d_print_tran_strvec(2*nb+2*ng, ws->res_m, 0);
#endif

	int kk;
	for(kk=0; kk<cws->iter_max & cws->mu>cws->mu_max; kk++)
		{

		// fact and solve kkt
		d_fact_solve_kkt_step_hard_dense_qp(qp, ws);

		// alpha
		d_compute_alpha_hard_qp(cws);
		cws->stat[5*kk+1] = cws->alpha;

		//
		d_update_var_hard_qp(cws);

#if 0
		printf("\nv\n");
		d_print_tran_strvec(nv, ws->v, 0);
		printf("\npi\n");
		d_print_tran_strvec(ne, ws->pi, 0);
		printf("\nlam\n");
		d_print_tran_strvec(2*nb+2*ng, ws->lam, 0);
		printf("\nt\n");
		d_print_tran_strvec(2*nb+2*ng, ws->t, 0);
#endif

		// compute residuals
		d_compute_res_hard_dense_qp(qp, qp_sol, ws);
		cws->mu = ws->res_mu;
		cws->stat[5*kk+2] = ws->res_mu;
#if 0

		printf("\nres_g\n");
		d_print_e_tran_strvec(nv, ws->res_g, 0);
		printf("\nres_b\n");
		d_print_e_tran_strvec(ne, ws->res_b, 0);
		printf("\nres_d\n");
		d_print_e_tran_strvec(2*nb+2*ng, ws->res_d, 0);
		printf("\nres_m\n");
		d_print_e_tran_strvec(2*nb+2*ng, ws->res_m, 0);
#endif

		}

	ws->iter = kk;

	return;

	}
