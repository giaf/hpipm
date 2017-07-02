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
#include <blasfeo_s_aux.h>

#include "../include/hpipm_d_ocp_qp.h"
#include "../include/hpipm_s_ocp_qp.h"
#include "../include/hpipm_d_ocp_qp_sol.h"
#include "../include/hpipm_d_ocp_qp_ipm_hard.h"
#include "../include/hpipm_m_ocp_qp_ipm_hard.h"
#include "../include/hpipm_d_core_qp_ipm_hard.h"
#include "../include/hpipm_d_core_qp_ipm_hard_aux.h"
#include "../include/hpipm_d_ocp_qp_kkt.h"
#include "../include/hpipm_m_ocp_qp_kkt.h"



int m_memsize_ipm_hard_ocp_qp(struct d_ocp_qp *qp, struct s_ocp_qp *s_qp, struct m_ipm_hard_ocp_qp_arg *arg)
	{

	// loop index
	int ii;

	// extract ocp qp size
	int N = qp->N;
	int *nx = qp->nx;
	int *nu = qp->nu;
	int *nb = qp->nb;
	int *ng = qp->ng;

	// compute core qp size and max size
	int nvt = 0;
	int net = 0;
	int nbt = 0;
	int ngt = 0;
	int nxM = 0;
	int nuM = 0;
	int nbM = 0;
	int ngM = 0;

	for(ii=0; ii<N; ii++)
		{
		nvt += nx[ii]+nu[ii];
		net += nx[ii+1];
		nbt += nb[ii];
		ngt += ng[ii];
		nxM = nx[ii]>nxM ? nx[ii] : nxM;
		nuM = nu[ii]>nuM ? nu[ii] : nuM;
		nbM = nb[ii]>nbM ? nb[ii] : nbM;
		ngM = ng[ii]>ngM ? ng[ii] : ngM;
		}
	ii = N;
	nvt += nx[ii]+nu[ii];
	nbt += nb[ii];
	ngt += ng[ii];
	nxM = nx[ii]>nxM ? nx[ii] : nxM;
	nuM = nu[ii]>nuM ? nu[ii] : nuM;
	nbM = nb[ii]>nbM ? nb[ii] : nbM;
	ngM = ng[ii]>ngM ? ng[ii] : ngM;

	int size = 0;

	size += (4+(N+1)*18)*sizeof(struct d_strvec); // dux dpi dt_lb dt_lg res_g res_b res_d res_d_lb res_d_ub res_d_lg res_d_ug res_m res_m_lb res_m_ub res_m_lg res_m_ug Qx_lb Qx_lg qx_lb qx_lg tmp_nbM tmp_ngM
	size += (1+(N+1)*9)*sizeof(struct s_strvec); // sdux sdpi sres_g sres_b sQx_lb sQx_lg, sqx_lb, sqx_lg tmp_nxM Pb
	size += (1+(N+1)*1)*sizeof(struct s_strmat); // L AL

	size += 1*d_size_strvec(nbM); // tmp_nbM
	size += 1*s_size_strvec(nxM); // tmp_nxM
	size += 2*d_size_strvec(nxM); // tmp_ngM
	for(ii=0; ii<=N; ii++) size += 2*s_size_strvec(nu[ii]+nx[ii]); // sdux sres_g
	for(ii=0; ii<N; ii++) size += 3*s_size_strvec(nx[ii+1]); // sdpi sres_b Pb
	for(ii=0; ii<=N; ii++) size += 2*s_size_strvec(nb[ii]); // sQx_lb sqx_lb
	for(ii=0; ii<=N; ii++) size += 2*s_size_strvec(ng[ii]); // sQx_lg sqx_lg
	for(ii=0; ii<=N; ii++) size += 1*s_size_strmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii]); // L
	size += 2*s_size_strmat(nuM+nxM+1, nxM+ngM); // AL

	size += 1*sizeof(struct d_ipm_hard_core_qp_workspace);
	size += 1*d_memsize_ipm_hard_core_qp(nvt, net, nbt, ngt, arg->iter_max);

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

	}



void m_create_ipm_hard_ocp_qp(struct d_ocp_qp *qp, struct s_ocp_qp *s_qp, struct m_ipm_hard_ocp_qp_arg *arg, struct m_ipm_hard_ocp_qp_workspace *workspace, void *mem)
	{

	// loop index
	int ii;

	// extract ocp qp size
	int N = qp->N;
	int *nx = qp->nx;
	int *nu = qp->nu;
	int *nb = qp->nb;
	int *ng = qp->ng;

	// compute core qp size and max size
	int nvt = 0;
	int net = 0;
	int nbt = 0;
	int ngt = 0;
	int nxM = 0;
	int nuM = 0;
	int nbM = 0;
	int ngM = 0;

	for(ii=0; ii<N; ii++)
		{
		nvt += nx[ii]+nu[ii];
		net += nx[ii+1];
		nbt += nb[ii];
		ngt += ng[ii];
		nxM = nx[ii]>nxM ? nx[ii] : nxM;
		nuM = nu[ii]>nuM ? nu[ii] : nuM;
		nbM = nb[ii]>nbM ? nb[ii] : nbM;
		ngM = ng[ii]>ngM ? ng[ii] : ngM;
		}
	ii = N;
	nvt += nx[ii]+nu[ii];
	nbt += nb[ii];
	ngt += ng[ii];
	nxM = nx[ii]>nxM ? nx[ii] : nxM;
	nuM = nu[ii]>nuM ? nu[ii] : nuM;
	nbM = nb[ii]>nbM ? nb[ii] : nbM;
	ngM = ng[ii]>ngM ? ng[ii] : ngM;


	// core struct
	struct d_ipm_hard_core_qp_workspace *sr_ptr = mem;

	// core workspace
	workspace->core_workspace = sr_ptr;
	sr_ptr += 1;
	struct d_ipm_hard_core_qp_workspace *rwork = workspace->core_workspace;


	// s matrix struct
	struct s_strmat *sm_ptr = (struct s_strmat *) sr_ptr;

	workspace->L = sm_ptr;
	sm_ptr += N+1;
	workspace->AL = sm_ptr;
	sm_ptr += 2;


	// s vector struct
	struct s_strvec *ssv_ptr = (struct s_strvec *) sm_ptr;

	workspace->sdux = ssv_ptr;
	ssv_ptr += N+1;
	workspace->sdpi = ssv_ptr;
	ssv_ptr += N+1;
	workspace->sres_g = ssv_ptr;
	ssv_ptr += N+1;
	workspace->sres_b = ssv_ptr;
	ssv_ptr += N+1;
	workspace->sQx_lb = ssv_ptr;
	ssv_ptr += N+1;
	workspace->sQx_lg = ssv_ptr;
	ssv_ptr += N+1;
	workspace->sqx_lb = ssv_ptr;
	ssv_ptr += N+1;
	workspace->sqx_lg = ssv_ptr;
	ssv_ptr += N+1;
	workspace->Pb = ssv_ptr;
	ssv_ptr += N+1;
	workspace->tmp_nxM = ssv_ptr;
	ssv_ptr += 1;


	// d vector struct
	struct d_strvec *sv_ptr = (struct d_strvec *) ssv_ptr;

	workspace->dux = sv_ptr;
	sv_ptr += N+1;
	workspace->dpi = sv_ptr;
	sv_ptr += N+1;
	workspace->dt_lb = sv_ptr;
	sv_ptr += N+1;
	workspace->dt_lg = sv_ptr;
	sv_ptr += N+1;
	workspace->res_g = sv_ptr;
	sv_ptr += N+1;
	workspace->res_b = sv_ptr;
	sv_ptr += N+1;
	workspace->res_d = sv_ptr;
	sv_ptr += 1;
	workspace->res_d_lb = sv_ptr;
	sv_ptr += N+1;
	workspace->res_d_ub = sv_ptr;
	sv_ptr += N+1;
	workspace->res_d_lg = sv_ptr;
	sv_ptr += N+1;
	workspace->res_d_ug = sv_ptr;
	sv_ptr += N+1;
	workspace->res_m = sv_ptr;
	sv_ptr += 1;
	workspace->res_m_lb = sv_ptr;
	sv_ptr += N+1;
	workspace->res_m_ub = sv_ptr;
	sv_ptr += N+1;
	workspace->res_m_lg = sv_ptr;
	sv_ptr += N+1;
	workspace->res_m_ug = sv_ptr;
	sv_ptr += N+1;
	workspace->Qx_lb = sv_ptr;
	sv_ptr += N+1;
	workspace->Qx_lg = sv_ptr;
	sv_ptr += N+1;
	workspace->qx_lb = sv_ptr;
	sv_ptr += N+1;
	workspace->qx_lg = sv_ptr;
	sv_ptr += N+1;
//	workspace->Pb = sv_ptr;
//	sv_ptr += N+1;
	workspace->tmp_nbM = sv_ptr;
	sv_ptr += 1;
//	workspace->tmp_nxM = sv_ptr;
//	sv_ptr += 1;
	workspace->tmp_ngM = sv_ptr;
	sv_ptr += 2;


	// align to typicl cache line size
	size_t s_ptr = (size_t) sv_ptr;
	s_ptr = (s_ptr+63)/64*64;


	// void stuf
	void *v_ptr = (void *) s_ptr;

	for(ii=0; ii<=N; ii++)
		{
		s_create_strmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], workspace->L+ii, v_ptr);
		v_ptr += (workspace->L+ii)->memory_size;
		}

	s_create_strmat(nuM+nxM+1, nxM+ngM, workspace->AL+0, v_ptr);
	v_ptr += (workspace->AL+0)->memory_size;

	s_create_strmat(nuM+nxM+1, nxM+ngM, workspace->AL+1, v_ptr);
	v_ptr += (workspace->AL+1)->memory_size;

	for(ii=0; ii<=N; ii++)
		{
		s_create_strvec(nu[ii]+nx[ii], workspace->sdux+ii, v_ptr);
		v_ptr += (workspace->sdux+ii)->memory_size;
		}

	for(ii=0; ii<N; ii++)
		{
		s_create_strvec(nx[ii+1], workspace->sdpi+ii, v_ptr);
		v_ptr += (workspace->sdpi+ii)->memory_size;
		}

	for(ii=0; ii<=N; ii++)
		{
		s_create_strvec(nu[ii]+nx[ii], workspace->sres_g+ii, v_ptr);
		v_ptr += (workspace->sdux+ii)->memory_size;
		}

	for(ii=0; ii<N; ii++)
		{
		s_create_strvec(nx[ii+1], workspace->sres_b+ii, v_ptr);
		v_ptr += (workspace->sdpi+ii)->memory_size;
		}

	for(ii=0; ii<N; ii++)
		{
		s_create_strvec(nx[ii+1], workspace->Pb+ii, v_ptr);
		v_ptr += (workspace->Pb+ii)->memory_size;
		}

	for(ii=0; ii<=N; ii++)
		{
		s_create_strvec(nb[ii], workspace->sQx_lb+ii, v_ptr);
		v_ptr += (workspace->sQx_lb+ii)->memory_size;
		}

	for(ii=0; ii<=N; ii++)
		{
		s_create_strvec(nb[ii], workspace->sqx_lb+ii, v_ptr);
		v_ptr += (workspace->sqx_lb+ii)->memory_size;
		}

	for(ii=0; ii<=N; ii++)
		{
		s_create_strvec(ng[ii], workspace->sQx_lg+ii, v_ptr);
		v_ptr += (workspace->sQx_lg+ii)->memory_size;
		}

	for(ii=0; ii<=N; ii++)
		{
		s_create_strvec(ng[ii], workspace->sqx_lg+ii, v_ptr);
		v_ptr += (workspace->sqx_lg+ii)->memory_size;
		}

	d_create_strvec(nbM, workspace->tmp_nbM, v_ptr);
	v_ptr += workspace->tmp_nbM->memory_size;

	s_create_strvec(nxM, workspace->tmp_nxM, v_ptr);
	v_ptr += workspace->tmp_nxM->memory_size;

	d_create_strvec(ngM, workspace->tmp_ngM+0, v_ptr);
	v_ptr += (workspace->tmp_ngM+0)->memory_size;

	d_create_strvec(ngM, workspace->tmp_ngM+1, v_ptr);
	v_ptr += (workspace->tmp_ngM+1)->memory_size;



	rwork->nv = nvt;
	rwork->ne = net;
	rwork->nb = nbt;
	rwork->ng = ngt;
	rwork->iter_max = arg->iter_max;
	d_create_ipm_hard_core_qp(rwork, v_ptr);
	v_ptr += workspace->core_workspace->memsize;

	rwork->alpha_min = arg->alpha_min;
	rwork->mu_max = arg->mu_max;
	rwork->mu0 = arg->mu0;
	rwork->nt_inv = 1.0/(2*nbt+2*ngt);


	// alias members of workspace and core_workspace
	v_ptr = rwork->dv;
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(nu[ii]+nx[ii], workspace->dux+ii, v_ptr);
		v_ptr += (nu[ii]+nx[ii])*sizeof(double);
		}
	v_ptr = rwork->dpi;
	for(ii=0; ii<N; ii++)
		{
		d_create_strvec(nx[ii+1], workspace->dpi+ii, v_ptr);
		v_ptr += (nx[ii+1])*sizeof(double);
		}
	v_ptr = rwork->dt_lb;
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(nb[ii], workspace->dt_lb+ii, v_ptr);
		v_ptr += (nb[ii])*sizeof(double);
		}
	v_ptr = rwork->dt_lg;
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(ng[ii], workspace->dt_lg+ii, v_ptr);
		v_ptr += (ng[ii])*sizeof(double);
		}
	v_ptr = rwork->res_g;
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(nu[ii]+nx[ii], workspace->res_g+ii, v_ptr);
		v_ptr += (nu[ii]+nx[ii])*sizeof(double);
		}
	v_ptr = rwork->res_b;
	for(ii=0; ii<N; ii++)
		{
		d_create_strvec(nx[ii+1], workspace->res_b+ii, v_ptr);
		v_ptr += (nx[ii+1])*sizeof(double);
		}
	v_ptr = rwork->res_d;
	d_create_strvec(2*nbt+2*ngt, workspace->res_d, v_ptr);
	v_ptr = rwork->res_d_lb;
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(nb[ii], workspace->res_d_lb+ii, v_ptr);
		v_ptr += (nb[ii])*sizeof(double);
		}
	v_ptr = rwork->res_d_ub;
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(nb[ii], workspace->res_d_ub+ii, v_ptr);
		v_ptr += (nb[ii])*sizeof(double);
		}
	v_ptr = rwork->res_d_lg;
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(ng[ii], workspace->res_d_lg+ii, v_ptr);
		v_ptr += (ng[ii])*sizeof(double);
		}
	v_ptr = rwork->res_d_ug;
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(ng[ii], workspace->res_d_ug+ii, v_ptr);
		v_ptr += (ng[ii])*sizeof(double);
		}
	v_ptr = rwork->res_m;
	d_create_strvec(2*nbt+2*ngt, workspace->res_m, v_ptr);
	v_ptr = rwork->res_m_lb;
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(nb[ii], workspace->res_m_lb+ii, v_ptr);
		v_ptr += (nb[ii])*sizeof(double);
		}
	v_ptr = rwork->res_m_ub;
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(nb[ii], workspace->res_m_ub+ii, v_ptr);
		v_ptr += (nb[ii])*sizeof(double);
		}
	v_ptr = rwork->res_m_lg;
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(ng[ii], workspace->res_m_lg+ii, v_ptr);
		v_ptr += (ng[ii])*sizeof(double);
		}
	v_ptr = rwork->res_m_ug;
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(ng[ii], workspace->res_m_ug+ii, v_ptr);
		v_ptr += (ng[ii])*sizeof(double);
		}
	v_ptr = rwork->Qx_lb;
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(nb[ii], workspace->Qx_lb+ii, v_ptr);
		v_ptr += (nb[ii])*sizeof(double);
		}
	v_ptr = rwork->Qx_lg;
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(ng[ii], workspace->Qx_lg+ii, v_ptr);
		v_ptr += (ng[ii])*sizeof(double);
		}
	v_ptr = rwork->qx_lb;
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(nb[ii], workspace->qx_lb+ii, v_ptr);
		v_ptr += (nb[ii])*sizeof(double);
		}
	v_ptr = rwork->qx_lg;
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(ng[ii], workspace->qx_lg+ii, v_ptr);
		v_ptr += (ng[ii])*sizeof(double);
		}
	workspace->stat = rwork->stat;

	return;

	}



void m_solve_ipm_hard_ocp_qp(struct d_ocp_qp *qp, struct s_ocp_qp *s_qp, struct d_ocp_qp_sol *qp_sol, struct m_ipm_hard_ocp_qp_workspace *ws)
	{

	struct d_ipm_hard_core_qp_workspace *cws = ws->core_workspace;

	// alias d_ocp_workspace to m_ocp_workspace
	struct d_ipm_hard_ocp_qp_workspace dws;
	dws.core_workspace = ws->core_workspace;
	dws.dux = ws->dux;
	dws.dpi = ws->dpi;
	dws.dt_lb = ws->dt_lb;
	dws.dt_lg = ws->dt_lg;
	dws.res_g = ws->res_g;
	dws.res_b = ws->res_b;
	dws.res_d = ws->res_d;
	dws.res_d_lb = ws->res_d_lb;
	dws.res_d_ub = ws->res_d_ub;
	dws.res_d_lg = ws->res_d_lg;
	dws.res_d_ug = ws->res_d_ug;
	dws.res_m = ws->res_m;
	dws.res_m_lb = ws->res_m_lb;
	dws.res_m_ub = ws->res_m_ub;
	dws.res_m_lg = ws->res_m_lg;
	dws.res_m_ug = ws->res_m_ug;
	dws.tmp_nbM = ws->tmp_nbM;
	dws.tmp_ngM = ws->tmp_ngM;

	// alias qp vectors into qp
	cws->d_lb = qp->d_lb->pa;
	cws->d_ub = qp->d_ub->pa;
	cws->d_lg = qp->d_lg->pa;
	cws->d_ug = qp->d_ug->pa;

	// alias qp vectors into qp_sol
	cws->v = qp_sol->ux->pa;
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
		//
		d_init_var_hard_ocp_qp(qp, qp_sol, &dws);
		//
		d_compute_res_hard_ocp_qp(qp, qp_sol, &dws);
		cws->mu = dws.res_mu;
		ws->res_mu = dws.res_mu;
		//
		m_fact_solve_kkt_step_hard_ocp_qp(qp, s_qp, ws);
		//
		cws->alpha = 1.0;
		d_update_var_hard_qp(cws);
		//
		d_compute_res_hard_ocp_qp(qp, qp_sol, &dws);
		cws->mu = dws.res_mu;
		ws->res_mu = dws.res_mu;
		//
		ws->compute_Pb = 1;
		m_solve_kkt_step_hard_ocp_qp(qp, s_qp, ws);
		//
		cws->alpha = 1.0;
		d_update_var_hard_qp(cws);
		//
		d_compute_res_hard_ocp_qp(qp, qp_sol, &dws);
		cws->mu = dws.res_mu;
		ws->res_mu = dws.res_mu;
		//
		ws->iter = 0;
		return;
		}

	// init solver
	d_init_var_hard_ocp_qp(qp, qp_sol, &dws);

	// compute residuals
	d_compute_res_hard_ocp_qp(qp, qp_sol, &dws);
	cws->mu = dws.res_mu;
	ws->res_mu = dws.res_mu;

	int kk;
	for(kk=0; kk<cws->iter_max & cws->mu>cws->mu_max; kk++)
		{

		// fact and solve kkt
		m_fact_solve_kkt_step_hard_ocp_qp(qp, s_qp, ws);

		// alpha
		d_compute_alpha_hard_qp(cws);
		cws->stat[5*kk+0] = cws->alpha;

		//
		d_update_var_hard_qp(cws);

		// compute residuals
		d_compute_res_hard_ocp_qp(qp, qp_sol, &dws);
		cws->mu = dws.res_mu;
		ws->res_mu = dws.res_mu;
		cws->stat[5*kk+1] = ws->res_mu;

		}
	
	ws->iter = kk;
	
	return;

	}



void m_solve_ipm2_hard_ocp_qp(struct d_ocp_qp *qp, struct s_ocp_qp *s_qp, struct d_ocp_qp_sol *qp_sol, struct m_ipm_hard_ocp_qp_workspace *ws)
	{

	struct d_ipm_hard_core_qp_workspace *cws = ws->core_workspace;

	// alias d_ocp_workspace to m_ocp_workspace
	struct d_ipm_hard_ocp_qp_workspace dws;
	dws.core_workspace = ws->core_workspace;
	dws.dux = ws->dux;
	dws.dpi = ws->dpi;
	dws.dt_lb = ws->dt_lb;
	dws.dt_lg = ws->dt_lg;
	dws.res_g = ws->res_g;
	dws.res_b = ws->res_b;
	dws.res_d = ws->res_d;
	dws.res_d_lb = ws->res_d_lb;
	dws.res_d_ub = ws->res_d_ub;
	dws.res_d_lg = ws->res_d_lg;
	dws.res_d_ug = ws->res_d_ug;
	dws.res_m = ws->res_m;
	dws.res_m_lb = ws->res_m_lb;
	dws.res_m_ub = ws->res_m_ub;
	dws.res_m_lg = ws->res_m_lg;
	dws.res_m_ug = ws->res_m_ug;
	dws.tmp_nbM = ws->tmp_nbM;
	dws.tmp_ngM = ws->tmp_ngM;

	// alias qp vectors into qp
	cws->d_lb = qp->d_lb->pa;
	cws->d_ub = qp->d_ub->pa;
	cws->d_lg = qp->d_lg->pa;
	cws->d_ug = qp->d_ug->pa;

	// alias qp vectors into qp_sol
	cws->v = qp_sol->ux->pa;
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

	double tmp;

	if(cws->nb+cws->ng==0)
		{
		//
		d_init_var_hard_ocp_qp(qp, qp_sol, &dws);
		//
		d_compute_res_hard_ocp_qp(qp, qp_sol, &dws);
		cws->mu = dws.res_mu;
		ws->res_mu = dws.res_mu;
		//
		m_fact_solve_kkt_step_hard_ocp_qp(qp, s_qp, ws);
		//
		cws->alpha = 1.0;
		d_update_var_hard_qp(cws);
		//
		d_compute_res_hard_ocp_qp(qp, qp_sol, &dws);
		cws->mu = dws.res_mu;
		ws->res_mu = dws.res_mu;
		//
		ws->compute_Pb = 1;
		m_solve_kkt_step_hard_ocp_qp(qp, s_qp, ws);
		//
		cws->alpha = 1.0;
		d_update_var_hard_qp(cws);
		//
		d_compute_res_hard_ocp_qp(qp, qp_sol, &dws);
		cws->mu = dws.res_mu;
		ws->res_mu = dws.res_mu;
		//
		ws->iter = 0;
		return;
		}

	// init solver
	d_init_var_hard_ocp_qp(qp, qp_sol, &dws);

	// compute residuals
	d_compute_res_hard_ocp_qp(qp, qp_sol, &dws);
	cws->mu = dws.res_mu;
	ws->res_mu = dws.res_mu;

#if 0
	printf("\nres_g\n");
	for(int ii=0; ii<=qp->N; ii++)
		{
		PRINT_E_TRAN_STRVEC(qp->nx[ii]+qp->nu[ii], ws->res_g+ii, 0);
		}
	printf("\nres_b\n");
	for(int ii=0; ii<qp->N; ii++)
		{
		PRINT_E_TRAN_STRVEC(qp->nx[ii+1], ws->res_b+ii, 0);
		}
	printf("\nres_d_lb\n");
	for(int ii=0; ii<qp->N; ii++)
		{
		PRINT_E_TRAN_STRVEC(qp->nb[ii], ws->res_d_lb+ii, 0);
		}
	printf("\nres_d_ub\n");
	for(int ii=0; ii<qp->N; ii++)
		{
		PRINT_E_TRAN_STRVEC(qp->nb[ii], ws->res_d_ub+ii, 0);
		}
	printf("\nres_d_lg\n");
	for(int ii=0; ii<qp->N; ii++)
		{
		PRINT_E_TRAN_STRVEC(qp->ng[ii], ws->res_d_lg+ii, 0);
		}
	printf("\nres_d_ug\n");
	for(int ii=0; ii<qp->N; ii++)
		{
		PRINT_E_TRAN_STRVEC(qp->ng[ii], ws->res_d_ug+ii, 0);
		}
	printf("\nres_m_lb\n");
	for(int ii=0; ii<qp->N; ii++)
		{
		PRINT_E_TRAN_STRVEC(qp->nb[ii], ws->res_m_lb+ii, 0);
		}
	printf("\nres_m_ub\n");
	for(int ii=0; ii<qp->N; ii++)
		{
		PRINT_E_TRAN_STRVEC(qp->nb[ii], ws->res_m_ub+ii, 0);
		}
	printf("\nres_m_lg\n");
	for(int ii=0; ii<qp->N; ii++)
		{
		PRINT_E_TRAN_STRVEC(qp->ng[ii], ws->res_m_lg+ii, 0);
		}
	printf("\nres_m_ug\n");
	for(int ii=0; ii<qp->N; ii++)
		{
		PRINT_E_TRAN_STRVEC(qp->ng[ii], ws->res_m_ug+ii, 0);
		}
	exit(1);
#endif

#if 0
	int ii;
	for(ii=0; ii<1; ii++)
		{
		cws->sigma = 1.0;
		cws->stat[5*kk+2] = cws->sigma;
		COMPUTE_CENTERING_CORRECTION_HARD_QP(cws);
		FACT_SOLVE_KKT_STEP_HARD_OCP_QP(qp, ws);
		COMPUTE_ALPHA_HARD_QP(cws);
		cws->stat[5*kk+3] = cws->alpha;
		UPDATE_VAR_HARD_QP(cws);
		COMPUTE_RES_HARD_OCP_QP(qp, qp_sol, ws);
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

		// fact and solve kkt
		m_fact_solve_kkt_step_hard_ocp_qp(qp, s_qp, ws);

#if 0
	printf("\ndux\n");
	for(int ii=0; ii<=qp->N; ii++)
		{
		PRINT_E_TRAN_STRVEC(qp->nx[ii]+qp->nu[ii], ws->dux+ii, 0);
		}
	printf("\ndpi\n");
	for(int ii=0; ii<qp->N; ii++)
		{
		PRINT_E_TRAN_STRVEC(qp->nx[ii+1], ws->dpi+ii, 0);
		}
	printf("\ndt\n");
	for(int ii=0; ii<qp->N; ii++)
		{
		PRINT_E_TRAN_STRVEC(2*qp->nb[ii]+2*qp->ng[ii], ws->dt_lb+ii, 0);
		}
	exit(1);
#endif
		// alpha
		d_compute_alpha_hard_qp(cws);
		cws->stat[5*kk+0] = cws->alpha;

		// mu_aff
		d_compute_mu_aff_hard_qp(cws);
		cws->stat[5*kk+1] = cws->mu_aff;

		tmp = cws->mu_aff/cws->mu;
		cws->sigma = tmp*tmp*tmp;
		cws->stat[5*kk+2] = cws->sigma;

		d_compute_centering_correction_hard_qp(cws);

		// fact and solve kkt
		ws->compute_Pb = 0;
		m_solve_kkt_step_hard_ocp_qp(qp, s_qp, ws);

#if 0
int ii;
for(ii=0; ii<=qp->N; ii++)
	d_print_tran_strvec(qp->nu[ii]+qp->nx[ii], ws->dux+ii, 0);
for(ii=0; ii<qp->N; ii++)
	d_print_tran_strvec(qp->nx[ii+1], ws->dpi+ii, 0);
exit(1);
#endif
		// alpha
		d_compute_alpha_hard_qp(cws);
		cws->stat[5*kk+3] = cws->alpha;

		//
		d_update_var_hard_qp(cws);

		// compute residuals
		d_compute_res_hard_ocp_qp(qp, qp_sol, &dws);
		cws->mu = dws.res_mu;
		ws->res_mu = dws.res_mu;
		cws->stat[5*kk+4] = ws->res_mu;

		}
	
	ws->iter = kk;
	
	return;

	}




