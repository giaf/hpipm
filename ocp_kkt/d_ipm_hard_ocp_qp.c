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

#include "../include/hpipm_d_ocp_qp.h"
#include "../include/hpipm_d_ipm_hard_ocp_qp.h"
#include "../include/hpipm_d_ipm_hard_core_qp.h"
#include "../include/hpipm_d_aux_ipm_hard.h"
#include "../include/hpipm_d_kkt_ocp_qp.h"



int d_memsize_ipm_hard_ocp_qp(struct d_ocp_qp *qp, struct d_ipm_hard_ocp_qp_arg *arg)
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

	size += (6+(N+1)*27)*sizeof(struct d_strvec); // ux pi lam lam_lb lam_ub lam_lg lam_ug t t_lb t_ub t_lg t_ug dux dpi dt_lb dt_lg res_g res_b res_d res_d_lb res_d_ub res_d_lg res_d_ug res_m res_m_lb res_m_ub res_m_lg res_m_ug Qx_lb qx_lb Pb tmp_nbM tmp_nxM
	size += (2+(N+1)*1)*sizeof(struct d_strmat); // L AL0 AL1

	size += 1*d_size_strvec(nbM); // tmp_nbM
	size += 1*d_size_strvec(nxM); // tmp_nxM
	for(ii=0; ii<N; ii++) size += 1*d_size_strvec(nx[ii+1]);
	for(ii=0; ii<=N; ii++) size += 1*d_size_strmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii]); // L
	size += 2*d_size_strmat(nuM+nxM+1, nxM+ngM); // AL0 AL1

	size += 1*sizeof(struct d_ipm_hard_core_qp_workspace);
	size += 1*d_memsize_ipm_hard_core_qp(nvt, net, nbt, ngt, arg->iter_max);

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

	}



void d_create_ipm_hard_ocp_qp(struct d_ocp_qp *qp, struct d_ipm_hard_ocp_qp_arg *arg, struct d_ipm_hard_ocp_qp_workspace *workspace, void *mem)
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


	// matrix struct
	struct d_strmat *sm_ptr = (struct d_strmat *) sr_ptr;

	workspace->L = sm_ptr;
	sm_ptr += N+1;
	workspace->AL0 = sm_ptr;
	sm_ptr += 1;
	workspace->AL1 = sm_ptr;
	sm_ptr += 1;


	// vector struct
	struct d_strvec *sv_ptr = (struct d_strvec *) sm_ptr;

	workspace->ux = sv_ptr;
	sv_ptr += N+1;
	workspace->pi = sv_ptr;
	sv_ptr += N+1;
	workspace->lam = sv_ptr;
	sv_ptr += 1;
	workspace->lam_lb = sv_ptr;
	sv_ptr += N+1;
	workspace->lam_ub = sv_ptr;
	sv_ptr += N+1;
	workspace->lam_lg = sv_ptr;
	sv_ptr += N+1;
	workspace->lam_ug = sv_ptr;
	sv_ptr += N+1;
	workspace->t = sv_ptr;
	sv_ptr += 1;
	workspace->t_lb = sv_ptr;
	sv_ptr += N+1;
	workspace->t_ub = sv_ptr;
	sv_ptr += N+1;
	workspace->t_lg = sv_ptr;
	sv_ptr += N+1;
	workspace->t_ug = sv_ptr;
	sv_ptr += N+1;
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
	workspace->qx_lb = sv_ptr;
	sv_ptr += N+1;
	workspace->Pb = sv_ptr;
	sv_ptr += N+1;
	workspace->tmp_nbM = sv_ptr;
	sv_ptr += 1;
	workspace->tmp_nxM = sv_ptr;
	sv_ptr += 1;


	// align to typicl cache line size
	size_t s_ptr = (size_t) sv_ptr;
	s_ptr = (s_ptr+63)/64*64;


	// void stuf
	void *v_ptr = (void *) s_ptr;

	for(ii=0; ii<=N; ii++)
		{
		d_create_strmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], workspace->L+ii, v_ptr);
		v_ptr += (workspace->L+ii)->memory_size;
		}

	d_create_strmat(nuM+nxM+1, nxM+ngM, workspace->AL0, v_ptr);
	v_ptr += workspace->AL0->memory_size;

	d_create_strmat(nuM+nxM+1, nxM+ngM, workspace->AL1, v_ptr);
	v_ptr += workspace->AL1->memory_size;

	for(ii=0; ii<N; ii++)
		{
		d_create_strvec(nx[ii+1], workspace->Pb+ii, v_ptr);
		v_ptr += (workspace->Pb+ii)->memory_size;
		}

	d_create_strvec(nbM, workspace->tmp_nbM, v_ptr);
	v_ptr += workspace->tmp_nbM->memory_size;

	d_create_strvec(nxM, workspace->tmp_nxM, v_ptr);
	v_ptr += workspace->tmp_nxM->memory_size;



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


	// alias members of workspace and core_workspace
	v_ptr = rwork->v;
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(nu[ii]+nx[ii], workspace->ux+ii, v_ptr);
		v_ptr += (nu[ii]+nx[ii])*sizeof(double);
		}
	v_ptr = rwork->pi;
	for(ii=0; ii<N; ii++)
		{
		d_create_strvec(nx[ii+1], workspace->pi+ii, v_ptr);
		v_ptr += (nx[ii+1])*sizeof(double);
		}
	v_ptr = rwork->lam;
	d_create_strvec(2*nbt+2*ngt, workspace->lam, v_ptr);
	v_ptr = rwork->lam_lb;
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(nb[ii], workspace->lam_lb+ii, v_ptr);
		v_ptr += (nb[ii])*sizeof(double);
		}
	v_ptr = rwork->lam_ub;
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(nb[ii], workspace->lam_ub+ii, v_ptr);
		v_ptr += (nb[ii])*sizeof(double);
		}
	v_ptr = rwork->lam_lg;
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(ng[ii], workspace->lam_lg+ii, v_ptr);
		v_ptr += (ng[ii])*sizeof(double);
		}
	v_ptr = rwork->lam_ug;
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(ng[ii], workspace->lam_ug+ii, v_ptr);
		v_ptr += (ng[ii])*sizeof(double);
		}
	v_ptr = rwork->t;
	d_create_strvec(2*nbt+2*ngt, workspace->t, v_ptr);
	v_ptr = rwork->t_lb;
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(nb[ii], workspace->t_lb+ii, v_ptr);
		v_ptr += (nb[ii])*sizeof(double);
		}
	v_ptr = rwork->t_ub;
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(nb[ii], workspace->t_ub+ii, v_ptr);
		v_ptr += (nb[ii])*sizeof(double);
		}
	v_ptr = rwork->t_lg;
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(ng[ii], workspace->t_lg+ii, v_ptr);
		v_ptr += (ng[ii])*sizeof(double);
		}
	v_ptr = rwork->t_ug;
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(ng[ii], workspace->t_ug+ii, v_ptr);
		v_ptr += (ng[ii])*sizeof(double);
		}
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
	v_ptr = rwork->qx_lb;
	for(ii=0; ii<=N; ii++)
		{
		d_create_strvec(nb[ii], workspace->qx_lb+ii, v_ptr);
		v_ptr += (nb[ii])*sizeof(double);
		}
	workspace->stat = rwork->stat;



	workspace->mu0 = arg->mu0;
	workspace->nt_inv = 1.0/(2*nbt+2*ngt);

	return;

	}



void d_solve_ipm_hard_ocp_qp(struct d_ocp_qp *qp, struct d_ipm_hard_ocp_qp_workspace *ws)
	{

	struct d_ipm_hard_core_qp_workspace *cws = ws->core_workspace;

	// alias qp vectors into core workspace

	// init solver
	d_init_var_hard_ocp_qp(qp, ws);

	// compute residuals
	d_compute_res_hard_ocp_qp(qp, ws);
	cws->mu = ws->res_mu;

	int kk;
	for(kk=0; kk<cws->iter_max & cws->mu>cws->mu_max; kk++)
		{

		// fact and solve kkt
		d_fact_solve_kkt_step_hard_ocp_qp(qp, ws);

		// alpha
		d_compute_alpha_hard_qp(cws);
		cws->stat[5*kk+1] = cws->alpha;

		//
		d_update_var_hard_qp(cws);

		// compute residuals
		d_compute_res_hard_ocp_qp(qp, ws);
		cws->mu = ws->res_mu;
		cws->stat[5*kk+2] = ws->res_mu;

//		break;
		}
	
	ws->iter = kk;
	
	return;

	}


