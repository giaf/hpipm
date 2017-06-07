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

	// compute core qp size
	int nvt = 0;
	int net = 0;
	int nbt = 0;
	int ngt = 0;

	for(ii=0; ii<N; ii++)
		{
		nvt += nx[ii]+nu[ii];
		net += nx[ii+1];
		nbt += nb[ii];
		ngt += ng[ii];
		}
	ii = N;
	nvt += nx[ii]+nu[ii];
	nbt += nb[ii];
	ngt += ng[ii];

	int size = 0;

	size += (N+1)*10*sizeof(struct d_strvec); // ux pi lam_lb lam_ub lam_lg lam_ug t_lb t_ub t_lg t_ug

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

	// compute core qp size
	int nvt = 0;
	int net = 0;
	int nbt = 0;
	int ngt = 0;

	for(ii=0; ii<N; ii++)
		{
		nvt += nx[ii]+nu[ii];
		net += nx[ii+1];
		nbt += nb[ii];
		ngt += ng[ii];
		}
	ii = N;
	nvt += nx[ii]+nu[ii];
	nbt += nb[ii];
	ngt += ng[ii];


	// core struct
	struct d_ipm_hard_core_qp_workspace *sr_ptr = mem;

	// core workspace
	workspace->core_workspace = sr_ptr;
	sr_ptr += 1;
	struct d_ipm_hard_core_qp_workspace *rwork = workspace->core_workspace;


	// matrix struct
	struct d_strmat *sm_ptr = (struct d_strmat *) sr_ptr;


	// vector struct
	struct d_strvec *sv_ptr = (struct d_strvec *) sm_ptr;

	workspace->ux = sv_ptr;
	sv_ptr += N+1;
	workspace->pi = sv_ptr;
	sv_ptr += N+1;
	workspace->lam_lb = sv_ptr;
	sv_ptr += N+1;
	workspace->lam_ub = sv_ptr;
	sv_ptr += N+1;
	workspace->lam_lg = sv_ptr;
	sv_ptr += N+1;
	workspace->lam_ug = sv_ptr;
	sv_ptr += N+1;
	workspace->t_lb = sv_ptr;
	sv_ptr += N+1;
	workspace->t_ub = sv_ptr;
	sv_ptr += N+1;
	workspace->t_lg = sv_ptr;
	sv_ptr += N+1;
	workspace->t_ug = sv_ptr;
	sv_ptr += N+1;


	// align to typicl cache line size
	size_t s_ptr = (size_t) sv_ptr;
	s_ptr = (s_ptr+63)/64*64;


	// void stuf
	void *v_ptr = (void *) s_ptr;





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




	workspace->mu0 = arg->mu0;
	workspace->nt_inv = 1.0/(2*nbt+2*ngt);

	return;

	}


void d_solve_ipm_hard_ocp_qp(struct d_ocp_qp *qp, struct d_ipm_hard_ocp_qp_workspace *ws)
	{

	struct d_ipm_hard_core_qp_workspace *rws = ws->core_workspace;

	// alias qp vectors into core workspace

	// init solver
	d_init_var_hard_ocp_qp(qp, ws);

	return;

	}


