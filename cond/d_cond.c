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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_d_blas.h>
#include <blasfeo_d_aux.h>

#include "../include/hpipm_d_ocp_qp.h"
#include "../include/hpipm_d_ocp_qp_sol.h"
#include "../include/hpipm_d_dense_qp.h"
#include "../include/hpipm_d_dense_qp_sol.h"
#include "../include/hpipm_d_cond.h"
#include "../include/hpipm_d_cond_aux.h"



void d_compute_qp_size_ocp2dense(int N, int *nx, int *nu, int *nb, int **idxb, int *ng, int *nvc, int *nec, int *nbc, int *ngc)
	{

	int ii, jj;

	nvc[0] = 0;
	nec[0] = 0;
	nbc[0] = 0;
	ngc[0] = 0;

	// first stage
	nvc[0] += nx[0]+nu[0];
	nbc[0] += nb[0];
	ngc[0] += ng[0];
	// remaining stages
	for(ii=1; ii<=N; ii++)
		{
		nvc[0] += nu[ii];
		for(jj=0; jj<nb[ii]; jj++)
			{
			if(idxb[ii][jj]<nu[ii]) // input constraint
				{
				nbc[0]++;
				}
			else // state constraint
				{
				ngc[0]++;
				}
			}
		ngc[0] += ng[ii];
		}

	return;

	}



int d_memsize_cond_qp_ocp2dense(struct d_ocp_qp *ocp_qp, struct d_dense_qp *dense_qp) // XXX + args for algorithm type ???
	{

	int ii;

	int N = ocp_qp->N;
	int *nx = ocp_qp->nx;
	int *nu = ocp_qp->nu;
	int *nb = ocp_qp->nb;
	int *ng = ocp_qp->ng;

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

	size += (2+2*(N+1))*sizeof(struct d_strmat); // Gamma L Lx AL
	size += (2+1*(N+1))*sizeof(struct d_strvec); // Gammab tmp_ngM tmp_nuxM

	int nu_tmp = 0;
	for(ii=0; ii<N; ii++)
		{
		nu_tmp += nu[ii];
		size += d_size_strmat(nu_tmp+nx[0]+1, nx[ii+1]); // Gamma
		}
	for(ii=0; ii<=N; ii++) 
		size += d_size_strmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii]); // L
	size += d_size_strmat(nxM+1, nxM); // Lx
	size += d_size_strmat(nuM+nxM+1, nxM); // AL
	for(ii=0; ii<N; ii++) 
		size += 1*d_size_strvec(nx[ii+1]); // Gammab
	size += d_size_strvec(ngM); // tmp_ngM
	size += 1*d_size_strvec(nuM+nxM); // tmp_nuxM
	size += 1*d_size_strvec(ngM); // tmp_ngM

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

	}



void d_create_cond_qp_ocp2dense(struct d_ocp_qp *ocp_qp, struct d_dense_qp *dense_qp, struct d_cond_qp_ocp2dense_workspace *cond_ws, void *mem)
	{

	int ii;

	int N = ocp_qp->N;
	int *nx = ocp_qp->nx;
	int *nu = ocp_qp->nu;
	int *nb = ocp_qp->nb;
	int *ng = ocp_qp->ng;

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


	// matrix struct
	struct d_strmat *sm_ptr = (struct d_strmat *) mem;

	cond_ws->Gamma = sm_ptr;
	sm_ptr += N+1;
	cond_ws->L = sm_ptr;
	sm_ptr += N+1;
	cond_ws->Lx = sm_ptr;
	sm_ptr += 1;
	cond_ws->AL = sm_ptr;
	sm_ptr += 1;


	// vector struct
	struct d_strvec *sv_ptr = (struct d_strvec *) sm_ptr;

	cond_ws->Gammab = sv_ptr;
	sv_ptr += N+1;
	cond_ws->tmp_ngM = sv_ptr;
	sv_ptr += 1;
	cond_ws->tmp_nuxM = sv_ptr;
	sv_ptr += 1;


	// align to typicl cache line size
	size_t s_ptr = (size_t) sv_ptr;
	s_ptr = (s_ptr+63)/64*64;


	// void stuf
	char *c_ptr = (char *) s_ptr;
	char *c_tmp;

	int nu_tmp = 0;
	for(ii=0; ii<N; ii++)
		{
		nu_tmp += nu[ii];
		d_create_strmat(nu_tmp+nx[0]+1, nx[ii+1], cond_ws->Gamma+ii, c_ptr);
		c_ptr += (cond_ws->Gamma+ii)->memory_size;
		}
	for(ii=0; ii<=N; ii++)
		{
		d_create_strmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], cond_ws->L+ii, c_ptr);
		c_ptr += (cond_ws->L+ii)->memory_size;
		}
	d_create_strmat(nxM+1, nxM, cond_ws->Lx, c_ptr);
	c_ptr += cond_ws->Lx->memory_size;
	d_create_strmat(nuM+nxM+1, nxM, cond_ws->AL, c_ptr);
	c_ptr += cond_ws->AL->memory_size;
	for(ii=0; ii<N; ii++)
		{
		d_create_strvec(nx[ii+1], cond_ws->Gammab+ii, c_ptr);
		c_ptr += (cond_ws->Gammab+ii)->memory_size;
		}
	d_create_strvec(ngM, cond_ws->tmp_ngM, c_ptr);
	c_ptr += cond_ws->tmp_ngM->memory_size;
	c_tmp = c_ptr;
	d_create_strvec(nuM+nxM, cond_ws->tmp_nuxM, c_ptr);
	c_ptr += cond_ws->tmp_nuxM->memory_size;

	cond_ws->memsize = d_memsize_cond_qp_ocp2dense(ocp_qp, dense_qp);

	return;

	}

	

void d_cond_qp_ocp2dense(struct d_ocp_qp *ocp_qp, struct d_dense_qp *dense_qp, struct d_cond_qp_ocp2dense_workspace *cond_ws)
	{

	d_compute_Gamma(ocp_qp, cond_ws);

	d_cond_RSQrq_N2nx3(ocp_qp, dense_qp->Hg, dense_qp->g, cond_ws);

	d_cond_DCtd(ocp_qp, dense_qp->idxb, dense_qp->d_lb, dense_qp->d_ub, dense_qp->Ct, dense_qp->d_lg, dense_qp->d_ug, cond_ws);

	return;

	}



void d_expand_sol_dense2ocp(struct d_ocp_qp *ocp_qp, struct d_dense_qp_sol *dense_qp_sol, struct d_ocp_qp_sol *ocp_qp_sol, struct d_cond_qp_ocp2dense_workspace *cond_ws)
	{

	d_expand_sol(ocp_qp, dense_qp_sol, ocp_qp_sol->ux, ocp_qp_sol->pi, ocp_qp_sol->lam_lb, ocp_qp_sol->lam_ub, ocp_qp_sol->lam_lg, ocp_qp_sol->lam_ug, ocp_qp_sol->t_lb, ocp_qp_sol->t_ub, ocp_qp_sol->t_lg, ocp_qp_sol->t_ug, cond_ws->tmp_nuxM, cond_ws->tmp_ngM);

	return;

	}
