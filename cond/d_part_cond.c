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
#include "../include/hpipm_d_part_cond.h"
#include "../include/hpipm_d_cond_aux.h"


void d_compute_qp_size_ocp2ocp(int N, int *nx, int *nu, int *nb, int **idxb, int *ng, int N2, int *nx2, int *nu2, int *nb2, int *ng2)
	{

	int ii, jj, kk;

	int N1 = N/N2; // (floor) horizon of small blocks
	int R1 = N - N2*N1; // the first R1 blocks have horizon N1+1
	int M1 = R1>0 ? N1+1 : N1; // (ceil) horizon of large blocks
	int T1; // horizon of current block

	int N_tmp = 0; // temporary sum of horizons
	int nbb; // box constr that remain box constr
	int nbg; // box constr that becomes general constr
	for(ii=0; ii<N2; ii++)
		{
		T1 = ii<R1 ? M1 : N1;
		nx2[ii] = nx[N_tmp+0];
		nu2[ii] = nu[N_tmp+0];
		nb2[ii] = nb[N_tmp+0];
		ng2[ii] = ng[N_tmp+0];
		for(jj=1; jj<T1; jj++)
			{
			nbb = 0;
			nbg = 0;
			for(kk=0; kk<nb[N_tmp+jj]; kk++)
				if(idxb[N_tmp+jj][kk]<nu[N_tmp+jj])
					nbb++;
				else
					nbg++;
			nx2[ii] += 0;
			nu2[ii] += nu[N_tmp+jj];
			nb2[ii] += nbb;
			ng2[ii] += ng[N_tmp+jj] + nbg;
			}
		N_tmp += T1;
		}
	nx2[N2] = nx[N];
	nu2[N2] = nu[N];
	nb2[N2] = nb[N];
	ng2[N2] = ng[N];

	return;

	}



int d_memsize_cond_qp_ocp2ocp(struct d_ocp_qp *ocp_qp, struct d_ocp_qp *part_dense_qp)
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
		size += d_size_strmat(nu_tmp+nx[0]+1, nx[ii+1]); // Gamma TODO
		}
	for(ii=0; ii<=N; ii++) 
		size += d_size_strmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii]); // L TODO
	size += d_size_strmat(nxM+1, nxM); // Lx
	size += d_size_strmat(nuM+nxM+1, nxM); // AL
	for(ii=0; ii<N; ii++) 
		size += 1*d_size_strvec(nx[ii+1]); // Gammab TODO
	size += d_size_strvec(ngM); // tmp_ngM
	size += 1*d_size_strvec(nuM+nxM); // tmp_nuxM
	size += 1*d_size_strvec(ngM); // tmp_ngM

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

	}



void d_create_cond_qp_ocp2ocp(struct d_ocp_qp *ocp_qp, struct d_ocp_qp *part_dense_qp, struct d_cond_qp_ocp2ocp_workspace *cond_ws, void *mem)
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

	cond_ws->memsize = d_memsize_cond_qp_ocp2ocp(ocp_qp, part_dense_qp);

	return;

	}

	

void d_cond_qp_ocp2ocp(struct d_ocp_qp *ocp_qp, struct d_ocp_qp *part_dense_qp, struct d_cond_qp_ocp2ocp_workspace *part_cond_ws)
	{

	struct d_cond_qp_ocp2dense_workspace cond_ws;
	struct d_ocp_qp tmp_ocp_qp;

	int ii;

	int N = ocp_qp->N;
	int N2 = part_dense_qp->N;
	int N1 = N/N2; // (floor) horizon of small blocks
	int R1 = N - N2*N1; // the first R1 blocks have horizon N1+1
	int M1 = R1>0 ? N1+1 : N1; // (ceil) horizon of large blocks
	int T1; // horizon of current block

	int N_tmp = 0; // temporary sum of horizons
	for(ii=0; ii<N2; ii++)
		{

		T1 = ii<R1 ? M1 : N1;

		// alias ocp_qp
		tmp_ocp_qp.N = T1;
		tmp_ocp_qp.nx = ocp_qp->nx+N_tmp;
		tmp_ocp_qp.nu = ocp_qp->nu+N_tmp;
		tmp_ocp_qp.nb = ocp_qp->nb+N_tmp;
		tmp_ocp_qp.ng = ocp_qp->ng+N_tmp;
		tmp_ocp_qp.idxb = ocp_qp->idxb+N_tmp;
		tmp_ocp_qp.BAbt = ocp_qp->BAbt+N_tmp;
		tmp_ocp_qp.b = ocp_qp->b+N_tmp;
		tmp_ocp_qp.RSQrq = ocp_qp->RSQrq+N_tmp;
		tmp_ocp_qp.rq = ocp_qp->rq+N_tmp;
		tmp_ocp_qp.DCt = ocp_qp->DCt+N_tmp;
		tmp_ocp_qp.d_lb = ocp_qp->d_lb+N_tmp;
		tmp_ocp_qp.d_ub = ocp_qp->d_ub+N_tmp;
		tmp_ocp_qp.d_lg = ocp_qp->d_lg+N_tmp;
		tmp_ocp_qp.d_ug = ocp_qp->d_ug+N_tmp;

		// alias part_cond_ws
		cond_ws.Gamma = part_cond_ws->Gamma+N_tmp;
		cond_ws.L = part_cond_ws->L+N_tmp;
		cond_ws.Gammab = part_cond_ws->Gammab+N_tmp;
		cond_ws.Lx = part_cond_ws->Lx;
		cond_ws.AL = part_cond_ws->AL;
		cond_ws.tmp_ngM = part_cond_ws->tmp_ngM;
		cond_ws.tmp_nuxM = part_cond_ws->tmp_nuxM;

		d_compute_Gamma(&tmp_ocp_qp, &cond_ws);

//		d_cond_RSQrq_N2nx3(&tmp_ocp_qp, part_dense_qp->RSQrq+ii, part_dense_qp->rq+ii, &cond_ws);

//		d_cond_DCtd(&tmp_ocp_qp, part_dense_qp->idxb, dense_qp->d_lb, dense_qp->d_ub, dense_qp->Ct, dense_qp->d_lg, dense_qp->d_ug, cond_ws);

		N_tmp += T1;

		}

	return;

	}




