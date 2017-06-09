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
#include "../include/hpipm_d_dense_qp.h"
#include "../include/hpipm_d_cond.h"
#include "../include/hpipm_d_aux_cond.h"



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

	int N = ocp_qp->N;
	int *nx = ocp_qp->nx;
	int *nu = ocp_qp->nu;

	int ii;

	int size = 0;

	size += (0+1*(N+1))*sizeof(struct d_strmat); // Gamma

	int nu_tmp = 0;
	for(ii=0; ii<N; ii++)
		{
		nu_tmp += nu[ii];
		size += d_size_strmat(nu_tmp+nx[0]+1, nx[ii+1]); // Gamma
		}

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

	}



void d_create_cond_qp_ocp2dense(struct d_ocp_qp *ocp_qp, struct d_dense_qp *dense_qp, struct d_cond_qp_ocp2dense_workspace *cond_ws, void *mem)
	{

	int N = ocp_qp->N;
	int *nx = ocp_qp->nx;
	int *nu = ocp_qp->nu;

	int ii;


	// matrix struct
	struct d_strmat *sm_ptr = (struct d_strmat *) mem;

	cond_ws->Gamma = sm_ptr;
	sm_ptr += N+1;


	// vector struct
	struct d_strvec *sv_ptr = (struct d_strvec *) sm_ptr;


	// align to typicl cache line size
	size_t s_ptr = (size_t) sv_ptr;
	s_ptr = (s_ptr+63)/64*64;


	// void stuf
	void *v_ptr = (void *) s_ptr;

	int nu_tmp = 0;
	for(ii=0; ii<N; ii++)
		{
		nu_tmp += nu[ii];
		d_create_strmat(nu_tmp+nx[0]+1, nx[ii+1], cond_ws->Gamma+ii, v_ptr);
		v_ptr += (cond_ws->Gamma+ii)->memory_size;
		}


	return;

	}

	

void d_cond_qp_ocp2dense(struct d_ocp_qp *ocp_qp, struct d_dense_qp *dense_qp, struct d_cond_qp_ocp2dense_workspace *cond_ws)
	{

	d_compute_Gamma(ocp_qp, cond_ws);

	return;

	}
