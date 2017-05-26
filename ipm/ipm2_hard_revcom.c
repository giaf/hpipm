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



#include "../include/hpipm_ipm2_hard_revcom.h"



void d_ipm2_hard_revcom(struct d_ipm2_hard_revcom_workspace *workspace)
	{

	int status, entry;
	int nv, ne, nb, ng, nv0, ne0, nb0, ng0;
	int it0;
	double *dptr;
	int *iptr;

	status = workspace->status;

	switch(status)
		{
		case IPMCORE_MEMSIZE:	goto memsize;
		case IPMCORE_MEMPART:	goto mempart;
		case IPMCORE_INIT	:	goto init;
		case IPMCORE_WAITING:	goto waiting;
		}
	workspace->status = IPMCORE_ERROR;
	return;



	/* compute memory size of work space */
memsize:
	
	nb = workspace->nb;
	nv0 = workspace->nv;
	ne0 = workspace->ne;
	nb0 = workspace->nb;
	ng0 = workspace->ng;
// if target avx
// nv0 = ...

	it0 = 0;

	it0 += 3*nv0*sizeof(double); // v dv res_q
	it0 += 3*ne0*sizeof(double); // pi dpi res_b
	it0 += 6*(2*nb0+2*ng0)*sizeof(double); // lam t dlam dt res_d res_m
	it0 += 2*nb0*sizeof(double); // Qx qx
	it0 += ng0*sizeof(double); // Dv
	it0 += 5*workspace->iter_max*sizeof(double); // conv_stat
	it0 += nb*sizeof(int); // idxb

// make multiple of cache line size?

	workspace->memsize = it0;

	return;



	/* partition workspace memory */
mempart:
	
	nb = workspace->nb;
	nv0 = workspace->nv;
	ne0 = workspace->ne;
	nb0 = workspace->nb;
	ng0 = workspace->ng;
// if target avx
// nv0 = ...

	dptr = (double *) workspace->mem;
	workspace->v = dptr; // v
	dptr += nv0;
	workspace->pi = dptr; // pi
	dptr += ne0;
	workspace->lam = dptr; // lam
	dptr += 2*nb0+2*ng0;
	workspace->t = dptr; // t
	dptr += 2*nb0+2*ng0;
	workspace->dv = dptr; // dv
	dptr += nv0;
	workspace->dpi = dptr; // dpi
	dptr += ne0;
	workspace->dlam = dptr; // dlam
	dptr += 2*nb0+2*ng0;
	workspace->dt = dptr; // dt
	dptr += 2*nb0+2*ng0;
	workspace->res_q = dptr; // res_q
	dptr += nv0;
	workspace->res_b = dptr; // res_b
	dptr += ne0;
	workspace->res_d = dptr; // res_d
	dptr += 2*nb0+2*ng0;
	workspace->res_m = dptr; // res_m
	dptr += 2*nb0+2*ng0;
	workspace->Qx = dptr; // Qx
	dptr += nb0;
	workspace->qx = dptr; // qx
	dptr += nb0;
	workspace->Dv = dptr; // Dv
	dptr += ng0;
	workspace->conv_stat = dptr; // conv_stat
	dptr += 5*workspace->iter_max;
	iptr = (int *) dptr;
	workspace->idxb = iptr; // idxb
	iptr += nb;

	return;



init:
	
	workspace->sigma = 0.0;

	return;



waiting:
	
	entry = workspace->entry;

	switch(entry)
		{
		case ITER_RES:	goto iter_res;
		case ALPHA_RES:	goto alpha_res;
		case EXIT_RES:	goto exit_res;
		}

	return;



	/* new iteration (ipm res) */
iter_res:
	
	d_update_hessian_gradient_res_hard(workspace);

	workspace->status = IPMCORE_FACT_SOLVE_KKT_COMP_DV;
	workspace->entry = ALPHA_RES;

	return;



	/* compute step lenght alpha (ipm res)*/
alpha_res:
	
	workspace->alpha = 1.0;

	d_compute_alpha_res_hard(workspace);

	workspace->conv_stat[5*(workspace->iter)+0] = workspace->sigma;
	workspace->conv_stat[5*(workspace->iter)+3] = workspace->alpha;
	workspace->alpha *= 0.995;

	d_update_var_res_hard(workspace);

	workspace->status = IPMCORE_COMP_RES;
	workspace->entry = EXIT_RES;
	
	return;



exit_res:
	
	workspace->iter++;

	if(workspace->iter < workspace->iter_max & \
	   workspace->mu > workspace->mu_max & \
	   workspace->alpha >= workspace->alpha_min)
		goto iter_res;
	
	workspace->status = IPMCORE_STOP;
	
	return;

	}
