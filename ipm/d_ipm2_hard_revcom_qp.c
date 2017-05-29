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



#include "../include/hpipm_d_ipm2_hard_revcom_qp.h"
#include "../include/hpipm_d_aux_ipm_hard.h"



int d_memsize_ipm2_hard_revcom_qp(int nv, int ne, int nb, int ng, int iter_max)
	{

	int size;

	int nv0 = nv;
	int ne0 = ne;
	int nb0 = nb;
	int ng0 = ng;
// if target avx
// nv0 = ...

	size = 0;

	size += 3*nv0*sizeof(double); // v dv res_g
	size += 3*ne0*sizeof(double); // pi dpi res_b
	size += 7*(2*nb0+2*ng0)*sizeof(double); // lam t dlam dt res_d res_m t_inv
	size += 2*nb0*sizeof(double); // Qx qx
	size += ng0*sizeof(double); // Dv
	size += 5*iter_max*sizeof(double); // conv_stat
	size += nb*sizeof(int); // idxb

	size = (size+63)/64*64; // make multiple of cache line size

	return size;

	}



void d_create_ipm2_hard_revcom_qp(struct d_ipm2_hard_revcom_qp_workspace *workspace, void *mem)
	{

	int nv = workspace->nv;
	int ne = workspace->ne;
	int nb = workspace->nb;
	int ng = workspace->ng;

	int nv0 = nv;
	int ne0 = ne;
	int nb0 = nb;
	int ng0 = ng;
// if target avx NO!!!!
// nv0 = ...

	workspace->memsize = d_memsize_ipm2_hard_revcom_qp(nv, ne, nb, ng, workspace->iter_max);

	double *d_ptr = (double *) mem;

	workspace->v = d_ptr; // v
	d_ptr += nv0;

	workspace->pi = d_ptr; // pi
	d_ptr += ne0;

	workspace->lam = d_ptr; // lam
	workspace->lam_lb = d_ptr; // lam_lb
	workspace->lam_lg = d_ptr+nb0; // lam_lg
	workspace->lam_ub = d_ptr+nb0+ng0; // lam_ub
	workspace->lam_ug = d_ptr+2*nb0+ng0; // lam_ug
	d_ptr += 2*nb0+2*ng0;

	workspace->t = d_ptr; // t
	workspace->t_lb = d_ptr; // t_lb
	workspace->t_lg = d_ptr+nb0; // t_lg
	workspace->t_ub = d_ptr+nb0+ng0; // t_ub
	workspace->t_ug = d_ptr+2*nb0+ng0; // t_ug
	d_ptr += 2*nb0+2*ng0;

	workspace->dv = d_ptr; // dv
	d_ptr += nv0;

	workspace->dpi = d_ptr; // dpi
	d_ptr += ne0;

	workspace->dlam = d_ptr; // dlam
	d_ptr += 2*nb0+2*ng0;

	workspace->dt = d_ptr; // dt
	d_ptr += 2*nb0+2*ng0;

	workspace->res_g = d_ptr; // res_g
	d_ptr += nv0;

	workspace->res_b = d_ptr; // res_b
	d_ptr += ne0;

	workspace->res_d = d_ptr; // res_d
	workspace->res_d_lb = d_ptr;
	workspace->res_d_lg = d_ptr+nb0;
	workspace->res_d_ub = d_ptr+nb0+ng0;
	workspace->res_d_ug = d_ptr+2*nb0+ng0;
	d_ptr += 2*nb0+2*ng0;

	workspace->res_m = d_ptr; // res_m
	d_ptr += 2*nb0+2*ng0;

	workspace->t_inv = d_ptr; // t_inv
	d_ptr += 2*nb0+2*ng0;

	workspace->Qx = d_ptr; // Qx
	d_ptr += nb0;

	workspace->qx = d_ptr; // qx
	d_ptr += nb0;

	workspace->Dv = d_ptr; // Dv
	d_ptr += ng0;

	workspace->conv_stat = d_ptr; // conv_stat
	d_ptr += 5*workspace->iter_max;

	int *i_ptr = (int *) d_ptr;

	workspace->idxb = i_ptr; // idxb
	i_ptr += nb;

	return;

	}



void d_ipm2_hard_revcom_qp(struct d_ipm2_hard_revcom_qp_workspace *workspace)
	{

	int status, entry;
	int nv, ne, nb, ng, nv0, ne0, nb0, ng0;
	int it0;
	double *dptr;
	int *iptr;



	workspace->status = IPMCORE_WAITING;



	entry = workspace->entry;

	switch(entry)
		{
		case INIT_RES:	goto init_res;
		case ITER_RES:	goto iter_res;
		case ALPHA_RES:	goto alpha_res;
		case EXIT_RES:	goto exit_res;
		}
	workspace->status = IPMCORE_ERROR;

	return;



	/* initialize */
init_res:
	
	d_init_var_hard(workspace);

//	workspace->sigma = 10.0; // XXX

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
