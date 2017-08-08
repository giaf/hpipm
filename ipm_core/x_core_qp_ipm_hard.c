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



int MEMSIZE_IPM_HARD_CORE_QP(int nv, int ne, int nb, int ng, int iter_max)
	{

	int size;

	int nv0 = nv;
	int ne0 = ne;
	int nb0 = nb;
	int ng0 = ng;
// if target avx
// nv0 = ...

	size = 0;

	size += 2*nv0*sizeof(REAL); // dv res_g
	size += 2*ne0*sizeof(REAL); // dpi res_b
	size += 5*(2*nb0+2*ng0)*sizeof(REAL); // dlam dt res_d res_m t_inv
	size += 2*(nb0+ng0)*sizeof(REAL); // Qx qx
	size += 5*iter_max*sizeof(REAL); // stat

	size = (size+63)/64*64; // make multiple of cache line size

	return size;

	}



void CREATE_IPM_HARD_CORE_QP(struct IPM_HARD_CORE_QP_WORKSPACE *workspace, void *mem)
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

	workspace->memsize = MEMSIZE_IPM_HARD_CORE_QP(nv, ne, nb, ng, workspace->iter_max);

	REAL *d_ptr = (REAL *) mem;

	workspace->t_inv = d_ptr; // t_inv
	workspace->t_inv_lb = d_ptr;
	workspace->t_inv_lg = d_ptr+nb0;
	workspace->t_inv_ub = d_ptr+nb0+ng0;
	workspace->t_inv_ug = d_ptr+2*nb0+ng0;
	d_ptr += 2*nb0+2*ng0;

	workspace->dv = d_ptr; // dv
	d_ptr += nv0;

	workspace->dpi = d_ptr; // dpi
	d_ptr += ne0;

	workspace->dlam = d_ptr; // dlam
	workspace->dlam_lb = d_ptr;
	workspace->dlam_lg = d_ptr+nb0;
	workspace->dlam_ub = d_ptr+nb0+ng0;
	workspace->dlam_ug = d_ptr+2*nb0+ng0;
	d_ptr += 2*nb0+2*ng0;

	workspace->dt = d_ptr; // dt
	workspace->dt_lb = d_ptr;
	workspace->dt_lg = d_ptr+nb0;
	workspace->dt_ub = d_ptr+nb0+ng0;
	workspace->dt_ug = d_ptr+2*nb0+ng0;
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
	workspace->res_m_lb = d_ptr;
	workspace->res_m_lg = d_ptr+nb0;
	workspace->res_m_ub = d_ptr+nb0+ng0;
	workspace->res_m_ug = d_ptr+2*nb0+ng0;
	d_ptr += 2*nb0+2*ng0;

	workspace->Qx = d_ptr; // Qx
	workspace->Qx_lb = d_ptr;
	workspace->Qx_lg = d_ptr+nb0;
	d_ptr += nb0+ng0;

	workspace->qx = d_ptr; // qx
	workspace->qx_lb = d_ptr;
	workspace->qx_lg = d_ptr+nb0;
	d_ptr += nb0+ng0;

	workspace->stat = d_ptr; // stat
	d_ptr += 5*workspace->iter_max;

	int *i_ptr = (int *) d_ptr;

	return;

	}


