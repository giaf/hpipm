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



int MEMSIZE_CORE_QP_IPM(int nv, int ne, int nc, int stat_max)
	{

	int size;

	int nv0 = nv;
	int ne0 = ne;
	int nc0 = nc;
// if target avx
// nv0 = ...

	size = 0;

	size += 2*nv0*sizeof(REAL); // dv res_g
	size += 2*ne0*sizeof(REAL); // dpi res_b
	size += 7*(2*nc0)*sizeof(REAL); // dlam dt res_d res_m t_inv Gamma gamma
	size += 2*(nc0)*sizeof(REAL); // Qx qx
	size += 5*stat_max*sizeof(REAL); // stat

	size = (size+63)/64*64; // make multiple of cache line size

	return size;

	}



void CREATE_CORE_QP_IPM(struct CORE_QP_IPM_WORKSPACE *workspace, void *mem)
	{

	int nv = workspace->nv;
	int ne = workspace->ne;
	int nc = workspace->nc;

	int nv0 = nv;
	int ne0 = ne;
	int nc0 = nc;
// if target avx NO!!!!
// nv0 = ...

	workspace->memsize = MEMSIZE_CORE_QP_IPM(nv, ne, nc, workspace->stat_max);

	REAL *d_ptr = (REAL *) mem;

	workspace->t_inv = d_ptr; // t_inv
	d_ptr += 2*nc0;

	workspace->dv = d_ptr; // dv
	d_ptr += nv0;

	workspace->dpi = d_ptr; // dpi
	d_ptr += ne0;

	workspace->dlam = d_ptr; // dlam
	d_ptr += 2*nc0;

	workspace->dt = d_ptr; // dt
	d_ptr += 2*nc0;

	workspace->res_g = d_ptr; // res_g
	d_ptr += nv0;

	workspace->res_b = d_ptr; // res_b
	d_ptr += ne0;

	workspace->res_d = d_ptr; // res_d
	d_ptr += 2*nc0;

	workspace->res_m = d_ptr; // res_m
	d_ptr += 2*nc0;

	workspace->Gamma = d_ptr; // Gamma
	d_ptr += 2*nc0;

	workspace->gamma = d_ptr; // gamma
	d_ptr += 2*nc0;

	workspace->Qx = d_ptr; // Qx
	d_ptr += nc0;

	workspace->qx = d_ptr; // qx
	d_ptr += nc0;

	workspace->stat = d_ptr; // stat
	d_ptr += 5*workspace->stat_max;

	int *i_ptr = (int *) d_ptr;

	return;

	}


