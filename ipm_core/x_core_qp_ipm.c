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



int MEMSIZE_CORE_QP_IPM(int nv, int ne, int nc)
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
	size += 8*nc0*sizeof(REAL); // dlam dt res_d res_m res_m_bkp t_inv Gamma gamma

	size = (size+63)/64*64; // make multiple of cache line size

	return size;

	}



void CREATE_CORE_QP_IPM(int nv, int ne, int nc, struct CORE_QP_IPM_WORKSPACE *workspace, void *mem)
	{

	workspace->nv = nv;
	workspace->ne = ne;
	workspace->nc = nc;

	int nv0 = nv;
	int ne0 = ne;
	int nc0 = nc;
// if target avx NO!!!!
// nv0 = ...

	REAL *d_ptr = (REAL *) mem;

	workspace->t_inv = d_ptr; // t_inv
	d_ptr += nc0;

	workspace->dv = d_ptr; // dv
	d_ptr += nv0;

	workspace->dpi = d_ptr; // dpi
	d_ptr += ne0;

	workspace->dlam = d_ptr; // dlam
	d_ptr += nc0;

	workspace->dt = d_ptr; // dt
	d_ptr += nc0;

	workspace->res_g = d_ptr; // res_g
	d_ptr += nv0;

	workspace->res_b = d_ptr; // res_b
	d_ptr += ne0;

	workspace->res_d = d_ptr; // res_d
	d_ptr += nc0;

	workspace->res_m = d_ptr; // res_m
	d_ptr += nc0;

	workspace->res_m_bkp = d_ptr; // res_m_bkp
	d_ptr += nc0;

	workspace->Gamma = d_ptr; // Gamma
	d_ptr += nc0;

	workspace->gamma = d_ptr; // gamma
	d_ptr += nc0;

	if(nc>0)
		workspace->nc_inv = 1.0/nc;
	else
		workspace->nc_inv = 0.0;


	workspace->lam_min = 0.0;
	workspace->t_min = 0.0;


	workspace->memsize = MEMSIZE_CORE_QP_IPM(nv, ne, nc);


	char * c_ptr = (char *) d_ptr;


	#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + workspace->memsize)
		{
		printf("\nCreate_core_qp_ipm: outsize memory bounds!\n\n");
		exit(1);
		}
#endif


return;

	}


