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



int MEMSIZE_OCP_QP_DIM(int N)
	{

	int size = 0;

	size += 10*(N+1)*sizeof(int);

	size = (size+8-1)/8*8;

	return size;

	}



void CREATE_OCP_QP_DIM(int N, struct OCP_QP_DIM *dim, void *memory)
	{

	// loop index
	int ii;

	char *c_ptr = memory;

	// nx
	dim->nx = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	// nu
	dim->nu = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	// nb
	dim->nb = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	// nbx
	dim->nbx = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	// nbu
	dim->nbu = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	// ng
	dim->ng = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	// ns
	dim->ns = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	// nsbx
	dim->nsbx = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	// nsbu
	dim->nsbu = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	// nsg
	dim->nsg = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);

	// N
	dim->N = N;

	dim->memsize = MEMSIZE_OCP_QP_DIM(N);

	return;

	}


void CVT_INT_TO_OCP_QP_DIM(int N, int *nx, int *nu, int *nbx, int *nbu, int *ng, int *nsbx, int *nsbu, int *nsg, struct OCP_QP_DIM *dim)
	{

	// loop index
	int ii;

	// N
	dim->N = N;

	// copy qp dim
	for(ii=0; ii<=N; ii++)
		dim->nx[ii] = nx[ii];
	for(ii=0; ii<=N; ii++)
		dim->nu[ii] = nu[ii];
	for(ii=0; ii<=N; ii++)
		dim->nb[ii] = nbx[ii]+nbu[ii];
	for(ii=0; ii<=N; ii++)
		dim->nbx[ii] = nbx[ii];
	for(ii=0; ii<=N; ii++)
		dim->nbu[ii] = nbu[ii];
	for(ii=0; ii<=N; ii++)
		dim->ng[ii] = ng[ii];
	for(ii=0; ii<=N; ii++)
		dim->ns[ii] = nsbx[ii]+nsbu[ii]+nsg[ii];
	for(ii=0; ii<=N; ii++)
		dim->nsbx[ii] = nsbx[ii];
	for(ii=0; ii<=N; ii++)
		dim->nsbu[ii] = nsbu[ii];
	for(ii=0; ii<=N; ii++)
		dim->nsg[ii] = nsg[ii];

	return;

	}
