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



int MEMSIZE_OCP_QP_SIZE(int N)
	{

	int size = 0;

	size += 7*(N+1)*sizeof(int);

	size = (size+8-1)/8*8;

	return size;

	}



void CREATE_OCP_QP_SIZE(int N, int *nx, int *nu, int *nbx, int *nbu, int *ng, int *ns, struct OCP_QP_SIZE *size, void *memory)
	{

	// loop index
	int ii;

	char *c_ptr = memory;

	// N
	size->N = N;

	// nx
	size->nx = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	// nu
	size->nu = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	// nb
	size->nb = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	// nbx
	size->nbx = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	// nbu
	size->nbu = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	// ng
	size->ng = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	// ns
	size->ns = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);

	// copy qp size
	for(ii=0; ii<=N; ii++)
		size->nx[ii] = nx[ii];
	for(ii=0; ii<=N; ii++)
		size->nu[ii] = nu[ii];
	for(ii=0; ii<=N; ii++)
		size->nb[ii] = nbx[ii]+nbu[ii];
	for(ii=0; ii<=N; ii++)
		size->nbx[ii] = nbx[ii];
	for(ii=0; ii<=N; ii++)
		size->nbu[ii] = nbu[ii];
	for(ii=0; ii<=N; ii++)
		size->ng[ii] = ng[ii];
	for(ii=0; ii<=N; ii++)
		size->ns[ii] = ns[ii];

	size->memsize = MEMSIZE_OCP_QP_SIZE(N);

	return;

	}

