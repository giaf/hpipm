/**************************************************************************************************
*                                                                                                 *
* This file is part of HPIPM.                                                                     *
*                                                                                                 *
* HPIPM -- High-Performance Interior Point Method.                                                *
* Copyright (C) 2017-2018 by Gianluca Frison.                                                     *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
*                                                                                                 *
* This program is free software: you can redistribute it and/or modify                            *
* it under the terms of the GNU General Public License as published by                            *
* the Free Software Foundation, either version 3 of the License, or                               *
* (at your option) any later version                                                              *.
*                                                                                                 *
* This program is distributed in the hope that it will be useful,                                 *
* but WITHOUT ANY WARRANTY; without even the implied warranty of                                  *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                   *
* GNU General Public License for more details.                                                    *
*                                                                                                 *
* You should have received a copy of the GNU General Public License                               *
* along with this program.  If not, see <https://www.gnu.org/licenses/>.                          *
*                                                                                                 *
* The authors designate this particular file as subject to the "Classpath" exception              *
* as provided by the authors in the LICENSE file that accompained this code.                      *
*                                                                                                 *
* Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             *
*                                                                                                 *
**************************************************************************************************/



int SIZEOF_OCP_QP_DIM()
	{
	return sizeof(struct OCP_QP_DIM);
	}



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
	for(ii=0; ii<=N; ii++)
		dim->nx[ii] = 0;
	// nu
	dim->nu = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	for(ii=0; ii<=N; ii++)
		dim->nu[ii] = 0;
	// nb
	dim->nb = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	for(ii=0; ii<=N; ii++)
		dim->nb[ii] = 0;
	// nbx
	dim->nbx = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	for(ii=0; ii<=N; ii++)
		dim->nbx[ii] = 0;
	// nbu
	dim->nbu = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	for(ii=0; ii<=N; ii++)
		dim->nbu[ii] = 0;
	// ng
	dim->ng = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	for(ii=0; ii<=N; ii++)
		dim->ng[ii] = 0;
	// ns
	dim->ns = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	for(ii=0; ii<=N; ii++)
		dim->ns[ii] = 0;
	// nsbx
	dim->nsbx = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	for(ii=0; ii<=N; ii++)
		dim->nsbx[ii] = 0;
	// nsbu
	dim->nsbu = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	for(ii=0; ii<=N; ii++)
		dim->nsbu[ii] = 0;
	// nsg
	dim->nsg = (int *) c_ptr;
	c_ptr += (N+1)*sizeof(int);
	for(ii=0; ii<=N; ii++)
		dim->ng[ii] = 0;

	// N
	dim->N = N;

	dim->memsize = MEMSIZE_OCP_QP_DIM(N);

	return;

	}


void CVT_INT_TO_OCP_QP_DIM(int N, int *nx, int *nu, int *nbx, int *nbu, int *ng, int *ns, struct OCP_QP_DIM *dim)
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
		dim->ns[ii] = ns[ii];//nsbx[ii]+nsbu[ii]+nsg[ii];
	for(ii=0; ii<=N; ii++)
		dim->nsbx[ii] = 0;//-1;//nsbx[ii];
	for(ii=0; ii<=N; ii++)
		dim->nsbu[ii] = 0;//-1;//nsbu[ii];
	for(ii=0; ii<=N; ii++)
		dim->nsg[ii] = 0;//-1;//nsg[ii];

	return;

	}



void SET_OCP_QP_DIM_NX(int stage, int nx, struct OCP_QP_DIM *dim)
	{
	dim->nx[stage] = nx;
	return;
	}



void SET_OCP_QP_DIM_NU(int stage, int nu, struct OCP_QP_DIM *dim)
	{
	dim->nu[stage] = nu;
	return;
	}



void SET_OCP_QP_DIM_NBX(int stage, int nbx, struct OCP_QP_DIM *dim)
	{
	dim->nbx[stage] = nbx;
	dim->nb[stage] = dim->nbx[stage] + dim->nbu[stage];
	return;
	}



void SET_OCP_QP_DIM_NBU(int stage, int nbu, struct OCP_QP_DIM *dim)
	{
	dim->nbu[stage] = nbu;
	dim->nb[stage] = dim->nbx[stage] + dim->nbu[stage];
	return;
	}



void SET_OCP_QP_DIM_NG(int stage, int ng, struct OCP_QP_DIM *dim)
	{
	dim->ng[stage] = ng;
	return;
	}



void SET_OCP_QP_DIM_NS(int stage, int ns, struct OCP_QP_DIM *dim)
	{
	dim->ns[stage] = ns;
	return;
	}
