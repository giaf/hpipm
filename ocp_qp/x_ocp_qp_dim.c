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



int OCP_QP_DIM_STRSIZE()
	{
	return sizeof(struct OCP_QP_DIM);
	}



int OCP_QP_DIM_MEMSIZE(int N)
	{

	int size = 0;

	size += 10*(N+1)*sizeof(int);

	size = (size+8-1)/8*8;

	return size;

	}



void OCP_QP_DIM_CREATE(int N, struct OCP_QP_DIM *dim, void *memory)
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
		dim->nsg[ii] = 0;

	// N
	dim->N = N;

	dim->memsize = OCP_QP_DIM_MEMSIZE(N);

	return;

	}


void OCP_QP_DIM_SET_ALL(int *nx, int *nu, int *nbx, int *nbu, int *ng, int *nsbx, int *nsbu, int *nsg, struct OCP_QP_DIM *dim)
	{

	// loop index
	int ii;

	// N
	int N = dim->N;

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



void OCP_QP_DIM_SET(char *field_name, int stage, int value, struct OCP_QP_DIM *dim)
	{
	if(hpipm_strcmp(field_name, "nx"))
		{ 
		OCP_QP_DIM_SET_NX(stage, value, dim);
		}
	else if(hpipm_strcmp(field_name, "nu"))
		{ 
		OCP_QP_DIM_SET_NU(stage, value, dim);
		}
	else if(hpipm_strcmp(field_name, "nbx"))
		{
		OCP_QP_DIM_SET_NBX(stage, value, dim);
		}
	else if(hpipm_strcmp(field_name, "nbu"))
		{
		OCP_QP_DIM_SET_NBU(stage, value, dim);
		}
	else if(hpipm_strcmp(field_name, "ng"))
		{
		OCP_QP_DIM_SET_NG(stage, value, dim);
		}
	else if(hpipm_strcmp(field_name, "nsbx"))
		{
		OCP_QP_DIM_SET_NSBX(stage, value, dim);
		}
	else if(hpipm_strcmp(field_name, "nsbu"))
		{
		OCP_QP_DIM_SET_NSBU(stage, value, dim);
		}
	else if(hpipm_strcmp(field_name, "nsg"))
		{
		OCP_QP_DIM_SET_NSG(stage, value, dim);
		}
	else 
		{
		printf("error [OCP_QP_DIM_SET]: unknown field name '%s'. Exiting.\n", field_name);
		exit(1);
		}
	return;
	}



void OCP_QP_DIM_SET_NX(int stage, int value, struct OCP_QP_DIM *dim)
	{
	dim->nx[stage] = value;
	return;
	}



void OCP_QP_DIM_SET_NU(int stage, int value, struct OCP_QP_DIM *dim)
	{
	dim->nu[stage] = value;
	return;
	}



void OCP_QP_DIM_SET_NBX(int stage, int value, struct OCP_QP_DIM *dim)
	{
	dim->nbx[stage] = value;
	dim->nb[stage] = dim->nbx[stage] + dim->nbu[stage];
	return;
	}



void OCP_QP_DIM_SET_NBU(int stage, int value, struct OCP_QP_DIM *dim)
	{
	dim->nbu[stage] = value;
	dim->nb[stage] = dim->nbx[stage] + dim->nbu[stage];
	return;
	}



void OCP_QP_DIM_SET_NG(int stage, int value, struct OCP_QP_DIM *dim)
	{
	dim->ng[stage] = value;
	return;
	}



void OCP_QP_DIM_SET_NSBX(int stage, int value, struct OCP_QP_DIM *dim)
	{
	dim->nsbx[stage] = value;
	dim->ns[stage] = dim->nsbx[stage] + dim->nsbu[stage] + dim->nsg[stage];
	return;
	}



void OCP_QP_DIM_SET_NSBU(int stage, int value, struct OCP_QP_DIM *dim)
	{
	dim->nsbu[stage] = value;
	dim->ns[stage] = dim->nsbx[stage] + dim->nsbu[stage] + dim->nsg[stage];
	return;
	}



void OCP_QP_DIM_SET_NSG(int stage, int value, struct OCP_QP_DIM *dim)
	{
	dim->nsg[stage] = value;
	dim->ns[stage] = dim->nsbx[stage] + dim->nsbu[stage] + dim->nsg[stage];
	return;
	}
