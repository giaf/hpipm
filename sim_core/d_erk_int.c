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



#if defined(RUNTIME_CHECKS)
#include <stdlib.h>
#include <stdio.h>
#endif

#include "../include/hpipm_d_rk_int.h"
#include "../include/hpipm_d_erk_int.h"



int d_memsize_erk_int(struct d_rk_data *rk_data, int nx, int nf, int np)
	{

	int ns = rk_data->ns;

	int nX = nx*(1+nf);

	int size = 0;

	size += 1*ns*nX*sizeof(double); // K
	size += 2*nX*sizeof(double); // x xt
	size += 1*np*sizeof(double); // p

	return size;

	}



void d_create_erk_int(struct d_rk_data *rk_data, int nx, int nf, int np, struct d_erk_workspace *ws, void *mem)
	{

	ws->rk_data = rk_data;
	ws->nx = nx;
	ws->nf = nf;
	ws->np = np;

	int ns = rk_data->ns;

	int nX = nx*(1+nf);

	double *d_ptr = mem;

	//
	ws->K = d_ptr;
	d_ptr += ns*nX;
	//
	ws->x = d_ptr;
	d_ptr += nX;
	//
	ws->xt = d_ptr;
	d_ptr += nX;
	//
	ws->p = d_ptr;
	d_ptr += np;


	ws->memsize = d_memsize_erk_int(rk_data, nx, nf, np);


	char *c_ptr = (char *) d_ptr;


#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + ws->memsize)
		{
		printf("\nCreate_erk_int: outsize memory bounds!\n\n");
		exit(1);
		}
#endif


	return;

	}



void d_init_erk_int(double *x0, double *fs0, double *p0, void (*ode)(int t, double *x, double *p, void *ode_args, double *xdot), void *ode_args, struct d_erk_workspace *ws)
	{

	int ii;

	int nx = ws->nx;
	int nf = ws->nf;
	int np = ws->np;

	int nX = nx*(1+nf);

	double *x = ws->x;
	double *p = ws->p;

	for(ii=0; ii<nx*nf; ii++)
		x[ii] = fs0[ii];

	for(ii=0; ii<nx; ii++)
		x[nx*nf+ii] = x0[ii];

	for(ii=0; ii<np; ii++)
		p[ii] = p0[ii];
	
	ws->ode = ode;
	ws->ode_args = ode_args;

//	d_print_mat(1, nx*nf, x, 1);
//	d_print_mat(1, np, p, 1);
//	printf("\n%p %p\n", ode, ode_args);

	return;

	}



void d_update_p_erk_int(double *p0, struct d_erk_workspace *ws)
	{

	int ii;

	int np = ws->np;

	double *p = ws->p;

	for(ii=0; ii<np; ii++)
		p[ii] = p0[ii];
	
	return;

	}



void d_erk_int(struct d_erk_args *erk_args, struct d_erk_workspace *ws)
	{

	int steps = erk_args->steps;
	double h = erk_args->h;

	struct d_rk_data *rk_data = ws->rk_data;
	int nx = ws->nx;
	int nf = ws->nf;
	double *K = ws->K;
	double *x = ws->x;
	double *p = ws->p;
	double *xt = ws->xt;

	int ns = rk_data->ns;
	double *A_rk = rk_data->A_rk;
	double *B_rk = rk_data->B_rk;
	double *C_rk = rk_data->C_rk;

	int ii, jj, step, ss;
	double t, a, b;

	int nX = nx*(1+nf);

//d_print_mat(nx, nf+1, x, nx);

	t = 0.0;
	for(step=0; step<steps; step++)
		{
		for(ss=0; ss<ns; ss++)
			{
			for(ii=0; ii<nX; ii++)
				xt[ii] = x[ii];
			for(ii=0; ii<ss; ii++)
				{
				a = A_rk[ss+ns*ii];
				if(a!=0)
					{
					a *= h;
					for(jj=0; jj<nX; jj++)
						xt[jj] += a*K[jj+ii*(nX)];
					}
				}
			ws->ode(t+h*C_rk[ss], xt, p, ws->ode_args, K+ss*(nX));
			}
//d_print_mat(nx, ns*(nf+1), K, nx);
		for(ss=0; ss<ns; ss++)
			{
			b = h*B_rk[ss];
			for(ii=0; ii<nX; ii++)
				x[ii] += b*K[ii+ss*(nX)];
			}
//d_print_mat(nx, nf+1, x, nx);
		t += h;
		}

//d_print_mat(nx, nf+1, x, nx);
	return;

	}




