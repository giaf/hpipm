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

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_d_aux.h>
#include <blasfeo_d_blas.h>

#include "../include/hpipm_d_rk_int.h"
#include "../include/hpipm_d_erk_int.h"



int d_memsize_erk_int(struct d_rk_data *rk_data, struct d_erk_arg *erk_arg, int nx, int nf, int np)
	{

	int ns = rk_data->ns;

	int nX = nx*(1+nf);

	int steps = erk_arg->steps;

	int size = 0;

	size += 1*nX*sizeof(double); // xt
	size += 1*np*sizeof(double); // p
	if(erk_arg->adj_sens==0)
		{
		size += 1*nX*sizeof(double); // x
		size += 1*ns*nX*sizeof(double); // K
		}
	else // (checkpoint)
		{
		size += 1*nX*(steps+1)*sizeof(double); // x
		size += 1*ns*nX*steps*sizeof(double); // K
		}

	return size;

	}



void d_create_erk_int(struct d_rk_data *rk_data, struct d_erk_arg *erk_arg, int nx, int nf, int np, struct d_erk_workspace *ws, void *mem)
	{

	ws->rk_data = rk_data;
	ws->erk_arg = erk_arg;
	ws->nx = nx;
	ws->nf = nf;
	ws->np = np;

	int ns = rk_data->ns;

	int nX = nx*(1+nf);

	int steps = erk_arg->steps;

	double *d_ptr = mem;

	if(erk_arg->adj_sens==0)
		{
		//
		ws->K = d_ptr;
		d_ptr += ns*nX;
		//
		ws->x = d_ptr;
		d_ptr += nX;
		}
	else
		{
		//
		ws->K = d_ptr;
		d_ptr += ns*nX*steps;
		//
		ws->x = d_ptr;
		d_ptr += nX*(steps+1);
		}
	//
	ws->xt = d_ptr;
	d_ptr += nX;
	//
	ws->p = d_ptr;
	d_ptr += np;


	ws->memsize = d_memsize_erk_int(rk_data, erk_arg, nx, nf, np);


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

	for(ii=0; ii<nx; ii++)
		x[ii] = x0[ii];

	for(ii=0; ii<nx*nf; ii++)
		x[nx+ii] = fs0[ii];

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



void d_erk_int(struct d_erk_workspace *ws)
	{

	int steps = ws->erk_arg->steps;
	double h = ws->erk_arg->h;
	int adj_sens = ws->erk_arg->adj_sens;

	struct d_rk_data *rk_data = ws->rk_data;
	int nx = ws->nx;
	int nf = ws->nf;
	double *K0 = ws->K;
	double *x0 = ws->x;
	double *x1 = ws->x;
	double *p = ws->p;
	double *xt = ws->xt;

	int ns = rk_data->ns;
	double *A_rk = rk_data->A_rk;
	double *B_rk = rk_data->B_rk;
	double *C_rk = rk_data->C_rk;

	struct d_strvec sxt; // XXX
	struct d_strvec sK; // XXX
	sxt.pa = xt; // XXX

	int ii, jj, step, ss;
	double t, a, b;

	int nX = nx*(1+nf);

	// forward sweep

	t = 0.0;
	for(step=0; step<steps; step++)
		{
		if(adj_sens!=0)
			{
			x0 = ws->x + step*nX;
			x1 = ws->x + (step+1)*nX;
			for(ii=0; ii<nX; ii++)
				x1[ii] = x0[ii];
			K0 = ws->K + ns*step*nX;
			}
		for(ss=0; ss<ns; ss++)
			{
			for(ii=0; ii<nX; ii++)
				xt[ii] = x0[ii];
			for(ii=0; ii<ss; ii++)
				{
				a = A_rk[ss+ns*ii];
				if(a!=0)
					{
					sK.pa = K0+ii*nX; // XXX
					a *= h;
#if 0
					daxpy_libstr(nX, a, &sK, 0, &sxt, 0, &sxt, 0); // XXX
#else
					for(jj=0; jj<nX; jj++)
						xt[jj] += a*K0[jj+ii*(nX)];
#endif
					}
				}
			ws->ode(t+h*C_rk[ss], xt, p, ws->ode_args, K0+ss*(nX));
			}
		for(ss=0; ss<ns; ss++)
			{
			b = h*B_rk[ss];
			for(ii=0; ii<nX; ii++)
				x1[ii] += b*K0[ii+ss*(nX)];
			}
		t += h;
		}
	
	// adjoint sweep
	if(adj_sens!=0)
		{
		for(step=steps-1; step>=0; step--)
			{
//			for(ss=ns-1; ss>=0; ss--)
//				{
//				}
			}
		}

	return;

	}




