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



int d_memsize_erk_int(struct d_erk_arg *erk_arg, int nx, int np, int nf_max, int na_max)
	{

	int ns = erk_arg->rk_data->ns;

	int nX = nx*(1+nf_max);

	int steps = erk_arg->steps;

	int size = 0;

	size += 1*np*sizeof(double); // p
	if(na_max>0)
		{
		size += 1*nX*(steps+1)*sizeof(double); // x
		size += 1*ns*nX*steps*sizeof(double); // K
		size += 1*nX*ns*sizeof(double); // xt
		size += 1*nf_max*(steps+1)*sizeof(double); // l // XXX *na_max ???
		size += 1*(nx+nf_max)*sizeof(double); // adj_in // XXX *na_max ???
		size += 1*nf_max*ns*sizeof(double); // adj_tmp // XXX *na_max ???
		}
	else
		{
		size += 1*nX*sizeof(double); // x
		size += 1*ns*nX*sizeof(double); // K
		size += 1*nX*sizeof(double); // xt
		}

	return size;

	}



void d_create_erk_int(struct d_erk_arg *erk_arg, int nx, int np, int nf_max, int na_max, struct d_erk_workspace *ws, void *mem)
	{

	ws->erk_arg = erk_arg;
	ws->nx = nx;
	ws->np = np;
	ws->nf_max = nf_max;
	ws->na_max = na_max;

	int ns = erk_arg->rk_data->ns;

	int nX = nx*(1+nf_max);

	int steps = erk_arg->steps;

	double *d_ptr = mem;

	//
	ws->p = d_ptr;
	d_ptr += np;
	//
	if(na_max>0)
		{
		//
		ws->x = d_ptr;
		d_ptr += nX*(steps+1);
		//
		ws->K = d_ptr;
		d_ptr += ns*nX*steps;
		//
		ws->xt = d_ptr;
		d_ptr += nX*ns;
		//
		ws->l = d_ptr;
		d_ptr += nf_max*(steps+1);
		//
		ws->adj_in = d_ptr;
		d_ptr += nx+nf_max;
		//
		ws->adj_tmp = d_ptr;
		d_ptr += nf_max*ns;
		}
	else
		{
		//
		ws->K = d_ptr;
		d_ptr += ns*nX;
		//
		ws->x = d_ptr;
		d_ptr += nX;
		//
		ws->xt = d_ptr;
		d_ptr += nX;
		}


	ws->memsize = d_memsize_erk_int(erk_arg, nx, np, nf_max, na_max);


	char *c_ptr = (char *) d_ptr;


#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + ws->memsize)
		{
		printf("\nCreate_erk_int: outsize memory bounds! %p %p\n\n", c_ptr, ((char *) mem) + ws->memsize);
		exit(1);
		}
#endif


	return;

	}



void d_init_erk_int(int nf, int na, double *x0, double *p0, double *fs0, double *bs0, void (*vde_for)(int t, double *x, double *p, void *ode_args, double *xdot), void (*vde_adj)(int t, double *adj_in, void *ode_args, double *adj_out), void *ode_args, struct d_erk_workspace *ws)
	{

	int ii;

	ws->nf = nf;
	ws->na = na;

	int nx = ws->nx;
	int np = ws->np;

	int nX = nx*(1+nf);
	int nA = np+nx; // XXX

	int steps = ws->erk_arg->steps;

	double *x = ws->x;
	double *p = ws->p;
	double *l = ws->l;

	for(ii=0; ii<nx; ii++)
		x[ii] = x0[ii];

	for(ii=0; ii<nx*nf; ii++)
		x[nx+ii] = fs0[ii];

	for(ii=0; ii<np; ii++)
		p[ii] = p0[ii];
	
	if(na>0) // TODO what if na>1 !!!
		{
		for(ii=0; ii<np; ii++)
			l[nA*steps+ii] = 0.0;
		for(ii=0; ii<nx; ii++)
			l[nA*steps+np+ii] = bs0[ii];
		}
	
//	ws->ode = ode;
	ws->vde_for = vde_for;
	ws->vde_adj = vde_adj;
	ws->ode_args = ode_args;

//	d_print_mat(1, nx*nf, x, 1);
//	d_print_mat(1, np, p, 1);
//	printf("\n%p %p\n", ode, ode_args);

	return;

	}



#if 0
void d_update_p_erk_int(double *p0, struct d_erk_workspace *ws)
	{

	int ii;

	int np = ws->np;

	double *p = ws->p;

	for(ii=0; ii<np; ii++)
		p[ii] = p0[ii];
	
	return;

	}
#endif



void d_erk_int(struct d_erk_workspace *ws)
	{

	int steps = ws->erk_arg->steps;
	double h = ws->erk_arg->h;

	struct d_rk_data *rk_data = ws->erk_arg->rk_data;
	int nx = ws->nx;
	int np = ws->np;
	int nf = ws->nf;
	int na = ws->na;
	double *K0 = ws->K;
	double *x0 = ws->x;
	double *x1 = ws->x;
	double *p = ws->p;
	double *xt = ws->xt;
	double *adj_in = ws->adj_in;
	double *adj_tmp = ws->adj_tmp;

	double *l0, *l1;

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
	int nA = nx+np; // XXX

//printf("\nnf %d na %d nX %d nA %d\n", nf, na, nX, nA);
	// forward sweep

	// TODO no need to save the entire [x Su Sx] & sens, but only [x] & sens !!!

	t = 0.0; // TODO plus time of multiple-shooting stage !!!
	for(step=0; step<steps; step++)
		{
		if(na>0)
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
					a *= h;
#if 0
					sK.pa = K0+ii*nX; // XXX
					daxpy_libstr(nX, a, &sK, 0, &sxt, 0, &sxt, 0); // XXX
#else
					for(jj=0; jj<nX; jj++)
						xt[jj] += a*K0[jj+ii*(nX)];
#endif
					}
				}
			ws->vde_for(t+h*C_rk[ss], xt, p, ws->ode_args, K0+ss*(nX));
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

	if(na>0)
		{
		t = steps*h; // TODO plus time of multiple-shooting stage !!!
		for(step=steps-1; step>=0; step--)
			{
			l0 = ws->l + step*nA;
			l1 = ws->l + (step+1)*nA;
			x0 = ws->x + step*nX;
			K0 = ws->K + ns*step*nX; // XXX save all x insead !!!
			// TODO save all x instead of K !!!
			for(ss=ns-1; ss>=0; ss--)
				{
				// x
				for(ii=0; ii<nx; ii++)
					adj_in[0+ii] = x0[ii];
				for(ii=0; ii<ss; ii++)
					{
					a = A_rk[ss+ns*ii];
					if(a!=0)
						{
						a *= h;
						for(jj=0; jj<nx; jj++)
							adj_in[0+jj] += a*K0[jj+ii*(nX)];
						}
					}
				// l
				b = h*B_rk[ss];
				for(ii=0; ii<nx; ii++)
					adj_in[nx+ii] = b*l1[np+ii];
				for(ii=ss+1; ii<ns; ii++)
					{
					a = A_rk[ii+ns*ss];
					if(a!=0)
						{
						a *= h;
						for(jj=0; jj<nx; jj++)
							adj_in[nx+jj] += a*adj_tmp[np+jj+ii*nA];
						}
					}
				// p
				for(ii=0; ii<np; ii++)
					adj_in[nx+nx+ii] = p[ii];
				// adj_vde
				ws->vde_adj(t+h*C_rk[ss], adj_in, ws->ode_args, adj_tmp+ss*nA);
				}
			// erk step
			for(ii=0; ii<nA; ii++) // TODO move in the erk step !!!
				l0[ii] = l1[ii];
			for(ss=0; ss<ns; ss++)
				for(ii=0; ii<nA; ii++)
					l0[ii] += adj_tmp[ii+ss*nA];
			t -= h;
			}
		}

	return;

	}




