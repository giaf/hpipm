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



#include <stdlib.h>
#include <stdio.h>

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_d_aux.h>
#include <blasfeo_d_blas.h>

#include <hpipm_d_sim_rk.h>
#include <hpipm_d_sim_erk.h>
#include <hpipm_aux_mem.h>



int d_sim_erk_arg_memsize()
	{
	return 0;
	}



void d_sim_erk_arg_create(struct d_sim_erk_arg *erk_arg, void *mem)
	{

	// zero memory (to avoid corrupted memory like e.g. NaN)
	int memsize = d_sim_erk_arg_memsize();
	hpipm_zero_memset(memsize, mem);

	erk_arg->memsize = memsize;

	return;

	}



void d_sim_erk_arg_set_all(struct d_sim_rk_data *rk_data, double h, int steps, struct d_sim_erk_arg *erk_arg)
	{

	erk_arg->rk_data = rk_data;
	erk_arg->h = h;
	erk_arg->steps = steps;
//	erk_arg->for_sens = for_sens;
//	erk_arg->adj_sens = adj_sens;

	return;

	}



int d_sim_erk_ws_memsize(struct d_sim_erk_arg *erk_arg, int nx, int np, int nf_max, int na_max)
	{

	int ns = erk_arg->rk_data->ns;

	int nX = nx*(1+nf_max);

	int steps = erk_arg->steps;

	int size = 0;

	size += 1*np*sizeof(double); // p
	size += 1*nX*sizeof(double); // x_for
	size += 1*ns*nX*sizeof(double); // K
	size += 1*nX*sizeof(double); // x_tmp
	if(na_max>0)
		{
//		size += 1*nX*(steps+1)*sizeof(double); // x_traj XXX
		size += 1*nx*(ns*steps+1)*sizeof(double); // x_traj
//		size += 1*ns*nX*steps*sizeof(double); // K
//		size += 1*nX*ns*sizeof(double); // x_tmp
		size += 1*nf_max*(steps+1)*sizeof(double); // l // XXX *na_max ???
		size += 1*(nx+nf_max)*sizeof(double); // adj_in // XXX *na_max ???
		size += 1*nf_max*ns*sizeof(double); // adj_tmp // XXX *na_max ???
		}

	return size;

	}



void d_sim_erk_ws_create(struct d_sim_erk_arg *erk_arg, int nx, int np, int nf_max, int na_max, struct d_sim_erk_ws *work, void *mem)
	{

	// zero memory (to avoid corrupted memory like e.g. NaN)
	int memsize = d_sim_erk_ws_memsize(erk_arg, nx, np, nf_max, na_max);
	hpipm_zero_memset(memsize, mem);

	work->erk_arg = erk_arg;
	work->nx = nx;
	work->np = np;
	work->nf_max = nf_max;
	work->na_max = na_max;

	int ns = erk_arg->rk_data->ns;

	int nX = nx*(1+nf_max);

	int steps = erk_arg->steps;

	double *d_ptr = mem;

	//
	work->p = d_ptr;
	d_ptr += np;
	//
	work->x_for = d_ptr;
	d_ptr += nX;
	//
	work->K = d_ptr;
	d_ptr += ns*nX;
	//
	work->x_tmp = d_ptr;
	d_ptr += nX;
	//
	if(na_max>0)
		{
		//
//		work->x_for = d_ptr;
//		d_ptr += nX*(steps+1);
		//
		work->x_traj = d_ptr;
		d_ptr += nx*(ns*steps+1);
		//
//		work->K = d_ptr;
//		d_ptr += ns*nX*steps;
		//
//		work->x_tmp = d_ptr;
//		d_ptr += nX*ns;
		//
		work->l = d_ptr;
		d_ptr += nf_max*(steps+1);
		//
		work->adj_in = d_ptr;
		d_ptr += nx+nf_max;
		//
		work->adj_tmp = d_ptr;
		d_ptr += nf_max*ns;
		}
	
	// default init
	work->nf = 0;
	work->na = 0;


	work->memsize = memsize;


	char *c_ptr = (char *) d_ptr;


#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + work->memsize)
		{
		printf("\nSIM_ERK_WS_CREATE: outsize memory bounds! %p %p\n\n", c_ptr, ((char *) mem) + work->memsize);
		exit(1);
		}
#endif


	return;

	}



void d_sim_erk_ws_set_all(int nf, int na, double *x, double *fs, double *bs, double *p, void (*ode)(int t, double *x, double *p, void *ode_args, double *xdot), void (*vde_for)(int t, double *x, double *p, void *ode_args, double *xdot), void (*vde_adj)(int t, double *adj_in, void *ode_args, double *adj_out), void *ode_args, struct d_sim_erk_ws *work)
	{

	int ii;

	// TODO check against nf_max and na_max
	if(nf>work->nf_max)
		{
		printf("\nerror: SIM_ERK_WS_SET_ALL: nf>nf_max: %d > %d\n", nf, work->nf_max);
		exit(1);
		}
	if(nf>work->nf_max | na>work->na_max)
		{
		printf("\nerror: SIM_ERK_WS_SET_ALL: na>na_max: %d > %d\n", na, work->na_max);
		exit(1);
		}

	work->nf = nf;
	work->na = na;

	int nx = work->nx;
	int np = work->np;

	int nX = nx*(1+nf);
	int nA = np+nx; // XXX

	int steps = work->erk_arg->steps;

	double *x_for = work->x_for;
	double *p_ws = work->p;
	double *l = work->l;

	for(ii=0; ii<nx; ii++)
		x_for[ii] = x[ii];

	for(ii=0; ii<nx*nf; ii++)
		x_for[nx+ii] = fs[ii];

	for(ii=0; ii<np; ii++)
		p_ws[ii] = p[ii];
	
	if(na>0) // TODO what if na>1 !!!
		{
		for(ii=0; ii<np; ii++)
			l[nA*steps+ii] = 0.0;
		for(ii=0; ii<nx; ii++)
			l[nA*steps+np+ii] = bs[ii];
		}
	
	work->ode = ode;
	work->vde_for = vde_for;
	work->vde_adj = vde_adj;

	work->ode_args = ode_args;

//	d_print_mat(1, nx*nf, x, 1);
//	d_print_mat(1, np, p, 1);
//	printf("\n%p %p\n", ode, ode_args);

	return;

	}



void d_sim_erk_ws_set_nf(int *nf, struct d_sim_erk_ws *work)
	{

	work->nf = *nf;
	
	return;

	}



// state
void d_sim_erk_ws_set_x(double *x, struct d_sim_erk_ws *work)
	{

	int ii;

	int nx = work->nx;

	double *x_ws = work->x_for;

	for(ii=0; ii<nx; ii++)
		x_ws[ii] = x[ii];
	
	return;

	}



// forward sensitivities
void d_sim_erk_ws_set_fs(double *fs, struct d_sim_erk_ws *work)
	{

	int ii;

	int nx = work->nx;
	int nf = work->nf;

	double *fs_ws = work->x_for+nx;

	for(ii=0; ii<nx*nf; ii++)
		fs_ws[ii] = fs[ii];
	
	return;

	}



// state
void d_sim_erk_ws_get_x(struct d_sim_erk_ws *work, double *x)
	{

	int ii;

	int nx = work->nx;

	double *x_ws = work->x_for;

	for(ii=0; ii<nx; ii++)
		x[ii] = x_ws[ii];
	
	return;

	}



void d_sim_erk_ws_set_p(double *p, struct d_sim_erk_ws *work)
	{

	int ii;

	int np = work->np;

	double *p_ws = work->p;

	for(ii=0; ii<np; ii++)
		p_ws[ii] = p[ii];
	
	return;

	}



void d_sim_erk_ws_set_ode(void (*ode)(int t, double *x, double *p, void *ode_args, double *xdot), struct d_sim_erk_ws *work)
	{

	work->ode = ode;
	
	return;

	}



void d_sim_erk_ws_set_vde_for(void (*vde_for)(int t, double *x, double *p, void *ode_args, double *xdot), struct d_sim_erk_ws *work)
	{

	work->vde_for = vde_for;
	
	return;

	}



void d_sim_erk_ws_set_ode_args(void *ode_args, struct d_sim_erk_ws *work)
	{

	work->ode_args = ode_args;

	return;

	}



// forward sensitivities
void d_sim_erk_ws_get_fs(struct d_sim_erk_ws *work, double *fs)
	{

	int ii;

	int nx = work->nx;
	int nf = work->nf;

	double *fs_ws = work->x_for+nx;

	for(ii=0; ii<nx*nf; ii++)
		fs[ii] = fs_ws[ii];
	
	return;

	}



void d_sim_erk_solve(struct d_sim_erk_arg *erk_arg, struct d_sim_erk_ws *work)
	{

	int steps = work->erk_arg->steps;
	double h = work->erk_arg->h;

	struct d_sim_rk_data *rk_data = erk_arg->rk_data;
	int nx = work->nx;
	int np = work->np;
	int nf = work->nf;
	int na = work->na;
	double *K0 = work->K;
	double *x0 = work->x_for;
	double *x1 = work->x_for;
	double *x_traj = work->x_traj;
	double *p = work->p;
	double *x_tmp = work->x_tmp;
	double *adj_in = work->adj_in;
	double *adj_tmp = work->adj_tmp;

	double *l0, *l1;

	int ns = rk_data->ns;
	double *A_rk = rk_data->A_rk;
	double *B_rk = rk_data->B_rk;
	double *C_rk = rk_data->C_rk;

	struct blasfeo_dvec sxt; // XXX
	struct blasfeo_dvec sK; // XXX
	sxt.pa = x_tmp; // XXX

	int ii, jj, step, ss;
	double t, a, b;

	int nX = nx*(1+nf);
	int nA = nx+np; // XXX

	if(na>0)
		{
		printf("\nSIM_ERK_SOLVE: na>0 not implemented yet!\n");
		exit(1);
		}

//printf("\nnf %d na %d nX %d nA %d\n", nf, na, nX, nA);
	// forward sweep

	// TODO no need to save the entire [x Su Sx] & sens, but only [x] & sens !!!

	t = 0.0; // TODO plus time of multiple-shooting stage !!!
	if(na>0)
		{
		x_traj = work->x_traj;
		for(ii=0; ii<nx; ii++)
			x_traj[ii] = x0[ii];
		x_traj += nx;
		}
	for(step=0; step<steps; step++)
		{
//		if(na>0)
//			{
//			x0 = work->x_for + step*nX;
//			x1 = work->x_for + (step+1)*nX;
//			for(ii=0; ii<nX; ii++)
//				x1[ii] = x0[ii];
//			K0 = work->K + ns*step*nX;
//			}
		for(ss=0; ss<ns; ss++)
			{
			for(ii=0; ii<nX; ii++)
				x_tmp[ii] = x0[ii];
			for(ii=0; ii<ss; ii++)
				{
				a = A_rk[ss+ns*ii];
				if(a!=0)
					{
					a *= h;
#if 0
					sK.pa = K0+ii*nX; // XXX
					blasfeo_daxpy(nX, a, &sK, 0, &sxt, 0, &sxt, 0); // XXX
#else
					for(jj=0; jj<nX; jj++)
						x_tmp[jj] += a*K0[jj+ii*(nX)];
#endif
					}
				}
			if(na>0)
				{
				for(ii=0; ii<nx; ii++)
					x_traj[ii] = x_tmp[ii];
				x_traj += nx;
				}
			if(nf>0)
				{
				work->vde_for(t+h*C_rk[ss], x_tmp, p, work->ode_args, K0+ss*(nX));
				}
			else
				{
				work->ode(t+h*C_rk[ss], x_tmp, p, work->ode_args, K0+ss*(nX));
				}
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
		x_traj = work->x_traj + nx*ns*steps;
		t = steps*h; // TODO plus time of multiple-shooting stage !!!
		for(step=steps-1; step>=0; step--)
			{
			l0 = work->l + step*nA;
			l1 = work->l + (step+1)*nA;
			x0 = work->x_for + step*nX;
			K0 = work->K + ns*step*nX; // XXX save all x insead !!!
			// TODO save all x instead of K !!!
			for(ss=ns-1; ss>=0; ss--)
				{
				// x
				for(ii=0; ii<nx; ii++)
					adj_in[0+ii] = x_traj[ii];
				x_traj -= nx;
//				for(ii=0; ii<nx; ii++)
//					adj_in[0+ii] = x0[ii];
//				for(ii=0; ii<ss; ii++)
//					{
//					a = A_rk[ss+ns*ii];
//					if(a!=0)
//						{
//						a *= h;
//						for(jj=0; jj<nx; jj++)
//							adj_in[0+jj] += a*K0[jj+ii*(nX)];
//						}
//					}
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
				work->vde_adj(t+h*C_rk[ss], adj_in, work->ode_args, adj_tmp+ss*nA);
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
