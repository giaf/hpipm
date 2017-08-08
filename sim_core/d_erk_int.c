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



#include "../include/hpipm_d_erk_int.h"



int d_memsize_erk_data(int ns)
	{

	int size = 0;

	size += 1*ns*ns*sizeof(double); // A
	size += 2*ns*sizeof(double); // B C

	return size;

	}



void d_create_erk_data(int ns, struct d_erk_data *erk_data, void *memory)
	{

	erk_data->ns = ns;

	double *d_ptr = memory;

	//
	erk_data->A_rk = d_ptr;
	d_ptr += ns*ns;
	//
	erk_data->B_rk = d_ptr;
	d_ptr += ns;
	//
	erk_data->C_rk = d_ptr;
	d_ptr += ns;

	return;

	}



void d_cvt_colmaj_to_erk_data(double *A, double *B, double *C, struct d_erk_data *erk_data)
	{

	int ii, jj;

	int ns = erk_data->ns;
	double *A_rk = erk_data->A_rk;
	double *B_rk = erk_data->B_rk;
	double *C_rk = erk_data->C_rk;

	for(jj=0; jj<ns; jj++)
		for(ii=0; ii<ns; ii++)
			A_rk[ii+ns*jj] = A[ii+ns*jj];
	for(ii=0; ii<ns; ii++)
		B_rk[ii] = B[ii];
	for(ii=0; ii<ns; ii++)
		C_rk[ii] = C[ii];

	return;

	}



void d_cvt_rowmaj_to_erk_data(double *A, double *B, double *C, struct d_erk_data *erk_data)
	{

	int ii, jj;

	int ns = erk_data->ns;
	double *A_rk = erk_data->A_rk;
	double *B_rk = erk_data->B_rk;
	double *C_rk = erk_data->C_rk;

	for(jj=0; jj<ns; jj++)
		for(ii=0; ii<ns; ii++)
			A_rk[jj+ns*ii] = A[ii+ns*jj];
	for(ii=0; ii<ns; ii++)
		B_rk[ii] = B[ii];
	for(ii=0; ii<ns; ii++)
		C_rk[ii] = C[ii];

	return;

	}



int d_memsize_erk_int(struct d_erk_data *erk_data, int nx)
	{

	int ns = erk_data->ns;

	int size = 0;

	size += ns*nx*sizeof(double); // K
	size += nx*sizeof(double); // xt

	return size;

	}



void d_create_erk_int(struct d_erk_data *erk_data, int nx, struct d_erk_workspace *workspace, void *memory)
	{

	workspace->erk_data = erk_data;
	workspace->nx = nx;

	int ns = erk_data->ns;

	double *d_ptr = memory;

	//
	workspace->K = d_ptr;
	d_ptr += ns*nx;
	//
	workspace->xt = d_ptr;
	d_ptr += nx;

	return;

	}



void d_erk_int(double *x0, double *p, double *xe, void (*ode)(int t, double *x, double *p, void *ode_args, double *xdot), void *ode_args, struct d_erk_args *erk_args, struct d_erk_workspace *workspace)
	{

	int steps = erk_args->steps;
	double h = erk_args->h;

	struct d_erk_data *erk_data = workspace->erk_data;
	int nx = workspace->nx;
	double *K = workspace->K;
	double *xt = workspace->xt;

	int ns = erk_data->ns;
	double *A_rk = erk_data->A_rk;
	double *B_rk = erk_data->B_rk;
	double *C_rk = erk_data->C_rk;

	int ii, jj, step, ss;
	double t, a, b;

	for(ii=0; ii<nx; ii++)
		xe[ii] = x0[ii];

	t = 0.0;
	for(step=0; step<steps; step++)
		{
		for(ss=0; ss<ns; ss++)
			{
			for(ii=0; ii<nx; ii++)
				xt[ii] = xe[ii];
			for(ii=0; ii<ss; ii++)
				{
				a = A_rk[ss+ns*ii];
				if(a!=0)
					{
					a *= h;
					for(jj=0; jj<nx; jj++)
						xt[jj] += a*K[jj+ii*nx];
					}
				}
			ode(t+h*C_rk[ss], xt, p, ode_args, K+ss*nx);
			}
		for(ss=0; ss<ns; ss++)
			{
			b = h*B_rk[ss];
			for(ii=0; ii<nx; ii++)
				xe[ii] += b*K[ii+ss*nx];
			}
		t += h;
		}

	return;

	}




