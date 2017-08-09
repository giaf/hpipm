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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_v_aux_ext_dep.h>
#include <blasfeo_d_aux_ext_dep.h>
#include <blasfeo_i_aux_ext_dep.h>
#include <blasfeo_d_aux.h>

#include "../include/hpipm_d_erk_int.h"

#include "d_tools.h"



#if ! defined(EXT_DEP)
/* creates a zero matrix */
void d_zeros(double **pA, int row, int col)
	{
	*pA = malloc((row*col)*sizeof(double));
	double *A = *pA;
	int i;
	for(i=0; i<row*col; i++) A[i] = 0.0;
	}
/* frees matrix */
void d_free(double *pA)
	{
	free( pA );
	}
/* prints a matrix in column-major format */
void d_print_mat(int m, int n, double *A, int lda)
	{
	int i, j;
	for(i=0; i<m; i++)
		{
		for(j=0; j<n; j++)
			{
			printf("%9.5f ", A[i+lda*j]);
			}
		printf("\n");
		}
	printf("\n");
	}	
/* prints the transposed of a matrix in column-major format */
void d_print_tran_mat(int row, int col, double *A, int lda)
	{
	int i, j;
	for(j=0; j<col; j++)
		{
		for(i=0; i<row; i++)
			{
			printf("%9.5f ", A[i+lda*j]);
			}
		printf("\n");
		}
	printf("\n");
	}	
/* prints a matrix in column-major format (exponential notation) */
void d_print_e_mat(int m, int n, double *A, int lda)
	{
	int i, j;
	for(i=0; i<m; i++)
		{
		for(j=0; j<n; j++)
			{
			printf("%e\t", A[i+lda*j]);
			}
		printf("\n");
		}
	printf("\n");
	}	
/* prints the transposed of a matrix in column-major format (exponential notation) */
void d_print_e_tran_mat(int row, int col, double *A, int lda)
	{
	int i, j;
	for(j=0; j<col; j++)
		{
		for(i=0; i<row; i++)
			{
			printf("%e\t", A[i+lda*j]);
			}
		printf("\n");
		}
	printf("\n");
	}	
/* creates a zero matrix aligned */
void int_zeros(int **pA, int row, int col)
	{
	void *temp = malloc((row*col)*sizeof(int));
	*pA = temp;
	int *A = *pA;
	int i;
	for(i=0; i<row*col; i++) A[i] = 0;
	}
/* frees matrix */
void int_free(int *pA)
	{
	free( pA );
	}
/* prints a matrix in column-major format */
void int_print_mat(int row, int col, int *A, int lda)
	{
	int i, j;
	for(i=0; i<row; i++)
		{
		for(j=0; j<col; j++)
			{
			printf("%d ", A[i+lda*j]);
			}
		printf("\n");
		}
	printf("\n");
	}	
#endif



/************************************************ 
Mass-spring system: nx/2 masses connected each other with springs (in a row), and the first and the last one to walls. nu (<=nx) controls act on the first nu masses. The system is sampled with sampling time Ts. 
************************************************/
void mass_spring_system(int nx, int nu, double *Ac, double *Bc)
	{

	int nx2 = nx*nx;

	int pp = nx/2; // number of masses
	
	double *T; d_zeros(&T, pp, pp);
	int ii;
	for(ii=0; ii<pp; ii++) T[ii*(pp+1)] = -2;
	for(ii=0; ii<pp-1; ii++) T[ii*(pp+1)+1] = 1;
	for(ii=1; ii<pp; ii++) T[ii*(pp+1)-1] = 1;

	double *Z; d_zeros(&Z, pp, pp);
	double *I; d_zeros(&I, pp, pp); for(ii=0; ii<pp; ii++) I[ii*(pp+1)]=1.0; // = eye(pp);
	dmcopy(pp, pp, Z, pp, Ac, nx);
	dmcopy(pp, pp, T, pp, Ac+pp, nx);
	dmcopy(pp, pp, I, pp, Ac+pp*nx, nx);
	dmcopy(pp, pp, Z, pp, Ac+pp*(nx+1), nx); 
	free(T);
	free(Z);
	free(I);
	
	d_zeros(&I, nu, nu); for(ii=0; ii<nu; ii++) I[ii*(nu+1)]=1.0; //I = eye(nu);
	dmcopy(nu, nu, I, nu, Bc+pp, nx);
	free(I);

	return;

	}



struct d_linear_system
	{
	double *Ac;
	double *Bc;
	int nx;
	int nu;
	int memsize;
	};



int d_memsize_linear_system(int nx, int nu)
	{
	int size = 0;
	size += (nx*nx+nx*nu)*sizeof(double);
	return size;
	}



void d_create_linear_system(int nx, int nu, struct d_linear_system *ls, void *memory)
	{
	ls->nx = nx;
	ls->nu = nu;
	char * c_ptr = (char *) memory;
	ls->Ac = (double *) c_ptr;
	c_ptr += nx*nx*sizeof(double);
	ls->Bc = (double *) c_ptr;
	c_ptr += nx*nu*sizeof(double);
	return;
	}



void d_cvt_colmaj_to_linear_system(double *A,  double *B, struct d_linear_system *ls)
	{
	int ii;
	int nx = ls->nx;
	int nu = ls->nu;
	double *Ac = ls->Ac;
	double *Bc = ls->Bc;
	for(ii=0; ii<nx*nx; ii++)
		Ac[ii] = A[ii];
	for(ii=0; ii<nx*nu; ii++)
		Bc[ii] = B[ii];
	return;
	}
		


void d_linear_ode(int t, double *x, double *u, void *ode_args, double *xdot)
	{
	struct d_linear_system *ls = ode_args;
	int ii, jj;
	int nx = ls->nx;
	int nu = ls->nu;
	double *Ac = ls->Ac;
	double *Bc = ls->Bc;
	for(ii=0; ii<nx; ii++)
		xdot[ii] = 0.0;
	for(jj=0; jj<nx; jj++)
		for(ii=0; ii<nx; ii++)
			xdot[ii] += Ac[ii+nx*jj] * x[jj];
	for(jj=0; jj<nu; jj++)
		for(ii=0; ii<nx; ii++)
			xdot[ii] += Bc[ii+nx*jj] * u[jj];
	return;
	}



void d_linear_vde(int t, double *x, double *u, void *ode_args, double *xdot)
	{
	struct d_linear_system *ls = ode_args;
	int ii, jj, kk;
	int nx = ls->nx;
	int nu = ls->nu;
	double *Ac = ls->Ac;
	double *Bc = ls->Bc;
	for(ii=0; ii<nx*(1+nx+nu); ii++)
		xdot[ii] = 0.0;
	for(kk=0; kk<1+nx+nu; kk++)
		for(jj=0; jj<nx; jj++)
			for(ii=0; ii<nx; ii++)
				xdot[ii+nx*kk] += Ac[ii+nx*jj] * x[jj+nx*kk];
	for(jj=0; jj<nu; jj++)
		for(ii=0; ii<nx; ii++)
			xdot[ii] += Bc[ii+nx*jj] * u[jj];
	double *tmp = xdot+nx;
	for(jj=0; jj<nu; jj++)
		for(ii=0; ii<nx; ii++)
			tmp[ii+nx*jj] += Bc[ii+nx*jj];
	return;
	}



int main()
	{

	int ii, jj;

/************************************************
* problem size
************************************************/	
	
	int nx = 8;
	int nu = 3;

/************************************************
* (continuous time) mass sprint system
************************************************/	
	
	int ls_memsize = d_memsize_linear_system(nx, nu);
	printf("\nls memsize = %d\n", ls_memsize);
	void *ls_memory = malloc(ls_memsize);

	struct d_linear_system ls;
	d_create_linear_system(nx, nu, &ls, ls_memory);

	double *Ac; d_zeros(&Ac, nx, nx);
	double *Bc; d_zeros(&Bc, nx, nu);

	mass_spring_system(nx, nu, Ac, Bc);

	d_print_mat(nx, nx, Ac, nx);
	d_print_mat(nx, nu, Bc, nx);

	d_cvt_colmaj_to_linear_system(Ac, Bc, &ls);

	d_print_mat(nx, nx, ls.Ac, nx);
	d_print_mat(nx, nu, ls.Bc, nx);

	double *x0; d_zeros(&x0, nx, 1);
	x0[0] = 2.5;
	x0[1] = 2.5;

	double *u; d_zeros(&u, nu, 1);
	u[0] = 1.0;

	double *xdot; d_zeros(&xdot, nx, 1);
	d_linear_ode(0, x0, u, &ls, xdot);

	d_print_mat(1, nx, x0, 1);
	d_print_mat(1, nx, xdot, 1);

/************************************************
* (discrete time) mass sprint system
************************************************/	

	double Ts = 0.5;

	double *A = malloc(nx*nx*sizeof(double));
	double *B = malloc(nx*nu*sizeof(double));
	double *T = malloc(nx*nx*sizeof(double));

	for(ii=0; ii<nx*nx; ii++)
		A[ii] = Ts*Ac[ii];
	expm(nx, A);

	for(ii=0; ii<nx*nx; ii++) T[ii] = A[ii];
	for(ii=0; ii<nx; ii++) T[ii*(nx+1)] -= 1.0;
	dgemm_nn_3l(nx, nu, nx, T, nx, Bc, nx, B, nx);

	int info;
	int *ipiv = malloc(nx*sizeof(int));
	dgesv_3l(nx, nu, Ac, nx, ipiv, B, nx, &info);

	d_print_mat(nx, nx, A, nx);
	d_print_mat(nx, nu, B, nx);

/************************************************
* simulation using analytic solution
************************************************/	

	double *xref; d_zeros(&xref, nx, 1);

	for(ii=0; ii<nx; ii++)
		xref[ii] = 0.0;
	for(jj=0; jj<nx; jj++)
		for(ii=0; ii<nx; ii++)
			xref[ii] += A[ii+nx*jj] * x0[jj];
	for(jj=0; jj<nu; jj++)
		for(ii=0; ii<nx; ii++)
			xref[ii] += B[ii+nx*jj] * u[jj];

	printf("\nx analytic\n");
	d_print_mat(1, nx, xref, 1);

/************************************************
* explicit rk4 integrator
************************************************/	
	
	int steps = 10;
	double h = Ts/steps;
	
#if 1
	// rk4
	int ns = 4; // number of stages
	double A_rk[] = {0.0, 0.0, 0.0, 0.0,
	                 0.5, 0.0, 0.0, 0.0,
	                 0.0, 0.5, 0.0, 0.0,
	                 0.0, 0.0, 1.0, 0.0};
	double B_rk[] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
	double C_rk[] = {0.0, 0.5, 0.5, 0.0};
#elif 1
	// midpoint rule
	int ns = 2; // number of stages
	double A_rk[] = {0.0, 0.0,
	                 0.5, 0.0};
	double B_rk[] = {0.0, 1.0};
	double C_rk[] = {0.0, 0.5};
#else
	// explicit euler
	int ns = 1; // number of stages
	double A_rk[] = {0.0};
	double B_rk[] = {1.0};
	double C_rk[] = {0.0};
#endif

	int memsize_erk_data = d_memsize_erk_data(ns);
	printf("\nmemsize erk data %d\n", memsize_erk_data);
	void *memory_erk_data = malloc(memsize_erk_data);

	struct d_erk_data erk_data;
	d_create_erk_data(ns, &erk_data, memory_erk_data);

	d_cvt_rowmaj_to_erk_data(A_rk, B_rk, C_rk, &erk_data);

#if 0
	double *memory_erk = malloc((nx+nx*ns)*sizeof(double));
	
	double *x_erk; d_zeros(&x_erk, nx, 1);
	double *ex_erk; d_zeros(&ex_erk, nx, 1);

	d_erk_int(&erk_data, steps, h, nx, x0, x_erk, &d_linear_ode, &ls, memory_erk);

	for(ii=0; ii<nx; ii++)
		ex_erk[ii] = x_erk[ii] - xref[ii];

	printf("\nx erk\n");
	d_print_mat(1, nx, x_erk, 1);
	printf("\nerror erk\n");
	d_print_e_mat(1, nx, ex_erk, 1);
#else

	int nf = nx+nu;
	int np = nu;

	int memsize_erk_int = d_memsize_erk_int(&erk_data, nx, nf, np);
	printf("\nmemsize erk int %d\n", memsize_erk_int);
	void *memory_erk = malloc(memsize_erk_int);

	struct d_erk_workspace erk_workspace;
	d_create_erk_int(&erk_data, nx, nf, np, &erk_workspace, memory_erk);

	struct d_erk_args erk_args;
	erk_args.steps = steps;
	erk_args.h = h;
	
	double *fs0; d_zeros(&fs0, nx*nf, 1);
	for(ii=0; ii<nx; ii++) fs0[nu*nx+ii*(nx+1)] = 1.0;

	double *ex_erk; d_zeros(&ex_erk, nx*(1+nx+nu), 1);

	d_init_erk_int(x0, fs0, u, &d_linear_vde, &ls, &erk_workspace);

	d_update_p_erk_int(u, &erk_workspace);

	d_erk_int(&erk_args, &erk_workspace);

	double *x_erk = erk_workspace.x;

	for(ii=0; ii<nx; ii++)
		ex_erk[ii] = x_erk[ii] - xref[ii];
	for(ii=0; ii<nx*nu; ii++)
		ex_erk[nx+ii] = x_erk[nx+ii] - B[ii];
	for(ii=0; ii<nx*nx; ii++)
		ex_erk[nx+nx*nu+ii] = x_erk[nx+nx*nu+ii] - A[ii];

	printf("\nx erk\n");
	d_print_mat(1, nx, x_erk, 1);
	d_print_mat(nx, nx, x_erk+nx+nx*nu, nx);
	d_print_mat(nx, nu, x_erk+nx, nx);
	printf("\nerror erk\n");
	d_print_e_mat(1, nx, ex_erk, 1);
	d_print_e_mat(nx, nx, ex_erk+nx+nx*nu, nx);
	d_print_e_mat(nx, nu, ex_erk+nx, nx);
#endif

/************************************************
* free memory
************************************************/	
	
	free(ls_memory);
	free(Ac);
	free(Bc);
	free(x0);
	free(u);
	free(xdot);
	free(A);
	free(B);
	free(T);
	free(ipiv);
	free(xref);
	free(fs0);
	free(memory_erk);
	free(ex_erk);

	return 0;

	}


