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
#include <math.h>
#include <sys/time.h>

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_v_aux_ext_dep.h>
#include <blasfeo_d_aux_ext_dep.h>
#include <blasfeo_i_aux_ext_dep.h>
#include <blasfeo_d_aux.h>

#include <hpipm_d_sim_rk.h>
#include <hpipm_d_sim_erk.h>
//#include <hpipm_d_sim_irk.h>
#include <hpipm_d_ocp_qp.h>
//#include <hpipm_d_ocp_qp_sim.h>

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
void d_print_exp_mat(int m, int n, double *A, int lda)
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
void d_print_exp_tran_mat(int row, int col, double *A, int lda)
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

	int ii;

	int nx2 = nx*nx;

	int pp = nx/2; // number of masses
	
	// Ac
	for(ii=0; ii<nx*nx; ii++)
		Ac[ii] = 0.0;
	for(ii=0; ii<pp; ii++)
		Ac[ii+nx*(pp+ii)] = 1.0;
	for(ii=0; ii<pp; ii++)
		Ac[pp+ii+nx*(ii)] = -2.0;
	for(ii=0; ii<pp-1; ii++)
		Ac[pp+ii+nx*(1+ii)] = 1.0;
	for(ii=0; ii<pp; ii++)
		Ac[1+pp+ii+nx*(ii)] = 1.0;

	// Bc
	for(ii=0; ii<nx*nu; ii++)
		Bc[ii] = 0.0;
	for(ii=0; ii<nu; ii++)
		Bc[pp+ii+nx*ii] = 1.0;

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



int d_linear_system_memsize(int nx, int nu)
	{
	int size = 0;
	size += (nx*nx+nx*nu)*sizeof(double);
	return size;
	}



void d_linear_system_create(int nx, int nu, struct d_linear_system *ls, void *memory)
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



void d_linear_system_set_all(double *A,  double *B, struct d_linear_system *ls)
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
		


void d_expl_linear_ode(int t, double *x, double *u, void *ode_args, double *xdot)
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



void d_expl_linear_vde(int t, double *x, double *u, void *ode_args, double *xdot)
	{
	struct d_linear_system *ls = ode_args;
	int ii, jj, kk;
	int nx = ls->nx;
	int nu = ls->nu;
	double *Ac = ls->Ac;
	double *Bc = ls->Bc;
	double *tmp;
//	for(ii=0; ii<nx*(1+nx+nu); ii++)
	for(ii=0; ii<nx*(1+nu+nx); ii++)
		xdot[ii] = 0.0;
	for(kk=0; kk<1+nx+nu; kk++)
		for(jj=0; jj<nx; jj++)
			for(ii=0; ii<nx; ii++)
				xdot[ii+nx*kk] += Ac[ii+nx*jj] * x[jj+nx*kk];
//	tmp = xdot+nx*(nx+nu);
	tmp = xdot;
	for(jj=0; jj<nu; jj++)
		for(ii=0; ii<nx; ii++)
			tmp[ii] += Bc[ii+nx*jj] * u[jj];
//	tmp = xdot;
	tmp = xdot+nx;
	for(jj=0; jj<nu; jj++)
		for(ii=0; ii<nx; ii++)
			tmp[ii+nx*jj] += Bc[ii+nx*jj];
	return;
	}



#if 0
void d_res_impl_linear_ode(int t, double *xdot, double *x, double *u, void *ode_args, double *res)
	{
	struct d_linear_system *ls = ode_args;
	int ii, jj;
	int nx = ls->nx;
	int nu = ls->nu;
	double *Ac = ls->Ac;
	double *Bc = ls->Bc;
	for(ii=0; ii<nx; ii++)
		res[ii] = - xdot[ii];
	for(jj=0; jj<nx; jj++)
		for(ii=0; ii<nx; ii++)
			res[ii] += Ac[ii+nx*jj] * x[jj];
	for(jj=0; jj<nu; jj++)
		for(ii=0; ii<nx; ii++)
			res[ii] += Bc[ii+nx*jj] * u[jj];
	return;
	}
#endif



#if 0
void d_res_impl_linear_vde(int t, double *xdot, double *x, double *u, void *ode_args, double *res)
	{
	struct d_linear_system *ls = ode_args;
	int ii, jj, kk;
	int nx = ls->nx;
	int nu = ls->nu;
	double *Ac = ls->Ac;
	double *Bc = ls->Bc;
	double *tmp;
	for(ii=0; ii<nx*(1+nx+nu); ii++)
		res[ii] = - xdot[ii];
	for(kk=0; kk<1+nx+nu; kk++)
		for(jj=0; jj<nx; jj++)
			for(ii=0; ii<nx; ii++)
				res[ii+nx*kk] += Ac[ii+nx*jj] * x[jj+nx*kk];
	tmp = res+nx*(nx+nu);
	for(jj=0; jj<nu; jj++)
		for(ii=0; ii<nx; ii++)
			tmp[ii] += Bc[ii+nx*jj] * u[jj];
	tmp = res;
	for(jj=0; jj<nu; jj++)
		for(ii=0; ii<nx; ii++)
			tmp[ii+nx*jj] += Bc[ii+nx*jj];
	return;
	}
#endif



#if 0
void d_jac_impl_linear_ode(int t, double *xdot, double *x, double *u, void *ode_args, double *jac)
	{
	struct d_linear_system *ls = ode_args;
	int ii, jj;
	int nx = ls->nx;
	double *Ac = ls->Ac;
	for(jj=0; jj<nx; jj++)
		for(ii=0; ii<nx; ii++)
			jac[ii+nx*jj] = Ac[ii+nx*jj];
	return;
	}
#endif



int main()
	{

	int ii, jj;

/************************************************
* problem size
************************************************/	
	
	int nx = 4;
	int nu = 1;

/************************************************
* (continuous time) mass sprint system
************************************************/	
	
	double *Ac; d_zeros(&Ac, nx, nx);
	double *Bc; d_zeros(&Bc, nx, nu);

	mass_spring_system(nx, nu, Ac, Bc);

	d_print_mat(nx, nx, Ac, nx);
	d_print_mat(nx, nu, Bc, nx);

	int ls_memsize = d_linear_system_memsize(nx, nu);
	printf("\nls memsize = %d\n", ls_memsize);
	void *ls_memory = malloc(ls_memsize);

	struct d_linear_system ls;
	d_linear_system_create(nx, nu, &ls, ls_memory);

	d_linear_system_set_all(Ac, Bc, &ls);

	d_print_mat(nx, nx, ls.Ac, nx);
	d_print_mat(nx, nu, ls.Bc, nx);

	double *x0; d_zeros(&x0, nx, 1);
	x0[0] = 2.5;
	x0[1] = 2.5;

	double *u; d_zeros(&u, nu, 1);
	u[0] = 0.0;

	double *xdot; d_zeros(&xdot, nx, 1);
	d_expl_linear_ode(0, x0, u, &ls, xdot);

	d_print_mat(1, nx, x0, 1);
	d_print_mat(1, nx, xdot, 1);

/************************************************
* (discrete time) mass sprint system
************************************************/	

	double Ts = 0.5;

	double *A = malloc(nx*nx*sizeof(double));
	double *B = malloc(nx*nu*sizeof(double));
	double *T = malloc(nx*nx*sizeof(double));
	int *ipiv = malloc(nx*sizeof(int));

	for(ii=0; ii<nx*nx; ii++)
		A[ii] = Ts*Ac[ii];
	expm(nx, A);

	for(ii=0; ii<nx*nx; ii++) T[ii] = A[ii];
	for(ii=0; ii<nx; ii++) T[ii*(nx+1)] -= 1.0;
	dgemm_nn_3l(nx, nu, nx, T, nx, Bc, nx, B, nx);

	int info = 0;
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

	int steps = 30;
	double h = Ts/steps;

	printf("\nTs=%e, steps=%d, h=%e\n", Ts, steps, h);

#if 1
	// rk4
	int ns = 4; // number of stages
	char rk_method[] = "ERK4";
//	double A_rk[] = {0.0, 0.0, 0.0, 0.0,
//	                 0.5, 0.0, 0.0, 0.0,
//	                 0.0, 0.5, 0.0, 0.0,
//	                 0.0, 0.0, 1.0, 0.0};
//	double B_rk[] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
//	double C_rk[] = {0.0, 0.5, 0.5, 0.0};
#elif 0
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
//	printf("\n%d %s\n", ns, rk_method);

	int size_rk_data = d_sim_rk_data_memsize(ns);
	printf("\nmemsize rk data %d\n", size_rk_data);
	void *mem_rk_data = malloc(size_rk_data);

	struct d_sim_rk_data rk_data;
	d_sim_rk_data_create(ns, &rk_data, mem_rk_data);

	d_sim_rk_data_init_default(rk_method, &rk_data);

//	d_print_mat(ns, ns, rk_data.A_rk, ns);

	
	int size_erk_arg = d_sim_erk_arg_memsize();
	printf("\nmemsize erk arg %d\n", size_erk_arg);
	void *mem_erk_arg = malloc(size_erk_arg);

	struct d_sim_erk_arg erk_arg;
	d_sim_erk_arg_create(&erk_arg, mem_erk_arg);

	d_sim_erk_arg_set_all(&rk_data, h, steps, &erk_arg);


	int nf_max = nu+nx;
	int na_max = 0;
	int size_erk_ws = d_sim_erk_ws_memsize(&erk_arg, nx, nu, nf_max, na_max);
	printf("\nmemsize erk ws %d\n", size_erk_ws);
	void *mem_erk_ws = malloc(size_erk_ws);

	struct d_sim_erk_ws erk_ws;
	d_sim_erk_ws_create(&erk_arg, nx, nu, nf_max, na_max, &erk_ws, mem_erk_ws);

	int nf = nu+nx;
	int na = 0;

	double *fs0; d_zeros(&fs0, nx*nf, 1);
	for(ii=0; ii<nx*nf; ii++)
		fs0[ii] = 0.0;
	for(ii=0; ii<nx; ii++)
		fs0[nu*nx+ii*(nx+1)] = 1.0;

	printf("\nforward seeds\n");
	d_print_mat(nx, nu, fs0, nx);
	d_print_mat(nx, nx, fs0+nu*nx, nx);

	d_sim_erk_ws_set_all(nf, na, x0, fs0, NULL, u, &d_expl_linear_ode, &d_expl_linear_vde, NULL, &ls, &erk_ws);

	d_sim_erk_solve(&erk_arg, &erk_ws);

	double *x_end; d_zeros(&x_end, nx, 1);

	d_sim_erk_ws_get_x(&erk_ws, x_end);

	printf("\nx erk %s\n", rk_method);
	d_print_mat(1, nx, x_end, 1);
	if(nf!=0)
		{
		d_print_mat(nx, nu, erk_ws.x_for+nx, nx);
		d_print_mat(nx, nx, erk_ws.x_for+nx+nu*nx, nx);

		}

#if 0
/************************************************
* explicit integrator
************************************************/	

#if 0
	int nf = nx+nu;
	int np = nu;

	int memsize_erk_int = d_memsize_erk_int(&rk_data, nx, nf, np);
	printf("\nmemsize erk int %d\n", memsize_erk_int);
	void *memory_erk = malloc(memsize_erk_int);

	struct d_erk_workspace erk_workspace;
	d_create_erk_int(&rk_data, nx, nf, np, &erk_workspace, memory_erk);

	struct d_erk_args erk_args;
	erk_args.steps = steps;
	erk_args.h = h;
	
	double *fs0; d_zeros(&fs0, nx*nf, 1);
	for(ii=0; ii<nx; ii++) fs0[nu*nx+ii*(nx+1)] = 1.0;

	double *ex_erk; d_zeros(&ex_erk, nx*(1+nx+nu), 1);

	d_init_erk_int(x0, fs0, u, &d_expl_linear_vde, &ls, &erk_workspace);

	d_update_p_erk_int(u, &erk_workspace);

	d_erk_int(&erk_args, &erk_workspace);

	double *x_erk = erk_workspace.x;

	for(ii=0; ii<nx; ii++)
		ex_erk[ii] = x_erk[nx*nf+ii] - xref[ii];
	for(ii=0; ii<nx*nu; ii++)
		ex_erk[nx+ii] = x_erk[ii] - B[ii];
	for(ii=0; ii<nx*nx; ii++)
		ex_erk[nx+nx*nu+ii] = x_erk[nx*nu+ii] - A[ii];

	struct blasfeo_dmat sBAbt;
	blasfeo_allocate_dmat(nx+nu+1, nx, &sBAbt);
	struct blasfeo_dvec sb;
	blasfeo_allocate_dvec(nx, &sb);

	int nxx[2] = {nx, nx};
	int nuu[1] = {nu};

	struct d_ocp_qp qp;
	qp.BAbt = &sBAbt;
	qp.b = &sb;
	qp.nx = nxx;
	qp.nu = nuu;

	d_print_mat(1, nx, xref, 1);
	d_cvt_erk_int_to_ocp_qp(0, &erk_workspace, xref, &qp);
	d_print_strmat(nx+nu+1, nx, &sBAbt, 0, 0);

	printf("\nx erk\n");
	d_print_mat(nx, nu+nx+1, x_erk, nx); // B A x_next
	d_print_mat(1, nx, x_erk+nx*nf, 1); // x_next
	d_print_mat(nx, nx, x_erk+nx*nu, nx); // A
	d_print_mat(nx, nu, x_erk, nx); // B
	printf("\nerror erk\n");
	d_print_exp_mat(1, nx, ex_erk, 1);
	d_print_exp_mat(nx, nx, ex_erk+nx+nx*nu, nx);
	d_print_exp_mat(nx, nu, ex_erk+nx, nx);
#endif

/************************************************
* implicit integrator
************************************************/	
	
	int memsize_irk_int = d_memsize_irk_int(&rk_data, nx, nf, np);
	printf("\nmemsize irk int %d\n", memsize_irk_int);
	void *memory_irk = malloc(memsize_irk_int);

	struct d_irk_workspace irk_workspace;
	d_create_irk_int(&rk_data, nx, nf, np, &irk_workspace, memory_irk);

	struct d_irk_args irk_args;
	irk_args.steps = steps;
	irk_args.h = h;
	irk_args.newton_iter = 2;

	d_init_irk_int(x0, fs0, u, &d_res_impl_linear_vde, &d_jac_impl_linear_ode, &ls, &irk_workspace);

//	d_update_p_irk_int(u, &irk_workspace);

	d_irk_int(&irk_args, &irk_workspace);

//	d_print_strmat(ns*nx, ns*nx, irk_workspace.JG, 0, 0);
//	d_print_strmat(ns*nx, nu+nx+1, irk_workspace.rG, 0, 0);
//	d_print_strmat(ns*nx, nu+nx+1, irk_workspace.K, 0, 0);

	d_print_mat(nx, nf+1, irk_workspace.x, nx);

	double *x_irk = irk_workspace.x;
	double *ex_irk; d_zeros(&ex_irk, nx*(1+nx+nu), 1);

	for(ii=0; ii<nx; ii++)
		ex_irk[ii] = x_irk[nx*nf+ii] - xref[ii];
	for(ii=0; ii<nx*nu; ii++)
		ex_irk[nx+ii] = x_irk[ii] - B[ii];
	for(ii=0; ii<nx*nx; ii++)
		ex_irk[nx+nx*nu+ii] = x_irk[nx*nu+ii] - A[ii];

	printf("\nx irk\n");
	d_print_mat(nx, nu+nx+1, x_irk, nx); // B A x_next
	d_print_mat(1, nx, x_irk+nx*nf, 1); // x_next
	d_print_mat(nx, nx, x_irk+nx*nu, nx); // A
	d_print_mat(nx, nu, x_irk, nx); // B
	printf("\nerror irk\n");
	d_print_exp_mat(1, nx, ex_irk, 1);
	d_print_exp_mat(nx, nx, ex_irk+nx+nx*nu, nx);
	d_print_exp_mat(nx, nu, ex_irk+nx, nx);

//	blasfeo_dgese(nu+nx+1, nx, 0.0, &sBAbt, 0, 0);
//	d_cvt_irk_int_to_ocp_qp(0, &irk_workspace, xref, &qp);
//	d_print_strmat(nx+nu+1, nx, &sBAbt, 0, 0);
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
	free(mem_rk_data);
	free(mem_erk_arg);
	free(mem_erk_ws);
	free(fs0);
	free(x_end);
//	free(memory_erk);
//	free(ex_erk);
//	free(memory_irk);
//	free(ex_irk);

//	blasfeo_free_dmat(&sBAbt);
//	blasfeo_free_dvec(&sb);

	return 0;

	}


