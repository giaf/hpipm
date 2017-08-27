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
#include <blasfeo_d_blas.h>

#include "../include/hpipm_d_rk_int.h"
#include "../include/hpipm_d_erk_int.h"
#include "../include/hpipm_d_ocp_qp.h"
#include "../include/hpipm_d_ocp_qp_sol.h"
#include "../include/hpipm_d_ocp_qp_ipm.h"
#include "../include/hpipm_d_ocp_qp_sim.h"

#include "d_tools.h"



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



void d_linear_vde0(int t, double *x, double *u, void *ode_args, double *xdot)
	{
	struct d_linear_system *ls = ode_args;
	int ii, jj, kk;
	int nx = ls->nx;
	int nu = ls->nu;
	double *Ac = ls->Ac;
	double *Bc = ls->Bc;
	double *tmp;
	for(ii=0; ii<nx*(nu+1); ii++)
		xdot[ii] = 0.0;
	for(kk=0; kk<nu+1; kk++)
		for(jj=0; jj<nx; jj++)
			for(ii=0; ii<nx; ii++)
				xdot[ii+nx*kk] += Ac[ii+nx*jj] * x[jj+nx*kk];
	tmp = xdot+nx*(nu);
	for(jj=0; jj<nu; jj++)
		for(ii=0; ii<nx; ii++)
			tmp[ii] += Bc[ii+nx*jj] * u[jj];
	tmp = xdot;
	for(jj=0; jj<nu; jj++)
		for(ii=0; ii<nx; ii++)
			tmp[ii+nx*jj] += Bc[ii+nx*jj];
	return;
	}



void d_linear_vde1(int t, double *x, double *u, void *ode_args, double *xdot)
	{
	struct d_linear_system *ls = ode_args;
	int ii, jj, kk;
	int nx = ls->nx;
	int nu = ls->nu;
	double *Ac = ls->Ac;
	double *Bc = ls->Bc;
	double *tmp;
	for(ii=0; ii<nx*(nu+nx+1); ii++)
		xdot[ii] = 0.0;
	for(kk=0; kk<nu+nx+1; kk++)
		for(jj=0; jj<nx; jj++)
			for(ii=0; ii<nx; ii++)
				xdot[ii+nx*kk] += Ac[ii+nx*jj] * x[jj+nx*kk];
	tmp = xdot+nx*(nu+nx);
	for(jj=0; jj<nu; jj++)
		for(ii=0; ii<nx; ii++)
			tmp[ii] += Bc[ii+nx*jj] * u[jj];
	tmp = xdot;
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
	
	int nx_ = 8;
	int nu_ = 3;

/************************************************
* (continuous time) mass sprint system
************************************************/	
	
	double *Ac; d_zeros(&Ac, nx_, nx_);
	double *Bc; d_zeros(&Bc, nx_, nu_);

	mass_spring_system(nx_, nu_, Ac, Bc);

	d_print_mat(nx_, nx_, Ac, nx_);
	d_print_mat(nx_, nu_, Bc, nx_);

	int ls_memsize = d_memsize_linear_system(nx_, nu_);
	printf("\nls memsize = %d\n", ls_memsize);
	void *ls_memory = malloc(ls_memsize);

	struct d_linear_system ls;
	d_create_linear_system(nx_, nu_, &ls, ls_memory);

	d_cvt_colmaj_to_linear_system(Ac, Bc, &ls);

	d_print_mat(nx_, nx_, ls.Ac, nx_);
	d_print_mat(nx_, nu_, ls.Bc, nx_);

	double *x0; d_zeros(&x0, nx_, 1);
	x0[0] = 2.5;
	x0[1] = 2.5;

	d_print_mat(1, nx_, x0, 1);

/************************************************
* (discrete time) mass sprint system
************************************************/	

	double Ts = 0.5;

	double *A = malloc(nx_*nx_*sizeof(double));
	double *B = malloc(nx_*nu_*sizeof(double));
	double *T = malloc(nx_*nx_*sizeof(double));
	int *ipiv = malloc(nx_*sizeof(int));

	for(ii=0; ii<nx_*nx_; ii++)
		A[ii] = Ts*Ac[ii];
	expm(nx_, A);

	for(ii=0; ii<nx_*nx_; ii++) T[ii] = A[ii];
	for(ii=0; ii<nx_; ii++) T[ii*(nx_+1)] -= 1.0;
	dgemm_nn_3l(nx_, nu_, nx_, T, nx_, Bc, nx_, B, nx_);

	int info = 0;
	dgesv_3l(nx_, nu_, Ac, nx_, ipiv, B, nx_, &info);

	double *b; d_zeros(&b, nx_, 1);

	d_print_mat(nx_, nx_, A, nx_);
	d_print_mat(nx_, nu_, B, nx_);
	d_print_mat(1, nx_, b, 1);

/************************************************
* quadratic cost function
************************************************/	
	
	double *Q; d_zeros(&Q, nx_, nx_);
	for(ii=0; ii<nx_; ii++) Q[ii*(nx_+1)] = 1.0;

	double *R; d_zeros(&R, nu_, nu_);
	for(ii=0; ii<nu_; ii++) R[ii*(nu_+1)] = 2.0;

	double *S; d_zeros(&S, nu_, nx_);

	double *q; d_zeros(&q, nx_, 1);
	for(ii=0; ii<nx_; ii++) q[ii] = 0.0;

	double *r; d_zeros(&r, nu_, 1);
	for(ii=0; ii<nu_; ii++) r[ii] = 0.0;

	d_print_mat(nx_, nx_, Q, nx_);
	d_print_mat(nu_, nu_, R, nu_);
	d_print_mat(nu_, nx_, S, nu_);
	d_print_mat(1, nx_, q, 1);
	d_print_mat(1, nu_, r, 1);

/************************************************
* integrator type and args
************************************************/	
	
#if 1
	// rk4
	int nsta = 4; // number of stages
	double A_rk[] = {0.0, 0.0, 0.0, 0.0,
	                 0.5, 0.0, 0.0, 0.0,
	                 0.0, 0.5, 0.0, 0.0,
	                 0.0, 0.0, 1.0, 0.0};
	double B_rk[] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
	double C_rk[] = {0.0, 0.5, 0.5, 0.0};
#elif 1
	// midpoint rule
	int nsta = 2; // number of stages
	double A_rk[] = {0.0, 0.0,
	                 0.5, 0.0};
	double B_rk[] = {0.0, 1.0};
	double C_rk[] = {0.0, 0.5};
#else
	// explicit euler
	int nsta = 1; // number of stages
	double A_rk[] = {0.0};
	double B_rk[] = {1.0};
	double C_rk[] = {0.0};
#endif

	// erk data structure
	int memsize_rk_data = d_memsize_rk_data(nsta);
	printf("\nmemsize rk data %d\n", memsize_rk_data);
	void *memory_rk_data = malloc(memsize_rk_data);

	struct d_rk_data rk_data;
	d_create_rk_data(nsta, &rk_data, memory_rk_data);

	d_cvt_rowmaj_to_rk_data(A_rk, B_rk, C_rk, &rk_data);

	// erk args structure
	struct d_erk_args erk_args;
	erk_args.steps = 30;
	erk_args.h = Ts/erk_args.steps;

/************************************************
* integrator workspace
************************************************/	
	
	// first stage

	int nf0 = nu_; // number of forward sensitivities
	int np0 = nu_; // numner of parameters

	// erk workspace structure
	int memsize_erk_int0 = d_memsize_erk_int(&rk_data, nx_, nf0, np0);
	printf("\nmemsize erk int0 %d\n", memsize_erk_int0);
	void *memory_erk0 = malloc(memsize_erk_int0);

	struct d_erk_workspace erk_workspace0;
	d_create_erk_int(&rk_data, nx_, nf0, np0, &erk_workspace0, memory_erk0);

	// forward sensitivities seeds
	double *fs0; d_zeros(&fs0, nx_*nf0, 1);

	// first stage

	int nf1 = nx_+nu_; // number of forward sensitivities
	int np1 = nu_; // numner of parameters

	// erk workspace structure
	int memsize_erk_int1 = d_memsize_erk_int(&rk_data, nx_, nf1, np1);
	printf("\nmemsize erk int1 %d\n", memsize_erk_int1);
	void *memory_erk1 = malloc(memsize_erk_int1);

	struct d_erk_workspace erk_workspace1;
	d_create_erk_int(&rk_data, nx_, nf1, np1, &erk_workspace1, memory_erk1);

	// forward sensitivities seeds
	double *fs1; d_zeros(&fs1, nx_*nf1, 1);
	for(ii=0; ii<nx_; ii++) fs1[nu_*nx_+ii*(nx_+1)] = 1.0;

/************************************************
* ocp qp
************************************************/	

	int N = 4;

	int nx[N+1];
	int nu[N+1];
	int nb[N+1];
	int ng[N+1];
	int ns[N+1];

	nx[0] = 0;
	nu[0] = nu_;
	nb[0] = 0;
	ng[0] = 0;
	ns[0] = 0;
	for(ii=1; ii<N; ii++)
		{
		nx[ii] = nx_;
		nu[ii] = nu_;
		nb[ii] = 0;
		ng[ii] = 0;
		ns[ii] = 0;
		}
	nx[N] = nx_;
	nu[N] = 0;
	nb[N] = 0;
	ng[N] = 0;
	ns[N] = 0;

	int qp_size = d_memsize_ocp_qp(N, nx, nu, nb, ng, ns);
	printf("\nqp size = %d\n", qp_size);
	void *qp_mem = malloc(qp_size);

	struct d_ocp_qp qp;
	d_create_ocp_qp(N, nx, nu, nb, ng, ns, &qp, qp_mem);

	// copy problem size
	qp.N = N;
	for(ii=0; ii<=N; ii++)
		{
		qp.nx[ii] = nx[ii];
		qp.nu[ii] = nu[ii];
		qp.nb[ii] = nb[ii];
		qp.ng[ii] = ng[ii];
		qp.ns[ii] = ns[ii];
		}

	// copy Hessian and gradient
	for(ii=0; ii<=N; ii++)
		{
		d_cvt_mat2strmat(nu[ii], nu[ii], R, nu_, qp.RSQrq+ii, 0, 0);
		d_cvt_tran_mat2strmat(nu[ii], nx[ii], S, nu_, qp.RSQrq+ii, nu[ii], 0);
		d_cvt_mat2strmat(nx[ii], nx[ii], Q, nx_, qp.RSQrq+ii, nu[ii], nu[ii]);
		d_cvt_tran_mat2strmat(nu[ii], 1, r, nu_, qp.RSQrq+ii, nu[ii]+nx[ii], 0);
		d_cvt_tran_mat2strmat(nx[ii], 1, q, nu_, qp.RSQrq+ii, nu[ii]+nx[ii], nu[ii]);
		d_cvt_vec2strvec(nu[ii], r, qp.rq+ii, 0);
		d_cvt_vec2strvec(nx[ii], q, qp.rq+ii, nu[ii]);
		}
	
//	for(ii=0; ii<=N; ii++)
//		d_print_strmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], qp.RSQrq+ii, 0, 0);

/************************************************
* ocp qp sol
************************************************/	
	
	int qp_sol_size = d_memsize_ocp_qp_sol(N, nx, nu, nb, ng, ns);
	printf("\nqp sol size = %d\n", qp_sol_size);
	void *qp_sol_mem = malloc(qp_sol_size);

	struct d_ocp_qp_sol qp_sol;
	d_create_ocp_qp_sol(N, nx, nu, nb, ng, ns, &qp_sol, qp_sol_mem);

/************************************************
* ipm
************************************************/	

	struct d_ipm_ocp_qp_arg arg;
	arg.alpha_min = 1e-8;
	arg.mu_max = 1e-12;
	arg.iter_max = 20;
	arg.mu0 = 2.0;

	int ipm_size = d_memsize_ipm_ocp_qp(&qp, &arg);
	printf("\nipm size = %d\n", ipm_size);
	void *ipm_mem = malloc(ipm_size);

	struct d_ipm_ocp_qp_workspace workspace;
	d_create_ipm_ocp_qp(&qp, &arg, &workspace, ipm_mem);

/************************************************
* SQP loop
************************************************/	
	
	int ss, nn;

	struct d_strmat BAbt0;
	d_allocate_strmat(nu_+nx_+1, nx_, &BAbt0);

	// solutoni
	double *u[N+1]; for(ii=0; ii<=N; ii++) d_zeros(u+ii, nu[ii], 1);
	double *x[N+1]; for(ii=0; ii<=N; ii++) d_zeros(x+ii, nx[ii], 1);
	
	// step
	double *du[N+1]; for(ii=0; ii<=N; ii++) d_zeros(du+ii, nu[ii], 1);
	double *dx[N+1]; for(ii=0; ii<=N; ii++) d_zeros(dx+ii, nx[ii], 1);
	double **dls;
	double **dus;
	double *dpi[N]; for(ii=0; ii<N; ii++) d_zeros(dpi+ii, nx[ii+1], 1);
	double *dlam_lb[N+1]; for(ii=0; ii<=N; ii++) d_zeros(dlam_lb+ii, nb[ii], 1);
	double *dlam_ub[N+1]; for(ii=0; ii<=N; ii++) d_zeros(dlam_ub+ii, nb[ii], 1);
	double *dlam_lg[N+1]; for(ii=0; ii<=N; ii++) d_zeros(dlam_lg+ii, ng[ii], 1);
	double *dlam_ug[N+1]; for(ii=0; ii<=N; ii++) d_zeros(dlam_ug+ii, ng[ii], 1);
	double **dlam_ls;
	double **dlam_us;

	// initialize solution to zero
	for(nn=0; nn<=N; nn++)
		for(ii=0; ii<nu[nn]; ii++)
			u[nn][ii] = 0.0;
	for(nn=0; nn<=N; nn++)
		for(ii=0; ii<nx[nn]; ii++)
			x[nn][ii] = 0.0;

	int sqp_steps = 1;
	for(ss=0; ss<sqp_steps; ss++)	
		{

		// initial stage
		// XXX x0 in the QP is zero since x0 in the nlp is initialized to x0 !!!
		nn = 0;
		// XXX it does not need the sensitivities wrt x here
		d_init_erk_int(x0, fs0, u[nn], &d_linear_vde0, &ls, &erk_workspace0);
		d_erk_int(&erk_args, &erk_workspace0);
		d_cvt_erk_int_to_ocp_qp(nn, &erk_workspace0, x[nn+1], &qp);

		// other stages
		for(nn=1; nn<N; nn++)
			{
			d_init_erk_int(x[nn], fs1, u[nn], &d_linear_vde1, &ls, &erk_workspace1);
			d_erk_int(&erk_args, &erk_workspace1);
			d_cvt_erk_int_to_ocp_qp(nn, &erk_workspace1, x[nn+1], &qp);
			}

		for(nn=0; nn<N; nn++)
			d_print_strmat(nu[nn]+nx[nn]+1, nx[nn+1], qp.BAbt+nn, 0, 0);

		d_solve_ipm2_ocp_qp(&qp, &qp_sol, &workspace);

		d_cvt_ocp_qp_sol_to_colmaj(&qp, &qp_sol, du, dx, dls, dus, dpi, dlam_lb, dlam_ub, dlam_lg, dlam_ug, dlam_ls, dlam_us);

		for(nn=0; nn<=N; nn++)
			for(ii=0; ii<nu[nn]; ii++)
				u[nn][ii] += du[nn][ii];
		for(nn=0; nn<=N; nn++)
			for(ii=0; ii<nx[nn]; ii++)
				x[nn][ii] += dx[nn][ii];

		}
	
	printf("\nu = \n");
	for(nn=0; nn<=N; nn++)
		d_print_mat(1, nu[nn], u[nn], 1);
	printf("\nx = \n");
	for(nn=0; nn<=N; nn++)
		d_print_mat(1, nx[nn], x[nn], 1);

/************************************************
* free memory
************************************************/	
	
	free(ls_memory);
	free(Ac);
	free(Bc);
	free(x0);
	free(A);
	free(B);
	free(T);
	free(ipiv);
	free(b);
	free(Q);
	free(S);
	free(R);
	free(q);
	free(r);
	free(memory_rk_data);
	free(memory_erk0);
	free(memory_erk1);
	free(qp_mem);
	free(qp_sol_mem);
	free(ipm_mem);

	for(ii=0; ii<N; ii++)
		{
		d_free(u[ii]);
		d_free(x[ii]);
		d_free(du[ii]);
		d_free(dx[ii]);
		d_free(dpi[ii]);
		d_free(dlam_lb[ii]);
		d_free(dlam_ub[ii]);
		d_free(dlam_lg[ii]);
		d_free(dlam_ug[ii]);
		}
	d_free(u[ii]);
	d_free(x[ii]);
	d_free(du[ii]);
	d_free(dx[ii]);
	d_free(dlam_lb[ii]);
	d_free(dlam_ub[ii]);
	d_free(dlam_lg[ii]);
	d_free(dlam_ug[ii]);

	d_free_strmat(&BAbt0);

	return 0;

	}


