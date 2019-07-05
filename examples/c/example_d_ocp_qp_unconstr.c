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
#include <sys/time.h>

#include <blasfeo_d_aux_ext_dep.h>

#include "../../include/hpipm_d_ocp_qp_ipm.h"
#include "../../include/hpipm_d_ocp_qp_dim.h"
#include "../../include/hpipm_d_ocp_qp.h"
#include "../../include/hpipm_d_ocp_qp_sol.h"



// qp data as global data
extern int N;
extern int *nx;
extern int *nu;
extern int *nbu;
extern int *nbx;
extern int *ng;
extern int *nsbx;
extern int *nsbu;
extern int *nsg;
extern double **hA;
extern double **hB;
extern double **hb;
extern double **hQ;
extern double **hR;
extern double **hS;
extern double **hq;
extern double **hr;
extern int **hidxbx;
extern double **hlbx;
extern double **hubx;
extern int **hidxbu;
extern double **hlbu;
extern double **hubu;
extern double **hC;
extern double **hD;
extern double **hlg;
extern double **hug;
extern double **hZl;
extern double **hZu;
extern double **hzl;
extern double **hzu;
extern int **hidxs;
extern double **hlls;
extern double **hlus;
//extern double **hu_guess;
//extern double **hx_guess;
//extern double **hsl_guess;
//extern double **hsu_guess;



// main
int main()
	{

	int ii, jj;

	int hpipm_return;

	int rep, nrep=10;

	struct timeval tv0, tv1;

/************************************************
* new first stage data
************************************************/

	double *x0 = malloc(nx[0]*sizeof(double));

	for(ii=0; ii<nx[0]; ii++)
		x0[ii] = hlbx[0][ii];

	d_print_mat(1, nx[0], x0, 1);

	int *idxbx0 = malloc(nx[0]*sizeof(int));

	for(ii=0; ii<nx[0]; ii++)
		idxbx0[ii] = ii;

	double *b0 = malloc(nx[1]*sizeof(double));

	for(ii=0; ii<nx[1]; ii++)
		b0[ii] = hb[0][ii];
	
	for(ii=0; ii<nx[1]; ii++)
		for(jj=0; jj<nx[0]; jj++)
			b0[ii] += hA[0][ii+nx[1]*jj]*x0[jj];
	
	double *r0 = malloc(nu[0]*sizeof(double));

	for(ii=0; ii<nu[0]; ii++)
		r0[ii] = hr[0][ii];
	
	for(ii=0; ii<nu[0]; ii++)
		for(jj=0; jj<nx[0]; jj++)
			r0[ii] += hS[0][ii+nu[0]*jj]*x0[jj];

/************************************************
* remove constraints
************************************************/

#if 0
	// keep x0 as an optimization variable, set its value using equal upper-lower bounds
	for(ii=1; ii<=N; ii++)
		nbx[ii] = 0;

	hlb[0] = x0;
	hub[0] = x0;
	hidxbx[0] = idxbx0;
#else
	// remove x0 from the optimization variables
	nx[0] = 0;

	for(ii=0; ii<=N; ii++)
		nbx[ii] = 0;

	hb[0] = b0;
	hr[0] = r0;
#endif

	for(ii=0; ii<=N; ii++)
		nbu[ii] = 0;

	for(ii=0; ii<=N; ii++)
		ng[ii] = 0;

	for(ii=0; ii<=N; ii++)
		nsbx[ii] = 0;

	for(ii=0; ii<=N; ii++)
		nsbu[ii] = 0;

	for(ii=0; ii<=N; ii++)
		nsg[ii] = 0;

/************************************************
* ocp qp dim
************************************************/

	int dim_size = d_ocp_qp_dim_memsize(N);
	void *dim_mem = malloc(dim_size);

	struct d_ocp_qp_dim dim;
	d_ocp_qp_dim_create(N, &dim, dim_mem);

	d_ocp_qp_dim_set_all(nx, nu, nbx, nbu, ng, nsbx, nsbu, nsg, &dim);

/************************************************
* ocp qp
************************************************/

	int qp_size = d_ocp_qp_memsize(&dim);
	void *qp_mem = malloc(qp_size);

	struct d_ocp_qp qp;
	d_ocp_qp_create(&dim, &qp, qp_mem);

	d_ocp_qp_set_all(hA, hB, hb, hQ, hS, hR, hq, hr, hidxbx, hlbx, hubx, hidxbu, hlbu, hubu, hC, hD, hlg, hug, hZl, hZu, hzl, hzu, hidxs, hlls, hlus, &qp);

/************************************************
* ocp qp sol
************************************************/

	int qp_sol_size = d_ocp_qp_sol_memsize(&dim);
	void *qp_sol_mem = malloc(qp_sol_size);

	struct d_ocp_qp_sol qp_sol;
	d_ocp_qp_sol_create(&dim, &qp_sol, qp_sol_mem);

/************************************************
* ipm arg
************************************************/

	int ipm_arg_size = d_ocp_qp_ipm_arg_memsize(&dim);
	void *ipm_arg_mem = malloc(ipm_arg_size);

	struct d_ocp_qp_ipm_arg arg;
	d_ocp_qp_ipm_arg_create(&dim, &arg, ipm_arg_mem);

//	enum hpipm_mode mode = SPEED_ABS;
	enum hpipm_mode mode = SPEED;
//	enum hpipm_mode mode = BALANCE;
//	enum hpipm_mode mode = ROBUST;
	d_ocp_qp_ipm_arg_set_default(mode, &arg);

	double mu0 = 1e4;
	int iter_max = 30;
	double tol_stat = 1e-4;
	double tol_eq = 1e-5;
	double tol_ineq = 1e-5;
	double tol_comp = 1e-5;
	double reg_prim = 1e-12;
	int warm_start = 0;
	int ric_alg = 0;

	d_ocp_qp_ipm_arg_set_mu0(&mu0, &arg);
	d_ocp_qp_ipm_arg_set_iter_max(&iter_max, &arg);
	d_ocp_qp_ipm_arg_set_tol_stat(&tol_stat, &arg);
	d_ocp_qp_ipm_arg_set_tol_eq(&tol_eq, &arg);
	d_ocp_qp_ipm_arg_set_tol_ineq(&tol_ineq, &arg);
	d_ocp_qp_ipm_arg_set_tol_comp(&tol_comp, &arg);
	d_ocp_qp_ipm_arg_set_reg_prim(&reg_prim, &arg);
	d_ocp_qp_ipm_arg_set_warm_start(&warm_start, &arg);
//	d_ocp_qp_ipm_arg_set_ric_alg(&ric_alg, &arg);

/************************************************
* ipm workspace
************************************************/

	int ipm_size = d_memsize_ocp_qp_ipm(&dim, &arg);
	void *ipm_mem = malloc(ipm_size);

	struct d_ocp_qp_ipm_workspace workspace;
	d_create_ocp_qp_ipm(&dim, &arg, &workspace, ipm_mem);

/************************************************
* ipm solver
************************************************/

	gettimeofday(&tv0, NULL); // start

	for(rep=0; rep<nrep; rep++)
		{
		// solution guess
//		for(ii=0; ii<=N; ii++)
//			d_cvt_colmaj_to_ocp_qp_sol_u(ii, hu_guess[ii], &qp_sol);
//		for(ii=0; ii<=N; ii++)
//			d_cvt_colmaj_to_ocp_qp_sol_x(ii, hx_guess[ii], &qp_sol);
//		for(ii=0; ii<=N; ii++)
//			d_cvt_colmaj_to_ocp_qp_sol_sl(ii, hsl_guess[ii], &qp_sol);
//		for(ii=0; ii<=N; ii++)
//			d_cvt_colmaj_to_ocp_qp_sol_su(ii, hsu_guess[ii], &qp_sol);

		// call solver
		hpipm_return = d_solve_ocp_qp_ipm(&qp, &qp_sol, &arg, &workspace);
		}

	gettimeofday(&tv1, NULL); // stop

	double time_ipm = (tv1.tv_sec-tv0.tv_sec)/(nrep+0.0)+(tv1.tv_usec-tv0.tv_usec)/(nrep*1e6);

/************************************************
* print solution info
************************************************/

    printf("\nHPIPM returned with flag %i.\n", hpipm_return);
    if(hpipm_return == 0)
		{
        printf("\n -> QP solved!\n");
		}
	else if(hpipm_return==1)
		{
        printf("\n -> Solver failed! Maximum number of iterations reached\n");
		}
	else if(hpipm_return==2)
		{
        printf("\n -> Solver failed! Minimum step lenght reached\n");
		}
	else if(hpipm_return==2)
		{
        printf("\n -> Solver failed! NaN in computations\n");
		}
	else
		{
        printf("\n -> Solver failed! Unknown return flag\n");
		}
    printf("\nAverage solution time over %i runs: %e [s]\n", nrep, time_ipm);
	printf("\n\n");

/************************************************
* extract and print solution
************************************************/

	// u

	int nu_max = nu[0];
	for(ii=1; ii<=N; ii++)
		if(nu[ii]>nu_max)
			nu_max = nu[ii];

	double *u = malloc(nu_max*sizeof(double));

	printf("\nu = \n");
	for(ii=0; ii<=N; ii++)
		{
		d_ocp_qp_sol_get("u", ii, &qp_sol, u);
		d_print_mat(1, nu[ii], u, 1);
		}

	// x

	int nx_max = nx[0];
	for(ii=1; ii<=N; ii++)
		if(nx[ii]>nx_max)
			nx_max = nx[ii];

	double *x = malloc(nx_max*sizeof(double));

	printf("\nx = \n");
	for(ii=0; ii<=N; ii++)
		{
		d_ocp_qp_sol_get("x", ii, &qp_sol, x);
		d_print_mat(1, nx[ii], x, 1);
		}

	// pi
	double *pi = malloc(nx_max*sizeof(double));

	printf("\npi = \n");
	for(ii=0; ii<N; ii++)
		{
		d_ocp_qp_sol_get("pi", ii, &qp_sol, pi);
		d_print_mat(1, nx[ii+1], pi, 1);
		}

/************************************************
* print ipm statistics
************************************************/

	printf("\nipm return = %d\n", hpipm_return);
	double res_stat = d_get_ocp_qp_ipm_res_stat(&workspace);
	double res_eq = d_get_ocp_qp_ipm_res_eq(&workspace);
	double res_ineq = d_get_ocp_qp_ipm_res_ineq(&workspace);
	double res_comp = d_get_ocp_qp_ipm_res_comp(&workspace);
	printf("\nipm residuals max: res_g = %e, res_b = %e, res_d = %e, res_m = %e\n", res_stat, res_eq, res_ineq, res_comp);

	int iter = d_get_ocp_qp_ipm_iter(&workspace);
	printf("\nipm iter = %d\n", iter);
	double *stat = d_get_ocp_qp_ipm_stat(&workspace);
	printf("\nalpha_aff\tmu_aff\t\tsigma\t\talpha\t\tmu\n");
	d_print_exp_tran_mat(5, iter, stat, 5);

	printf("\nocp ipm time = %e [s]\n\n", time_ipm);

/************************************************
* free memory and return
************************************************/

	free(x0);
	free(idxbx0);
	free(b0);
	free(r0);

    free(dim_mem);
    free(qp_mem);
	free(qp_sol_mem);
	free(ipm_arg_mem);
	free(ipm_mem);

	free(u);
	free(x);
	free(pi);

	return 0;

	}



