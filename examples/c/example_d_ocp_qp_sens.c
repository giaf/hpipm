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
#include <sys/time.h>

#include <blasfeo_common.h>
#include <blasfeo_d_aux_ext_dep.h>

#include "../../include/hpipm_d_ocp_qp_ipm.h"
#include "../../include/hpipm_d_ocp_qp_dim.h"
#include "../../include/hpipm_d_ocp_qp.h"
#include "../../include/hpipm_d_ocp_qp_sol.h"
#include "../../include/hpipm_d_ocp_qp_utils.h"
#include "../../include/hpipm_timing.h"



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
// arg
extern int mode;
extern int iter_max;
extern double alpha_min;
extern double mu0;
extern double tol_stat;
extern double tol_eq;
extern double tol_ineq;
extern double tol_comp;
extern double reg_prim;
extern int warm_start;
extern int pred_corr;
extern int ric_alg;



// main
int main()
	{

	int ii, jj;

	int hpipm_status;

	int rep, nrep=10;

	struct timeval tv0, tv1;

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

	d_ocp_qp_ipm_arg_set_default(mode, &arg);

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

	int ipm_size = d_ocp_qp_ipm_ws_memsize(&dim, &arg);
	void *ipm_mem = malloc(ipm_size);

	struct d_ocp_qp_ipm_ws workspace;
	d_ocp_qp_ipm_ws_create(&dim, &arg, &workspace, ipm_mem);

/************************************************
* ipm solver
************************************************/

//	gettimeofday(&tv0, NULL); // start
	hpipm_timer timer;
	hpipm_tic(&timer);

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
		d_ocp_qp_ipm_solve(&qp, &qp_sol, &arg, &workspace);
		d_ocp_qp_ipm_get_status(&workspace, &hpipm_status);
		}

//	gettimeofday(&tv1, NULL); // stop
//	double time_ipm = (tv1.tv_sec-tv0.tv_sec)/(nrep+0.0)+(tv1.tv_usec-tv0.tv_usec)/(nrep*1e6);
	double time_ipm = hpipm_toc(&timer) / nrep;

// XXX
//exit(1);

/************************************************
* print solution info
************************************************/

    printf("\nHPIPM returned with flag %i.\n", hpipm_status);
    if(hpipm_status == 0)
		{
        printf("\n -> QP solved!\n");
		}
	else if(hpipm_status==1)
		{
        printf("\n -> Solver failed! Maximum number of iterations reached\n");
		}
	else if(hpipm_status==2)
		{
        printf("\n -> Solver failed! Minimum step lenght reached\n");
		}
	else if(hpipm_status==2)
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
		d_ocp_qp_sol_get_u(ii, &qp_sol, u);
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
		d_ocp_qp_sol_get_x(ii, &qp_sol, x);
		d_print_mat(1, nx[ii], x, 1);
		}

	// pi
	double *pi = malloc(nx_max*sizeof(double));

	printf("\npi = \n");
	for(ii=0; ii<N; ii++)
		{
		d_ocp_qp_sol_get_pi(ii, &qp_sol, pi);
		d_print_mat(1, nx[ii+1], pi, 1);
		}

/************************************************
* print ipm statistics
************************************************/

	int iter; d_ocp_qp_ipm_get_iter(&workspace, &iter);
	double res_stat; d_ocp_qp_ipm_get_max_res_stat(&workspace, &res_stat);
	double res_eq; d_ocp_qp_ipm_get_max_res_eq(&workspace, &res_eq);
	double res_ineq; d_ocp_qp_ipm_get_max_res_ineq(&workspace, &res_ineq);
	double res_comp; d_ocp_qp_ipm_get_max_res_comp(&workspace, &res_comp);
	double *stat; d_ocp_qp_ipm_get_stat(&workspace, &stat);
	int stat_m; d_ocp_qp_ipm_get_stat_m(&workspace, &stat_m);

	printf("\nipm return = %d\n", hpipm_status);
	printf("\nipm residuals max: res_g = %e, res_b = %e, res_d = %e, res_m = %e\n", res_stat, res_eq, res_ineq, res_comp);

	printf("\nipm iter = %d\n", iter);
	printf("\nalpha_aff\tmu_aff\t\tsigma\t\talpha\t\tmu\t\tres_stat\tres_eq\t\tres_ineq\tres_comp\n");
	d_print_exp_tran_mat(stat_m, iter+1, stat, stat_m);

	printf("\nocp ipm time = %e [s]\n\n", time_ipm);

/************************************************
* predict solution of QP with new RHS
************************************************/

	void *qp1_mem = malloc(qp_size);
	struct d_ocp_qp qp1;
	d_ocp_qp_create(&dim, &qp1, qp1_mem);

	void *qp_sol1_mem = malloc(qp_sol_size);
	struct d_ocp_qp_sol qp_sol1;
	d_ocp_qp_sol_create(&dim, &qp_sol1, qp_sol1_mem);

	// slightly modify RHS of QP
	d_ocp_qp_copy_all(&qp, &qp1);

	double *lbx0_tmp = malloc(nx[0]*sizeof(double));
	double *ubx0_tmp = malloc(nx[0]*sizeof(double));
	//
	for(ii=0; ii<nx[0]; ii++)
		lbx0_tmp[ii] = hlbx[0][ii];
	for(ii=0; ii<nx[0]; ii++)
		ubx0_tmp[ii] = hubx[0][ii];
	lbx0_tmp[0] = 1.1*hlbx[0][0];
	ubx0_tmp[0] = 1.1*hubx[0][0];
	//
//	lbx0_tmp[1] = 0.95*hlbx[0][1];
//	ubx0_tmp[1] = 0.95*hubx[0][1];

	d_ocp_qp_set_lbx(0, lbx0_tmp, &qp1);
	d_ocp_qp_set_ubx(0, ubx0_tmp, &qp1);

	hpipm_tic(&timer);

	for(rep=0; rep<nrep; rep++)
		{
		d_ocp_qp_ipm_predict(&qp1, &qp_sol1, &arg, &workspace);
		}

	double time_pred = hpipm_toc(&timer) / nrep;

    printf("\nAverage prediction time over %i runs: %e [s]\n", nrep, time_pred);
	printf("\n\n");

	// predicted solution

	// u
	printf("\nu_pred = \n");
	for(ii=0; ii<=N; ii++)
		{
		d_ocp_qp_sol_get_u(ii, &qp_sol1, u);
		d_print_mat(1, nu[ii], u, 1);
		}

	// x
	printf("\nx_pred = \n");
	for(ii=0; ii<=N; ii++)
		{
		d_ocp_qp_sol_get_x(ii, &qp_sol1, x);
		d_print_mat(1, nx[ii], x, 1);
		}

	// pi
	printf("\npi_pred = \n");
	for(ii=0; ii<N; ii++)
		{
		d_ocp_qp_sol_get_pi(ii, &qp_sol1, pi);
		d_print_mat(1, nx[ii+1], pi, 1);
		}

	d_ocp_qp_ipm_get_max_res_stat(&workspace, &res_stat);
	d_ocp_qp_ipm_get_max_res_eq(&workspace, &res_eq);
	d_ocp_qp_ipm_get_max_res_ineq(&workspace, &res_ineq);
	d_ocp_qp_ipm_get_max_res_comp(&workspace, &res_comp);
	printf("\nprediction residuals max: res_g = %e, res_b = %e, res_d = %e, res_m = %e\n", res_stat, res_eq, res_ineq, res_comp);

/************************************************
* sensitivity of solution of QP
************************************************/

	void *qp2_mem = malloc(qp_size);
	struct d_ocp_qp qp2;
	d_ocp_qp_create(&dim, &qp2, qp1_mem);

	// new sol struct
	void *qp_sol2_mem = malloc(qp_sol_size);
	struct d_ocp_qp_sol qp_sol2;
	d_ocp_qp_sol_create(&dim, &qp_sol2, qp_sol2_mem);

	// set I to param at RHS
	d_ocp_qp_copy_all(&qp, &qp2);

	d_ocp_qp_set_rhs_zero(&qp2);

	double one = 1.0;
	int index = 0;
//	d_ocp_qp_set_el_lbx(0, index, &one, &qp2);
//	d_ocp_qp_set_el_ubx(0, index, &one, &qp2);
	d_ocp_qp_set_el("lbx", 0, index, &one, &qp2);
	d_ocp_qp_set_el("ubx", 0, index, &one, &qp2);

//	d_ocp_qp_print(&dim, &qp2);
//	exit(1);

	// sensitivity solution
//	int comp_res_pred = 0;
//	d_ocp_qp_ipm_arg_set_comp_res_pred(&comp_res_pred, &arg);

	d_ocp_qp_ipm_sens(&qp2, &qp_sol2, &arg, &workspace);

	// u
	printf("\nu_sens = \n");
	for(ii=0; ii<=N; ii++)
		{
		d_ocp_qp_sol_get_u(ii, &qp_sol2, u);
		d_print_mat(1, nu[ii], u, 1);
		}

	// x
	printf("\nx_sens = \n");
	for(ii=0; ii<=N; ii++)
		{
		d_ocp_qp_sol_get_x(ii, &qp_sol2, x);
		d_print_mat(1, nx[ii], x, 1);
		}

	// pi
	printf("\npi_sens = \n");
	for(ii=0; ii<N; ii++)
		{
		d_ocp_qp_sol_get_pi(ii, &qp_sol2, pi);
		d_print_mat(1, nx[ii+1], pi, 1);
		}

/************************************************
* free memory and return
************************************************/

    free(dim_mem);
    free(qp_mem);
    free(qp1_mem);
    free(qp2_mem);
	free(qp_sol_mem);
	free(qp_sol1_mem);
	free(qp_sol2_mem);
	free(ipm_arg_mem);
	free(ipm_mem);

	free(lbx0_tmp);
	free(ubx0_tmp);

	free(u);
	free(x);
	free(pi);

	return 0;

	}


