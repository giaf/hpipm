/**************************************************************************************************
*                                                                                                 *
* This file is part of HPIPM.                                                                     *
*                                                                                                 *
* HPIPM -- High-Performance Interior Point Method.                                                *
* Copyright (C) 2019 by Gianluca Frison.                                                          *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
*                                                                                                 *
* The 2-Clause BSD License                                                                        *
*                                                                                                 *
* Redistribution and use in source and binary forms, with or without                              *
* modification, are permitted provided that the following conditions are met:                     *
*                                                                                                 *
* 1. Redistributions of source code must retain the above copyright notice, this                  *
*    list of conditions and the following disclaimer.                                             *
* 2. Redistributions in binary form must reproduce the above copyright notice,                    *
*    this list of conditions and the following disclaimer in the documentation                    *
*    and/or other materials provided with the distribution.                                       *
*                                                                                                 *
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND                 *
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED                   *
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE                          *
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR                 *
* ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES                  *
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;                    *
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND                     *
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT                      *
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS                   *
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                                    *
*                                                                                                 *
* Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             *
*                                                                                                 *
**************************************************************************************************/

#include <stdlib.h>
#include <stdio.h>

#include <blasfeo_d_aux_ext_dep.h>

#include <hpipm_d_dense_qp_ipm.h>
#include <hpipm_d_dense_qp_dim.h>
#include <hpipm_d_dense_qp.h>
#include <hpipm_d_dense_qp_sol.h>
#include <hpipm_d_ocp_qp_dim.h>
#include <hpipm_d_ocp_qp.h>
#include <hpipm_d_ocp_qp_sol.h>
#include <hpipm_d_ocp_qp_res.h>
#include <hpipm_d_ocp_qp_seed.h>
#include <hpipm_d_ocp_qp_utils.h>
#include <hpipm_d_cond.h>
#include <hpipm_timing.h>



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
extern int *nbue;
extern int *nbxe;
extern int *nge;
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
extern double **hlbx_mask;
extern double **hubx;
extern double **hubx_mask;
extern int **hidxbu;
extern double **hlbu;
extern double **hlbu_mask;
extern double **hubu;
extern double **hubu_mask;
extern double **hC;
extern double **hD;
extern double **hlg;
extern double **hlg_mask;
extern double **hug;
extern double **hug_mask;
extern double **hZl;
extern double **hZu;
extern double **hzl;
extern double **hzu;
extern int **hidxs;
extern double **hlls;
extern double **hlls_mask;
extern double **hlus;
extern double **hlus_mask;
extern int **hidxe;
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

	int ii;

	int hpipm_status;

	int rep, nrep=10;

	hpipm_timer timer;

/************************************************
* ocp qp dim
************************************************/

	hpipm_size_t dim_size = d_ocp_qp_dim_memsize(N);
	void *dim_mem = malloc(dim_size);

	struct d_ocp_qp_dim dim;
	d_ocp_qp_dim_create(N, &dim, dim_mem);

	// unified setter
	d_ocp_qp_dim_set_all(nx, nu, nbx, nbu, ng, nsbx, nsbu, nsg, &dim);

	// additional single setters

	// set number of inequality constr to be considered as equality constr
	// (ignored in this example)
	//for(ii=0; ii<=N; ii++)
	//	{
	//	d_ocp_qp_dim_set_nbxe(ii, nbxe[ii], &dim);
	//	d_ocp_qp_dim_set_nbue(ii, nbue[ii], &dim);
	//	d_ocp_qp_dim_set_nge(ii, nge[ii], &dim);
	//	}

/************************************************
* dense qp dim cond
************************************************/

	hpipm_size_t dim2_size = d_dense_qp_dim_memsize();
	void *dim2_mem = malloc(dim2_size);

	struct d_dense_qp_dim dim2;
	d_dense_qp_dim_create(&dim2, dim2_mem);

	d_cond_qp_compute_dim(&dim, &dim2);

/************************************************
* ocp qp
************************************************/

	hpipm_size_t qp_size = d_ocp_qp_memsize(&dim);
	void *qp_mem = malloc(qp_size);

	struct d_ocp_qp qp;
	d_ocp_qp_create(&dim, &qp, qp_mem);

	// unified setter
	d_ocp_qp_set_all(hA, hB, hb, hQ, hS, hR, hq, hr, hidxbx, hlbx, hubx, hidxbu, hlbu, hubu, hC, hD, hlg, hug, hZl, hZu, hzl, hzu, hidxs, hlls, hlus, &qp);

	// additional single setters

	// mark the inequality constr to be considered as equality constr
	// (ignored in this example)
	//for(ii=0; ii<=N; ii++)
	//	{
	//	d_ocp_qp_set_idxe(ii, hidxe[ii], &qp);
	//	}

	// set inequality constraints mask
	for(ii=0; ii<=N; ii++)
		{
		d_ocp_qp_set_lbu_mask(ii, hlbu_mask[ii], &qp);
		d_ocp_qp_set_ubu_mask(ii, hubu_mask[ii], &qp);
		d_ocp_qp_set_lbx_mask(ii, hlbx_mask[ii], &qp);
		d_ocp_qp_set_ubx_mask(ii, hubx_mask[ii], &qp);
		d_ocp_qp_set_lg_mask(ii, hlg_mask[ii], &qp);
		d_ocp_qp_set_ug_mask(ii, hug_mask[ii], &qp);
		d_ocp_qp_set_lls_mask(ii, hlls_mask[ii], &qp);
		d_ocp_qp_set_lus_mask(ii, hlus_mask[ii], &qp);
		}

/************************************************
* dense qp cond
************************************************/

	hpipm_size_t qp2_size = d_dense_qp_memsize(&dim2);
	void *qp2_mem = malloc(qp2_size);

	struct d_dense_qp qp2;
	d_dense_qp_create(&dim2, &qp2, qp2_mem);

/************************************************
* ocp qp sol
************************************************/

	hpipm_size_t qp_sol_size = d_ocp_qp_sol_memsize(&dim);
	void *qp_sol_mem = malloc(qp_sol_size);

	struct d_ocp_qp_sol qp_sol;
	d_ocp_qp_sol_create(&dim, &qp_sol, qp_sol_mem);

/************************************************
* dense qp sol cond
************************************************/

	hpipm_size_t qp_sol2_size = d_dense_qp_sol_memsize(&dim2);
	void *qp_sol2_mem = malloc(qp_sol2_size);

	struct d_dense_qp_sol qp_sol2;
	d_dense_qp_sol_create(&dim2, &qp_sol2, qp_sol2_mem);

/************************************************
* cond arg
************************************************/

	hpipm_size_t cond_arg_size = d_cond_qp_arg_memsize();
	void *cond_arg_mem = malloc(cond_arg_size);

	struct d_cond_qp_arg cond_arg;
	d_cond_qp_arg_create(&cond_arg, cond_arg_mem);

	d_cond_qp_arg_set_default(&cond_arg);

//	d_cond_qp_arg_set_ric_alg(0, &part_cond_arg);
//	d_cond_qp_arg_set_comp_dual_sol_eq(0, &part_cond_arg);

/************************************************
* ipm arg
************************************************/

	hpipm_size_t ipm_arg_size = d_dense_qp_ipm_arg_memsize(&dim2);
	void *ipm_arg_mem = malloc(ipm_arg_size);

	struct d_dense_qp_ipm_arg arg;
	d_dense_qp_ipm_arg_create(&dim2, &arg, ipm_arg_mem);

	d_dense_qp_ipm_arg_set_default(mode, &arg);

	d_dense_qp_ipm_arg_set_mu0(&mu0, &arg);
	d_dense_qp_ipm_arg_set_iter_max(&iter_max, &arg);
	d_dense_qp_ipm_arg_set_tol_stat(&tol_stat, &arg);
	d_dense_qp_ipm_arg_set_tol_eq(&tol_eq, &arg);
	d_dense_qp_ipm_arg_set_tol_ineq(&tol_ineq, &arg);
	d_dense_qp_ipm_arg_set_tol_comp(&tol_comp, &arg);
	d_dense_qp_ipm_arg_set_reg_prim(&reg_prim, &arg);
	d_dense_qp_ipm_arg_set_warm_start(&warm_start, &arg);

/************************************************
* cond workspace
************************************************/

	hpipm_size_t cond_size = d_cond_qp_ws_memsize(&dim, &cond_arg);
	void *cond_mem = malloc(cond_size);

	struct d_cond_qp_ws cond_ws;
	d_cond_qp_ws_create(&dim, &cond_arg, &cond_ws, cond_mem);

/************************************************
* ipm workspace
************************************************/

	hpipm_size_t ipm_size = d_dense_qp_ipm_ws_memsize(&dim2, &arg);
	void *ipm_mem = malloc(ipm_size);

	struct d_dense_qp_ipm_ws workspace;
	d_dense_qp_ipm_ws_create(&dim2, &arg, &workspace, ipm_mem);

/************************************************
* cond
************************************************/

	hpipm_tic(&timer);

	for(rep=0; rep<nrep; rep++)
		{
		d_cond_qp_cond(&qp, &qp2, &cond_arg, &cond_ws);
		}

	double time_cond = hpipm_toc(&timer) / nrep;

/************************************************
* ipm solver
************************************************/

	hpipm_tic(&timer);

	for(rep=0; rep<nrep; rep++)
		{
		d_dense_qp_ipm_solve(&qp2, &qp_sol2, &arg, &workspace);
		d_dense_qp_ipm_get_status(&workspace, &hpipm_status);
		}

	double time_ipm = hpipm_toc(&timer) / nrep;

/************************************************
* part expand
************************************************/

	hpipm_tic(&timer);

	for(rep=0; rep<nrep; rep++)
		{
		d_cond_qp_expand_sol(&qp, &qp_sol2, &qp_sol, &cond_arg, &cond_ws);
		}

	double time_expa = hpipm_toc(&timer) / nrep;

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

	int iter; d_dense_qp_ipm_get_iter(&workspace, &iter);
	double res_stat; d_dense_qp_ipm_get_max_res_stat(&workspace, &res_stat);
	double res_eq; d_dense_qp_ipm_get_max_res_eq(&workspace, &res_eq);
	double res_ineq; d_dense_qp_ipm_get_max_res_ineq(&workspace, &res_ineq);
	double res_comp; d_dense_qp_ipm_get_max_res_comp(&workspace, &res_comp);
	double *stat; d_dense_qp_ipm_get_stat(&workspace, &stat);
	int stat_m; d_dense_qp_ipm_get_stat_m(&workspace, &stat_m);

	printf("\nipm return = %d\n", hpipm_status);
	printf("\nipm residuals max: res_g = %e, res_b = %e, res_d = %e, res_m = %e\n", res_stat, res_eq, res_ineq, res_comp);

	printf("\nipm iter = %d\n", iter);
	printf("\nalpha_aff\tmu_aff\t\tsigma\t\talpha_prim\talpha_dual\tmu\t\tres_stat\tres_eq\t\tres_ineq\tres_comp\tdual_gap\tobj\t\tlq fact\t\titref pred\titref corr\tlin res stat\tlin res eq\tlin res ineq\tlin res comp\n");
	d_print_exp_tran_mat(stat_m, iter+1, stat, stat_m);

	printf("\npart cond time = %e [s]\n\n", time_cond);
	printf("\ndense ipm time = %e [s]\n\n", time_ipm);
	printf("\npart expa time = %e [s]\n\n", time_expa);
	printf("\ntotal time = %e [s]\n\n", time_cond+time_ipm+time_expa);

/************************************************
* residuals of original QP
************************************************/

	hpipm_size_t res_size = d_ocp_qp_res_memsize(&dim);
	void *res_mem = malloc(res_size);
	struct d_ocp_qp_res res;
	d_ocp_qp_res_create(&dim, &res, res_mem);

	hpipm_size_t res_ws_size = d_ocp_qp_res_ws_memsize(&dim);
	void *res_ws_mem = malloc(res_ws_size);
	struct d_ocp_qp_res_ws res_ws;
	d_ocp_qp_res_ws_create(&dim, &res_ws, res_ws_mem);

	d_ocp_qp_res_compute(&qp, &qp_sol, &res, &res_ws);
	d_ocp_qp_res_compute_inf_norm(&res);

	d_ocp_qp_res_get_max_res_stat(&res, &res_stat);
	d_ocp_qp_res_get_max_res_eq(&res, &res_eq);
	d_ocp_qp_res_get_max_res_ineq(&res, &res_ineq);
	d_ocp_qp_res_get_max_res_comp(&res, &res_comp);

	printf("\noriginal qp residuals max: res_g = %e, res_b = %e, res_d = %e, res_m = %e\n", res_stat, res_eq, res_ineq, res_comp);

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

	// cond RHS
	d_cond_qp_cond_rhs(&qp1, &qp2, &cond_arg, &cond_ws);

	// solve
	d_dense_qp_ipm_predict(&qp2, &qp_sol2, &arg, &workspace);

	// expand sol
	d_cond_qp_expand_sol(&qp1, &qp_sol2, &qp_sol1, &cond_arg, &cond_ws);

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

	d_ocp_qp_res_compute(&qp1, &qp_sol1, &res, &res_ws);
	d_ocp_qp_res_compute_inf_norm(&res);

	d_ocp_qp_res_get_max_res_stat(&res, &res_stat);
	d_ocp_qp_res_get_max_res_eq(&res, &res_eq);
	d_ocp_qp_res_get_max_res_ineq(&res, &res_ineq);
	d_ocp_qp_res_get_max_res_comp(&res, &res_comp);

	printf("\nprediction residuals max: res_g = %e, res_b = %e, res_d = %e, res_m = %e\n", res_stat, res_eq, res_ineq, res_comp);

/************************************************
* sensitivity of solution of QP
************************************************/

	// seed struct
	hpipm_size_t seed_size = d_ocp_qp_seed_memsize(&dim);
	void *seed_mem = malloc(seed_size);
	struct d_ocp_qp_seed seed;
	d_ocp_qp_seed_create(&dim, &seed, seed_mem);

	hpipm_size_t seed2_size = d_dense_qp_seed_memsize(&dim2);
	void *seed2_mem = malloc(seed2_size);
	struct d_dense_qp_seed seed2;
	d_dense_qp_seed_create(&dim2, &seed2, seed2_mem);

	// new sol struct
	void *sens_mem = malloc(qp_sol_size);
	struct d_ocp_qp_sol sens;
	d_ocp_qp_sol_create(&dim, &sens, sens_mem);

	void *sens2_mem = malloc(qp_sol2_size);
	struct d_dense_qp_sol sens2;
	d_dense_qp_sol_create(&dim2, &sens2, sens2_mem);

	// set seeds to zero
	d_ocp_qp_seed_set_zero(&seed);

	// set I to param
	double *seed_x0 = malloc(nx[0]*sizeof(double));
	for(ii=0; ii<nx[0]; ii++)
		seed_x0[ii] = 0.0;
	int index = 0;
	seed_x0[index] = 1.0;
	int stage = 0;
	d_ocp_qp_seed_set_seed_lbx(stage, seed_x0, &seed);
	d_ocp_qp_seed_set_seed_ubx(stage, seed_x0, &seed);

	// print seeds
	//d_ocp_qp_seed_print(seed.dim, &seed);

	// cond RHS
	d_cond_qp_cond_seed(&qp, &seed, &seed2, &cond_arg, &cond_ws);

	// forward sensitivity of solution
	d_dense_qp_ipm_sens_frw(&qp2, &seed2, &sens2, &arg, &workspace);
	// forward sensitivity of solution
	//d_dense_qp_ipm_sens_adj(&qp2, &seed2, &sens2, &arg, &workspace);

	// expand sens
	d_cond_qp_expand_sol(&qp, &sens2, &sens, &cond_arg, &cond_ws);

	// u
	printf("\nu_sens = \n");
	for(ii=0; ii<=N; ii++)
		{
		d_ocp_qp_sol_get_u(ii, &sens, u);
		d_print_mat(1, nu[ii], u, 1);
		}

	// x
	printf("\nx_sens = \n");
	for(ii=0; ii<=N; ii++)
		{
		d_ocp_qp_sol_get_x(ii, &sens, x);
		d_print_mat(1, nx[ii], x, 1);
		}

	// pi
	printf("\npi_sens = \n");
	for(ii=0; ii<N; ii++)
		{
		d_ocp_qp_sol_get_pi(ii, &sens, pi);
		d_print_mat(1, nx[ii+1], pi, 1);
		}

	// print solution sensitivities
	//d_ocp_qp_sol_print(sens.dim, &sens);

/************************************************
* free memory and return
************************************************/

    free(dim_mem);
    free(dim2_mem);
    free(qp_mem);
    free(qp1_mem);
    free(qp2_mem);
	free(qp_sol_mem);
	free(qp_sol2_mem);
	free(seed_mem);
	free(seed2_mem);
	free(sens_mem);
	free(sens2_mem);
	free(cond_arg_mem);
	free(ipm_arg_mem);
	free(cond_mem);
	free(ipm_mem);
	free(res_mem);
	free(res_ws_mem);

	free(u);
	free(x);
	free(pi);

	return 0;

	}



