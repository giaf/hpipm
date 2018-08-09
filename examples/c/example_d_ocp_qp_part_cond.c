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

#include <blasfeo_d_aux_ext_dep.h>

#include "../../include/hpipm_d_ocp_qp_ipm.h"
#include "../../include/hpipm_d_ocp_qp_dim.h"
#include "../../include/hpipm_d_ocp_qp.h"
#include "../../include/hpipm_d_ocp_qp_sol.h"
#include "../include/hpipm_d_part_cond.h"



// qp data as global data
extern int N;
extern int *nx;
extern int *nu;
extern int *nbu;
extern int *nbx;
extern int *ng;
extern int *ns;
extern double **hA;
extern double **hB;
extern double **hb;
extern double **hQ;
extern double **hR;
extern double **hS;
extern double **hq;
extern double **hr;
extern int **hidxb;
extern double **hlb;
extern double **hub;
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



// main
int main()
	{

	int ii;

	int hpipm_return;

	int rep, nrep=10;

	struct timeval tv0, tv1;

/************************************************
* ocp qp dim
************************************************/

	int dim_size = d_memsize_ocp_qp_dim(N);
	void *dim_mem = malloc(dim_size);

	struct d_ocp_qp_dim dim;
	d_create_ocp_qp_dim(N, &dim, dim_mem);

	d_cvt_int_to_ocp_qp_dim(N, nx, nu, nbx, nbu, ng, ns, &dim);

/************************************************
* ocp qp dim part cond
************************************************/

	// horizon length of partially condensed OCP QP
	int N2 = 5;

	int dim_size2 = d_memsize_ocp_qp_dim(N2);
	void *dim_mem2 = malloc(dim_size2);

	struct d_ocp_qp_dim dim2;
	d_create_ocp_qp_dim(N2, &dim2, dim_mem2);

	int *block_size = malloc((N+1)*sizeof(int));
	d_compute_block_size_cond_qp_ocp2ocp(N, N2, block_size);
//	block_size[0] = 1;
//	block_size[1] = 1;
//	printf("\nblock_size\n");
//	int_print_mat(1, N2+1, block_size, 1);

	d_compute_qp_dim_ocp2ocp(&dim, block_size, &dim2);

/************************************************
* ocp qp
************************************************/

	int qp_size = d_memsize_ocp_qp(&dim);
	void *qp_mem = malloc(qp_size);

	struct d_ocp_qp qp;
	d_create_ocp_qp(&dim, &qp, qp_mem);

	d_cvt_colmaj_to_ocp_qp(hA, hB, hb, hQ, hS, hR, hq, hr, hidxb, hlb, hub, hC, hD, hlg, hug, hZl, hZu, hzl, hzu, hidxs, hlls, hlus, &qp);

/************************************************
* ocp qp part cond
************************************************/

	int qp_size2 = d_memsize_ocp_qp(&dim2);
	void *qp_mem2 = malloc(qp_size2);

	struct d_ocp_qp qp2;
	d_create_ocp_qp(&dim2, &qp2, qp_mem2);

/************************************************
* ocp qp sol
************************************************/

	int qp_sol_size = d_memsize_ocp_qp_sol(&dim);
	void *qp_sol_mem = malloc(qp_sol_size);

	struct d_ocp_qp_sol qp_sol;
	d_create_ocp_qp_sol(&dim, &qp_sol, qp_sol_mem);

/************************************************
* ocp qp sol part cond
************************************************/

	int qp_sol_size2 = d_memsize_ocp_qp_sol(&dim2);
	void *qp_sol_mem2 = malloc(qp_sol_size2);

	struct d_ocp_qp_sol qp_sol2;
	d_create_ocp_qp_sol(&dim2, &qp_sol2, qp_sol_mem2);

/************************************************
* part cond arg
************************************************/

	int part_cond_arg_size = d_memsize_cond_qp_ocp2ocp_arg(dim2.N);
	void *part_cond_arg_mem = malloc(part_cond_arg_size);

	struct d_cond_qp_ocp2ocp_arg part_cond_arg;
	d_create_cond_qp_ocp2ocp_arg(dim2.N, &part_cond_arg, part_cond_arg_mem);

	d_set_default_cond_qp_ocp2ocp_arg(dim2.N, &part_cond_arg);

/************************************************
* ipm arg
************************************************/

	int ipm_arg_size = d_memsize_ocp_qp_ipm_arg(&dim2);
	void *ipm_arg_mem = malloc(ipm_arg_size);

	struct d_ocp_qp_ipm_arg arg;
	d_create_ocp_qp_ipm_arg(&dim2, &arg, ipm_arg_mem);

//	enum hpipm_mode mode = SPEED_ABS;
	enum hpipm_mode mode = SPEED;
//	enum hpipm_mode mode = BALANCE;
//	enum hpipm_mode mode = ROBUST;
	d_set_default_ocp_qp_ipm_arg(mode, &arg);

	d_set_ocp_qp_ipm_arg_mu0(1e4, &arg);
	d_set_ocp_qp_ipm_arg_iter_max(30, &arg);
	d_set_ocp_qp_ipm_arg_tol_stat(1e-4, &arg);
	d_set_ocp_qp_ipm_arg_tol_eq(1e-5, &arg);
	d_set_ocp_qp_ipm_arg_tol_ineq(1e-5, &arg);
	d_set_ocp_qp_ipm_arg_tol_comp(1e-5, &arg);
	d_set_ocp_qp_ipm_arg_reg_prim(1e-12, &arg);

/************************************************
* part cond workspace
************************************************/

	int part_cond_size = d_memsize_cond_qp_ocp2ocp(&dim, block_size, &dim2, &part_cond_arg);
	void *part_cond_mem = malloc(part_cond_size);

	struct d_cond_qp_ocp2ocp_workspace part_cond_ws;
	d_create_cond_qp_ocp2ocp(&dim, block_size, &dim2, &part_cond_arg, &part_cond_ws, part_cond_mem);

/************************************************
* ipm workspace
************************************************/

	int ipm_size = d_memsize_ocp_qp_ipm(&dim2, &arg);
	void *ipm_mem = malloc(ipm_size);

	struct d_ocp_qp_ipm_workspace workspace;
	d_create_ocp_qp_ipm(&dim2, &arg, &workspace, ipm_mem);

/************************************************
* part cond
************************************************/

	gettimeofday(&tv0, NULL); // start

	for(rep=0; rep<nrep; rep++)
		{
		d_cond_qp_ocp2ocp(&qp, &qp2, &part_cond_arg, &part_cond_ws);
		}

	gettimeofday(&tv1, NULL); // stop

	double time_cond = (tv1.tv_sec-tv0.tv_sec)/(nrep+0.0)+(tv1.tv_usec-tv0.tv_usec)/(nrep*1e6);

/************************************************
* ipm solver
************************************************/

	gettimeofday(&tv0, NULL); // start

	for(rep=0; rep<nrep; rep++)
		{
		hpipm_return = d_solve_ocp_qp_ipm(&qp2, &qp_sol2, &arg, &workspace);
		}

	gettimeofday(&tv1, NULL); // stop

	double time_ipm = (tv1.tv_sec-tv0.tv_sec)/(nrep+0.0)+(tv1.tv_usec-tv0.tv_usec)/(nrep*1e6);

/************************************************
* part expand
************************************************/

	gettimeofday(&tv0, NULL); // start

	for(rep=0; rep<nrep; rep++)
		{
		d_expand_sol_ocp2ocp(&qp, &qp2, &qp_sol2, &qp_sol, &part_cond_arg, &part_cond_ws);
		}

	gettimeofday(&tv1, NULL); // stop

	double time_expa = (tv1.tv_sec-tv0.tv_sec)/(nrep+0.0)+(tv1.tv_usec-tv0.tv_usec)/(nrep*1e6);

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
		d_cvt_ocp_qp_sol_to_colmaj_u(ii, &qp_sol, u);
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
		d_cvt_ocp_qp_sol_to_colmaj_x(ii, &qp_sol, x);
		d_print_mat(1, nx[ii], x, 1);
		}

/************************************************
* print ipm statistics
************************************************/

	int iter = d_get_ocp_qp_ipm_iter(&workspace);
	double res_stat = d_get_ocp_qp_ipm_res_stat(&workspace);
	double res_eq = d_get_ocp_qp_ipm_res_eq(&workspace);
	double res_ineq = d_get_ocp_qp_ipm_res_ineq(&workspace);
	double res_comp = d_get_ocp_qp_ipm_res_comp(&workspace);
	double *stat = d_get_ocp_qp_ipm_stat(&workspace);

	printf("\nipm return = %d\n", hpipm_return);
	printf("\nipm residuals max: res_g = %e, res_b = %e, res_d = %e, res_m = %e\n", res_stat, res_eq, res_ineq, res_comp);

	printf("\nipm iter = %d\n", iter);
	printf("\nalpha_aff\tmu_aff\t\tsigma\t\talpha\t\tmu\n");
	d_print_exp_tran_mat(5, iter, stat, 5);

	printf("\npart cond time = %e [s]\n\n", time_cond);
	printf("\nocp ipm time = %e [s]\n\n", time_ipm);
	printf("\npart expa time = %e [s]\n\n", time_expa);
	printf("\ntotal time = %e [s]\n\n", time_cond+time_ipm+time_expa);

/************************************************
* free memory and return
************************************************/

    free(dim_mem);
    free(dim_mem2);
    free(block_size);
    free(qp_mem);
    free(qp_mem2);
	free(qp_sol_mem);
	free(qp_sol_mem2);
	free(part_cond_arg_mem);
	free(ipm_arg_mem);
	free(part_cond_mem);
	free(ipm_mem);

	free(u);
	free(x);

	return 0;

	}


