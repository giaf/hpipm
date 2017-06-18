/**************************************************************************************************
*                                                                                                 *
* This file is part of HPIPM.                                                                     *
*                                                                                                 *
* HPIPM -- High Performance Interior Point Method.                                                *
* Copyright (C) 2017 by Gianluca Frison.                                                          *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
*                                                                                                 *
* HPMPC is free software; you can redistribute it and/or                                          *
* modify it under the terms of the GNU Lesser General Public                                      *
* License as published by the Free Software Foundation; either                                    *
* version 2.1 of the License, or (at your option) any later version.                              *
*                                                                                                 *
* HPMPC is distributed in the hope that it will be useful,                                        *
* but WITHOUT ANY WARRANTY; without even the implied warranty of                                  *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                                            *
* See the GNU Lesser General Public License for more details.                                     *
*                                                                                                 *
* You should have received a copy of the GNU Lesser General Public                                *
* License along with HPMPC; if not, write to the Free Software                                    *
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

#include "../include/hpipm_d_dense_qp.h"
#include "../include/hpipm_d_dense_qp_sol.h"
#include "../include/hpipm_d_dense_qp_ipm_hard.h"



#define PRINT 1



#if ! defined(EXT_DEP)
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
#endif



int main()
	{

	int ii;

/************************************************
* qp dimension and data
************************************************/	

	int nv = 2;
	int ne = 1;
	int nb = 2;
	int ng = 0;

	double H[] = {4.0, 1.0, 1.0, 2.0};
	double g[] = {1.0, 1.0};
	double A[] = {1.0, 1.0};
	double b[] = {1.0};
//	double d_lb[] = {0.0, 0.0};
//	double d_ub[] = {INFINITY, INFINITY};
	double d_lb[] = {-1.0, -1.0};
	double d_ub[] = {1.5, 0.5};
	int idxb[] = {0, 1};
	double C[] = {1.0, 0.0, 0.0, 1.0};
	double d_lg[] = {-1.0, -1.0};
	double d_ug[] = {1.5, 0.5};

/************************************************
* dense qp
************************************************/	

	int qp_size = d_memsize_dense_qp(nv, ne, nb, ng);
	printf("\nqp size = %d\n", qp_size);
	void *qp_mem = malloc(qp_size);

	struct d_dense_qp qp;
	d_create_dense_qp(nv, ne, nb, ng, &qp, qp_mem);
	d_cvt_colmaj_to_dense_qp(H, g, A, b, idxb, d_lb, d_ub, C, d_lg, d_ug, &qp);

#if 0
	d_print_strmat(nv+1, nv, qp.Hg, 0, 0);
	d_print_strmat(ne, nv, qp.A, 0, 0);
	d_print_strmat(nv, ng, qp.Ct, 0, 0);
	d_print_strvec(nv, qp.g, 0);
	d_print_strvec(ne, qp.b, 0);
	d_print_strvec(2*nb+2*ng, qp.d, 0);
	d_print_strvec(nb, qp.d_lb, 0);
	d_print_strvec(nb, qp.d_ub, 0);
	d_print_strvec(ng, qp.d_lg, 0);
	d_print_strvec(ng, qp.d_ug, 0);
#endif

/************************************************
* dense qp sol
************************************************/	

	int qp_sol_size = d_memsize_dense_qp_sol(nv, ne, nb, ng);
	printf("\nqp sol size = %d\n", qp_sol_size);
	void *qp_sol_mem = malloc(qp_sol_size);

	struct d_dense_qp_sol qp_sol;
	d_create_dense_qp_sol(nv, ne, nb, ng, &qp_sol, qp_sol_mem);

/************************************************
* ipm
************************************************/	

	struct d_ipm_hard_dense_qp_arg arg;
	arg.alpha_min = 1e-8;
	arg.mu_max = 1e-12;
	arg.iter_max = 20;
	arg.mu0 = 1.0;

	int ipm_size = d_memsize_ipm_hard_dense_qp(&qp, &arg);
	printf("\nipm size = %d\n", ipm_size);
	void *ipm_mem = malloc(ipm_size);

	struct d_ipm_hard_dense_qp_workspace workspace;
	d_create_ipm_hard_dense_qp(&qp, &arg, &workspace, ipm_mem);

	int rep, nrep=1000;

	struct timeval tv0, tv1;

	gettimeofday(&tv0, NULL); // start

	for(rep=0; rep<nrep; rep++)
		{
//		d_solve_ipm_hard_dense_qp(&qp, &qp_sol, &workspace);
		d_solve_ipm2_hard_dense_qp(&qp, &qp_sol, &workspace);
		}

	gettimeofday(&tv1, NULL); // stop

	double time_dense_ipm = (tv1.tv_sec-tv0.tv_sec)/(nrep+0.0)+(tv1.tv_usec-tv0.tv_usec)/(nrep*1e6);

	printf("\nsolution\n\n");
	printf("\nv\n");
	d_print_mat(1, nv, qp_sol.v->pa, 1);
	printf("\npi\n");
	d_print_mat(1, ne, qp_sol.pi->pa, 1);
	printf("\nlam_lb\n");
	d_print_mat(1, nb, qp_sol.lam_lb->pa, 1);
	printf("\nlam_ub\n");
	d_print_mat(1, nb, qp_sol.lam_ub->pa, 1);
	printf("\nlam_lg\n");
	d_print_mat(1, ng, qp_sol.lam_lg->pa, 1);
	printf("\nlam_ug\n");
	d_print_mat(1, ng, qp_sol.lam_ug->pa, 1);
	printf("\nt_lb\n");
	d_print_mat(1, nb, qp_sol.t_lb->pa, 1);
	printf("\nt_ub\n");
	d_print_mat(1, nb, qp_sol.t_ub->pa, 1);
	printf("\nt_lg\n");
	d_print_mat(1, ng, qp_sol.t_lg->pa, 1);
	printf("\nt_ug\n");
	d_print_mat(1, ng, qp_sol.t_ug->pa, 1);

	printf("\nresiduals\n\n");
	printf("\nres_g\n");
	d_print_e_mat(1, nv, workspace.res_g->pa, 1);
	printf("\nres_b\n");
	d_print_e_mat(1, ne, workspace.res_b->pa, 1);
	printf("\nres_d\n");
	d_print_e_mat(1, 2*nb+2*ng, workspace.res_d->pa, 1);
	printf("\nres_m\n");
	d_print_e_mat(1, 2*nb+2*ng, workspace.res_m->pa, 1);
	printf("\nres_mu\n");
	printf("\n%e\n\n", workspace.res_mu);

	printf("\nipm iter = %d\n", workspace.iter);
	printf("\nalpha_aff\tmu_aff\t\tsigma\t\talpha\t\tmu\n");
	d_print_e_tran_mat(5, workspace.iter, workspace.stat, 5);

	printf("\ndense ipm time = %e [s]\n\n", time_dense_ipm);

/************************************************
* free memory
************************************************/	

	free(qp_mem);
	free(qp_sol_mem);
	free(ipm_mem);

/************************************************
* return
************************************************/	

	return 0;

	}

	
