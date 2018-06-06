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
#include "../include/hpipm_d_dense_qp_res.h"
#include "../include/hpipm_d_dense_qp_ipm.h"



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

#if 1
	int nv = 2;
	int ne = 1;
	int nb = 2;
	int ng = 0;
	int ns = 2;
	int nsb = 2;
	int nsg = 0;

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
	double Zl[] = {1e3, 1e3};
	double Zu[] = {1e3, 1e3};
	double zl[] = {1e2, 1e2};
	double zu[] = {1e2, 1e2};
	int idxs[] = {0, 1};
	double d_ls[] = {0, 0};
	double d_us[] = {0, 0};
#else
	int nv = 10;
	int ne = 0;
	int nb = nv;
	int ng = 0;
	int ns = 0;

	double H[nv*nv]; for(ii=0; ii<nv*nv; ii++) H[ii] = 0.0; for(ii=0; ii<nv; ii++) H[ii*(nv+1)] = 1.0;
	double g[nv]; for(ii=0; ii<nv; ii++) g[ii] = 10.0;
	double A[ne*nv]; for(ii=0; ii<ne*nv; ii++) A[ii] = 0.0;
	double b[ne]; for(ii=0; ii<ne; ii++) b[ii] = 0.0;
	double d_lb[nb]; for(ii=0; ii<nb; ii++) d_lb[ii] = -1.0;
	double d_ub[nb]; for(ii=0; ii<nb; ii++) d_ub[ii] =  1.0;
	int idxb[nb]; for(ii=0; ii<nb; ii++) idxb[ii] = ii;
	double C[0];
	double d_lg[0];
	double d_ug[0];
	double Zl[0];
	double Zu[0];
	double zl[0];
	double zu[0];
	int idxs[0];
	double d_ls[] = {0, 0};
	double d_us[] = {0, 0};
#endif

/************************************************
* dense qp dim
************************************************/

	int dense_qp_dim_size = d_memsize_dense_qp_dim();
	printf("\nqp dim size = %d\n", dense_qp_dim_size);
	void *dense_qp_dim_mem = malloc(dense_qp_dim_size);

	struct d_dense_qp_dim qp_dim;
	d_create_dense_qp_dim(&qp_dim, dense_qp_dim_mem);

	d_cvt_int_to_dense_qp_dim(nv, ne, nb, ng, nsb, nsg, &qp_dim);

/************************************************
* dense qp
************************************************/

	int qp_size = d_memsize_dense_qp(&qp_dim);
	printf("\nqp size = %d\n", qp_size);
	void *qp_mem = malloc(qp_size);

	struct d_dense_qp qp;
	d_create_dense_qp(&qp_dim, &qp, qp_mem);
	d_cvt_colmaj_to_dense_qp(H, g, A, b, idxb, d_lb, d_ub, C, d_lg, d_ug, Zl, Zu, zl, zu, idxs, d_ls, d_us, &qp);

#if 1
	printf("\nH = \n");
	blasfeo_print_dmat(nv, nv, qp.Hv, 0, 0);
	printf("\nA = \n");
	blasfeo_print_dmat(ne, nv, qp.A, 0, 0);
	printf("\nCt = \n");
	blasfeo_print_dmat(nv, ng, qp.Ct, 0, 0);
	printf("\ng = \n");
	blasfeo_print_dvec(nv, qp.gz, 0);
	printf("\nb = \n");
	blasfeo_print_dvec(ne, qp.b, 0);
	printf("\nd = \n");
	blasfeo_print_dvec(2*nb+2*ng, qp.d, 0);
#endif

/************************************************
* dense qp sol
************************************************/

	int qp_sol_size = d_memsize_dense_qp_sol(&qp_dim);
	printf("\nqp sol size = %d\n", qp_sol_size);
	void *qp_sol_mem = malloc(qp_sol_size);

	struct d_dense_qp_sol qp_sol;
	d_create_dense_qp_sol(&qp_dim, &qp_sol, qp_sol_mem);

/************************************************
* ipm arg
************************************************/

	int ipm_arg_size = d_memsize_dense_qp_ipm_arg(&qp_dim);
	printf("\nipm arg size = %d\n", ipm_arg_size);
	void *ipm_arg_mem = malloc(ipm_arg_size);

	struct d_dense_qp_ipm_arg arg;
	d_create_dense_qp_ipm_arg(&qp_dim, &arg, ipm_arg_mem);
//	enum dense_qp_ipm_mode mode = SPEED_ABS;
	enum dense_qp_ipm_mode mode = SPEED;
//	enum dense_qp_ipm_mode mode = BALANCE;
//	enum dense_qp_ipm_mode mode = ROBUST;
	d_set_default_dense_qp_ipm_arg(mode, &arg);

//	arg.alpha_min = 1e-8;
//	arg.res_g_max = 1e-8;
//	arg.res_b_max = 1e-8;
//	arg.res_d_max = 1e-12;
//	arg.res_m_max = 1e-12;
//	arg.mu0 = 10.0;
//	arg.iter_max = 10;
//	arg.stat_max = 10;
//	arg.pred_corr = 1;
//	arg.scale = 1;

/************************************************
* ipm
************************************************/

	int ipm_size = d_memsize_dense_qp_ipm(&qp_dim, &arg);
	printf("\nipm size = %d\n", ipm_size);
	void *ipm_mem = malloc(ipm_size);

	struct d_dense_qp_ipm_workspace workspace;
	d_create_dense_qp_ipm(&qp_dim, &arg, &workspace, ipm_mem);

	int rep, nrep=1000;

	int hpipm_return; // 0 normal; 1 max iter

	struct timeval tv0, tv1;

	gettimeofday(&tv0, NULL); // start

	for(rep=0; rep<nrep; rep++)
		{
		hpipm_return = d_solve_dense_qp_ipm(&qp, &qp_sol, &arg, &workspace);
		}

	gettimeofday(&tv1, NULL); // stop

	double time_dense_ipm = (tv1.tv_sec-tv0.tv_sec)/(nrep+0.0)+(tv1.tv_usec-tv0.tv_usec)/(nrep*1e6);

	printf("\nsolution\n\n");
	printf("\nv\n");
	d_print_mat(1, nv, qp_sol.v->pa, 1);
	printf("\nls\n");
	d_print_mat(1, ns, qp_sol.v->pa+nv, 1);
	printf("\nus\n");
	d_print_mat(1, ns, qp_sol.v->pa+nv+ns, 1);
	printf("\npi\n");
	d_print_mat(1, ne, qp_sol.pi->pa, 1);
	printf("\nlam_lb\n");
	d_print_mat(1, nb, qp_sol.lam->pa+0, 1);
	printf("\nlam_ub\n");
	d_print_mat(1, nb, qp_sol.lam->pa+nb+ng, 1);
	printf("\nlam_lg\n");
	d_print_mat(1, ng, qp_sol.lam->pa+nb, 1);
	printf("\nlam_ug\n");
	d_print_mat(1, ng, qp_sol.lam->pa+2*nb+ng, 1);
	printf("\nlam_ls\n");
	d_print_mat(1, ns, qp_sol.lam->pa+2*nb+2*ng, 1);
	printf("\nlam_us\n");
	d_print_mat(1, ns, qp_sol.lam->pa+2*nb+2*ng+ns, 1);
	printf("\nt_lb\n");
	d_print_mat(1, nb, qp_sol.t->pa+0, 1);
	printf("\nt_ub\n");
	d_print_mat(1, nb, qp_sol.t->pa+nb+ng, 1);
	printf("\nt_lg\n");
	d_print_mat(1, ng, qp_sol.t->pa+nb, 1);
	printf("\nt_ug\n");
	d_print_mat(1, ng, qp_sol.t->pa+2*nb+ng, 1);
	printf("\nt_ls\n");
	d_print_mat(1, ns, qp_sol.t->pa+2*nb+2*ng, 1);
	printf("\nt_us\n");
	d_print_mat(1, ns, qp_sol.t->pa+2*nb+2*ng+ns, 1);

	printf("\nresiduals\n\n");
	printf("\nres_g\n");
	d_print_e_mat(1, nv+2*ns, workspace.res->res_g->pa, 1);
	printf("\nres_b\n");
	d_print_e_mat(1, ne, workspace.res->res_b->pa, 1);
	printf("\nres_d\n");
	d_print_e_mat(1, 2*nb+2*ng+2*ns, workspace.res->res_d->pa, 1);
	printf("\nres_m\n");
	d_print_e_mat(1, 2*nb+2*ng+2*ns, workspace.res->res_m->pa, 1);
	printf("\nres_mu\n");
	printf("\n%e\n\n", workspace.res->res_mu);

/************************************************
* print ipm statistics
************************************************/

	printf("\nipm return = %d\n", hpipm_return);
	printf("\nipm residuals max: res_g = %e, res_b = %e, res_d = %e, res_m = %e\n", workspace.qp_res[0], workspace.qp_res[1], workspace.qp_res[2], workspace.qp_res[3]);

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


