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
#include "../include/hpipm_d_ipm_hard_dense_qp.h"

#include "tools.h"



#define PRINT 1



int main()
	{

	int ii;

/************************************************
* problem dimension and data
************************************************/	

	int nv = 2;
	int ne = 0;
	int nb = 0;
	int ng = 2;

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
* problem struct
************************************************/	

	int qp_size = d_memsize_dense_qp(nv, ne, nb, ng);
	printf("\nqp size = %d\n", qp_size);
	void *qp_mem = malloc(qp_size);

	struct d_dense_qp qp;
	d_create_dense_qp(nv, ne, nb, ng, &qp, qp_mem);
	d_cvt_colmaj_to_dense_qp(H, g, A, b, idxb, d_lb, d_ub, C, d_lg, d_ug, &qp);

	d_print_strmat(nv, nv, qp.H, 0, 0);
	d_print_strmat(ne, nv, qp.A, 0, 0);
	d_print_strmat(nv, ng, qp.Ct, 0, 0);
	d_print_strvec(nv, qp.g, 0);
	d_print_strvec(ne, qp.b, 0);
	d_print_strvec(2*nb+2*ng, qp.d, 0);
	d_print_strvec(nb, qp.d_lb, 0);
	d_print_strvec(nb, qp.d_ub, 0);
	d_print_strvec(ng, qp.d_lg, 0);
	d_print_strvec(ng, qp.d_ug, 0);

/************************************************
* ipm
************************************************/	

	struct d_ipm_hard_dense_qp_arg arg;
	arg.alpha_min = 1e-8;
	arg.mu_max = 1e-12;
	arg.iter_max = 10;
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
		d_solve_ipm_hard_dense_qp(&qp, &workspace);
		}

	gettimeofday(&tv1, NULL); // stop

	double time0 = (tv1.tv_sec-tv0.tv_sec)/(nrep+0.0)+(tv1.tv_usec-tv0.tv_usec)/(nrep*1e6);

	gettimeofday(&tv0, NULL); // start

	for(rep=0; rep<nrep; rep++)
		{
		d_solve_ipm_hard_dense_qp(&qp, &workspace);
		}

	gettimeofday(&tv1, NULL); // stop

	double time1 = (tv1.tv_sec-tv0.tv_sec)/(nrep+0.0)+(tv1.tv_usec-tv0.tv_usec)/(nrep*1e6);

	printf("\nsol time = %e %e [s]\n\n", time0, time1);

	printf("\nsolution\n\n");
	printf("\nv\n");
	d_print_tran_strvec(nv, workspace.v, 0);
	printf("\npi\n");
	d_print_tran_strvec(ne, workspace.pi, 0);
	printf("\nlam\n");
	d_print_tran_strvec(2*nb+2*ng, workspace.lam, 0);
	printf("\nt\n");
	d_print_tran_strvec(2*nb+2*ng, workspace.t, 0);

	printf("\nresiduals\n\n");
	printf("\nres_g\n");
	d_print_e_tran_strvec(nv, workspace.res_g, 0);
	printf("\nres_b\n");
	d_print_e_tran_strvec(ne, workspace.res_b, 0);
	printf("\nres_d\n");
	d_print_e_tran_strvec(2*nb+2*ng, workspace.res_d, 0);
	printf("\nres_m\n");
	d_print_e_tran_strvec(2*nb+2*ng, workspace.res_m, 0);
	printf("\nres_mu\n");
	printf("\n%e\n\n", workspace.res_mu);

/************************************************
* free memory
************************************************/	

	free(qp_mem);
	free(ipm_mem);

/************************************************
* return
************************************************/	

	return 0;

	}

	
