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
#include "../include/hpipm_d_ipm2_hard_dense_qp.h"
#include "../include/hpipm_d_ipm2_hard_revcom_qp.h"

#include "tools.h"



#define PRINT 1



int main()
	{

	int ii;

/************************************************
* problem dimension and data
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
	double d_ub[] = {2.0, 2.0};
	int idxb[] = {0, 1};
	double C[] = {};
	double d_lg[] = {};
	double d_ug[] = {};

/************************************************
* data struct
************************************************/	

	struct d_strmat sH;
	d_allocate_strmat(nv, nv, &sH);
	d_cvt_mat2strmat(nv, nv, H, 2, &sH, 0, 0);

	struct d_strmat sA;
	d_allocate_strmat(ne, nv, &sA);
	d_cvt_mat2strmat(ne, nv, A, 1, &sA, 0, 0);

	struct d_strmat sC;
	d_allocate_strmat(nv, ng, &sC);
	d_cvt_mat2strmat(nv, ng, C, 0, &sC, 0, 0);

	struct d_strvec sg;
	d_allocate_strvec(nv, &sg);
	d_cvt_vec2strvec(nv, g, &sg, 0);

	struct d_strvec sb;
	d_allocate_strvec(ne, &sb);
	d_cvt_vec2strvec(ne, b, &sb, 0);

	struct d_strvec sd_lb;
	d_allocate_strvec(nb, &sd_lb);
	d_cvt_vec2strvec(nb, d_lb, &sd_lb, 0);

	struct d_strvec sd_ub;
	d_allocate_strvec(nb, &sd_ub);
	d_cvt_vec2strvec(nb, d_ub, &sd_ub, 0);

	struct d_strvec sd_lg;
	d_allocate_strvec(ng, &sd_lg);
	d_cvt_vec2strvec(ng, d_lg, &sd_lg, 0);

	struct d_strvec sd_ug;
	d_allocate_strvec(ng, &sd_ug);
	d_cvt_vec2strvec(ng, d_ug, &sd_ug, 0);

/************************************************
* problem struct
************************************************/	

	int qp_size = d_memsize_dense_qp(nv, ne, nb, ng);
	printf("\nqp size = %d\n", qp_size);
	void *qp_mem = malloc(qp_size);

	struct d_dense_qp qp;
	d_create_dense_qp(nv, ne, nb, ng, &qp, qp_mem);
	d_init_dense_qp(&sH, &sA, &sC, &sg, &sb, &sd_lb, &sd_ub, &sd_lg, &sd_ug, idxb, &qp);

	d_print_strmat(nv, nv, qp.H, 0, 0);
	d_print_strmat(ne, nv, qp.A, 0, 0);
	d_print_strmat(nv, ng, qp.Ct, 0, 0);
	d_print_strvec(nv, qp.g, 0);
	d_print_strvec(ne, qp.b, 0);
	d_print_strvec(nb, qp.d_lb, 0);
	d_print_strvec(nb, qp.d_ub, 0);
	d_print_strvec(ng, qp.d_lg, 0);
	d_print_strvec(ng, qp.d_ug, 0);

/************************************************
* ipm
************************************************/	

	struct d_ipm2_hard_dense_qp_arg arg;
	arg.alpha_min = 1e-8;
	arg.mu_max = 1e-12;
	arg.iter_max = 10;
	arg.mu0 = 1.0;

	int ipm_size = d_memsize_ipm2_hard_dense_qp(&qp, &arg);
	printf("\nipm size = %d\n", ipm_size);
	void *ipm_mem = malloc(ipm_size);

	struct d_ipm2_hard_dense_qp_workspace workspace;
	d_create_ipm2_hard_dense_qp(&qp, &arg, &workspace, ipm_mem);

	d_solve_ipm2_hard_dense_qp(&qp, &workspace);

	printf("\n%f\n\n", workspace.revcom_workspace->sigma);
	d_print_strvec(nv, workspace.v, 0);
	d_print_strvec(ne, workspace.pi, 0);
	d_print_strvec(2*nb+2*ng, workspace.lam, 0);
	d_print_strvec(2*nb+2*ng, workspace.t, 0);

/************************************************
* free memory
************************************************/	

	d_free_strmat(&sH);
	d_free_strmat(&sA);
	d_free_strmat(&sC);
	d_free_strvec(&sg);
	d_free_strvec(&sb);
	d_free_strvec(&sd_lb);
	d_free_strvec(&sd_ub);
	d_free_strvec(&sd_lg);
	d_free_strvec(&sd_ug);
	free(qp_mem);
	free(ipm_mem);

/************************************************
* return
************************************************/	

	return 0;

	}

	
