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

#include "../include/hpipm_d_dense_kkt.h"

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
	int nc = 0;

	double H[] = {4, 1, 1, 2};
	double g[] = {1, 1};
	double A[] = {1, 1};
	double be[] = {1};
	double lb[] = {0, 0};
	double ub[] = {INFINITY, INFINITY};
	double Ct[] = {};
	double lc[] = {};
	double uc[] = {};

/************************************************
* data struct
************************************************/	

	struct d_strmat sH;
	d_allocate_strmat(nv, nv, &sH);
	d_cvt_mat2strmat(nv, nv, H, 2, &sH, 0, 0);

	struct d_strmat sA;
	d_allocate_strmat(ne, nv, &sA);
	d_cvt_mat2strmat(ne, nv, A, 1, &sA, 0, 0);

	struct d_strmat sCt;
	d_allocate_strmat(nv, nc, &sCt);
	d_cvt_mat2strmat(nv, nc, Ct, 0, &sCt, 0, 0);

	struct d_strvec sg;
	d_allocate_strvec(nv, &sg);
	d_cvt_vec2strvec(nv, g, &sg, 0);

	struct d_strvec sbe;
	d_allocate_strvec(ne, &sbe);
	d_cvt_vec2strvec(ne, be, &sbe, 0);

	struct d_strvec slb;
	d_allocate_strvec(nb, &slb);
	d_cvt_vec2strvec(nb, lb, &slb, 0);

	struct d_strvec sub;
	d_allocate_strvec(nb, &sub);
	d_cvt_vec2strvec(nb, ub, &sub, 0);

	struct d_strvec slc;
	d_allocate_strvec(nc, &slc);
	d_cvt_vec2strvec(nc, lc, &slc, 0);

	struct d_strvec suc;
	d_allocate_strvec(nc, &suc);
	d_cvt_vec2strvec(nc, uc, &suc, 0);

/************************************************
* problem struct
************************************************/	

	struct d_dense_qp_dim qp_dim;
	d_cast_dense_qp_dim(nv, ne, nb, nc, &qp_dim);
	printf("\n%d %d %d %d\n", qp_dim.nv, qp_dim.ne, qp_dim.nb, qp_dim.nc);

	int qp_size = d_size_dense_qp(&qp_dim);
	printf("\nqp size = %d\n", qp_size);
	void *qp_mem = malloc(qp_size);

#if 0
	void *ptr = qp_mem;

	struct d_dense_qp_vec qp_vec;
	d_create_dense_qp_vec(&qp_dim, &qp_vec, ptr);
	ptr += d_size_dense_qp_vec(&qp_dim);
	d_init_dense_qp_vec(&qp_dim, &sg, &sbe, &slb, &sub, &slc, &suc, &qp_vec);
	d_print_strvec(nv, qp_vec.g, 0);
	d_print_strvec(ne, qp_vec.be, 0);
	d_print_strvec(nb, qp_vec.lb, 0);
	d_print_strvec(nb, qp_vec.ub, 0);
	d_print_strvec(nc, qp_vec.lc, 0);
	d_print_strvec(nc, qp_vec.uc, 0);

	struct d_dense_qp_mat qp_mat;
	d_create_dense_qp_mat(&qp_dim, &qp_mat, ptr);
	ptr += d_size_dense_qp_mat(&qp_dim);
	d_init_dense_qp_mat(&qp_dim, &sH, &sA, &sCt, &qp_mat);
	d_print_strmat(nv, nv, qp_mat.H, 0, 0);
	d_print_strmat(ne, nv, qp_mat.A, 0, 0);
	d_print_strmat(nv, nc, qp_mat.Ct, 0, 0);
#else
	struct d_dense_qp qp;
	d_create_dense_qp(&qp_dim, &qp, qp_mem);
	d_init_dense_qp(&qp_dim, &sg, &sbe, &slb, &sub, &slc, &suc, &sH, &sA, &sCt, &qp);

	d_print_strmat(nv, nv, qp.mat->H, 0, 0);
	d_print_strmat(ne, nv, qp.mat->A, 0, 0);
	d_print_strmat(nv, nc, qp.mat->Ct, 0, 0);
	d_print_strvec(nv, qp.vec->g, 0);
	d_print_strvec(ne, qp.vec->be, 0);
	d_print_strvec(nb, qp.vec->lb, 0);
	d_print_strvec(nb, qp.vec->ub, 0);
	d_print_strvec(nc, qp.vec->lc, 0);
	d_print_strvec(nc, qp.vec->uc, 0);
#endif

/************************************************
* free memory
************************************************/	

	free(qp_mem);

/************************************************
* return
************************************************/	

	return 0;

	}

	
