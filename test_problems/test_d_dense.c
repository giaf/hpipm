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
	double be[] = {1.0};
	double lb[] = {0.0, 0.0};
	double ub[] = {INFINITY, INFINITY};
	int idxb[] = {0, 1};
	double C[] = {};
	double lg[] = {};
	double ug[] = {};

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

	struct d_strvec sbe;
	d_allocate_strvec(ne, &sbe);
	d_cvt_vec2strvec(ne, be, &sbe, 0);

	struct d_strvec slb;
	d_allocate_strvec(nb, &slb);
	d_cvt_vec2strvec(nb, lb, &slb, 0);

	struct d_strvec sub;
	d_allocate_strvec(nb, &sub);
	d_cvt_vec2strvec(nb, ub, &sub, 0);

	struct d_strvec slg;
	d_allocate_strvec(ng, &slg);
	d_cvt_vec2strvec(ng, lg, &slg, 0);

	struct d_strvec sug;
	d_allocate_strvec(ng, &sug);
	d_cvt_vec2strvec(ng, ug, &sug, 0);

/************************************************
* problem struct
************************************************/	

	int qp_size = d_size_dense_qp(nv, ne, nb, ng);
	printf("\nqp size = %d\n", qp_size);
	void *qp_mem = malloc(qp_size);

	struct d_dense_qp qp;
	d_create_dense_qp(nv, ne, nb, ng, &qp, qp_mem);
	d_init_dense_qp(&sH, &sA, &sC, &sg, &sbe, &slb, &sub, &slg, &sug, idxb, &qp);

	d_print_strmat(nv, nv, qp.H, 0, 0);
	d_print_strmat(ne, nv, qp.A, 0, 0);
	d_print_strmat(nv, ng, qp.Ct, 0, 0);
	d_print_strvec(nv, qp.g, 0);
	d_print_strvec(ne, qp.be, 0);
	d_print_strvec(nb, qp.lb, 0);
	d_print_strvec(nb, qp.ub, 0);
	d_print_strvec(ng, qp.lg, 0);
	d_print_strvec(ng, qp.ug, 0);

/************************************************
* free memory
************************************************/	

	free(qp_mem);

/************************************************
* return
************************************************/	

	return 0;

	}

	
