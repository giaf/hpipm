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
#include "../include/hpipm_d_dense_qp_utils.h"



int main()
	{

	int ii;

/************************************************
* qp dimension and data
************************************************/

	int nv = 2;
	int ne = 0;
	int nb = 0;
	int ng = 0;
	int nq = 1;
	int ns = 0;
	int nsb = 0;
	int nsg = 0;

	double H[] = {1.0, 0.0, 0.0, 1.0};
	double g[] = {0.0, 1.0};
	double Hq[] = {1.0, 0.0, 0.0, 1.0};
	double gq[] = {0.0, 0.0};
	double uq[] = {1.0};

/************************************************
* dense qp dim
************************************************/

	int dense_qp_dim_size = d_dense_qp_dim_memsize();
	printf("\nqp dim size = %d\n", dense_qp_dim_size);
	void *qp_dim_mem = malloc(dense_qp_dim_size);

	struct d_dense_qp_dim qp_dim;
	d_dense_qp_dim_create(&qp_dim, qp_dim_mem);

	d_dense_qp_dim_set("nv", nv, &qp_dim);
	d_dense_qp_dim_set("nq", nq, &qp_dim);

	d_dense_qp_dim_print(&qp_dim);

/************************************************
* dense qp
************************************************/

	int qp_size = d_dense_qp_memsize(&qp_dim);
	printf("\nqp size = %d\n", qp_size);
	void *qp_mem = malloc(qp_size);

	struct d_dense_qp qp;
	d_dense_qp_create(&qp_dim, &qp, qp_mem);

	// test setters

	d_dense_qp_set_H(H, &qp);
	d_dense_qp_set_g(g, &qp);
	d_dense_qp_set_Hq(Hq, &qp);
	d_dense_qp_set_gq(gq, &qp);
	d_dense_qp_set_uq(uq, &qp);

	d_dense_qp_print(&qp_dim, &qp);

/************************************************
* dense qp sol
************************************************/

	int qp_sol_size = d_dense_qp_sol_memsize(&qp_dim);
	printf("\nqp sol size = %d\n", qp_sol_size);
	void *qp_sol_mem = malloc(qp_sol_size);

	struct d_dense_qp_sol qp_sol;
	d_dense_qp_sol_create(&qp_dim, &qp_sol, qp_sol_mem);

	d_dense_qp_sol_print(&qp_dim, &qp_sol);

/************************************************
* free memory
************************************************/

	free(qp_dim_mem);
	free(qp_mem);
	free(qp_sol_mem);
//	free(ipm_mem);

/************************************************
* return
************************************************/

	return 0;

	}
