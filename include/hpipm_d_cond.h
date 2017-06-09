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



struct d_cond_qp_ocp2dense_workspace
	{
	struct d_strmat *Gamma;
	struct d_strmat *L;
	struct d_strmat *Lx;
	struct d_strmat *AL;
	};



struct d_cond_mem
	{
	struct d_strmat *Gamma;
	struct d_strmat *L;
	};



struct d_cond_work
	{
	int *cond_RSQrq_N2nx3_sizes;
	int NN; // hotizon lenght
	int mem_size; // memory size in bytes
	};



//
void d_compute_qp_size_ocp2dense(int N, int *nx, int *nu, int *nb, int **idxb, int *ng, int *nvc, int *nec, int *nbc, int *ngc);
//
int d_memsize_cond_qp_ocp2dense(struct d_ocp_qp *ocp_qp, struct d_dense_qp *dense_qp); // XXX + args for algorithm type ???
//
void d_create_cond_qp_ocp2dense(struct d_ocp_qp *ocp_qp, struct d_dense_qp *dense_qp, struct d_cond_qp_ocp2dense_workspace *cond_ws, void *mem);
//
void d_cond_qp_ocp2dense(struct d_ocp_qp *ocp_qp, struct d_dense_qp *dense_qp, struct d_cond_qp_ocp2dense_workspace *cond_ws);
