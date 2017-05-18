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



#include <blasfeo_target.h>
#include <blasfeo_common.h>



struct s_ocp_qp
	{
	int NN; // hotizon lenght
	int *nx; // number of states
	int *nu; // number of inputs
	int *nb; // number of box constraints
	int **idxb; // index of box constraints
	int *ng; // number of general constraints
	struct s_strmat *sBAbt;
	struct s_strvec *sb;
	struct s_strmat *sRSQrq;
	struct s_strvec *srq;
	struct s_strmat *sDCt;
	struct s_strvec *slb;
	struct s_strvec *sub;
	struct s_strvec *slg;
	struct s_strvec *sug;
	};



struct s_cond_mem
	{
	struct s_strmat *sGamma;
	struct s_strmat *sL;
	};



struct s_cond_work
	{
	int *cond_RSQrq_N2nx3_sizes;
	};



//
int s_size_ocp_qp(int N, int *nx, int *nu, int *nb, int *ng);
//
void s_create_ocp_qp(int N, int *nx, int *nu, int *nb, int *ng, struct s_ocp_qp *str_out, void *memory);
//
void s_cast_ocp_qp(int N, int *nx, int *nu, int *nb, int **idxb, int *ng, struct s_strmat *sBAbt, struct s_strvec *sb, struct s_strmat *sRSQrq, struct s_strvec *srq, struct s_strmat *sDCt, struct s_strvec *slb, struct s_strvec *sub, struct s_strvec *slg, struct s_strvec *sug, struct s_ocp_qp *str_out);
//
void s_copy_ocp_qp(struct s_ocp_qp *str_in, struct s_ocp_qp *str_out);
