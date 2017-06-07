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




struct d_ipm_hard_ocp_qp_workspace
	{
	struct d_ipm_hard_core_qp_workspace *core_workspace;
	struct d_strvec *ux;
	struct d_strvec *pi;
	struct d_strvec *lam_lb;
	struct d_strvec *lam_ub;
	struct d_strvec *lam_lg;
	struct d_strvec *lam_ug;
	struct d_strvec *t_lb;
	struct d_strvec *t_ub;
	struct d_strvec *t_lg;
	struct d_strvec *t_ug;
	double nt_inv; // 1.0/nt, where nt is the total number of constraints
	double mu0;
	};



struct d_ipm_hard_ocp_qp_arg
	{
	double alpha_min; // exit cond on step length
	double mu_max; // exit cond on duality measure
	double mu0; // initial value for duality measure
	int iter_max; // exit cond in iter number
	};



//
int d_memsize_ipm_hard_ocp_qp(struct d_ocp_qp *qp, struct d_ipm_hard_ocp_qp_arg *arg);

