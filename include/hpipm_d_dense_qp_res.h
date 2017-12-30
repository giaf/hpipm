/**************************************************************************************************
*                                                                                                 *
* This file is part of HPIPM.                                                                     *
*                                                                                                 *
* HPIPM -- High Performance Interior Point Method.                                                *
* Copyright (C) 2017 by Gianluca Frison.                                                          *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
*                                                                                                 *
* HPIPM is free software; you can redistribute it and/or                                          *
* modify it under the terms of the GNU Lesser General Public                                      *
* License as published by the Free Software Foundation; either                                    *
* version 2.1 of the License, or (at your option) any later version.                              *
*                                                                                                 *
* HPIPM is distributed in the hope that it will be useful,                                        *
* but WITHOUT ANY WARRANTY; without even the implied warranty of                                  *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                                            *
* See the GNU Lesser General Public License for more details.                                     *
*                                                                                                 *
* You should have received a copy of the GNU Lesser General Public                                *
* License along with HPIPM; if not, write to the Free Software                                    *
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA                  *
*                                                                                                 *
* Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             *
*                                                                                                 *
**************************************************************************************************/

#ifndef HPIPM_D_DENSE_QP_RES_H_
#define HPIPM_D_DENSE_QP_RES_H_



#include <blasfeo_target.h>
#include <blasfeo_common.h>



#ifdef __cplusplus
extern "C" {
#endif



struct d_dense_qp_res
	{
	struct d_dense_qp_dim *dim;
	struct blasfeo_dvec *res_g; // q-residuals
	struct blasfeo_dvec *res_b; // b-residuals
	struct blasfeo_dvec *res_d; // d-residuals
	struct blasfeo_dvec *res_m; // m-residuals
	double res_mu; // mu-residual
	int memsize;
	};



struct d_dense_qp_res_workspace
	{
	struct blasfeo_dvec *tmp_nbg; // work space of size nbM+ngM
	struct blasfeo_dvec *tmp_ns; // work space of size nsM
	int memsize;
	};



//
int d_memsize_dense_qp_res(struct d_dense_qp_dim *dim);
//
void d_create_dense_qp_res(struct d_dense_qp_dim *dim, struct d_dense_qp_res *res, void *mem);
//
int d_memsize_dense_qp_res_workspace(struct d_dense_qp_dim *dim);
//
void d_create_dense_qp_res_workspace(struct d_dense_qp_dim *dim, struct d_dense_qp_res_workspace *workspace, void *mem);
//
void d_cvt_dense_qp_res_to_colmaj(struct d_dense_qp_res *res, double *res_g, double *res_ls, double *res_us, double *res_b, double *res_d_lb, double *res_d_ub, double *res_d_lg, double *res_d_ug, double *res_d_ls, double *res_d_us, double *res_m_lb, double *res_m_ub, double *res_m_lg, double *res_m_ug, double *res_m_ls, double *res_m_us);
//
void d_cvt_dense_qp_res_to_rowmaj(struct d_dense_qp_res *res, double *res_g, double *res_ls, double *res_us, double *res_b, double *res_d_lb, double *res_d_ub, double *res_d_lg, double *res_d_ug, double *res_d_ls, double *res_d_us, double *res_m_lb, double *res_m_ub, double *res_m_lg, double *res_m_ug, double *res_m_ls, double *res_m_us);



#ifdef __cplusplus
}	// #extern "C"
#endif


#endif // HPIPM_D_DENSE_QP_RES_H_


