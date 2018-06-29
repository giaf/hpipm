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

#ifndef HPIPM_S_DENSE_QP_RES_H_
#define HPIPM_S_DENSE_QP_RES_H_



#include <blasfeo_target.h>
#include <blasfeo_common.h>



#ifdef __cplusplus
extern "C" {
#endif



struct s_dense_qp_res
	{
	struct s_dense_qp_dim *dim;
	struct blasfeo_svec *res_g; // q-residuals
	struct blasfeo_svec *res_b; // b-residuals
	struct blasfeo_svec *res_d; // d-residuals
	struct blasfeo_svec *res_m; // m-residuals
	float res_mu; // mu-residual
	int memsize;
	};



struct s_dense_qp_res_workspace
	{
	struct blasfeo_svec *tmp_nbg; // work space of size nbM+ngM
	struct blasfeo_svec *tmp_ns; // work space of size nsM
	int memsize;
	};



//
int s_memsize_dense_qp_res(struct s_dense_qp_dim *dim);
//
void s_create_dense_qp_res(struct s_dense_qp_dim *dim, struct s_dense_qp_res *res, void *mem);
//
int s_memsize_dense_qp_res_workspace(struct s_dense_qp_dim *dim);
//
void s_create_dense_qp_res_workspace(struct s_dense_qp_dim *dim, struct s_dense_qp_res_workspace *workspace, void *mem);
//
void s_cvt_dense_qp_res_to_colmaj(struct s_dense_qp_res *res, float *res_g, float *res_ls, float *res_us, float *res_b, float *res_d_lb, float *res_d_ub, float *res_d_lg, float *res_d_ug, float *res_d_ls, float *res_d_us, float *res_m_lb, float *res_m_ub, float *res_m_lg, float *res_m_ug, float *res_m_ls, float *res_m_us);
//
void s_cvt_dense_qp_res_to_rowmaj(struct s_dense_qp_res *res, float *res_g, float *res_ls, float *res_us, float *res_b, float *res_d_lb, float *res_d_ub, float *res_d_lg, float *res_d_ug, float *res_d_ls, float *res_d_us, float *res_m_lb, float *res_m_ub, float *res_m_lg, float *res_m_ug, float *res_m_ls, float *res_m_us);



#ifdef __cplusplus
}	// #extern "C"
#endif


#endif // HPIPM_D_DENSE_QP_RES_H_



