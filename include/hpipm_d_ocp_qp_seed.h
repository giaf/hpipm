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

#ifndef HPIPM_D_OCP_QP_SEED_H_
#define HPIPM_D_OCP_QP_SEED_H_



#include <blasfeo_target.h>
#include <blasfeo_common.h>

#include <hpipm_common.h>
#include <hpipm_d_ocp_qp_dim.h>



#ifdef __cplusplus
extern "C" {
#endif



struct d_ocp_qp_seed
	{
	struct d_ocp_qp_dim *dim;
	struct blasfeo_dvec *seed_g;
	struct blasfeo_dvec *seed_b;
	struct blasfeo_dvec *seed_d;
	struct blasfeo_dvec *seed_m;
	hpipm_size_t memsize;
	};



//
hpipm_size_t d_ocp_qp_seed_strsize();
//
hpipm_size_t d_ocp_qp_seed_memsize(struct d_ocp_qp_dim *ocp_dim);
//
void d_ocp_qp_seed_create(struct d_ocp_qp_dim *ocp_dim, struct d_ocp_qp_seed *seed, void *mem);
//
void d_ocp_qp_seed_set_zero(struct d_ocp_qp_seed *seed);
//
void d_ocp_qp_seed_set(char *field, int stage, double *vec, struct d_ocp_qp_seed *qp_seed);
//
void d_ocp_qp_seed_set_seed_r(int stage, double *vec, struct d_ocp_qp_seed *qp_seed);
//
void d_ocp_qp_seed_set_seed_q(int stage, double *vec, struct d_ocp_qp_seed *qp_seed);
//
void d_ocp_qp_seed_set_seed_zl(int stage, double *vec, struct d_ocp_qp_seed *qp_seed);
//
void d_ocp_qp_seed_set_seed_zu(int stage, double *vec, struct d_ocp_qp_seed *qp_seed);
//
void d_ocp_qp_seed_set_seed_b(int stage, double *vec, struct d_ocp_qp_seed *qp_seed);
//
void d_ocp_qp_seed_set_seed_lb(int stage, double *vec, struct d_ocp_qp_seed *qp_seed);
//
void d_ocp_qp_seed_set_seed_lbu(int stage, double *vec, struct d_ocp_qp_seed *qp_seed);
//
void d_ocp_qp_seed_set_seed_lbx(int stage, double *vec, struct d_ocp_qp_seed *qp_seed);
//
void d_ocp_qp_seed_set_seed_ub(int stage, double *vec, struct d_ocp_qp_seed *qp_seed);
//
void d_ocp_qp_seed_set_seed_ubu(int stage, double *vec, struct d_ocp_qp_seed *qp_seed);
//
void d_ocp_qp_seed_set_seed_ubx(int stage, double *vec, struct d_ocp_qp_seed *qp_seed);
//
void d_ocp_qp_seed_set_seed_lg(int stage, double *vec, struct d_ocp_qp_seed *qp_seed);
//
void d_ocp_qp_seed_set_seed_ug(int stage, double *vec, struct d_ocp_qp_seed *qp_seed);
//
void d_ocp_qp_seed_set_seed_ls(int stage, double *vec, struct d_ocp_qp_seed *qp_seed);
//
void d_ocp_qp_seed_set_seed_us(int stage, double *vec, struct d_ocp_qp_seed *qp_seed);




#ifdef __cplusplus
}	// #extern "C"
#endif


#endif // HPIPM_D_OCP_QP_SEED_H_

