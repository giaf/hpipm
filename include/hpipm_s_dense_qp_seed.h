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

#ifndef HPIPM_S_DENSE_QP_SEED_H_
#define HPIPM_S_DENSE_QP_SEED_H_



#include <blasfeo_target.h>
#include <blasfeo_common.h>

#include <hpipm_s_dense_qp_dim.h>



#ifdef __cplusplus
extern "C" {
#endif



struct s_dense_qp_seed
	{
	struct s_dense_qp_dim *dim;
	struct blasfeo_svec *seed_g;
	struct blasfeo_svec *seed_b;
	struct blasfeo_svec *seed_d;
	struct blasfeo_svec *seed_m;
	hpipm_size_t memsize;
	};



//
hpipm_size_t s_dense_qp_seed_memsize(struct s_dense_qp_dim *dim);
//
void s_dense_qp_seed_create(struct s_dense_qp_dim *dim, struct s_dense_qp_seed *seed, void *mem);
//
void s_dense_qp_seed_set(char *field, void *value, struct s_dense_qp_seed *qp_seed);
//
void s_dense_qp_seed_set_seed_g(float *seed_g, struct s_dense_qp_seed *seed);
//
void s_dense_qp_seed_set_seed_zl(float *seed_zl, struct s_dense_qp_seed *seed);
//
void s_dense_qp_seed_set_seed_zu(float *seed_zu, struct s_dense_qp_seed *seed);
//
void s_dense_qp_seed_set_seed_b(float *seed_b, struct s_dense_qp_seed *seed);
//
void s_dense_qp_seed_set_seed_lb(float *seed_lb, struct s_dense_qp_seed *seed);
//
void s_dense_qp_seed_set_seed_ub(float *seed_ub, struct s_dense_qp_seed *seed);
//
void s_dense_qp_seed_set_seed_lg(float *seed_lg, struct s_dense_qp_seed *seed);
//
void s_dense_qp_seed_set_seed_ug(float *seed_ug, struct s_dense_qp_seed *seed);
//
void s_dense_qp_seed_set_seed_ls(float *seed_ls, struct s_dense_qp_seed *seed);
//
void s_dense_qp_seed_set_seed_us(float *seed_us, struct s_dense_qp_seed *seed);




#ifdef __cplusplus
}	// #extern "C"
#endif


#endif // HPIPM_S_DENSE_QP_SEED_H_




