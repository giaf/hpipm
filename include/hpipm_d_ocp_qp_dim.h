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

#ifndef HPIPM_D_OCP_QP_DIM_H_
#define HPIPM_D_OCP_QP_DIM_H_

#include "hpipm_common.h"

#ifdef __cplusplus
extern "C" {
#endif



struct d_ocp_qp_dim
	{
	int *nx; // number of states
	int *nu; // number of inputs
	int *nb; // number of (two-sided) box constraints
	int *nbx; // number of (two-sided) state box constraints
	int *nbu; // number of (two-sided) input box constraints
	int *ng; // number of (two-sided) general constraints
	int *ns; // number of soft constraints
	int *nsbx; // number of (two-sided) soft state box constraints
	int *nsbu; // number of (two-sided) soft input box constraints
	int *nsg; // number of (two-sided) soft general constraints
	int *nbxe; // number of state box constraints which are equality
	int *nbue; // number of input box constraints which are equality
	int *nge; // number of general constraints which are equality
	int N; // horizon length
	hpipm_size_t memsize;
	};



//
hpipm_size_t d_ocp_qp_dim_strsize();
//
hpipm_size_t d_ocp_qp_dim_memsize(int N);
//
void d_ocp_qp_dim_create(int N, struct d_ocp_qp_dim *qp_dim, void *memory);
//
void d_ocp_qp_dim_copy_all(struct d_ocp_qp_dim *dim_orig, struct d_ocp_qp_dim *dim_dest);
//
void d_ocp_qp_dim_set_all(int *nx, int *nu, int *nbx, int *nbu, int *ng, int *nsbx, int *nsbu, int *nsg, struct d_ocp_qp_dim *dim);
//
void d_ocp_qp_dim_set(char *field, int stage, int value, struct d_ocp_qp_dim *dim);
//
void d_ocp_qp_dim_set_nx(int stage, int value, struct d_ocp_qp_dim *dim);
//
void d_ocp_qp_dim_set_nu(int stage, int value, struct d_ocp_qp_dim *dim);
//
void d_ocp_qp_dim_set_nbx(int stage, int value, struct d_ocp_qp_dim *dim);
//
void d_ocp_qp_dim_set_nbu(int stage, int value, struct d_ocp_qp_dim *dim);
//
void d_ocp_qp_dim_set_ng(int stage, int value, struct d_ocp_qp_dim *dim);
//
void d_ocp_qp_dim_set_ns(int stage, int value, struct d_ocp_qp_dim *dim);
//
void d_ocp_qp_dim_set_nsbx(int stage, int value, struct d_ocp_qp_dim *dim);
//
void d_ocp_qp_dim_set_nsbu(int stage, int value, struct d_ocp_qp_dim *dim);
//
void d_ocp_qp_dim_set_nsg(int stage, int value, struct d_ocp_qp_dim *dim);
//
void d_ocp_qp_dim_set_nbxe(int stage, int value, struct d_ocp_qp_dim *dim);
//
void d_ocp_qp_dim_set_nbue(int stage, int value, struct d_ocp_qp_dim *dim);
//
void d_ocp_qp_dim_set_nge(int stage, int value, struct d_ocp_qp_dim *dim);
//
void d_ocp_qp_dim_deepcopy(struct d_ocp_qp_dim *dim_s, struct d_ocp_qp_dim *dim_d);
//
void d_ocp_qp_dim_get(struct d_ocp_qp_dim *dim, char *field, int stage, int *value);
//
void d_ocp_qp_dim_get_N(struct d_ocp_qp_dim *dim, int *value);
//
void d_ocp_qp_dim_get_nx(struct d_ocp_qp_dim *dim, int stage, int *value);
//
void d_ocp_qp_dim_get_nu(struct d_ocp_qp_dim *dim, int stage, int *value);
//
void d_ocp_qp_dim_get_nbx(struct d_ocp_qp_dim *dim, int stage, int *value);
//
void d_ocp_qp_dim_get_nbu(struct d_ocp_qp_dim *dim, int stage, int *value);
//
void d_ocp_qp_dim_get_ng(struct d_ocp_qp_dim *dim, int stage, int *value);
//
void d_ocp_qp_dim_get_ns(struct d_ocp_qp_dim *dim, int stage, int *value);
//
void d_ocp_qp_dim_get_nsbx(struct d_ocp_qp_dim *dim, int stage, int *value);
//
void d_ocp_qp_dim_get_nsbu(struct d_ocp_qp_dim *dim, int stage, int *value);
//
void d_ocp_qp_dim_get_nsg(struct d_ocp_qp_dim *dim, int stage, int *value);
//
void d_ocp_qp_dim_get_nbxe(struct d_ocp_qp_dim *dim, int stage, int *value);
//
void d_ocp_qp_dim_get_nbue(struct d_ocp_qp_dim *dim, int stage, int *value);
//
void d_ocp_qp_dim_get_nge(struct d_ocp_qp_dim *dim, int stage, int *value);




#ifdef __cplusplus
}	// #extern "C"
#endif



#endif // HPIPM_D_OCP_QP_DIM_H_
