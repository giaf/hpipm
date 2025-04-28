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

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_d_aux.h>

#include <hpipm_d_dense_qp_dim.h>
#include <hpipm_d_dense_qp_seed.h>
#include <hpipm_aux_string.h>
#include <hpipm_aux_mem.h>



#define BLASFEO_VECEL BLASFEO_DVECEL



#define CREATE_STRVEC blasfeo_create_dvec
#define DENSE_QP_DIM d_dense_qp_dim
#define DENSE_QP_SEED d_dense_qp_seed
#define PACK_VEC blasfeo_pack_dvec
#define REAL double
#define SIZE_STRVEC blasfeo_memsize_dvec
#define STRVEC blasfeo_dvec
#define VECSC blasfeo_dvecsc
#define VECSE blasfeo_dvecse



#define DENSE_QP_SEED_MEMSIZE d_dense_qp_seed_memsize
#define DENSE_QP_SEED_CREATE d_dense_qp_seed_create
#define DENSE_QP_SEED_SET d_dense_qp_seed_set
#define DENSE_QP_SEED_SET_SEED_G d_dense_qp_seed_set_seed_g
#define DENSE_QP_SEED_SET_SEED_ZL d_dense_qp_seed_set_seed_zl
#define DENSE_QP_SEED_SET_SEED_ZU d_dense_qp_seed_set_seed_zu
#define DENSE_QP_SEED_SET_SEED_B d_dense_qp_seed_set_seed_b
#define DENSE_QP_SEED_SET_SEED_LB d_dense_qp_seed_set_seed_lb
#define DENSE_QP_SEED_SET_SEED_UB d_dense_qp_seed_set_seed_ub
#define DENSE_QP_SEED_SET_SEED_LG d_dense_qp_seed_set_seed_lg
#define DENSE_QP_SEED_SET_SEED_UG d_dense_qp_seed_set_seed_ug
#define DENSE_QP_SEED_SET_SEED_LS d_dense_qp_seed_set_seed_ls
#define DENSE_QP_SEED_SET_SEED_US d_dense_qp_seed_set_seed_us



#include "x_dense_qp_seed.c"


