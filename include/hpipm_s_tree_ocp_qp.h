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

#ifndef HPIPM_S_TREE_OCP_QP_H_
#define HPIPM_S_TREE_OCP_QP_H_



#include <blasfeo_target.h>
#include <blasfeo_common.h>

#include "hpipm_s_tree_ocp_qp_dim.h"



#ifdef __cplusplus
extern "C" {
#endif



struct s_tree_ocp_qp
	{
	struct s_tree_ocp_qp_dim *dim;
	struct blasfeo_smat *BAbt; // Nn-1
	struct blasfeo_smat *RSQrq; // Nn
	struct blasfeo_smat *DCt; // Nn
	struct blasfeo_svec *b; // Nn-1
	struct blasfeo_svec *rqz; // Nn
	struct blasfeo_svec *d; // Nn
	struct blasfeo_svec *m; // Nn
	struct blasfeo_svec *Z; // Nn
	int **idxb; // index of box constraints // Nn
	int **idxs; // index of soft constraints
	int memsize; // memory size in bytes
	};



//
int s_memsize_tree_ocp_qp(struct s_tree_ocp_qp_dim *dim);
//
void s_create_tree_ocp_qp(struct s_tree_ocp_qp_dim *dim, struct s_tree_ocp_qp *qp, void *memory);
//
void s_cvt_colmaj_to_tree_ocp_qp(float **A, float **B, float **b, float **Q, float **S, float **R, float **q, float **r, int **idxb, float **d_lb, float **d_ub, float **C, float **D, float **d_lg, float **d_ug, float **Zl, float **Zu, float **zl, float **zu, int **idxs, float **d_ls, float **d_us, struct s_tree_ocp_qp *qp);



#ifdef __cplusplus
} /* extern "C" */
#endif



#endif // HPIPM_S_TREE_OCP_QP_H_
