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



#include <stdlib.h>
#include <stdio.h>

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_m_aux.h>

#include "../include/hpipm_d_dense_qp_dim.h"
#include "../include/hpipm_d_dense_qp.h"
#include "../include/hpipm_s_dense_qp_dim.h"
#include "../include/hpipm_s_dense_qp.h"



void cvt_d2s_dense_qp(struct d_dense_qp *qpd, struct s_dense_qp *qps)
	{

	int ii;

	int nv = qpd->dim->nv;
	int ne = qpd->dim->ne;
	int nb = qpd->dim->nb;
	int ng = qpd->dim->ng;
	int ns = qpd->dim->ns;

	blasfeo_cvt_d2s_mat(nv, nv, qpd->Hv, 0, 0, qps->Hv, 0, 0);
	blasfeo_cvt_d2s_vec(2*ns, qpd->Z, 0, qps->Z, 0);
	blasfeo_cvt_d2s_vec(nv+2*ns, qpd->gz, 0, qps->gz, 0);
	blasfeo_cvt_d2s_mat(ne, nv, qpd->A, 0, 0, qps->A, 0, 0);
	blasfeo_cvt_d2s_vec(ne, qpd->b, 0, qps->b, 0);
	blasfeo_cvt_d2s_mat(nv, ng, qpd->Ct, 0, 0, qps->Ct, 0, 0);
	blasfeo_cvt_d2s_vec(2*nb+2*ng+2*ns, qpd->d, 0, qps->d, 0);
	blasfeo_cvt_d2s_vec(2*nb+2*ng+2*ns, qpd->m, 0, qps->m, 0);
	for(ii=0; ii<nb; ii++) qps->idxb[ii] = qpd->idxb[ii];
	for(ii=0; ii<nb; ii++) qps->idxs[ii] = qpd->idxs[ii];

	return;

	}
