/**************************************************************************************************
*                                                                                                 *
* This file is part of HPIPM.                                                                     *
*                                                                                                 *
* HPIPM -- High-Performance Interior Point Method.                                                *
* Copyright (C) 2017-2018 by Gianluca Frison.                                                     *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
*                                                                                                 *
* This program is free software: you can redistribute it and/or modify                            *
* it under the terms of the GNU General Public License as published by                            *
* the Free Software Foundation, either version 3 of the License, or                               *
* (at your option) any later version                                                              *.
*                                                                                                 *
* This program is distributed in the hope that it will be useful,                                 *
* but WITHOUT ANY WARRANTY; without even the implied warranty of                                  *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                   *
* GNU General Public License for more details.                                                    *
*                                                                                                 *
* You should have received a copy of the GNU General Public License                               *
* along with this program.  If not, see <https://www.gnu.org/licenses/>.                          *
*                                                                                                 *
* The authors designate this particular file as subject to the "Classpath" exception              *
* as provided by the authors in the LICENSE file that accompained this code.                      *
*                                                                                                 *
* Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             *
*                                                                                                 *
**************************************************************************************************/



#if defined(RUNTIME_CHECKS)
#include <stdlib.h>
#endif

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_m_aux.h>

#include <hpipm_d_ocp_qp.h>
#include <hpipm_s_ocp_qp.h>



void m_cvt_d_ocp_qp_to_s_ocp_qp(struct d_ocp_qp *d_qp, struct s_ocp_qp *s_qp)
	{

	// TODO check that they have equal sizes !!!!!

	int N = d_qp->N;
	int *nx = d_qp->nx;
	int *nu = d_qp->nu;
	int *nb = d_qp->nb;
	int *ng = d_qp->ng;

	int ii, jj;

	for(ii=0; ii<N; ii++)
		{
		m_cvt_d2blasfeo_smat(nu[ii]+nx[ii]+1, nx[ii+1], d_qp->BAbt+ii, 0, 0, s_qp->BAbt+ii, 0, 0);
		m_cvt_d2blasfeo_smat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], d_qp->RSQrq+ii, 0, 0, s_qp->RSQrq+ii, 0, 0);
		m_cvt_d2blasfeo_smat(nu[ii]+nx[ii], ng[ii], d_qp->DCt+ii, 0, 0, s_qp->DCt+ii, 0, 0);
		m_cvt_d2blasfeo_svec(nx[ii+1], d_qp->b+ii, 0, s_qp->b+ii, 0);
		m_cvt_d2blasfeo_svec(nu[ii]+nx[ii], d_qp->rq+ii, 0, s_qp->rq+ii, 0);
		m_cvt_d2blasfeo_svec(nb[ii], d_qp->d_lb+ii, 0, s_qp->d_lb+ii, 0);
		m_cvt_d2blasfeo_svec(nb[ii], d_qp->d_ub+ii, 0, s_qp->d_ub+ii, 0);
		m_cvt_d2blasfeo_svec(ng[ii], d_qp->d_lg+ii, 0, s_qp->d_lg+ii, 0);
		m_cvt_d2blasfeo_svec(ng[ii], d_qp->d_ug+ii, 0, s_qp->d_ug+ii, 0);
		for(jj=0; jj<nb[ii]; jj++) s_qp->idxb[ii][jj] = d_qp->idxb[ii][jj];
		}
	ii = N;
	m_cvt_d2blasfeo_smat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], d_qp->RSQrq+ii, 0, 0, s_qp->RSQrq+ii, 0, 0);
	m_cvt_d2blasfeo_smat(nu[ii]+nx[ii], ng[ii], d_qp->DCt+ii, 0, 0, s_qp->DCt+ii, 0, 0);
	m_cvt_d2blasfeo_svec(nu[ii]+nx[ii], d_qp->rq+ii, 0, s_qp->rq+ii, 0);
	m_cvt_d2blasfeo_svec(nb[ii], d_qp->d_lb+ii, 0, s_qp->d_lb+ii, 0);
	m_cvt_d2blasfeo_svec(nb[ii], d_qp->d_ub+ii, 0, s_qp->d_ub+ii, 0);
	m_cvt_d2blasfeo_svec(ng[ii], d_qp->d_lg+ii, 0, s_qp->d_lg+ii, 0);
	m_cvt_d2blasfeo_svec(ng[ii], d_qp->d_ug+ii, 0, s_qp->d_ug+ii, 0);
	for(jj=0; jj<nb[ii]; jj++) s_qp->idxb[ii][jj] = d_qp->idxb[ii][jj];

	return;

	}
