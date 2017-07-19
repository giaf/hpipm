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



#include <math.h>

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_d_aux.h>
#include <blasfeo_s_aux.h>
#include <blasfeo_m_aux.h>
#include <blasfeo_d_blas.h>
#include <blasfeo_s_blas.h>

#include "../include/hpipm_d_ocp_qp.h"
#include "../include/hpipm_s_ocp_qp.h"
#include "../include/hpipm_d_ocp_qp_sol.h"
#include "../include/hpipm_d_ocp_qp_ipm_hard.h"
#include "../include/hpipm_m_ocp_qp_ipm_hard.h"
#include "../include/hpipm_d_core_qp_ipm_hard.h"
#include "../include/hpipm_d_core_qp_ipm_hard_aux.h"



// backward Riccati recursion
void m_fact_solve_kkt_step_hard_ocp_qp(struct d_ocp_qp *d_qp, struct s_ocp_qp *s_qp, struct m_ipm_hard_ocp_qp_workspace *ws)
	{

	int N = s_qp->N;
	int *nx = s_qp->nx;
	int *nu = s_qp->nu;
	int *nb = s_qp->nb;
	int *ng = s_qp->ng;

	struct s_strmat *BAbt = s_qp->BAbt;
	struct s_strmat *RSQrq = s_qp->RSQrq;
	struct s_strmat *DCt = s_qp->DCt;
	struct d_strmat *d_DCt = d_qp->DCt;
	int **idxb = s_qp->idxb;
	int **d_idxb = d_qp->idxb;

	struct s_strmat *L = ws->L;
	struct s_strmat *AL = ws->AL;
	struct d_strvec *d_res_b = ws->res_b;
	struct d_strvec *d_res_g = ws->res_g;
	struct s_strvec *res_b = ws->sres_b;
	struct s_strvec *res_g = ws->sres_g;
	struct d_strvec *d_dux = ws->dux;
	struct d_strvec *d_dpi = ws->dpi;
	struct s_strvec *dux = ws->sdux;
	struct s_strvec *dpi = ws->sdpi;
	struct d_strvec *d_dt_lb = ws->dt_lb;
	struct d_strvec *d_dt_lg = ws->dt_lg;
	struct d_strvec *d_Qx_lg = ws->Qx_lg;
	struct d_strvec *d_Qx_lb = ws->Qx_lb;
	struct d_strvec *d_qx_lg = ws->qx_lg;
	struct d_strvec *d_qx_lb = ws->qx_lb;
	struct s_strvec *Qx_lg = ws->sQx_lg;
	struct s_strvec *Qx_lb = ws->sQx_lb;
	struct s_strvec *qx_lg = ws->sqx_lg;
	struct s_strvec *qx_lb = ws->sqx_lb;
	struct s_strvec *Pb = ws->Pb;
	struct s_strvec *tmp_nxM = ws->tmp_nxM;
	struct s_strvec *Sx = ws->sSx;
	struct s_strvec *Si = ws->sSi;

	//
	int ii, jj;
	float d_tmp0, d_tmp1;
	float *d_ptr0, *d_ptr1;

	//
	struct d_ipm_hard_core_qp_workspace *cws = ws->core_workspace;

//	if(nb>0 | ng>0)
//		{
		d_compute_Qx_qx_hard_qp(cws);
//		}



	// cvt double => single
	for(ii=0; ii<N; ii++)
		{
		m_cvt_d2s_strvec(nu[ii]+nx[ii], d_res_g+ii, 0, res_g+ii, 0);
		m_cvt_d2s_strvec(nx[ii+1], d_res_b+ii, 0, res_b+ii, 0);
		m_cvt_d2s_strvec(nb[ii], d_Qx_lb+ii, 0, Qx_lb+ii, 0);
		m_cvt_d2s_strvec(nb[ii], d_qx_lb+ii, 0, qx_lb+ii, 0);
		m_cvt_d2s_strvec(ng[ii], d_Qx_lg+ii, 0, Qx_lg+ii, 0);
		m_cvt_d2s_strvec(ng[ii], d_qx_lg+ii, 0, qx_lg+ii, 0);
		}
	ii = N;
	m_cvt_d2s_strvec(nu[ii]+nx[ii], d_res_g+ii, 0, res_g+ii, 0);
	m_cvt_d2s_strvec(nb[ii], d_Qx_lb+ii, 0, Qx_lb+ii, 0);
	m_cvt_d2s_strvec(nb[ii], d_qx_lb+ii, 0, qx_lb+ii, 0);
	m_cvt_d2s_strvec(ng[ii], d_Qx_lg+ii, 0, Qx_lg+ii, 0);
	m_cvt_d2s_strvec(ng[ii], d_qx_lg+ii, 0, qx_lg+ii, 0);


#if 0
	for(ii=0; ii<=N; ii++)
		s_print_e_tran_strvec(nu[ii]+nx[ii], res_g+ii, 0);
	for(ii=0; ii<N; ii++)
		s_print_e_tran_strvec(nx[ii+1], res_b+ii, 0);
//	exit(1);
#endif


	// factorization and backward substitution

	// last stage
#if defined(DOUBLE_PRECISION)
	strcp_l_libstr(nu[N]+nx[N], RSQrq+N, 0, 0, L+N, 0, 0); // TODO dtrcp_l_libstr with m and n, for m>=n
#else
	sgecp_libstr(nu[N]+nx[N], nu[N]+nx[N], RSQrq+N, 0, 0, L+N, 0, 0); // TODO dtrcp_l_libstr with m and n, for m>=n
#endif
	srowin_libstr(nu[N]+nx[N], 1.0, res_g+N, 0, L+N, nu[N]+nx[N], 0);
	if(ws->scale)
		{
		if(nb[N]>0)
			{
			sdiaad_sp_libstr(nb[N], 1.0, Qx_lb+N, 0, idxb[N], L+N, 0, 0);
			srowad_sp_libstr(nb[N], 1.0, qx_lb+N, 0, idxb[N], L+N, nu[N]+nx[N], 0);
			}
		if(ng[N]>0)
			{
			sgemm_r_diag_libstr(nu[N]+nx[N], ng[N], 1.0, DCt+N, 0, 0, Qx_lg+N, 0, 0.0, AL+0, 0, 0, AL+0, 0, 0);
			srowin_libstr(ng[N], 1.0, qx_lg+N, 0, AL+0, nu[N]+nx[N], 0);
			ssyrk_ln_mn_libstr(nu[N]+nx[N]+1, nu[N]+nx[N], ng[N], 1.0, AL+0, 0, 0, DCt+N, 0, 0, 1.0, L+N, 0, 0, L+N, 0, 0);
			}
		sdiaex_libstr(nu[N]+nx[N], 1.0, L+N, 0, 0, Sx+N, 0);
		d_ptr0 = (Sx+N)->pa;
		d_ptr1 = (Si+N)->pa;
		for(jj=0; jj<nu[N]+nx[N]; jj++)
			{
			d_tmp0 = sqrt(d_ptr0[jj]);
//			d_tmp0 = sqrt(d_tmp0);
			d_tmp1 = sqrt(d_tmp0);
			d_tmp0 = sqrt(d_tmp0*d_tmp1);
			d_tmp0 = 1.0;
			d_ptr0[jj] = d_tmp0;
			d_ptr1[jj] = 1.0/d_tmp0;
			}
//		s_print_strmat(nu[N]+nx[N]+1, nu[N]+nx[N], L+N, 0, 0);
		sgemm_r_diag_libstr(nu[N]+nx[N]+1, nu[N]+nx[N], 1.0, L+N, 0, 0, Si+N, 0, 0.0, L+N, 0, 0, L+N, 0, 0);
//		s_print_strmat(nu[N]+nx[N]+1, nu[N]+nx[N], L+N, 0, 0);
		sgemm_l_diag_libstr(nu[N]+nx[N], nu[N]+nx[N], 1.0, Si+N, 0, L+N, 0, 0, 0.0, L+N, 0, 0, L+N, 0, 0);
//		s_print_strmat(nu[N]+nx[N]+1, nu[N]+nx[N], L+N, 0, 0);

		spotrf_l_mn_libstr(nu[N]+nx[N]+1, nu[N]+nx[N], L+N, 0, 0, L+N, 0, 0);
		}
	else
		{
		if(nb[N]>0)
			{
			sdiaad_sp_libstr(nb[N], 1.0, Qx_lb+N, 0, idxb[N], L+N, 0, 0);
			srowad_sp_libstr(nb[N], 1.0, qx_lb+N, 0, idxb[N], L+N, nu[N]+nx[N], 0);
			}
		if(ng[N]>0)
			{
			sgemm_r_diag_libstr(nu[N]+nx[N], ng[N], 1.0, DCt+N, 0, 0, Qx_lg+N, 0, 0.0, AL+0, 0, 0, AL+0, 0, 0);
			srowin_libstr(ng[N], 1.0, qx_lg+N, 0, AL+0, nu[N]+nx[N], 0);
			ssyrk_spotrf_ln_libstr(nu[N]+nx[N]+1, nu[N]+nx[N], ng[N], AL+0, 0, 0, DCt+N, 0, 0, L+N, 0, 0, L+N, 0, 0);
			}
		else
			{
			spotrf_l_mn_libstr(nu[N]+nx[N]+1, nu[N]+nx[N], L+N, 0, 0, L+N, 0, 0);
			}
		}

	// middle stages
	for(ii=0; ii<N; ii++)
		{
//		printf("\n%d\n", N-ii-1);
#if defined(DOUBLE_PRECISION)
		strcp_l_libstr(nu[N-ii-1]+nx[N-ii-1], RSQrq+(N-ii-1), 0, 0, L+(N-ii-1), 0, 0);
#else
		sgecp_libstr(nu[N-ii-1]+nx[N-ii-1], nu[N-ii-1]+nx[N-ii-1], RSQrq+(N-ii-1), 0, 0, L+(N-ii-1), 0, 0);
#endif
		srowin_libstr(nu[N-ii-1]+nx[N-ii-1], 1.0, res_g+(N-ii-1), 0, L+(N-ii-1), nu[N-ii-1]+nx[N-ii-1], 0);

		if(ws->scale)
			{
			if(nb[N-ii-1]>0)
				{
				sdiaad_sp_libstr(nb[N-ii-1], 1.0, Qx_lb+N-ii-1, 0, idxb[N-ii-1], L+N-ii-1, 0, 0);
				srowad_sp_libstr(nb[N-ii-1], 1.0, qx_lb+N-ii-1, 0, idxb[N-ii-1], L+N-ii-1, nu[N-ii-1]+nx[N-ii-1], 0);
				}
			if(ng[N-ii-1]>0)
				{
				sgemm_r_diag_libstr(nu[N-ii-1]+nx[N-ii-1], ng[N-ii-1], 1.0, DCt+N-ii-1, 0, 0, Qx_lg+N-ii-1, 0, 0.0, AL+1, 0, 0, AL+1, 0, 0);
				srowin_libstr(ng[N-ii-1], 1.0, qx_lg+N-ii-1, 0, AL+1, nu[N-ii-1]+nx[N-ii-1], 0);
				ssyrk_ln_mn_libstr(nu[N-ii-1]+nx[N-ii-1]+1, nu[N-ii-1]+nx[N-ii-1], ng[N-ii-1], 1.0, AL+1, 0, 0, DCt+N-ii-1, 0, 0, 1.0, L+N-ii-1, 0, 0, L+N-ii-1, 0, 0);
				}
			sdiaex_libstr(nu[N-ii-1]+nx[N-ii-1], 1.0, L+N-ii-1, 0, 0, Sx+N-ii-1, 0);
			d_ptr0 = (Sx+N-ii-1)->pa;
			d_ptr1 = (Si+N-ii-1)->pa;
			for(jj=0; jj<nu[N-ii-1]+nx[N-ii-1]; jj++)
				{
				d_tmp0 = sqrt(d_ptr0[jj]);
//				d_tmp0 = sqrt(d_tmp0);
				d_tmp1 = sqrt(d_tmp0);
				d_tmp0 = sqrt(d_tmp0*d_tmp1);
				d_tmp0 = 1.0;
				d_ptr0[jj] = d_tmp0;
				d_ptr1[jj] = 1.0/d_tmp0;
				}

			sgecp_libstr(nu[N-ii-1]+nx[N-ii-1], nx[N-ii], BAbt+(N-ii-1), 0, 0, AL, 0, 0);
			srowin_libstr(nx[N-ii], 1.0, res_b+(N-ii-1), 0, AL, nu[N-ii-1]+nx[N-ii-1], 0);
//			s_print_strmat(nu[N-ii-1]+nx[N-ii-1]+1, nx[N-ii], AL, 0, 0);
			sgemm_r_diag_libstr(nu[N-ii-1]+nx[N-ii-1]+1, nx[N-ii], 1.0, AL, 0, 0, Sx+N-ii, nu[N-ii], 0.0, AL, 0, 0, AL, 0, 0);
//			s_print_strmat(nu[N-ii-1]+nx[N-ii-1]+1, nx[N-ii], AL, 0, 0);
			sgemm_l_diag_libstr(nu[N-ii-1]+nx[N-ii-1], nx[N-ii], 1.0, Si+N-ii-1, 0, AL, 0, 0, 0.0, AL, 0, 0, AL, 0, 0);
//			s_print_strmat(nu[N-ii-1]+nx[N-ii-1]+1, nx[N-ii], AL, 0, 0);

			strmm_rlnn_libstr(nu[N-ii-1]+nx[N-ii-1]+1, nx[N-ii], 1.0, L+(N-ii), nu[N-ii], nu[N-ii], AL, 0, 0, AL, 0, 0);
			srowex_libstr(nx[N-ii], 1.0, AL, nu[N-ii-1]+nx[N-ii-1], 0, tmp_nxM, 0);
			strmv_lnn_libstr(nx[N-ii], nx[N-ii], L+(N-ii), nu[N-ii], nu[N-ii], tmp_nxM, 0, Pb+(N-ii-1), 0);
			sgead_libstr(1, nx[N-ii], 1.0, L+(N-ii), nu[N-ii]+nx[N-ii], nu[N-ii], AL, nu[N-ii-1]+nx[N-ii-1], 0);

			sgemm_r_diag_libstr(nu[N-ii-1]+nx[N-ii-1]+1, nu[N-ii-1]+nx[N-ii-1], 1.0, L+N-ii-1, 0, 0, Si+N-ii-1, 0, 0.0, L+N-ii-1, 0, 0, L+N-ii-1, 0, 0);
			sgemm_l_diag_libstr(nu[N-ii-1]+nx[N-ii-1], nu[N-ii-1]+nx[N-ii-1], 1.0, Si+N-ii-1, 0, L+N-ii-1, 0, 0, 0.0, L+N-ii-1, 0, 0, L+N-ii-1, 0, 0);

			ssyrk_spotrf_ln_libstr(nu[N-ii-1]+nx[N-ii-1]+1, nu[N-ii-1]+nx[N-ii-1], nx[N-ii], AL, 0, 0, AL, 0, 0, L+(N-ii-1), 0, 0, L+(N-ii-1), 0, 0);
//			ssyrk_ln_mn_libstr(nu[N-ii-1]+nx[N-ii-1]+1, nu[N-ii-1]+nx[N-ii-1], nx[N-ii], 1.0, AL, 0, 0, AL, 0, 0, 1.0, L+(N-ii-1), 0, 0, L+(N-ii-1), 0, 0);
//			spotrf_l_mn_libstr(nu[N-ii-1]+nx[N-ii-1]+1, nu[N-ii-1]+nx[N-ii-1], L+(N-ii-1), 0, 0, L+(N-ii-1), 0, 0);
//			s_print_e_strmat(nu[N-ii-1]+nx[N-ii-1]+1, nu[N-ii-1]+nx[N-ii-1], L+N-ii-1, 0, 0);
			}
		else
			{
			sgecp_libstr(nu[N-ii-1]+nx[N-ii-1], nx[N-ii], BAbt+(N-ii-1), 0, 0, AL, 0, 0);
			srowin_libstr(nx[N-ii], 1.0, res_b+(N-ii-1), 0, AL, nu[N-ii-1]+nx[N-ii-1], 0);
			strmm_rlnn_libstr(nu[N-ii-1]+nx[N-ii-1]+1, nx[N-ii], 1.0, L+(N-ii), nu[N-ii], nu[N-ii], AL, 0, 0, AL, 0, 0);
			srowex_libstr(nx[N-ii], 1.0, AL, nu[N-ii-1]+nx[N-ii-1], 0, tmp_nxM, 0);
			strmv_lnn_libstr(nx[N-ii], nx[N-ii], L+(N-ii), nu[N-ii], nu[N-ii], tmp_nxM, 0, Pb+(N-ii-1), 0);
			sgead_libstr(1, nx[N-ii], 1.0, L+(N-ii), nu[N-ii]+nx[N-ii], nu[N-ii], AL, nu[N-ii-1]+nx[N-ii-1], 0);

			if(nb[N-ii-1]>0)
				{
				sdiaad_sp_libstr(nb[N-ii-1], 1.0, Qx_lb+(N-ii-1), 0, idxb[N-ii-1], L+(N-ii-1), 0, 0);
				srowad_sp_libstr(nb[N-ii-1], 1.0, qx_lb+(N-ii-1), 0, idxb[N-ii-1], L+(N-ii-1), nu[N-ii-1]+nx[N-ii-1], 0);
				}

			if(ng[N-ii-1]>0)
				{
				sgemm_r_diag_libstr(nu[N-ii-1]+nx[N-ii-1], ng[N-ii-1], 1.0, DCt+N-ii-1, 0, 0, Qx_lg+N-ii-1, 0, 0.0, AL+0, 0, nx[N-ii], AL+0, 0, nx[N-ii]);
				srowin_libstr(ng[N-ii-1], 1.0, qx_lg+N-ii-1, 0, AL+0, nu[N-ii-1]+nx[N-ii-1], nx[N-ii]);
				sgecp_libstr(nu[N-ii-1]+nx[N-ii-1], nx[N-ii], AL+0, 0, 0, AL+1, 0, 0);
				sgecp_libstr(nu[N-ii-1]+nx[N-ii-1], ng[N-ii-1], DCt+N-ii-1, 0, 0, AL+1, 0, nx[N-ii]);
				ssyrk_spotrf_ln_libstr(nu[N-ii-1]+nx[N-ii-1]+1, nu[N-ii-1]+nx[N-ii-1], nx[N-ii]+ng[N-ii-1], AL+0, 0, 0, AL+1, 0, 0, L+N-ii-1, 0, 0, L+N-ii-1, 0, 0);
				}
			else
				{
				ssyrk_spotrf_ln_libstr(nu[N-ii-1]+nx[N-ii-1]+1, nu[N-ii-1]+nx[N-ii-1], nx[N-ii], AL, 0, 0, AL, 0, 0, L+(N-ii-1), 0, 0, L+(N-ii-1), 0, 0);
				}
			}

//		d_print_strmat(nu[N-ii-1]+nx[N-ii-1]+1, nu[N-ii-1]+nx[N-ii-1], L+(N-ii-1), 0, 0);
		}

	// forward substitution

	// first stage
	ii = 0;
	if(ws->scale)
		{
		srowex_libstr(nu[ii]+nx[ii], -1.0, L+(ii), nu[ii]+nx[ii], 0, dux+ii, 0);
		strsv_ltn_libstr(nu[ii]+nx[ii], L+ii, 0, 0, dux+ii, 0, dux+ii, 0);

		sgecp_libstr(nu[ii]+nx[ii], nx[ii+1], BAbt+(ii), 0, 0, AL, 0, 0);
		sgemm_r_diag_libstr(nu[ii]+nx[ii], nx[ii+1], 1.0, AL, 0, 0, Sx+ii+1, nu[ii+1], 0.0, AL, 0, 0, AL, 0, 0);
		sgemm_l_diag_libstr(nu[ii]+nx[ii], nx[ii+1], 1.0, Si+ii, 0, AL, 0, 0, 0.0, AL, 0, 0, AL, 0, 0);
		sveccp_libstr(nx[ii+1], res_b+ii, 0, dux+(ii+1), nu[ii+1]);
		d_ptr0 = (dux+ii+1)->pa+nu[ii+1];
		d_ptr1 = (Sx+ii+1)->pa+nu[ii+1];
		for(jj=0; jj<nx[ii+1]; jj++)
			d_ptr0[jj] *= d_ptr1[jj];
		sgemv_t_libstr(nu[ii]+nx[ii], nx[ii+1], 1.0, AL, 0, 0, dux+ii, 0, 1.0, dux+(ii+1), nu[ii+1], dux+(ii+1), nu[ii+1]);

		srowex_libstr(nx[ii+1], 1.0, L+(ii+1), nu[ii+1]+nx[ii+1], nu[ii+1], tmp_nxM, 0);
		strmv_ltn_libstr(nx[ii+1], nx[ii+1], L+(ii+1), nu[ii+1], nu[ii+1], dux+(ii+1), nu[ii+1], dpi+ii, 0);
		saxpy_libstr(nx[ii+1], 1.0, tmp_nxM, 0, dpi+ii, 0, dpi+ii, 0);
		strmv_lnn_libstr(nx[ii+1], nx[ii+1], L+(ii+1), nu[ii+1], nu[ii+1], dpi+ii, 0, dpi+ii, 0);
		}
	else
		{
		srowex_libstr(nu[ii]+nx[ii], -1.0, L+(ii), nu[ii]+nx[ii], 0, dux+ii, 0);
		strsv_ltn_libstr(nu[ii]+nx[ii], L+ii, 0, 0, dux+ii, 0, dux+ii, 0);
		sgemv_t_libstr(nu[ii]+nx[ii], nx[ii+1], 1.0, BAbt+ii, 0, 0, dux+ii, 0, 1.0, res_b+ii, 0, dux+(ii+1), nu[ii+1]);
		srowex_libstr(nx[ii+1], 1.0, L+(ii+1), nu[ii+1]+nx[ii+1], nu[ii+1], tmp_nxM, 0);
		strmv_ltn_libstr(nx[ii+1], nx[ii+1], L+(ii+1), nu[ii+1], nu[ii+1], dux+(ii+1), nu[ii+1], dpi+ii, 0);
		saxpy_libstr(nx[ii+1], 1.0, tmp_nxM, 0, dpi+ii, 0, dpi+ii, 0);
		strmv_lnn_libstr(nx[ii+1], nx[ii+1], L+(ii+1), nu[ii+1], nu[ii+1], dpi+ii, 0, dpi+ii, 0);
		}

//	d_print_tran_strvec(nu[ii]+nx[ii], dux+ii, 0);

	// middle stages
	for(ii=1; ii<N; ii++)
		{
		if(ws->scale)
			{
			srowex_libstr(nu[ii], -1.0, L+(ii), nu[ii]+nx[ii], 0, dux+ii, 0);
			strsv_ltn_mn_libstr(nu[ii]+nx[ii], nu[ii], L+ii, 0, 0, dux+ii, 0, dux+ii, 0);

			sgecp_libstr(nu[ii]+nx[ii], nx[ii+1], BAbt+(ii), 0, 0, AL, 0, 0);
			sgemm_r_diag_libstr(nu[ii]+nx[ii], nx[ii+1], 1.0, AL, 0, 0, Sx+ii+1, nu[ii+1], 0.0, AL, 0, 0, AL, 0, 0);
			sgemm_l_diag_libstr(nu[ii]+nx[ii], nx[ii+1], 1.0, Si+ii, 0, AL, 0, 0, 0.0, AL, 0, 0, AL, 0, 0);
			sveccp_libstr(nx[ii+1], res_b+ii, 0, dux+(ii+1), nu[ii+1]);
			d_ptr0 = (dux+ii+1)->pa+nu[ii+1];
			d_ptr1 = (Sx+ii+1)->pa+nu[ii+1];
			for(jj=0; jj<nx[ii+1]; jj++)
				d_ptr0[jj] *= d_ptr1[jj];
			sgemv_t_libstr(nu[ii]+nx[ii], nx[ii+1], 1.0, AL, 0, 0, dux+ii, 0, 1.0, dux+(ii+1), nu[ii+1], dux+(ii+1), nu[ii+1]);

			srowex_libstr(nx[ii+1], 1.0, L+(ii+1), nu[ii+1]+nx[ii+1], nu[ii+1], tmp_nxM, 0);
			strmv_ltn_libstr(nx[ii+1], nx[ii+1], L+(ii+1), nu[ii+1], nu[ii+1], dux+(ii+1), nu[ii+1], dpi+ii, 0);
			saxpy_libstr(nx[ii+1], 1.0, tmp_nxM, 0, dpi+ii, 0, dpi+ii, 0);
			strmv_lnn_libstr(nx[ii+1], nx[ii+1], L+(ii+1), nu[ii+1], nu[ii+1], dpi+ii, 0, dpi+ii, 0);
			}
		else
			{
			srowex_libstr(nu[ii], -1.0, L+(ii), nu[ii]+nx[ii], 0, dux+ii, 0);
			strsv_ltn_mn_libstr(nu[ii]+nx[ii], nu[ii], L+ii, 0, 0, dux+ii, 0, dux+ii, 0);
			sgemv_t_libstr(nu[ii]+nx[ii], nx[ii+1], 1.0, BAbt+ii, 0, 0, dux+ii, 0, 1.0, res_b+ii, 0, dux+(ii+1), nu[ii+1]);
			srowex_libstr(nx[ii+1], 1.0, L+(ii+1), nu[ii+1]+nx[ii+1], nu[ii+1], tmp_nxM, 0);
			strmv_ltn_libstr(nx[ii+1], nx[ii+1], L+(ii+1), nu[ii+1], nu[ii+1], dux+(ii+1), nu[ii+1], dpi+ii, 0);
			saxpy_libstr(nx[ii+1], 1.0, tmp_nxM, 0, dpi+ii, 0, dpi+ii, 0);
			strmv_lnn_libstr(nx[ii+1], nx[ii+1], L+(ii+1), nu[ii+1], nu[ii+1], dpi+ii, 0, dpi+ii, 0);
			}

//		d_print_tran_strvec(nu[ii]+nx[ii], dux+ii, 0);
		}
	
#if 1
	for(ii=0; ii<=N; ii++)
		s_print_e_tran_strvec(nu[ii]+nx[ii], dux+ii, 0);
	for(ii=0; ii<N; ii++)
		s_print_e_tran_strvec(nx[ii+1], dpi+ii, 0);
//	exit(1);
#endif
#if 0
	for(ii=0; ii<=N; ii++)
		s_print_e_strmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], L+ii, 0, 0);
#endif

	// scale solution back
	if(ws->scale)
		{
		for(ii=0; ii<=N; ii++)
			{
			d_ptr0 = (dux+ii)->pa;
			d_ptr1 = (Si+ii)->pa;
			for(jj=0; jj<nu[ii]+nx[ii]; jj++)
				d_ptr0[jj] *= d_ptr1[jj];
			}
		for(ii=0; ii<N; ii++)
			{
			d_ptr0 = (dpi+ii)->pa;
			d_ptr1 = (Sx+ii+1)->pa+nu[ii+1];
			for(jj=0; jj<nx[ii+1]; jj++)
				d_ptr0[jj] *= d_ptr1[jj];
			}
		}


	// cvt single => double
	for(ii=0; ii<N; ii++)
		{
		m_cvt_s2d_strvec(nu[ii]+nx[ii], dux+ii, 0, d_dux+ii, 0);
		m_cvt_s2d_strvec(nx[ii+1], dpi+ii, 0, d_dpi+ii, 0);
		}
	ii = N;
	m_cvt_s2d_strvec(nu[ii]+nx[ii], dux+ii, 0, d_dux+ii, 0);



//	if(nb>0)
//		{
		for(ii=0; ii<=N; ii++)
			dvecex_sp_libstr(nb[ii], 1.0, d_idxb[ii], d_dux+ii, 0, d_dt_lb+ii, 0);
//		}

//	if(ng>0)
//		{
		for(ii=0; ii<=N; ii++)
			dgemv_t_libstr(nu[ii]+nx[ii], ng[ii], 1.0, d_DCt+ii, 0, 0, d_dux+ii, 0, 0.0, d_dt_lg+ii, 0, d_dt_lg+ii, 0);
//		}

//	if(nb>0 | ng>0)
//		{
		d_compute_lam_t_hard_qp(cws);
//		}

	return;

	}



// backward Riccati recursion
void m_solve_kkt_step_hard_ocp_qp(struct d_ocp_qp *d_qp, struct s_ocp_qp *s_qp, struct m_ipm_hard_ocp_qp_workspace *ws)
	{

	int N = s_qp->N;
	int *nx = s_qp->nx;
	int *nu = s_qp->nu;
	int *nb = s_qp->nb;
	int *ng = s_qp->ng;

	struct s_strmat *BAbt = s_qp->BAbt;
	struct s_strmat *RSQrq = s_qp->RSQrq;
	struct s_strmat *DCt = s_qp->DCt;
	struct d_strmat *d_DCt = d_qp->DCt;
	int **idxb = s_qp->idxb;
	int **d_idxb = d_qp->idxb;

	struct s_strmat *L = ws->L;
	struct s_strmat *AL = ws->AL;
	struct d_strvec *d_res_b = ws->res_b;
	struct d_strvec *d_res_g = ws->res_g;
	struct s_strvec *res_b = ws->sres_b;
	struct s_strvec *res_g = ws->sres_g;
	struct d_strvec *d_dux = ws->dux;
	struct d_strvec *d_dpi = ws->dpi;
	struct s_strvec *dux = ws->sdux;
	struct s_strvec *dpi = ws->sdpi;
	struct d_strvec *d_dt_lb = ws->dt_lb;
	struct d_strvec *d_dt_lg = ws->dt_lg;
	struct d_strvec *d_qx_lg = ws->qx_lg;
	struct d_strvec *d_qx_lb = ws->qx_lb;
	struct s_strvec *qx_lg = ws->sqx_lg;
	struct s_strvec *qx_lb = ws->sqx_lb;
	struct s_strvec *Pb = ws->Pb;
	struct s_strvec *tmp_nxM = ws->tmp_nxM;

	//
	int ii;

	struct d_ipm_hard_core_qp_workspace *cws = ws->core_workspace;

//	if(nb>0 | ng>0)
//		{
		d_compute_qx_hard_qp(cws);
//		}



	// cvt double => single
	for(ii=0; ii<N; ii++)
		{
		m_cvt_d2s_strvec(nu[ii]+nx[ii], d_res_g+ii, 0, res_g+ii, 0);
		m_cvt_d2s_strvec(nx[ii+1], d_res_b+ii, 0, res_b+ii, 0);
		m_cvt_d2s_strvec(nb[ii], d_qx_lb+ii, 0, qx_lb+ii, 0);
		m_cvt_d2s_strvec(ng[ii], d_qx_lg+ii, 0, qx_lg+ii, 0);
		}
	ii = N;
	m_cvt_d2s_strvec(nu[ii]+nx[ii], d_res_g+ii, 0, res_g+ii, 0);
	m_cvt_d2s_strvec(nb[ii], d_qx_lb+ii, 0, qx_lb+ii, 0);
	m_cvt_d2s_strvec(ng[ii], d_qx_lg+ii, 0, qx_lg+ii, 0);



	// backward substitution

	// last stage
	sveccp_libstr(nu[N]+nx[N], res_g+N, 0, dux+N, 0);
	if(nb[N]>0)
		{
		svecad_sp_libstr(nb[N], 1.0, qx_lb+N, 0, idxb[N], dux+N, 0);
		}
	// general constraints
	if(ng[N]>0)
		{
		sgemv_n_libstr(nu[N]+nx[N], ng[N], 1.0, DCt+N, 0, 0, qx_lg+N, 0, 1.0, dux+N, 0, dux+N, 0);
		}

	// middle stages
	for(ii=0; ii<N-1; ii++)
		{
		sveccp_libstr(nu[N-ii-1]+nx[N-ii-1], res_g+N-ii-1, 0, dux+N-ii-1, 0);
		if(nb[N-ii-1]>0)
			{
			svecad_sp_libstr(nb[N-ii-1], 1.0, qx_lb+N-ii-1, 0, idxb[N-ii-1], dux+N-ii-1, 0);
			}
		if(ng[N-ii-1]>0)
			{
			sgemv_n_libstr(nu[N-ii-1]+nx[N-ii-1], ng[N-ii-1], 1.0, DCt+N-ii-1, 0, 0, qx_lg+N-ii-1, 0, 1.0, dux+N-ii-1, 0, dux+N-ii-1, 0);
			}
		if(ws->compute_Pb)
			{
			strmv_ltn_libstr(nx[N-ii], nx[N-ii], L+(N-ii), nu[N-ii], nu[N-ii], res_b+N-ii-1, 0, Pb+(N-ii-1), 0);
			strmv_lnn_libstr(nx[N-ii], nx[N-ii], L+(N-ii), nu[N-ii], nu[N-ii], Pb+(N-ii-1), 0, Pb+(N-ii-1), 0);
			}
		saxpy_libstr(nx[N-ii], 1.0, dux+N-ii, nu[N-ii], Pb+N-ii-1, 0, tmp_nxM, 0);
		sgemv_n_libstr(nu[N-ii-1]+nx[N-ii-1], nx[N-ii], 1.0, BAbt+N-ii-1, 0, 0, tmp_nxM, 0, 1.0, dux+N-ii-1, 0, dux+N-ii-1, 0);
		strsv_lnn_mn_libstr(nu[N-ii-1]+nx[N-ii-1], nu[N-ii-1], L+N-ii-1, 0, 0, dux+N-ii-1, 0, dux+N-ii-1, 0);
		}

	// first stage
	ii = N-1;
	sveccp_libstr(nu[N-ii-1]+nx[N-ii-1], res_g+N-ii-1, 0, dux+N-ii-1, 0);
	if(nb[N-ii-1]>0)
		{
		svecad_sp_libstr(nb[N-ii-1], 1.0, qx_lb+N-ii-1, 0, idxb[N-ii-1], dux+N-ii-1, 0);
		}
	if(ng[N-ii-1]>0)
		{
		sgemv_n_libstr(nu[N-ii-1]+nx[N-ii-1], ng[N-ii-1], 1.0, DCt+N-ii-1, 0, 0, qx_lg+N-ii-1, 0, 1.0, dux+N-ii-1, 0, dux+N-ii-1, 0);
		}
	if(ws->compute_Pb)
		{
		strmv_ltn_libstr(nx[N-ii], nx[N-ii], L+(N-ii), nu[N-ii], nu[N-ii], res_b+N-ii-1, 0, Pb+(N-ii-1), 0);
		strmv_lnn_libstr(nx[N-ii], nx[N-ii], L+(N-ii), nu[N-ii], nu[N-ii], Pb+(N-ii-1), 0, Pb+(N-ii-1), 0);
		}
	saxpy_libstr(nx[N-ii], 1.0, dux+N-ii, nu[N-ii], Pb+N-ii-1, 0, tmp_nxM, 0);
	sgemv_n_libstr(nu[N-ii-1]+nx[N-ii-1], nx[N-ii], 1.0, BAbt+N-ii-1, 0, 0, tmp_nxM, 0, 1.0, dux+N-ii-1, 0, dux+N-ii-1, 0);
	strsv_lnn_libstr(nu[N-ii-1]+nx[N-ii-1], L+N-ii-1, 0, 0, dux+N-ii-1, 0, dux+N-ii-1, 0);

	// first stage
	ii = 0;
	sveccp_libstr(nx[ii+1], dux+ii+1, nu[ii+1], dpi+ii, 0);
	svecsc_libstr(nu[ii]+nx[ii], -1.0, dux+ii, 0);
	strsv_ltn_libstr(nu[ii]+nx[ii], L+ii, 0, 0, dux+ii, 0, dux+ii, 0);
	sgemv_t_libstr(nu[ii]+nx[ii], nx[ii+1], 1.0, BAbt+ii, 0, 0, dux+ii, 0, 1.0, res_b+ii, 0, dux+ii+1, nu[ii+1]);
	sveccp_libstr(nx[ii+1], dux+ii+1, nu[ii+1], tmp_nxM, 0);
	strmv_ltn_libstr(nx[ii+1], nx[ii+1], L+ii+1, nu[ii+1], nu[ii+1], tmp_nxM, 0, tmp_nxM, 0);
	strmv_lnn_libstr(nx[ii+1], nx[ii+1], L+ii+1, nu[ii+1], nu[ii+1], tmp_nxM, 0, tmp_nxM, 0);
	saxpy_libstr(nx[ii+1], 1.0, tmp_nxM, 0, dpi+ii, 0, dpi+ii, 0);

	// middle stages
	for(ii=1; ii<N; ii++)
		{
		sveccp_libstr(nx[ii+1], dux+ii+1, nu[ii+1], dpi+ii, 0);
		svecsc_libstr(nu[ii], -1.0, dux+ii, 0);
		strsv_ltn_mn_libstr(nu[ii]+nx[ii], nu[ii], L+ii, 0, 0, dux+ii, 0, dux+ii, 0);
		sgemv_t_libstr(nu[ii]+nx[ii], nx[ii+1], 1.0, BAbt+ii, 0, 0, dux+ii, 0, 1.0, res_b+ii, 0, dux+ii+1, nu[ii+1]);
		sveccp_libstr(nx[ii+1], dux+ii+1, nu[ii+1], tmp_nxM, 0);
		strmv_ltn_libstr(nx[ii+1], nx[ii+1], L+ii+1, nu[ii+1], nu[ii+1], tmp_nxM, 0, tmp_nxM, 0);
		strmv_lnn_libstr(nx[ii+1], nx[ii+1], L+ii+1, nu[ii+1], nu[ii+1], tmp_nxM, 0, tmp_nxM, 0);
		saxpy_libstr(nx[ii+1], 1.0, tmp_nxM, 0, dpi+ii, 0, dpi+ii, 0);
		}



	// cvt single => double
	for(ii=0; ii<N; ii++)
		{
		m_cvt_s2d_strvec(nu[ii]+nx[ii], dux+ii, 0, d_dux+ii, 0);
		m_cvt_s2d_strvec(nx[ii+1], dpi+ii, 0, d_dpi+ii, 0);
		}
	ii = N;
	m_cvt_s2d_strvec(nu[ii]+nx[ii], dux+ii, 0, d_dux+ii, 0);



//	if(nb>0)
//		{
		for(ii=0; ii<=N; ii++)
			dvecex_sp_libstr(nb[ii], 1.0, d_idxb[ii], d_dux+ii, 0, d_dt_lb+ii, 0);
//		}

//	if(ng>0)
//		{
		for(ii=0; ii<=N; ii++)
			dgemv_t_libstr(nu[ii]+nx[ii], ng[ii], 1.0, d_DCt+ii, 0, 0, d_dux+ii, 0, 0.0, d_dt_lg+ii, 0, d_dt_lg+ii, 0);
//		}

//	if(nb>0 | ng>0)
//		{
		d_compute_lam_t_hard_qp(cws);
//		}

	return;

	}


