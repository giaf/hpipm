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

#include <mmintrin.h>
#include <xmmintrin.h>  // SSE
#include <emmintrin.h>  // SSE2
#include <pmmintrin.h>  // SSE3
#include <smmintrin.h>  // SSE4
#include <immintrin.h>  // AVX

#include "../include/hpipm_d_core_qp_ipm.h"



void d_compute_Gamma_gamma_qp(double *res_d, double *res_m, struct d_core_qp_ipm_workspace *cws)
	{

	int nc = cws->nc;

	double *lam = cws->lam;
	double *t = cws->t;
	double *t_inv = cws->t_inv;
	double *Gamma = cws->Gamma;
	double *gamma = cws->gamma;

	__m256d
		y_ones,
		y_tmp0, y_t_inv0, y_lam0,
		y_tmp1, y_t_inv1, y_lam1;

	y_ones = _mm256_set_pd( 1.0, 1.0, 1.0, 1.0 );

	// local variables
	int ii;

	ii = 0;
#if 0
	for(; ii<nc-7; ii+=8)
		{
		y_t_inv0 = _mm256_div_pd( y_ones, _mm256_loadu_pd( &t[ii+0] ) );
		y_t_inv1 = _mm256_div_pd( y_ones, _mm256_loadu_pd( &t[ii+4] ) );
		_mm256_storeu_pd( &t_inv[ii+0], y_t_inv0 );
		_mm256_storeu_pd( &t_inv[ii+4], y_t_inv1 );
		y_lam0 = _mm256_loadu_pd( &lam[ii+0] );
		y_lam1 = _mm256_loadu_pd( &lam[ii+4] );
		y_tmp0 = _mm256_mul_pd( y_t_inv0, y_lam0 );
		y_tmp1 = _mm256_mul_pd( y_t_inv1, y_lam1 );
		_mm256_storeu_pd( &Gamma[ii+0], y_tmp0 );
		_mm256_storeu_pd( &Gamma[ii+4], y_tmp1 );
		y_tmp0 = _mm256_mul_pd( y_lam0, _mm256_loadu_pd( &res_d[ii+0] ) );
		y_tmp1 = _mm256_mul_pd( y_lam1, _mm256_loadu_pd( &res_d[ii+4] ) );
		y_tmp0 = _mm256_sub_pd( _mm256_loadu_pd( &res_m[ii+0] ), y_tmp0 );
		y_tmp1 = _mm256_sub_pd( _mm256_loadu_pd( &res_m[ii+4] ), y_tmp1 );
		y_tmp0 = _mm256_mul_pd( y_t_inv0, y_tmp0 );
		y_tmp1 = _mm256_mul_pd( y_t_inv1, y_tmp1 );
		_mm256_storeu_pd( &gamma[ii+0], y_tmp0 );
		_mm256_storeu_pd( &gamma[ii+4], y_tmp1 );
		}
#endif
	for(; ii<nc-3; ii+=4)
		{
		y_t_inv0 = _mm256_div_pd( y_ones, _mm256_loadu_pd( &t[ii] ) );
		_mm256_storeu_pd( &t_inv[ii], y_t_inv0 );
		y_lam0 = _mm256_loadu_pd( &lam[ii] );
		y_tmp0 = _mm256_mul_pd( y_t_inv0, y_lam0 );
		_mm256_storeu_pd( &Gamma[ii], y_tmp0 );
		y_tmp0 = _mm256_mul_pd( y_lam0, _mm256_loadu_pd( &res_d[ii] ) );
		y_tmp0 = _mm256_sub_pd( _mm256_loadu_pd( &res_m[ii] ), y_tmp0 );
		y_tmp0 = _mm256_mul_pd( y_t_inv0, y_tmp0 );
		_mm256_storeu_pd( &gamma[ii], y_tmp0 );
		}
	for(; ii<nc; ii++)
		{
		t_inv[ii] = 1.0/t[ii];
		Gamma[ii] = t_inv[ii]*lam[ii];
		gamma[ii] = t_inv[ii]*(res_m[ii]-lam[ii]*res_d[ii]);
		}

	return;

	}



void d_compute_gamma_qp(double *res_d, double *res_m, struct d_core_qp_ipm_workspace *cws)
	{

	int nc = cws->nc;

	double *lam = cws->lam;
	double *t_inv = cws->t_inv;
	double *gamma = cws->gamma;

	__m256d
		y_tmp0, y_lam0,
		y_tmp1, y_lam1;

	// local variables
	int ii;

	ii = 0;
#if 0
	for(; ii<nc-7; ii+=8)
		{
		y_lam0 = _mm256_loadu_pd( &lam[ii+0] );
		y_lam1 = _mm256_loadu_pd( &lam[ii+4] );
		y_tmp0 = _mm256_mul_pd( y_lam0, _mm256_loadu_pd( &res_d[ii+0] ) );
		y_tmp1 = _mm256_mul_pd( y_lam1, _mm256_loadu_pd( &res_d[ii+4] ) );
		y_tmp0 = _mm256_sub_pd( _mm256_loadu_pd( &res_m[ii+0] ), y_tmp0 );
		y_tmp1 = _mm256_sub_pd( _mm256_loadu_pd( &res_m[ii+4] ), y_tmp1 );
		y_tmp0 = _mm256_mul_pd( _mm256_loadu_pd( &t_inv[ii+0] ), y_tmp0 );
		y_tmp1 = _mm256_mul_pd( _mm256_loadu_pd( &t_inv[ii+4] ), y_tmp1 );
		_mm256_storeu_pd( &gamma[ii+0], y_tmp0 );
		_mm256_storeu_pd( &gamma[ii+4], y_tmp1 );
		}
#endif
	for(; ii<nc-3; ii+=4)
		{
		y_lam0 = _mm256_loadu_pd( &lam[ii] );
		y_tmp0 = _mm256_mul_pd( y_lam0, _mm256_loadu_pd( &res_d[ii] ) );
		y_tmp0 = _mm256_sub_pd( _mm256_loadu_pd( &res_m[ii] ), y_tmp0 );
		y_tmp0 = _mm256_mul_pd( _mm256_loadu_pd( &t_inv[ii] ), y_tmp0 );
		_mm256_storeu_pd( &gamma[ii], y_tmp0 );
		}
	for(; ii<nc; ii++)
		{
		gamma[ii] = t_inv[ii]*(res_m[ii]-lam[ii]*res_d[ii]);
		}

	return;

	}



void d_compute_lam_t_qp(double *res_d, double *res_m, double *dlam, double *dt, struct d_core_qp_ipm_workspace *cws)
	{

	int nc = cws->nc;

	double *lam = cws->lam;
	double *t_inv = cws->t_inv;

	__m256d
		y_sign,
		y_tmp0, y_tmp2, y_dt0,
		y_tmp1, y_tmp3, y_dt1;

	long long long_sign = 0x8000000000000000;
	y_sign = _mm256_broadcast_sd( (double *) &long_sign );

	// local variables
	int ii;

	ii = 0;
#if 0
	for(; ii<nc-7; ii+=8)
		{
		y_dt0 = _mm256_loadu_pd( &dt[ii+0] );
		y_dt1 = _mm256_loadu_pd( &dt[ii+4] );
		y_dt0 = _mm256_sub_pd( y_dt0, _mm256_loadu_pd( &res_d[ii+0] ) );
		y_dt1 = _mm256_sub_pd( y_dt1, _mm256_loadu_pd( &res_d[ii+4] ) );
		_mm256_storeu_pd( &dt[ii+0], y_dt0 );
		_mm256_storeu_pd( &dt[ii+4], y_dt1 );
		y_tmp0 = _mm256_mul_pd( y_dt0, _mm256_loadu_pd( &lam[ii+0] ) );
		y_tmp1 = _mm256_mul_pd( y_dt1, _mm256_loadu_pd( &lam[ii+4] ) );
		y_tmp2 = _mm256_loadu_pd( &t_inv[ii+0] );
		y_tmp3 = _mm256_loadu_pd( &t_inv[ii+4] );
		y_tmp0 = _mm256_add_pd( y_tmp0, _mm256_loadu_pd( &res_m[ii+0] ) );
		y_tmp1 = _mm256_add_pd( y_tmp1, _mm256_loadu_pd( &res_m[ii+4] ) );
		y_tmp2 = _mm256_xor_pd( y_tmp2, y_sign );
		y_tmp3 = _mm256_xor_pd( y_tmp3, y_sign );
		y_tmp0 = _mm256_mul_pd( y_tmp0, y_tmp2 );
		y_tmp1 = _mm256_mul_pd( y_tmp1, y_tmp3 );
		_mm256_storeu_pd( &dlam[ii+0], y_tmp0 );
		_mm256_storeu_pd( &dlam[ii+4], y_tmp1 );
		}
#endif
	for(; ii<nc-3; ii+=4)
		{
		y_dt0 = _mm256_loadu_pd( &dt[ii+0] );
		y_dt0 = _mm256_sub_pd( y_dt0, _mm256_loadu_pd( &res_d[ii+0] ) );
		_mm256_storeu_pd( &dt[ii+0], y_dt0 );
		y_tmp0 = _mm256_mul_pd( y_dt0, _mm256_loadu_pd( &lam[ii+0] ) );
		y_tmp2 = _mm256_loadu_pd( &t_inv[ii+0] );
		y_tmp0 = _mm256_add_pd( y_tmp0, _mm256_loadu_pd( &res_m[ii+0] ) );
		y_tmp2 = _mm256_xor_pd( y_tmp2, y_sign );
		y_tmp0 = _mm256_mul_pd( y_tmp0, y_tmp2 );
		_mm256_storeu_pd( &dlam[ii+0], y_tmp0 );
		}
	for(; ii<nc; ii++)
		{
		dt[ii] -= res_d[ii];
		// TODO compute lamda alone ???
		dlam[ii] = - t_inv[ii] * (lam[ii]*dt[ii] + res_m[ii]);
		}
	
	return;

	}



void d_compute_alpha_qp(struct d_core_qp_ipm_workspace *cws)
	{
	
	// extract workspace members
	int nc = cws->nc;

	double *lam = cws->lam;
	double *t = cws->t;
	double *dlam = cws->dlam;
	double *dt = cws->dt;

	double alpha_prim = - 1.0;
	double alpha_dual = - 1.0;
	double alpha = - 1.0;

	__m256d
		y_alpha0, y_alpha1;

	__m128d
		x_alpha0, x_alpha1;

	__m128
		s_zeros, s_mones,
		s_tmp0, s_tmp2, s_mask0, s_mask2, s_alpha0, s_alpha2,
		s_tmp1, s_tmp3, s_mask1, s_mask3, s_alpha1, s_alpha3;
	
	s_mones  = _mm_set_ps( -1.0, -1.0, -1.0, -1.0 );
	s_zeros = _mm_setzero_ps( );

	s_alpha0 = _mm_set_ps( -1.0, -1.0, -1.0, -1.0 );
	s_alpha2 = _mm_set_ps( -1.0, -1.0, -1.0, -1.0 );

	// local variables
	int ii;

	ii = 0;
#if 0
	s_alpha1 = _mm_set_ps( -1.0, -1.0, -1.0, -1.0 );
	s_alpha3 = _mm_set_ps( -1.0, -1.0, -1.0, -1.0 );
	for(; ii<nc-7; ii+=8)
		{
		s_tmp0 = _mm256_cvtpd_ps( _mm256_loadu_pd( &dlam[ii+0] ) );
		s_tmp1 = _mm256_cvtpd_ps( _mm256_loadu_pd( &dlam[ii+4] ) );
		s_tmp2 = _mm256_cvtpd_ps( _mm256_loadu_pd( &dt[ii+0] ) );
		s_tmp3 = _mm256_cvtpd_ps( _mm256_loadu_pd( &dt[ii+4] ) );
		s_mask0 = _mm_cmp_ps( s_tmp0, s_zeros, 0x01 );
		s_mask1 = _mm_cmp_ps( s_tmp1, s_zeros, 0x01 );
		s_mask2 = _mm_cmp_ps( s_tmp2, s_zeros, 0x01 );
		s_mask3 = _mm_cmp_ps( s_tmp3, s_zeros, 0x01 );
		s_tmp0 = _mm_div_ps( _mm256_cvtpd_ps( _mm256_loadu_pd( &lam[ii+0] ) ), s_tmp0 );
		s_tmp1 = _mm_div_ps( _mm256_cvtpd_ps( _mm256_loadu_pd( &lam[ii+4] ) ), s_tmp1 );
		s_tmp2 = _mm_div_ps( _mm256_cvtpd_ps( _mm256_loadu_pd( &t[ii+0] ) ), s_tmp2 );
		s_tmp3 = _mm_div_ps( _mm256_cvtpd_ps( _mm256_loadu_pd( &t[ii+4] ) ), s_tmp3 );
		s_tmp0 = _mm_blendv_ps( s_mones, s_tmp0, s_mask0 );
		s_tmp1 = _mm_blendv_ps( s_mones, s_tmp1, s_mask1 );
		s_tmp2 = _mm_blendv_ps( s_mones, s_tmp2, s_mask2 );
		s_tmp3 = _mm_blendv_ps( s_mones, s_tmp3, s_mask3 );
		s_alpha0 = _mm_max_ps( s_alpha0, s_tmp0 );
		s_alpha1 = _mm_max_ps( s_alpha1, s_tmp1 );
		s_alpha2 = _mm_max_ps( s_alpha2, s_tmp2 );
		s_alpha3 = _mm_max_ps( s_alpha3, s_tmp3 );
		}
	s_alpha0 = _mm_max_ps( s_alpha0, s_alpha1 );
	s_alpha2 = _mm_max_ps( s_alpha2, s_alpha3 );
#endif
	for(; ii<nc-3; ii+=4)
		{
		s_tmp0 = _mm256_cvtpd_ps( _mm256_loadu_pd( &dlam[ii+0] ) );
		s_tmp2 = _mm256_cvtpd_ps( _mm256_loadu_pd( &dt[ii+0] ) );
		s_mask0 = _mm_cmp_ps( s_tmp0, s_zeros, 0x01 );
		s_mask2 = _mm_cmp_ps( s_tmp2, s_zeros, 0x01 );
		s_tmp0 = _mm_div_ps( _mm256_cvtpd_ps( _mm256_loadu_pd( &lam[ii+0] ) ), s_tmp0 );
		s_tmp2 = _mm_div_ps( _mm256_cvtpd_ps( _mm256_loadu_pd( &t[ii+0] ) ), s_tmp2 );
		s_tmp0 = _mm_blendv_ps( s_mones, s_tmp0, s_mask0 );
		s_tmp2 = _mm_blendv_ps( s_mones, s_tmp2, s_mask2 );
		s_alpha0 = _mm_max_ps( s_alpha0, s_tmp0 );
		s_alpha2 = _mm_max_ps( s_alpha2, s_tmp2 );
		}
	for(; ii<nc; ii++)
		{

		if( alpha_dual*dlam[ii]>lam[ii] )
			{
			alpha_dual = lam[ii] / dlam[ii];
			}
		if( alpha_prim*dt[ii]>t[ii] )
			{
			alpha_prim = t[ii] / dt[ii];
			}

		}

	y_alpha0 = _mm256_cvtps_pd( s_alpha0 );
	y_alpha1 = _mm256_cvtps_pd( s_alpha2 );
	x_alpha0 = _mm_max_pd( _mm256_extractf128_pd( y_alpha0, 0x1 ), _mm256_castpd256_pd128( y_alpha0 ) );
	x_alpha1 = _mm_max_pd( _mm256_extractf128_pd( y_alpha1, 0x1 ), _mm256_castpd256_pd128( y_alpha1 ) );
	x_alpha0 = _mm_max_sd( x_alpha0, _mm_permute_pd( x_alpha0, 0x1 ) );
	x_alpha1 = _mm_max_sd( x_alpha1, _mm_permute_pd( x_alpha1, 0x1 ) );
	x_alpha0 = _mm_max_sd( x_alpha0, _mm_load_sd( &alpha_dual ) );
	x_alpha1 = _mm_max_sd( x_alpha1, _mm_load_sd( &alpha_prim ) );
	_mm_store_sd( &alpha_dual, x_alpha0 );
	_mm_store_sd( &alpha_prim, x_alpha1 );

	alpha = alpha_prim>alpha_dual ? alpha_prim : alpha_dual;

	// store alpha
	cws->alpha_prim = - alpha_prim;
	cws->alpha_dual = - alpha_dual;
	cws->alpha = - alpha;

	return;

	}
	


void d_update_var_qp(struct d_core_qp_ipm_workspace *cws)
	{
	
	// extract workspace members
	int nv = cws->nv;
	int ne = cws->ne;
	int nc = cws->nc;

	double *v = cws->v;
	double *pi = cws->pi;
	double *lam = cws->lam;
	double *t = cws->t;
	double *dv = cws->dv;
	double *dpi = cws->dpi;
	double *dlam = cws->dlam;
	double *dt = cws->dt;
	double alpha = cws->alpha;
	double alpha_prim = cws->alpha_prim;
	double alpha_dual = cws->alpha_dual;

	__m256d
		y_tmp0, y_tmp1,
		y_alpha,
		y_lam_min, y_t_min, y_mask0, y_mask1;
	
#if 0
	if(alpha<1.0)
		alpha *= 0.995;
#else
//	alpha_prim = alpha_prim * ((1.0-alpha)*0.99 + alpha*0.9999999);
//	alpha_dual = alpha_dual * ((1.0-alpha)*0.99 + alpha*0.9999999);
	alpha_prim = alpha_prim * ((1.0-alpha_prim)*0.99 + alpha_prim*0.9999999);
	alpha_dual = alpha_dual * ((1.0-alpha_dual)*0.99 + alpha_dual*0.9999999);
	alpha = alpha * ((1.0-alpha)*0.99 + alpha*0.9999999);
#endif

	// local variables
	int ii;

#if 1

	y_alpha = _mm256_broadcast_sd( &alpha );

	// update v
	ii = 0;
	for(; ii<nv-3; ii+=4)
		{
		y_tmp0 = _mm256_mul_pd( y_alpha, _mm256_loadu_pd( &dv[ii] ) );
		y_tmp0 = _mm256_add_pd( y_tmp0, _mm256_loadu_pd( &v[ii] ) );
		_mm256_storeu_pd( &v[ii], y_tmp0 );
		}
	for(; ii<nv; ii++)
		{
		v[ii] += alpha * dv[ii];
		}

	// update pi
	ii = 0;
	for(; ii<ne-3; ii+=4)
		{
		y_tmp0 = _mm256_mul_pd( y_alpha, _mm256_loadu_pd( &dpi[ii] ) );
		y_tmp0 = _mm256_add_pd( y_tmp0, _mm256_loadu_pd( &pi[ii] ) );
		_mm256_storeu_pd( &pi[ii], y_tmp0 );
		}
	for(; ii<ne; ii++)
		{
		pi[ii] += alpha * dpi[ii];
		}

	// update lam and t
	y_lam_min = _mm256_broadcast_sd( &cws->lam_min );
	y_t_min = _mm256_broadcast_sd( &cws->t_min );
	ii = 0;
	for(; ii<nc-3; ii+=4)
		{
		y_tmp0 = _mm256_mul_pd( y_alpha, _mm256_loadu_pd( &dlam[ii] ) );
		y_tmp1 = _mm256_mul_pd( y_alpha, _mm256_loadu_pd( &dt[ii] ) );
		y_tmp0 = _mm256_add_pd( y_tmp0, _mm256_loadu_pd( &lam[ii] ) );
		y_tmp1 = _mm256_add_pd( y_tmp1, _mm256_loadu_pd( &t[ii] ) );
		// max does not preserve NaN !!!
//		y_tmp0 = _mm256_max_pd( y_tmp0, y_lam_min );
//		y_tmp1 = _mm256_max_pd( y_tmp1, y_t_min );
		y_mask0 = _mm256_cmp_pd( y_tmp0, y_lam_min, 2 );
		y_mask1 = _mm256_cmp_pd( y_tmp1, y_t_min, 2 );
		y_tmp0 = _mm256_blendv_pd( y_tmp0, y_lam_min, y_mask0 );
		y_tmp1 = _mm256_blendv_pd( y_tmp1, y_t_min, y_mask1 );
		_mm256_storeu_pd( &lam[ii], y_tmp0 );
		_mm256_storeu_pd( &t[ii], y_tmp1 );
		}
	for(; ii<nc; ii++)
		{
		lam[ii] += alpha * dlam[ii];
		t[ii] += alpha * dt[ii];
		lam[ii] = lam[ii]<=cws->lam_min ? cws->lam_min : lam[ii];
		t[ii] = t[ii]<=cws->t_min ? cws->t_min : t[ii];
		}

#else // split step

	// update v
	for(ii=0; ii<nv; ii++)
		{
		v[ii] += alpha_prim * dv[ii];
		}

	// update pi
	for(ii=0; ii<ne; ii++)
		{
//		pi[ii] += alpha_prim * dpi[ii];
		pi[ii] += alpha_dual * dpi[ii];
		}

	// update lam
	for(ii=0; ii<nc; ii++)
		{
		lam[ii] += alpha_dual * dlam[ii];
		lam[ii] = lam[ii]<=cws->lam_min ? cws->lam_min : lam[ii];
		}

	// update t
	for(ii=0; ii<nc; ii++)
		{
		t[ii] += alpha_prim * dt[ii];
		t[ii] = t[ii]<=cws->t_min ? cws->t_min : t[ii];
		}

#endif

	return;

	}



void d_compute_mu_aff_qp(struct d_core_qp_ipm_workspace *cws)
	{

	int ii;

	// extract workspace members
	int nc = cws->nc;

	double *lam = cws->lam;
	double *t = cws->t;
	double *dlam = cws->dlam;
	double *dt = cws->dt;
	double alpha = cws->alpha;
	// this affects the minimum value of signa !!!
//		alpha *= 0.99;

	__m256d
		y_tmp0, y_tmp1,
		y_alpha, y_mu;
	
	__m128d
		x_mu;

	double mu = 0;

	y_mu = _mm256_setzero_pd( );

	y_alpha = _mm256_broadcast_sd( &alpha );

	ii = 0;
	for(; ii<nc-3; ii+=4)
		{
		y_tmp0 = _mm256_mul_pd( y_alpha, _mm256_loadu_pd( &dlam[ii] ) );
		y_tmp1 = _mm256_mul_pd( y_alpha, _mm256_loadu_pd( &dt[ii] ) );
		y_tmp0 = _mm256_add_pd( y_tmp0, _mm256_loadu_pd( &lam[ii] ) );
		y_tmp1 = _mm256_add_pd( y_tmp1, _mm256_loadu_pd( &t[ii] ) );
		y_tmp0 = _mm256_mul_pd( y_tmp0, y_tmp1 );
		y_mu = _mm256_add_pd( y_mu, y_tmp0 );
		}
	for(; ii<nc; ii++)
		{
		mu += (lam[ii] + alpha*dlam[ii]) * (t[ii] + alpha*dt[ii]);
		}
	
	x_mu = _mm_add_pd( _mm256_castpd256_pd128( y_mu ), _mm256_extractf128_pd( y_mu, 0x1 ) );
	x_mu = _mm_hadd_pd( x_mu, x_mu );
	x_mu = _mm_add_sd( x_mu, _mm_load_sd( &mu ) );
	_mm_store_sd( &mu, x_mu );

	cws->mu_aff = mu*cws->nc_inv;

	return;

	}



void d_backup_res_m(struct d_core_qp_ipm_workspace *cws)
	{

	int ii;

	// extract workspace members
	int nc = cws->nc;

	double *res_m = cws->res_m;
	double *res_m_bkp = cws->res_m_bkp;

	__m256d
		y_tmp0;

	ii = 0;
	for(; ii<nc-3; ii+=4)
		{
		y_tmp0 = _mm256_loadu_pd( &res_m[ii] );
		_mm256_storeu_pd( &res_m_bkp[ii], y_tmp0 );
		}
	for(; ii<nc; ii++)
		{
		res_m_bkp[ii] = res_m[ii];
		}

	return;

	}



void d_compute_centering_correction_qp(struct d_core_qp_ipm_workspace *cws)
	{

	int ii;

	// extract workspace members
	int nc = cws->nc;

	double *dlam = cws->dlam;
	double *dt = cws->dt;
	double *res_m = cws->res_m;
	double *res_m_bkp = cws->res_m_bkp;

	__m256d
		y_tmp0,
		y_sigma_mu;

	double sigma_mu = cws->sigma*cws->mu;

	y_sigma_mu = _mm256_broadcast_sd( &sigma_mu );

	ii = 0;
	for(; ii<nc-3; ii+=4)
		{
		y_tmp0 = _mm256_mul_pd( _mm256_loadu_pd( &dt[ii] ), _mm256_loadu_pd( &dlam[ii] ) );
		y_tmp0 = _mm256_add_pd( y_tmp0, _mm256_loadu_pd( &res_m_bkp[ii] ) );
		y_tmp0 = _mm256_sub_pd( y_tmp0, y_sigma_mu );
		_mm256_storeu_pd( &res_m[ii], y_tmp0 );
		}
	for(; ii<nc; ii++)
		{
		res_m[ii] = res_m_bkp[ii] + dt[ii] * dlam[ii] - sigma_mu;
		}

	return;

	}



void d_compute_centering_qp(struct d_core_qp_ipm_workspace *cws)
	{

	int ii;

	// extract workspace members
	int nc = cws->nc;

	double *dlam = cws->dlam;
	double *dt = cws->dt;
	double *res_m = cws->res_m;
	double *res_m_bkp = cws->res_m_bkp;

	__m256d
		y_tmp0,
		y_sigma_mu;

	double sigma_mu = cws->sigma*cws->mu;

	y_sigma_mu = _mm256_broadcast_sd( &sigma_mu );

	ii = 0;
	for(; ii<nc-3; ii+=4)
		{
		y_tmp0 = _mm256_sub_pd( _mm256_loadu_pd( &res_m_bkp[ii] ), y_sigma_mu );
		_mm256_storeu_pd( &res_m[ii], y_tmp0 );
		}
	for(; ii<nc; ii++)
		{
		res_m[ii] = res_m_bkp[ii] - sigma_mu;
		}

	return;

	}




