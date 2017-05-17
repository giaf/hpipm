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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_v_aux_ext_dep.h>
#include <blasfeo_d_aux_ext_dep.h>
#include <blasfeo_i_aux_ext_dep.h>
#include <blasfeo_d_aux.h>
#include <blasfeo_d_blas.h>

#include "../include/hpipm_ocp_kkt.h"

#include "tools.h"



#define KEEP_X0 0

#define PRINT 1



/************************************************ 
Mass-spring system: nx/2 masses connected each other with springs (in a row), and the first and the last one to walls. nu (<=nx) controls act on the first nu masses. The system is sampled with sampling time Ts. 
************************************************/
void mass_spring_system(double Ts, int nx, int nu, int N, double *A, double *B, double *b, double *x0)
	{

	int nx2 = nx*nx;

	int info = 0;

	int pp = nx/2; // number of masses
	
/************************************************
* build the continuous time system 
************************************************/
	
	double *T; d_zeros(&T, pp, pp);
	int ii;
	for(ii=0; ii<pp; ii++) T[ii*(pp+1)] = -2;
	for(ii=0; ii<pp-1; ii++) T[ii*(pp+1)+1] = 1;
	for(ii=1; ii<pp; ii++) T[ii*(pp+1)-1] = 1;

	double *Z; d_zeros(&Z, pp, pp);
	double *I; d_zeros(&I, pp, pp); for(ii=0; ii<pp; ii++) I[ii*(pp+1)]=1.0; // = eye(pp);
	double *Ac; d_zeros(&Ac, nx, nx);
	dmcopy(pp, pp, Z, pp, Ac, nx);
	dmcopy(pp, pp, T, pp, Ac+pp, nx);
	dmcopy(pp, pp, I, pp, Ac+pp*nx, nx);
	dmcopy(pp, pp, Z, pp, Ac+pp*(nx+1), nx); 
	free(T);
	free(Z);
	free(I);
	
	d_zeros(&I, nu, nu); for(ii=0; ii<nu; ii++) I[ii*(nu+1)]=1.0; //I = eye(nu);
	double *Bc; d_zeros(&Bc, nx, nu);
	dmcopy(nu, nu, I, nu, Bc+pp, nx);
	free(I);
	
/************************************************
* compute the discrete time system 
************************************************/

	double *bb; d_zeros(&bb, nx, 1);
	dmcopy(nx, 1, bb, nx, b, nx);
		
	dmcopy(nx, nx, Ac, nx, A, nx);
	dscal_3l(nx2, Ts, A);
	expm(nx, A);
	
	d_zeros(&T, nx, nx);
	d_zeros(&I, nx, nx); for(ii=0; ii<nx; ii++) I[ii*(nx+1)]=1.0; //I = eye(nx);
	dmcopy(nx, nx, A, nx, T, nx);
	daxpy_3l(nx2, -1.0, I, T);
	dgemm_nn_3l(nx, nu, nx, T, nx, Bc, nx, B, nx);
	free(T);
	free(I);
	
	int *ipiv = (int *) malloc(nx*sizeof(int));
	dgesv_3l(nx, nu, Ac, nx, ipiv, B, nx, &info);
	free(ipiv);

	free(Ac);
	free(Bc);
	free(bb);
	
			
/************************************************
* initial state 
************************************************/
	
	if(nx==4)
		{
		x0[0] = 5;
		x0[1] = 10;
		x0[2] = 15;
		x0[3] = 20;
		}
	else
		{
		int jj;
		for(jj=0; jj<nx; jj++)
			x0[jj] = 1;
		}

	}



int main()
	{

	int ii;

	int N = 10;
	int nx_ = 8;
	int nu_ = 3;

	// stage-wise variant size
	int nx[N+1];
#if KEEP_X0
	nx[0] = nx_;
#else
	nx[0] = 0;
#endif
	for(ii=1; ii<=N; ii++)
		nx[ii] = nx_;

	int nu[N+1];
	for(ii=0; ii<N; ii++)
		nu[ii] = nu_;
	nu[N] = 0;

	int nb[N+1];
#if KEEP_X0
	nb[0] = nu[0]+nx[0]/2;
#else
	nb[0] = nu[0];
#endif
	for(ii=1; ii<N; ii++)
		nb[ii] = nu[1]+nx[1]/2;
	nb[N] = nx[N]/2;

	int ng[N+1];
	for(ii=0; ii<=N; ii++)
		ng[ii] = 0;
		

/************************************************
* dynamical system
************************************************/	

	double *A; d_zeros(&A, nx_, nx_); // states update matrix

	double *B; d_zeros(&B, nx_, nu_); // inputs matrix

	double *b; d_zeros_align(&b, nx_, 1); // states offset
	double *x0; d_zeros_align(&x0, nx_, 1); // initial state

	double Ts = 0.5; // sampling time
	mass_spring_system(Ts, nx_, nu_, N, A, B, b, x0);
	
	for(ii=0; ii<nx_; ii++)
		b[ii] = 0.1;
	
	for(ii=0; ii<nx_; ii++)
		x0[ii] = 0;
	x0[0] = 2.5;
	x0[1] = 2.5;

#if PRINT
	d_print_mat(nx_, nx_, A, nx_);
	d_print_mat(nx_, nu_, B, nu_);
	d_print_mat(1, nx_, b, 1);
	d_print_mat(1, nx_, x0, 1);
#endif

	struct d_strmat sA;
	d_allocate_strmat(nx_, nx_, &sA);
	d_cvt_mat2strmat(nx_, nx_, A, nx_, &sA, 0, 0);
#if PRINT
	d_print_strmat(nx_, nx_, &sA, 0, 0);
#endif

	struct d_strvec sx0;
	d_allocate_strvec(nx_, &sx0);
	d_cvt_vec2strvec(nx_, x0, &sx0, 0);
#if PRINT
	d_print_tran_strvec(nx_, &sx0, 0);
#endif

	struct d_strvec sb0;
	d_allocate_strvec(nx_, &sb0);
	d_cvt_vec2strvec(nx_, b, &sb0, 0);
#if PRINT
	d_print_tran_strvec(nx_, &sb0, 0);
#endif
#if ! KEEP_X0
	dgemv_n_libstr(nx_, nx_, 1.0, &sA, 0, 0, &sx0, 0, 1.0, &sb0, 0, &sb0, 0);
#endif
#if PRINT
	d_print_tran_strvec(nx_, &sb0, 0);
#endif

	struct d_strmat sBAbt0;
	d_allocate_strmat(nu[0]+nx[0]+1, nx[1], &sBAbt0);
	d_cvt_tran_mat2strmat(nx[1], nu[0], B, nx_, &sBAbt0, 0, 0);
	d_cvt_tran_mat2strmat(nx[1], nx[0], A, nx_, &sBAbt0, nu[0], 0);
	drowin_libstr(nx[1], 1.0, &sb0, 0, &sBAbt0, nu[0]+nx[0], 0);
#if PRINT
	d_print_strmat(nu[0]+nx[0]+1, nx[1], &sBAbt0, 0, 0);
#endif

	struct d_strmat sBAbt1;
	struct d_strvec sb1;
	if(N>1)
		{
		d_allocate_strmat(nu[1]+nx[1]+1, nx[2], &sBAbt1);
		d_cvt_tran_mat2strmat(nx[2], nu[1], B, nx_, &sBAbt1, 0, 0);
		d_cvt_tran_mat2strmat(nx[2], nx[1], A, nx_, &sBAbt1, nu[1], 0);
		d_cvt_tran_mat2strmat(nx[2], 1, b, nx_, &sBAbt1, nu[1]+nx[1], 0);
#if PRINT
		d_print_strmat(nu[1]+nx[1]+1, nx[2], &sBAbt1, 0, 0);
#endif
		d_allocate_strvec(nx_, &sb1);
		d_cvt_vec2strvec(nx_, b, &sb1, 0);
		}

/************************************************
* cost function
************************************************/	
	
	double *Q; d_zeros(&Q, nx_, nx_);
	for(ii=0; ii<nx_; ii++) Q[ii*(nx_+1)] = 1.0;

	double *R; d_zeros(&R, nu_, nu_);
	for(ii=0; ii<nu_; ii++) R[ii*(nu_+1)] = 2.0;

	double *S; d_zeros(&S, nu_, nx_); // S=0, so no need to update r0

	double *q; d_zeros(&q, nx_, 1);
	for(ii=0; ii<nx_; ii++) q[ii] = 0.1;

	double *r; d_zeros(&r, nu_, 1);
	for(ii=0; ii<nu_; ii++) r[ii] = 0.2;

	struct d_strmat sRSQrq0;
	struct d_strvec srq0;
	d_allocate_strmat(nu[0]+nx[0]+1, nu[0]+nx[0], &sRSQrq0);
	d_cvt_mat2strmat(nu[0], nu[0], R, nu_, &sRSQrq0, 0, 0);
	d_cvt_tran_mat2strmat(nu[0], nx[0], S, nu_, &sRSQrq0, nu[0], 0);
	d_cvt_mat2strmat(nx[0], nx[0], Q, nx_, &sRSQrq0, nu[0], nu[0]);
	d_cvt_tran_mat2strmat(nu[0], 1, r, nu_, &sRSQrq0, nu[0]+nx[0], 0);
	d_cvt_tran_mat2strmat(nx[0], 1, q, nx_, &sRSQrq0, nu[0]+nx[0], nu[0]);
#if PRINT
	d_print_strmat(nu[0]+nx[0]+1, nu[0]+nx[0], &sRSQrq0, 0, 0);
#endif
	d_allocate_strvec(nu[0]+nx[0], &srq0);
	d_cvt_vec2strvec(nu[0], r, &srq0, 0);
	d_cvt_vec2strvec(nx[0], q, &srq0, nu[0]);
#if PRINT
	d_print_tran_strvec(nu[0]+nx[0], &srq0, 0);
#endif

	struct d_strmat sRSQrq1;
	struct d_strvec srq1;
	if(N>1)
		{
		d_allocate_strmat(nu[1]+nx[1]+1, nu[1]+nx[1], &sRSQrq1);
		d_cvt_mat2strmat(nu[1], nu[1], R, nu_, &sRSQrq1, 0, 0);
		d_cvt_tran_mat2strmat(nu[1], nx[1], S, nu_, &sRSQrq1, nu[1], 0);
		d_cvt_mat2strmat(nx[1], nx[1], Q, nx_, &sRSQrq1, nu[1], nu[1]);
		d_cvt_tran_mat2strmat(nu[1], 1, r, nu_, &sRSQrq1, nu[1]+nx[1], 0);
		d_cvt_tran_mat2strmat(nx[1], 1, q, nx_, &sRSQrq1, nu[1]+nx[1], nu[1]);
#if PRINT
		d_print_strmat(nu[1]+nx[1]+1, nu[1]+nx[1], &sRSQrq1, 0, 0);
#endif
		d_allocate_strvec(nu[1]+nx[1], &srq1);
		d_cvt_vec2strvec(nu[1], r, &srq1, 0);
		d_cvt_vec2strvec(nx[1], q, &srq1, nu[1]);
#if PRINT
		d_print_tran_strvec(nu[1]+nx[1], &srq1, 0);
#endif
		}

	struct d_strmat sRSQrqN;
	struct d_strvec srqN;
	d_allocate_strmat(nu[N]+nx[N]+1, nu[N]+nx[N], &sRSQrqN);
	d_cvt_mat2strmat(nu[N], nu[N], R, nu_, &sRSQrqN, 0, 0);
	d_cvt_tran_mat2strmat(nu[N], nx[N], S, nu_, &sRSQrqN, nu[N], 0);
	d_cvt_mat2strmat(nx[N], nx[N], Q, nx_, &sRSQrqN, nu[N], nu[N]);
	d_cvt_tran_mat2strmat(nu[N], 1, r, nu_, &sRSQrqN, nu[N]+nx[N], 0);
	d_cvt_tran_mat2strmat(nx[N], 1, q, nx_, &sRSQrqN, nu[N]+nx[N], nu[N]);
#if PRINT
	d_print_strmat(nu[N]+nx[N]+1, nu[N]+nx[N], &sRSQrqN, 0, 0);
#endif
	d_allocate_strvec(nu[N]+nx[N], &srqN);
	d_cvt_vec2strvec(nu[N], r, &srqN, 0);
	d_cvt_vec2strvec(nx[N], q, &srqN, nu[N]);
#if PRINT
	d_print_tran_strvec(nu[N]+nx[N], &srqN, 0);
#endif

	// maximum element in cost functions
	double mu0 = 2.0;

/************************************************
* box & general constraints
************************************************/	

	int *idxb0; int_zeros(&idxb0, nb[0], 1);
	double *d0; d_zeros(&d0, 2*nb[0]+2*ng[0], 1);
	for(ii=0; ii<nb[0]; ii++)
		{
//		if(0) // input
		if(ii<nu[0]) // input
			{
			d0[ii]             = - 0.5; // umin
			d0[nb[0]+ng[0]+ii] =   0.5; // umax
			}
		else // state
			{
			d0[ii]             = - 4.0; // xmin
			d0[nb[0]+ng[0]+ii] =   4.0; // xmax
			}
//		idxb0[ii] = nu[0]+nx[0]/2+ii;
		idxb0[ii] = ii;
		}
	for(ii=0; ii<ng[0]; ii++)
		{
		if(ii<nu[0]) // input
			{
			d0[nb[0]+ii]             = - 0.5; // umin
			d0[nb[0]+ng[0]+nb[0]+ii] =   0.5; // umax
			}
		else // state
			{
			d0[nb[0]+ii]             = - 4.0; // xmin
			d0[nb[0]+ng[0]+nb[0]+ii] =   4.0; // xmax
			}
		}
#if PRINT
	int_print_mat(1, nb[0], idxb0, 1);
	d_print_mat(1, 2*nb[0]+2*ng[0], d0, 1);
#endif

	int *idxb1; int_zeros(&idxb1, nb[1], 1);
	double *d1; d_zeros(&d1, 2*nb[1]+2*ng[1], 1);
	for(ii=0; ii<nb[1]; ii++)
		{
//		if(0) // input
		if(ii<nu[1]) // input
			{
			d1[ii]             = - 0.5; // umin
			d1[nb[1]+ng[1]+ii] =   0.5; // umax
			}
		else // state
			{
			d1[ii]             = - 4.0; // xmin
			d1[nb[1]+ng[1]+ii] =   4.0; // xmax
			}
		idxb1[ii] = ii;
//		idxb1[ii] = nu[1]+nx[1]/2+ii;
		}
	for(ii=0; ii<ng[1]; ii++)
		{
		if(ii<nu[1]) // input
			{
			d1[nb[1]+ii]             = - 0.5; // umin
			d1[nb[1]+ng[1]+nb[1]+ii] =   0.5; // umax
			}
		else // state
			{
			d1[nb[1]+ii]             = - 4.0; // xmin
			d1[nb[1]+ng[1]+nb[1]+ii] =   4.0; // xmax
			}
		}
#if PRINT
	int_print_mat(1, nb[1], idxb1, 1);
	d_print_mat(1, 2*nb[1]+2*ng[1], d1, 1);
#endif

	int *idxbN; int_zeros(&idxbN, nb[N], 1);
	double *dN; d_zeros(&dN, 2*nb[N]+2*ng[N], 1);
	for(ii=0; ii<nb[N]; ii++)
		{
		dN[ii]             = - 4.0; // xmin
		dN[nb[N]+ng[N]+ii] =   4.0; // xmax
		idxbN[ii] = ii;
		}
	for(ii=0; ii<ng[N]; ii++)
		{
		dN[nb[N]+ii]             = - 4.0; // dmin
		dN[nb[N]+ng[N]+nb[N]+ii] =   4.0; // dmax
		}
#if PRINT
	int_print_mat(1, nb[N], idxbN, 1);
	d_print_mat(1, 2*nb[N]+2*ng[N], dN, 1);
#endif

	double *DC0; d_zeros(&DC0, ng[0], nu[0]+nx[0]);
	for(ii=0; ii<ng[0]; ii++)
		DC0[ii*(ng[0]+1)] = 1.0;
	double *DC1; d_zeros(&DC1, ng[1], nu[1]+nx[1]);
	for(ii=0; ii<ng[1]; ii++)
		DC1[ii*(ng[1]+1)] = 1.0;
	double *DCN; d_zeros(&DCN, ng[N], nx[N]);
	for(ii=0; ii<ng[N]; ii++)
		DCN[ii*(ng[N]+1)] = 1.0;

	struct d_strmat sDCt0;
	d_allocate_strmat(nu[0]+nx[0], ng[0], &sDCt0);
	d_cvt_tran_mat2strmat(ng[0], nu[0]+nx[0], DC0, ng[0], &sDCt0, 0, 0);
#if PRINT
	d_print_strmat(nu[0]+nx[0], ng[0], &sDCt0, 0, 0);
#endif
	struct d_strvec slb0;
	struct d_strvec sub0;
	struct d_strvec slg0;
	struct d_strvec sug0;
	d_allocate_strvec(nb[0], &slb0);
	d_allocate_strvec(nb[0], &sub0);
	d_allocate_strvec(ng[0], &slg0);
	d_allocate_strvec(ng[0], &sug0);
	d_cvt_vec2strvec(nb[0], d0+0, &slb0, 0);
	d_cvt_vec2strvec(nb[0], d0+nb[0], &sub0, 0);
	d_cvt_vec2strvec(ng[0], d0+2*nb[0], &slg0, 0);
	d_cvt_vec2strvec(ng[0], d0+2*nb[0]+ng[0], &sug0, 0);
#if PRINT
	d_print_tran_strvec(nb[0], &slb0, 0);
	d_print_tran_strvec(nb[0], &sub0, 0);
	d_print_tran_strvec(ng[0], &slg0, 0);
	d_print_tran_strvec(ng[0], &sug0, 0);
#endif

	struct d_strmat sDCt1;
	d_allocate_strmat(nu[1]+nx[1], ng[1], &sDCt1);
	d_cvt_tran_mat2strmat(ng[1], nu[1]+nx[1], DC1, ng[1], &sDCt1, 0, 0);
#if PRINT
	d_print_strmat(nu[1]+nx[1], ng[1], &sDCt1, 0, 0);
#endif
	struct d_strvec slb1;
	struct d_strvec sub1;
	struct d_strvec slg1;
	struct d_strvec sug1;
	d_allocate_strvec(nb[1], &slb1);
	d_allocate_strvec(nb[1], &sub1);
	d_allocate_strvec(ng[1], &slg1);
	d_allocate_strvec(ng[1], &sug1);
	d_cvt_vec2strvec(nb[1], d1+0, &slb1, 0);
	d_cvt_vec2strvec(nb[1], d1+nb[1], &sub1, 0);
	d_cvt_vec2strvec(ng[1], d1+2*nb[1], &slg1, 0);
	d_cvt_vec2strvec(ng[1], d1+2*nb[1]+ng[1], &sug1, 0);
#if PRINT
	d_print_tran_strvec(nb[1], &slb1, 0);
	d_print_tran_strvec(nb[1], &sub1, 0);
	d_print_tran_strvec(ng[1], &slg1, 0);
	d_print_tran_strvec(ng[1], &sug1, 0);
#endif

	struct d_strmat sDCtN;
	d_allocate_strmat(nx[N], ng[N], &sDCtN);
	d_cvt_tran_mat2strmat(ng[N], nx[N], DCN, ng[N], &sDCtN, 0, 0);
#if PRINT
	d_print_strmat(nx[N], ng[N], &sDCtN, 0, 0);
#endif
	struct d_strvec slbN;
	struct d_strvec subN;
	struct d_strvec slgN;
	struct d_strvec sugN;
	d_allocate_strvec(nb[N], &slbN);
	d_allocate_strvec(nb[N], &subN);
	d_allocate_strvec(ng[N], &slgN);
	d_allocate_strvec(ng[N], &sugN);
	d_cvt_vec2strvec(nb[N], dN+0, &slbN, 0);
	d_cvt_vec2strvec(nb[N], dN+nb[N], &subN, 0);
	d_cvt_vec2strvec(ng[N], dN+2*nb[N], &slgN, 0);
	d_cvt_vec2strvec(ng[N], dN+2*nb[N]+ng[N], &sugN, 0);
#if PRINT
	d_print_tran_strvec(nb[N], &slbN, 0);
	d_print_tran_strvec(nb[N], &subN, 0);
	d_print_tran_strvec(ng[N], &slgN, 0);
	d_print_tran_strvec(ng[N], &sugN, 0);
#endif

/************************************************
* libstr ip2 solver
************************************************/	

	struct d_strmat hsBAbt[N];
	struct d_strvec hsb[N];
	struct d_strmat hsRSQrq[N+1];
	struct d_strvec hsrq[N+1];
	struct d_strmat hsDCt[N+1];
	struct d_strvec hslb[N+1];
	struct d_strvec hsub[N+1];
	struct d_strvec hslg[N+1];
	struct d_strvec hsug[N+1];
	int *hidxb[N+1];
//	struct d_strvec hsux[N+1];
//	struct d_strvec hspi[N+1];
//	struct d_strvec hslam[N+1];
//	struct d_strvec hst[N+1];

	hsBAbt[0] = sBAbt0;
	hsb[0] = sb0;
	hsRSQrq[0] = sRSQrq0;
	hsrq[0] = srq0;
	hsDCt[0] = sDCt0;
	hslb[0] = slb0;
	hsub[0] = sub0;
	hslg[0] = slg0;
	hsug[0] = sug0;
	hidxb[0] = idxb0;
//	d_allocate_strvec(nu[0]+nx[0], &hsux[0]);
//	d_allocate_strvec(nx[1], &hspi[1]);
//	d_allocate_strvec(2*nb[0]+2*ng[0], &hslam[0]);
//	d_allocate_strvec(2*nb[0]+2*ng[0], &hst[0]);
	for(ii=1; ii<N; ii++)
		{
		hsBAbt[ii] = sBAbt1;
		hsb[ii] = sb1;
		hsRSQrq[ii] = sRSQrq1;
		hsrq[ii] = srq1;
		hsDCt[ii] = sDCt1;
		hslb[ii] = slb1;
		hsub[ii] = sub1;
		hslg[ii] = slg1;
		hsug[ii] = sug1;
		hidxb[ii] = idxb1;
//		d_allocate_strvec(nu[ii]+nx[ii], &hsux[ii]);
//		d_allocate_strvec(nx[ii+1], &hspi[ii+1]);
//		d_allocate_strvec(2*nb[ii]+2*ng[ii], &hslam[ii]);
//		d_allocate_strvec(2*nb[ii]+2*ng[ii], &hst[ii]);
		}
	hsRSQrq[N] = sRSQrqN;
	hsrq[N] = srqN;
	hsDCt[N] = sDCtN;
	hslb[N] = slbN;
	hsub[N] = subN;
	hslg[N] = slgN;
	hsug[N] = sugN;
	hidxb[N] = idxbN;
//	d_allocate_strvec(nu[N]+nx[N], &hsux[N]);
//	d_allocate_strvec(2*nb[N]+2*ng[N], &hslam[N]);
//	d_allocate_strvec(2*nb[N]+2*ng[N], &hst[N]);
	
/************************************************
* ocp qp structure
************************************************/	

	struct d_ocp_qp str_in;
	d_cast_ocp_qp(N, nx, nu, nb, hidxb, ng, hsBAbt, hsb, hsRSQrq, hsrq, hsDCt, hslb, hsub, hslg, hsug, &str_in);

	d_print_strmat(str_in.sBAbt[0].m, str_in.sBAbt[0].n, str_in.sBAbt+0, 0, 0);
	d_print_strmat(str_in.sBAbt[1].m, str_in.sBAbt[1].n, str_in.sBAbt+1, 0, 0);

	int size_out = d_size_ocp_qp(N, nx, nu, nb, ng);
	printf("\nocp qp size = %d\n\n", size_out);

	void *mem_out;
	v_zeros_align(&mem_out, size_out);

	struct d_ocp_qp str_out;
	d_create_ocp_qp(N, nx, nu, nb, ng, &str_out, mem_out);
	d_copy_ocp_qp(&str_in, &str_out);

	d_print_strmat(str_out.sBAbt[0].m, str_out.sBAbt[0].n, str_out.sBAbt+0, 0, 0);
	d_print_strmat(str_out.sBAbt[1].m, str_out.sBAbt[1].n, str_out.sBAbt+1, 0, 0);

/************************************************
* free memory
************************************************/	

	d_free(A);
	d_free(B);
	d_free(b);
	d_free(x0);
	int_free(idxb0);
	d_free(d0);
	int_free(idxb1);
	d_free(d1);
	int_free(idxbN);
	d_free(dN);

	d_free_strmat(&sA);
	d_free_strvec(&sx0);
	d_free_strmat(&sBAbt0);
	d_free_strvec(&sb0);
	d_free_strmat(&sRSQrq0);
	d_free_strvec(&srq0);
	d_free_strmat(&sRSQrqN);
	d_free_strvec(&srqN);
	d_free_strmat(&sDCt0);
	d_free_strvec(&slb0);
	d_free_strvec(&sub0);
	d_free_strvec(&slg0);
	d_free_strvec(&sug0);
	d_free_strmat(&sDCt1);
	d_free_strvec(&slb1);
	d_free_strvec(&sub1);
	d_free_strvec(&slg1);
	d_free_strvec(&sug1);
	d_free_strmat(&sDCtN);
	d_free_strvec(&slbN);
	d_free_strvec(&subN);
	d_free_strvec(&slgN);
	d_free_strvec(&sugN);
	if(N>1)
		{
		d_free_strmat(&sBAbt1);
		d_free_strvec(&sb1);
		d_free_strmat(&sRSQrq1);
		d_free_strvec(&srq1);
		}
//	d_free_strvec(&hsux[0]);
//	d_free_strvec(&hspi[1]);
//	d_free_strvec(&hslam[0]);
//	d_free_strvec(&hst[0]);
//	d_free_strvec(&hsrrq[0]);
//	d_free_strvec(&hsrb[0]);
	for(ii=1; ii<N; ii++)
		{
//		d_free_strvec(&hsux[ii]);
//		d_free_strvec(&hspi[ii+1]);
//		d_free_strvec(&hslam[ii]);
//		d_free_strvec(&hst[ii]);
//		d_free_strvec(&hsrrq[ii]);
//		d_free_strvec(&hsrb[ii]);
		}
//	d_free_strvec(&hsux[N]);
//	d_free_strvec(&hslam[N]);
//	d_free_strvec(&hst[N]);
//	d_free_strvec(&hsrrq[N]);

	v_free_align(mem_out);

/************************************************
* return
************************************************/	

	return 0;

	}

	
