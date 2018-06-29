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
#include <math.h>
#include <sys/time.h>

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_v_aux_ext_dep.h>
#include <blasfeo_d_aux_ext_dep.h>
#include <blasfeo_i_aux_ext_dep.h>
#include <blasfeo_d_aux.h>
#include <blasfeo_d_blas.h>

#include "../include/hpipm_d_ocp_qp_dim.h"
#include "../include/hpipm_d_ocp_qp.h"
#include "../include/hpipm_d_ocp_qp_sol.h"
#include "../include/hpipm_d_dense_qp.h"
#include "../include/hpipm_d_dense_qp_sol.h"
#include "../include/hpipm_d_dense_qp_res.h"
#include "../include/hpipm_d_dense_qp_ipm.h"
#include "../include/hpipm_d_cond.h"

#include "d_tools.h"



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


	// local variables

	int ii, jj;

	int rep, nrep=1000;

	struct timeval tv0, tv1;



	// problem size

	int nx_ = 8; // number of states (it has to be even for the mass-spring system test problem)
	int nu_ = 3; // number of inputs (controllers) (it has to be at least 1 and at most nx/2 for the mass-spring system test problem)
	int N  = 8; // horizon lenght



	// stage-wise variant size

	int nx[N+1];
#if KEEP_X0
	nx[0] = nx_;
#else
	nx[0] = 0;
#endif
	for(ii=1; ii<=N; ii++)
		nx[ii] = nx_;
//	nx[N] = 0;

	int nu[N+1];
	for(ii=0; ii<N; ii++)
		nu[ii] = nu_;
	nu[N] = 0;

	int nbu[N+1];
	for (ii=0; ii<=N; ii++)
		nbu[ii] = nu[ii];

#if 1
	int nbx[N+1];
#if KEEP_X0
	nbx[0] = nx[0]/2;
#else
	nbx[0] = 0;
#endif
	for(ii=1; ii<=N; ii++)
		nbx[ii] = nx[ii]/2;

	int nb[N+1];
	for (ii=0; ii<=N; ii++)
		nb[ii] = nbu[ii]+nbx[ii];

	int ng[N+1];
	ng[0] = 0;
	for(ii=1; ii<N; ii++)
		ng[ii] = 0;
	ng[N] = 0;

	int ns[N+1];
	ns[0] = 0;
	for(ii=1; ii<N; ii++)
		ns[ii] = nx[ii]/2;
	ns[N] = nx[N]/2;

	int nsbx[N+1];
	for(ii=0; ii<=N; ii++)
		nsbx[ii] = ns[ii];

	int nsg[N+1];
	int nsbu[N+1];
	for(ii=0; ii<=N; ii++)
	{
		nsg[ii] = 0;
		nsbu[ii] = 0;
	}
//	ns[N] = 0;
#elif 0
	int nb[N+1];
	nb[0] = 0;
	for(ii=1; ii<N; ii++)
		nb[ii] = 0;
	nb[N] = 0;

	int ng[N+1];
#if KEEP_X0
	ng[0] = nu[0]+nx[0]/2;
#else
	ng[0] = nu[0];
#endif
	for(ii=1; ii<N; ii++)
		ng[ii] = nu[1]+nx[1]/2;
	ng[N] = nx[N]/2;
#else
	int nb[N+1];
	nb[0] = nu[0] + nx[0]/2;
	for(ii=1; ii<N; ii++)
		nb[ii] = nu[ii] + nx[ii]/2;
	nb[N] = nu[N] + nx[N]/2;

	int ng[N+1];
#if KEEP_X0
	ng[0] = nx[0]/2;
#else
	ng[0] = 0;
#endif
	for(ii=1; ii<N; ii++)
		ng[ii] = nx[1]/2;
	ng[N] = nx[N]/2;
#endif

/************************************************
* dynamical system
************************************************/

	double *A; d_zeros(&A, nx_, nx_); // states update matrix

	double *B; d_zeros(&B, nx_, nu_); // inputs matrix

	double *b; d_zeros(&b, nx_, 1); // states offset
	double *x0; d_zeros(&x0, nx_, 1); // initial state

	double Ts = 0.5; // sampling time
	mass_spring_system(Ts, nx_, nu_, N, A, B, b, x0);

	for(jj=0; jj<nx_; jj++)
		b[jj] = 0.0;

	for(jj=0; jj<nx_; jj++)
		x0[jj] = 0;
	x0[0] = 2.5;
	x0[1] = 2.5;

	double *b0; d_zeros(&b0, nx_, 1);
	dgemv_n_3l(nx_, nx_, A, nx_, x0, b0);
	daxpy_3l(nx_, 1.0, b, b0);

#if PRINT
	d_print_mat(nx_, nx_, A, nx_);
	d_print_mat(nx_, nu_, B, nu_);
	d_print_mat(1, nx_, b, 1);
	d_print_mat(1, nx_, x0, 1);
	d_print_mat(1, nx_, b0, 1);
#endif

/************************************************
* cost function
************************************************/

	double *Q; d_zeros(&Q, nx_, nx_);
	for(ii=0; ii<nx_; ii++) Q[ii*(nx_+1)] = 0.0;

	double *R; d_zeros(&R, nu_, nu_);
	for(ii=0; ii<nu_; ii++) R[ii*(nu_+1)] = 2.0;

	double *S; d_zeros(&S, nu_, nx_);

	double *q; d_zeros(&q, nx_, 1);
	for(ii=0; ii<nx_; ii++) q[ii] = 0.0;

	double *r; d_zeros(&r, nu_, 1);
	for(ii=0; ii<nu_; ii++) r[ii] = 0.0;

	double *r0; d_zeros(&r0, nu_, 1);
	dgemv_n_3l(nu_, nx_, S, nu_, x0, r0);
	daxpy_3l(nu_, 1.0, r, r0);

#if PRINT
	d_print_mat(nx_, nx_, Q, nx_);
	d_print_mat(nu_, nu_, R, nu_);
	d_print_mat(nu_, nx_, S, nu_);
	d_print_mat(1, nx_, q, 1);
	d_print_mat(1, nu_, r, 1);
	d_print_mat(1, nu_, r0, 1);
#endif

	// maximum element in cost functions
	double mu0;
	if(ns[1]>0 | ns[N]>0)
		mu0 = 1000.0;
	else
		mu0 = 2.0;


/************************************************
* box & general constraints
************************************************/

	int *idxb0; int_zeros(&idxb0, nb[0], 1);
	double *d_lb0; d_zeros(&d_lb0, nb[0], 1);
	double *d_ub0; d_zeros(&d_ub0, nb[0], 1);
	double *d_lg0; d_zeros(&d_lg0, ng[0], 1);
	double *d_ug0; d_zeros(&d_ug0, ng[0], 1);
	for(ii=0; ii<nb[0]; ii++)
		{
		if(ii<nu[0]) // input
			{
			d_lb0[ii] = - 0.5; // umin
			d_ub0[ii] =   0.5; // umax
			}
		else // state
			{
			d_lb0[ii] = - 4.0; // xmin
			d_ub0[ii] =   4.0; // xmax
			}
		idxb0[ii] = ii;
		}
	for(ii=0; ii<ng[0]; ii++)
		{
		if(ii<nu[0]-nb[0]) // input
			{
			d_lg0[ii] = - 0.5; // umin
			d_ug0[ii] =   0.5; // umax
			}
		else // state
			{
			d_lg0[ii] = - 4.0; // xmin
			d_ug0[ii] =   4.0; // xmax
			}
		}

	int *idxb1; int_zeros(&idxb1, nb[1], 1);
	double *d_lb1; d_zeros(&d_lb1, nb[1], 1);
	double *d_ub1; d_zeros(&d_ub1, nb[1], 1);
	double *d_lg1; d_zeros(&d_lg1, ng[1], 1);
	double *d_ug1; d_zeros(&d_ug1, ng[1], 1);
	for(ii=0; ii<nb[1]; ii++)
		{
		if(ii<nu[1]) // input
			{
			d_lb1[ii] = - 0.5; // umin
			d_ub1[ii] =   0.5; // umax
			}
		else // state
			{
			d_lb1[ii] = - 1.0; // xmin
			d_ub1[ii] =   1.0; // xmax
			}
		idxb1[ii] = ii;
		}
	for(ii=0; ii<ng[1]; ii++)
		{
		if(ii<nu[1]-nb[1]) // input
			{
			d_lg1[ii] = - 0.5; // umin
			d_ug1[ii] =   0.5; // umax
			}
		else // state
			{
			d_lg1[ii] = - 4.0; // xmin
			d_ug1[ii] =   4.0; // xmax
			}
		}


	int *idxbN; int_zeros(&idxbN, nb[N], 1);
	double *d_lbN; d_zeros(&d_lbN, nb[N], 1);
	double *d_ubN; d_zeros(&d_ubN, nb[N], 1);
	double *d_lgN; d_zeros(&d_lgN, ng[N], 1);
	double *d_ugN; d_zeros(&d_ugN, ng[N], 1);
	for(ii=0; ii<nb[N]; ii++)
		{
		d_lbN[ii] = - 1.0; // xmin
		d_ubN[ii] =   1.0; // xmax
		idxbN[ii] = ii;
		}
	for(ii=0; ii<ng[N]; ii++)
		{
		d_lgN[ii] = - 4.0; // dmin
		d_ugN[ii] =   4.0; // dmax
		}

	double *C0; d_zeros(&C0, ng[0], nx[0]);
	double *D0; d_zeros(&D0, ng[0], nu[0]);
	for(ii=0; ii<nu[0]-nb[0] & ii<ng[0]; ii++)
		D0[ii+(nb[0]+ii)*ng[0]] = 1.0;
	for(; ii<ng[0]; ii++)
		C0[ii+(nb[0]+ii-nu[0])*ng[0]] = 1.0;

	double *C1; d_zeros(&C1, ng[1], nx[1]);
	double *D1; d_zeros(&D1, ng[1], nu[1]);
	for(ii=0; ii<nu[1]-nb[1] & ii<ng[1]; ii++)
		D1[ii+(nb[1]+ii)*ng[1]] = 1.0;
	for(; ii<ng[1]; ii++)
		C1[ii+(nb[1]+ii-nu[1])*ng[1]] = 1.0;

	double *CN; d_zeros(&CN, ng[N], nx[N]);
	double *DN; d_zeros(&DN, ng[N], nu[N]);
	for(ii=0; ii<nu[N]-nb[N] & ii<ng[N]; ii++)
		DN[ii+(nb[N]+ii)*ng[N]] = 1.0;
	for(; ii<ng[N]; ii++)
		CN[ii+(nb[N]+ii-nu[N])*ng[N]] = 1.0;

#if PRINT
	// box constraints
	int_print_mat(1, nb[0], idxb0, 1);
	d_print_mat(1, nb[0], d_lb0, 1);
	d_print_mat(1, nb[0], d_ub0, 1);
	int_print_mat(1, nb[1], idxb1, 1);
	d_print_mat(1, nb[1], d_lb1, 1);
	d_print_mat(1, nb[1], d_ub1, 1);
	int_print_mat(1, nb[N], idxbN, 1);
	d_print_mat(1, nb[N], d_lbN, 1);
	d_print_mat(1, nb[N], d_ubN, 1);
	// general constraints
	d_print_mat(1, ng[0], d_lg0, 1);
	d_print_mat(1, ng[0], d_ug0, 1);
	d_print_mat(ng[0], nu[0], D0, ng[0]);
	d_print_mat(ng[0], nx[0], C0, ng[0]);
	d_print_mat(1, ng[1], d_lg1, 1);
	d_print_mat(1, ng[1], d_ug1, 1);
	d_print_mat(ng[1], nu[1], D1, ng[1]);
	d_print_mat(ng[1], nx[1], C1, ng[1]);
	d_print_mat(1, ng[N], d_lgN, 1);
	d_print_mat(1, ng[N], d_ugN, 1);
	d_print_mat(ng[N], nu[N], DN, ng[N]);
	d_print_mat(ng[N], nx[N], CN, ng[N]);
#endif

/************************************************
* soft constraints
************************************************/

	double *Zl0; d_zeros(&Zl0, ns[0], 1);
	for(ii=0; ii<ns[0]; ii++)
		Zl0[ii] = 1e3;
	double *Zu0; d_zeros(&Zu0, ns[0], 1);
	for(ii=0; ii<ns[0]; ii++)
		Zu0[ii] = 1e3;
	double *zl0; d_zeros(&zl0, ns[0], 1);
	for(ii=0; ii<ns[0]; ii++)
		zl0[ii] = 0e2;
	double *zu0; d_zeros(&zu0, ns[0], 1);
	for(ii=0; ii<ns[0]; ii++)
		zu0[ii] = 1e2;
	int *idxs0; int_zeros(&idxs0, ns[0], 1);
	for(ii=0; ii<ns[0]; ii++)
		idxs0[ii] = nu[0]+ii;
	double *d_ls0; d_zeros(&d_ls0, ns[0], 1);
	for(ii=0; ii<ns[0]; ii++)
		d_ls0[ii] = -1.0;
	double *d_us0; d_zeros(&d_us0, ns[0], 1);
	for(ii=0; ii<ns[0]; ii++)
		d_us0[ii] = 0.0;

	double *Zl1; d_zeros(&Zl1, ns[1], 1);
	for(ii=0; ii<ns[1]; ii++)
		Zl1[ii] = 1e3;
	double *Zu1; d_zeros(&Zu1, ns[1], 1);
	for(ii=0; ii<ns[1]; ii++)
		Zu1[ii] = 1e3;
	double *zl1; d_zeros(&zl1, ns[1], 1);
	for(ii=0; ii<ns[1]; ii++)
		zl1[ii] = 0e2;
	double *zu1; d_zeros(&zu1, ns[1], 1);
	for(ii=0; ii<ns[1]; ii++)
		zu1[ii] = 1e2;
	int *idxs1; int_zeros(&idxs1, ns[1], 1);
	for(ii=0; ii<ns[1]; ii++)
		idxs1[ii] = nu[1]+ii;
	double *d_ls1; d_zeros(&d_ls1, ns[1], 1);
	for(ii=0; ii<ns[1]; ii++)
		d_ls1[ii] = -1.0;
	double *d_us1; d_zeros(&d_us1, ns[1], 1);
	for(ii=0; ii<ns[1]; ii++)
		d_us1[ii] = 0.0;

	double *ZlN; d_zeros(&ZlN, ns[N], 1);
	for(ii=0; ii<ns[N]; ii++)
		ZlN[ii] = 1e3;
	double *ZuN; d_zeros(&ZuN, ns[N], 1);
	for(ii=0; ii<ns[N]; ii++)
		ZuN[ii] = 1e3;
	double *zlN; d_zeros(&zlN, ns[N], 1);
	for(ii=0; ii<ns[N]; ii++)
		zlN[ii] = 0e2;
	double *zuN; d_zeros(&zuN, ns[N], 1);
	for(ii=0; ii<ns[N]; ii++)
		zuN[ii] = 1e2;
	int *idxsN; int_zeros(&idxsN, ns[N], 1);
	for(ii=0; ii<ns[N]; ii++)
		idxsN[ii] = nu[N]+ii;
	double *d_lsN; d_zeros(&d_lsN, ns[N], 1);
	for(ii=0; ii<ns[N]; ii++)
		d_lsN[ii] = -1.0;
	double *d_usN; d_zeros(&d_usN, ns[N], 1);
	for(ii=0; ii<ns[N]; ii++)
		d_usN[ii] = 0.0;

#if 1
	// soft constraints
	int_print_mat(1, ns[0], idxs0, 1);
	d_print_mat(1, ns[0], Zl0, 1);
	d_print_mat(1, ns[0], Zu0, 1);
	d_print_mat(1, ns[0], zl0, 1);
	d_print_mat(1, ns[0], zu0, 1);
	d_print_mat(1, ns[0], d_ls0, 1);
	d_print_mat(1, ns[0], d_us0, 1);
	int_print_mat(1, ns[1], idxs1, 1);
	d_print_mat(1, ns[1], Zl1, 1);
	d_print_mat(1, ns[1], Zu1, 1);
	d_print_mat(1, ns[1], zl1, 1);
	d_print_mat(1, ns[1], zu1, 1);
	d_print_mat(1, ns[1], d_ls1, 1);
	d_print_mat(1, ns[1], d_us1, 1);
	int_print_mat(1, ns[N], idxsN, 1);
	d_print_mat(1, ns[N], ZlN, 1);
	d_print_mat(1, ns[N], ZuN, 1);
	d_print_mat(1, ns[N], zlN, 1);
	d_print_mat(1, ns[N], zuN, 1);
	d_print_mat(1, ns[N], d_lsN, 1);
	d_print_mat(1, ns[N], d_usN, 1);
#endif

/************************************************
* array of matrices
************************************************/

	double *hA[N];
	double *hB[N];
	double *hb[N];
	double *hQ[N+1];
	double *hS[N+1];
	double *hR[N+1];
	double *hq[N+1];
	double *hr[N+1];
	int *hidxb[N+1];
	double *hd_lb[N+1];
	double *hd_ub[N+1];
	double *hC[N+1];
	double *hD[N+1];
	double *hd_lg[N+1];
	double *hd_ug[N+1];
	double *hZl[N+1];
	double *hZu[N+1];
	double *hzl[N+1];
	double *hzu[N+1];
	int *hidxs[N+1]; // XXX
	double *hd_ls[N+1];
	double *hd_us[N+1];

	hA[0] = A;
	hB[0] = B;
	hb[0] = b0;
	hQ[0] = Q;
	hS[0] = S;
	hR[0] = R;
	hq[0] = q;
	hr[0] = r0;
	hidxb[0] = idxb0;
	hd_lb[0] = d_lb0;
	hd_ub[0] = d_ub0;
	hC[0] = C0;
	hD[0] = D0;
	hd_lg[0] = d_lg0;
	hd_ug[0] = d_ug0;
	hZl[0] = Zl0;
	hZu[0] = Zu0;
	hzl[0] = zl0;
	hzu[0] = zu0;
	hidxs[0] = idxs0;
	hd_ls[0] = d_ls0;
	hd_us[0] = d_us0;
	for(ii=1; ii<N; ii++)
		{
		hA[ii] = A;
		hB[ii] = B;
		hb[ii] = b;
		hQ[ii] = Q;
		hS[ii] = S;
		hR[ii] = R;
		hq[ii] = q;
		hr[ii] = r;
		hidxb[ii] = idxb1;
		hd_lb[ii] = d_lb1;
		hd_ub[ii] = d_ub1;
		hd_lg[ii] = d_lg1;
		hd_ug[ii] = d_ug1;
		hC[ii] = C1;
		hD[ii] = D1;
		hZl[ii] = Zl1;
		hZu[ii] = Zu1;
		hzl[ii] = zl1;
		hzu[ii] = zu1;
		hidxs[ii] = idxs1;
		hd_ls[ii] = d_ls1;
		hd_us[ii] = d_us1;
		}
	hQ[N] = Q;
	hS[N] = S;
	hR[N] = R;
	hq[N] = q;
	hr[N] = r;
	hidxb[N] = idxbN;
	hd_lb[N] = d_lbN;
	hd_ub[N] = d_ubN;
	hd_lg[N] = d_lgN;
	hd_ug[N] = d_ugN;
	hC[N] = CN;
	hD[N] = DN;
	hZl[N] = ZlN;
	hZu[N] = ZuN;
	hzl[N] = zlN;
	hzu[N] = zuN;
	hidxs[N] = idxsN;
	hd_ls[N] = d_lsN;
	hd_us[N] = d_usN;

/************************************************
* ocp qp dim
************************************************/

	int dim_size = d_memsize_ocp_qp_dim(N);
	printf("\ndim size = %d\n", dim_size);
	void *dim_mem = malloc(dim_size);

	struct d_ocp_qp_dim dim;
	d_create_ocp_qp_dim(N, &dim, dim_mem);
	d_cvt_int_to_ocp_qp_dim(N, nx, nu, nbx, nbu, ng, nsbx, nsbu, nsg, &dim);

/************************************************
* ocp qp
************************************************/

	int ocp_qp_size = d_memsize_ocp_qp(&dim);
	printf("\nqp size = %d\n", ocp_qp_size);
	void *ocp_qp_mem = malloc(ocp_qp_size);

	struct d_ocp_qp ocp_qp;
	d_create_ocp_qp(&dim, &ocp_qp, ocp_qp_mem);
	d_cvt_colmaj_to_ocp_qp(hA, hB, hb, hQ, hS, hR, hq, hr, hidxb, hd_lb, hd_ub, hC, hD, hd_lg, hd_ug, hZl, hZu, hzl, hzu, hidxs, hd_ls, hd_us, &ocp_qp);

#if 0
	printf("\nN = %d\n", ocp_qp.N);
	for(ii=0; ii<N; ii++)
		blasfeo_print_dmat(ocp_qp.nu[ii]+ocp_qp.nx[ii]+1, ocp_qp.nx[ii+1], ocp_qp.BAbt+ii, 0, 0);
	for(ii=0; ii<N; ii++)
		blasfeo_print_tran_dvec(ocp_qp.nx[ii+1], ocp_qp.b+ii, 0);
	for(ii=0; ii<=N; ii++)
		blasfeo_print_dmat(ocp_qp.nu[ii]+ocp_qp.nx[ii]+1, ocp_qp.nu[ii]+ocp_qp.nx[ii], ocp_qp.RSQrq+ii, 0, 0);
	for(ii=0; ii<=N; ii++)
		blasfeo_print_tran_dvec(ocp_qp.nu[ii]+ocp_qp.nx[ii], ocp_qp.rq+ii, 0);
	for(ii=0; ii<=N; ii++)
		int_print_mat(1, nb[ii], ocp_qp.idxb[ii], 1);
	for(ii=0; ii<=N; ii++)
		blasfeo_print_tran_dvec(ocp_qp.nb[ii], ocp_qp.d_lb+ii, 0);
	for(ii=0; ii<=N; ii++)
		blasfeo_print_tran_dvec(ocp_qp.nb[ii], ocp_qp.d_ub+ii, 0);
	for(ii=0; ii<=N; ii++)
		blasfeo_print_dmat(ocp_qp.nu[ii]+ocp_qp.nx[ii], ocp_qp.ng[ii], ocp_qp.DCt+ii, 0, 0);
	for(ii=0; ii<=N; ii++)
		blasfeo_print_tran_dvec(ocp_qp.ng[ii], ocp_qp.d_lg+ii, 0);
	for(ii=0; ii<=N; ii++)
		blasfeo_print_tran_dvec(ocp_qp.ng[ii], ocp_qp.d_ug+ii, 0);
#endif

/************************************************
* dense qp dim
************************************************/

	int dense_qp_dim_size = d_memsize_dense_qp_dim();
	printf("\nqp dim size = %d\n", dense_qp_dim_size);
	void *dense_qp_dim_mem = malloc(dense_qp_dim_size);

	struct d_dense_qp_dim qp_dim;
	d_create_dense_qp_dim(&qp_dim, dense_qp_dim_mem);

	d_compute_qp_dim_ocp2dense(&dim, &qp_dim);

	int nvc = qp_dim.nv;
	int nec = qp_dim.ne;
	int nbc = qp_dim.nb;
	int ngc = qp_dim.ng;
	int nsc = qp_dim.ns;

	printf("\nnv = %d, ne = %d, nb = %d, ng = %d, ns = %d\n\n", nvc, nec, nbc, ngc, nsc);

/************************************************
* dense qp
************************************************/

	// qp
	int dense_qp_size = d_memsize_dense_qp(&qp_dim);
	printf("\nqp size = %d\n", dense_qp_size);
	void *dense_qp_mem = malloc(dense_qp_size);

	struct d_dense_qp dense_qp;
	d_create_dense_qp(&qp_dim, &dense_qp, dense_qp_mem);

	// arg
	int cond_arg_size = d_memsize_cond_qp_ocp2dense_arg();
	printf("\ncond_arg size = %d\n", cond_arg_size);
	void *cond_arg_mem = malloc(cond_arg_size);

	struct d_cond_qp_ocp2dense_arg cond_arg;
	d_create_cond_qp_ocp2dense_arg(&cond_arg, cond_arg_mem);
	d_set_default_cond_qp_ocp2dense_arg(&cond_arg);

	// ws
	int cond_size = d_memsize_cond_qp_ocp2dense(&dim, &cond_arg);
	printf("\ncond size = %d\n", cond_size);
	void *cond_mem = malloc(cond_size);

	struct d_cond_qp_ocp2dense_workspace cond_ws;
	d_create_cond_qp_ocp2dense(&dim, &cond_arg, &cond_ws, cond_mem);

	/* cond */

	gettimeofday(&tv0, NULL); // start

	for(rep=0; rep<nrep; rep++)
		{
		d_cond_qp_ocp2dense(&ocp_qp, &dense_qp, &cond_arg, &cond_ws);
		}

	gettimeofday(&tv1, NULL); // stop

	double time_cond = (tv1.tv_sec-tv0.tv_sec)/(nrep+0.0)+(tv1.tv_usec-tv0.tv_usec)/(nrep*1e6);

#if 1
	printf("\ncond data\n\n");
	blasfeo_print_dmat(nvc, nvc, dense_qp.Hv, 0, 0);
	blasfeo_print_dmat(nec, nvc, dense_qp.A, 0, 0);
	blasfeo_print_dmat(nvc, ngc, dense_qp.Ct, 0, 0);
	blasfeo_print_tran_dvec(nvc, dense_qp.gz, 0);
	blasfeo_print_tran_dvec(nec, dense_qp.b, 0);
	blasfeo_print_tran_dvec(2*nbc+2*ngc+2*nsc, dense_qp.d, 0);
	blasfeo_print_tran_dvec(nbc, dense_qp.d, 0);
	blasfeo_print_tran_dvec(nbc, dense_qp.d, nbc+ngc);
	blasfeo_print_tran_dvec(ngc, dense_qp.d, nbc);
	blasfeo_print_tran_dvec(ngc, dense_qp.d, 2*nbc+ngc);
	blasfeo_print_tran_dvec(nsc, dense_qp.d, 2*nbc+2*ngc);
	blasfeo_print_tran_dvec(nsc, dense_qp.d, 2*nbc+2*ngc+nsc);
	blasfeo_print_tran_dvec(nsc, dense_qp.Z, 0);
	blasfeo_print_tran_dvec(nsc, dense_qp.Z, nsc);
	blasfeo_print_tran_dvec(nsc, dense_qp.gz, nvc);
	blasfeo_print_tran_dvec(nsc, dense_qp.gz, nvc+nsc);
#endif

	/* update cond */

	// index of updated dynamics
	int idxc[N];
	for(ii=0; ii<N; ii++)
		idxc[ii] = 0;
	idxc[0] = 1;

//	blasfeo_dgese(nvc+1, nvc, 0.0, dense_qp.Hv, 0, 0);

	gettimeofday(&tv0, NULL); // start

	for(rep=0; rep<nrep; rep++)
		{
		d_update_cond_qp_ocp2dense(idxc, &ocp_qp, &dense_qp, &cond_arg, &cond_ws);
		}

	gettimeofday(&tv1, NULL); // stop

	double time_update_cond = (tv1.tv_sec-tv0.tv_sec)/(nrep+0.0)+(tv1.tv_usec-tv0.tv_usec)/(nrep*1e6);

#if 1
	printf("\nupdate cond data\n\n");
	blasfeo_print_dmat(nvc, nvc, dense_qp.Hv, 0, 0);
	blasfeo_print_dmat(nec, nvc, dense_qp.A, 0, 0);
	blasfeo_print_dmat(nvc, ngc, dense_qp.Ct, 0, 0);
	blasfeo_print_tran_dvec(nvc, dense_qp.gz, 0);
	blasfeo_print_tran_dvec(nec, dense_qp.b, 0);
	blasfeo_print_tran_dvec(2*nbc+2*ngc+2*nsc, dense_qp.d, 0);
	blasfeo_print_tran_dvec(nbc, dense_qp.d, 0);
	blasfeo_print_tran_dvec(nbc, dense_qp.d, nbc+ngc);
	blasfeo_print_tran_dvec(ngc, dense_qp.d, nbc);
	blasfeo_print_tran_dvec(ngc, dense_qp.d, 2*nbc+ngc);
	blasfeo_print_tran_dvec(nsc, dense_qp.d, 2*nbc+2*ngc);
	blasfeo_print_tran_dvec(nsc, dense_qp.d, 2*nbc+2*ngc+nsc);
	blasfeo_print_tran_dvec(nsc, dense_qp.Z, 0);
	blasfeo_print_tran_dvec(nsc, dense_qp.Z, nsc);
	blasfeo_print_tran_dvec(nsc, dense_qp.gz, nvc);
	blasfeo_print_tran_dvec(nsc, dense_qp.gz, nvc+nsc);
#endif
//return 0;

	/* cond rhs */

	gettimeofday(&tv0, NULL); // start

	for(rep=0; rep<nrep; rep++)
		{
		d_cond_rhs_qp_ocp2dense(&ocp_qp, &dense_qp, &cond_arg, &cond_ws);
		}

	gettimeofday(&tv1, NULL); // stop

	double time_cond_rhs = (tv1.tv_sec-tv0.tv_sec)/(nrep+0.0)+(tv1.tv_usec-tv0.tv_usec)/(nrep*1e6);

#if 1
	printf("\ncond rhs data\n\n");
	blasfeo_print_tran_dvec(nvc, dense_qp.gz, 0);
	blasfeo_print_tran_dvec(nec, dense_qp.b, 0);
	blasfeo_print_tran_dvec(2*nbc+2*ngc+2*nsc, dense_qp.d, 0);
	blasfeo_print_tran_dvec(nbc, dense_qp.d, 0);
	blasfeo_print_tran_dvec(nbc, dense_qp.d, nbc+ngc);
	blasfeo_print_tran_dvec(ngc, dense_qp.d, nbc);
	blasfeo_print_tran_dvec(ngc, dense_qp.d, 2*nbc+ngc);
	blasfeo_print_tran_dvec(nsc, dense_qp.d, 2*nbc+2*ngc);
	blasfeo_print_tran_dvec(nsc, dense_qp.d, 2*nbc+2*ngc+nsc);
	blasfeo_print_tran_dvec(nsc, dense_qp.Z, 0);
	blasfeo_print_tran_dvec(nsc, dense_qp.Z, nsc);
	blasfeo_print_tran_dvec(nsc, dense_qp.gz, nvc);
	blasfeo_print_tran_dvec(nsc, dense_qp.gz, nvc+nsc);
#endif

#if 0
	int nu_tmp = 0;
	for(ii=0; ii<N; ii++)
		{
		nu_tmp += nu[ii];
		blasfeo_print_dmat(nu_tmp+nx[0]+1, nx[ii+1], cond_ws.Gamma+ii, 0, 0);
		}
#endif

/************************************************
* dense qp sol
************************************************/

	int dense_qp_sol_size = d_memsize_dense_qp_sol(&qp_dim);
	printf("\ndense qp sol size = %d\n", dense_qp_sol_size);
	void *dense_qp_sol_mem = malloc(dense_qp_sol_size);

	struct d_dense_qp_sol dense_qp_sol;
	d_create_dense_qp_sol(&qp_dim, &dense_qp_sol, dense_qp_sol_mem);

/************************************************
* ipm arg
************************************************/

	int ipm_arg_size = d_memsize_dense_qp_ipm_arg(&qp_dim);
	printf("\nipm arg size = %d\n", ipm_arg_size);
	void *ipm_arg_mem = malloc(ipm_arg_size);

	struct d_dense_qp_ipm_arg dense_arg;
	d_create_dense_qp_ipm_arg(&qp_dim, &dense_arg, ipm_arg_mem);
//	enum dense_qp_ipm_mode mode = SPEED_ABS;
	enum dense_qp_ipm_mode mode = SPEED;
//	enum dense_qp_ipm_mode mode = BALANCE;
//	enum dense_qp_ipm_mode mode = ROBUST;
	d_set_default_dense_qp_ipm_arg(mode, &dense_arg);

//	arg.alpha_min = 1e-8;
//	arg.res_g_max = 1e-8;
//	arg.res_b_max = 1e-8;
//	arg.res_d_max = 1e-12;
//	arg.res_m_max = 1e-12;
//	dense_arg.mu0 = mu0;
//	arg.iter_max = 20;
//	arg.stat_max = 100;
//	arg.pred_corr = 1;

/************************************************
* ipm
************************************************/

	int dense_ipm_size = d_memsize_dense_qp_ipm(&qp_dim, &dense_arg);
	printf("\ndense ipm size = %d\n", dense_ipm_size);
	void *dense_ipm_mem = malloc(dense_ipm_size);

	struct d_dense_qp_ipm_workspace dense_workspace;
	d_create_dense_qp_ipm(&qp_dim, &dense_arg, &dense_workspace, dense_ipm_mem);

	int hpipm_return; // 0 normal; 1 max iter

	gettimeofday(&tv0, NULL); // start

	for(rep=0; rep<nrep; rep++)
		{
		hpipm_return = d_solve_dense_qp_ipm(&dense_qp, &dense_qp_sol, &dense_arg, &dense_workspace);
		}

	gettimeofday(&tv1, NULL); // stop

	double time_dense_ipm = (tv1.tv_sec-tv0.tv_sec)/(nrep+0.0)+(tv1.tv_usec-tv0.tv_usec)/(nrep*1e6);


	printf("\nsolution\n\n");
	printf("\nv\n");
	blasfeo_print_tran_dvec(nvc, dense_qp_sol.v, 0);
	printf("\nls\n");
	blasfeo_print_tran_dvec(nsc, dense_qp_sol.v, nvc);
	printf("\nus\n");
	blasfeo_print_tran_dvec(nsc, dense_qp_sol.v, nvc+nsc);
	printf("\npi\n");
	blasfeo_print_tran_dvec(nec, dense_qp_sol.pi, 0);
	printf("\nlam_lb\n");
	blasfeo_print_tran_dvec(nbc, dense_qp_sol.lam, 0);
	printf("\nlam_ub\n");
	blasfeo_print_tran_dvec(nbc, dense_qp_sol.lam, nbc+ngc);
	printf("\nlam_lg\n");
	blasfeo_print_tran_dvec(ngc, dense_qp_sol.lam, nbc);
	printf("\nlam_ug\n");
	blasfeo_print_tran_dvec(ngc, dense_qp_sol.lam, 2*nbc+ngc);
	printf("\nlam_ls\n");
	blasfeo_print_tran_dvec(nsc, dense_qp_sol.lam, 2*nbc+2*ngc);
	printf("\nlam_us\n");
	blasfeo_print_tran_dvec(nsc, dense_qp_sol.lam, 2*nbc+2*ngc+nsc);
	printf("\nt_lb\n");
	blasfeo_print_tran_dvec(nbc, dense_qp_sol.t, 0);
	printf("\nt_ub\n");
	blasfeo_print_tran_dvec(nbc, dense_qp_sol.t, nbc+ngc);
	printf("\nt_lg\n");
	blasfeo_print_tran_dvec(ngc, dense_qp_sol.t, nbc);
	printf("\nt_ug\n");
	blasfeo_print_tran_dvec(ngc, dense_qp_sol.t, 2*nbc+ngc);
	printf("\nt_ls\n");
	blasfeo_print_tran_dvec(nsc, dense_qp_sol.t, 2*nbc+2*ngc);
	printf("\nt_us\n");
	blasfeo_print_tran_dvec(nsc, dense_qp_sol.t, 2*nbc+2*ngc+nsc);

	printf("\nresiduals\n\n");
	printf("\nres_g\n");
	blasfeo_print_exp_tran_dvec(nvc+2*nsc, dense_workspace.res->res_g, 0);
	printf("\nres_b\n");
	blasfeo_print_exp_tran_dvec(nec, dense_workspace.res->res_b, 0);
	printf("\nres_d\n");
	blasfeo_print_exp_tran_dvec(2*nbc+2*ngc+2*nsc, dense_workspace.res->res_d, 0);
	printf("\nres_m\n");
	blasfeo_print_exp_tran_dvec(2*nbc+2*ngc+2*nsc, dense_workspace.res->res_m, 0);
	printf("\nres_mu\n");
	printf("\n%e\n\n", dense_workspace.res->res_mu);

/************************************************
* print ipm statistics
************************************************/

	printf("\nipm return = %d\n", hpipm_return);
	printf("\nipm residuals max: res_g = %e, res_b = %e, res_d = %e, res_m = %e\n", dense_workspace.qp_res[0], dense_workspace.qp_res[1], dense_workspace.qp_res[2], dense_workspace.qp_res[3]);

	printf("\nipm iter = %d\n", dense_workspace.iter);
	printf("\nalpha_aff\tmu_aff\t\tsigma\t\talpha\t\tmu\n");
	d_print_e_tran_mat(5, dense_workspace.iter, dense_workspace.stat, 5);

	printf("\ncond time           = %e [s]\n", time_cond);
	printf("\nupdate cond time    = %e [s]\n", time_update_cond);
	printf("\ncond rhs time       = %e [s]\n", time_cond_rhs);
	printf("\ncond dense ipm time = %e [s]\n\n", time_dense_ipm);

/************************************************
* ocp qp sol
************************************************/

	int ocp_qp_sol_size = d_memsize_ocp_qp_sol(&dim);
	printf("\nocp qp sol size = %d\n", ocp_qp_sol_size);
	void *ocp_qp_sol_mem = malloc(ocp_qp_sol_size);

	struct d_ocp_qp_sol ocp_qp_sol;
	d_create_ocp_qp_sol(&dim, &ocp_qp_sol, ocp_qp_sol_mem);

/************************************************
* expand solution
************************************************/

	d_expand_sol_dense2ocp(&ocp_qp, &dense_qp_sol, &ocp_qp_sol, &cond_arg, &cond_ws);

	double *u[N+1]; for(ii=0; ii<=N; ii++) d_zeros(u+ii, nu[ii], 1);
	double *x[N+1]; for(ii=0; ii<=N; ii++) d_zeros(x+ii, nx[ii], 1);
	double *ls[N+1]; for(ii=0; ii<=N; ii++) d_zeros(ls+ii, ns[ii], 1);
	double *us[N+1]; for(ii=0; ii<=N; ii++) d_zeros(us+ii, ns[ii], 1);
	double *pi[N]; for(ii=0; ii<N; ii++) d_zeros(pi+ii, nx[ii+1], 1);
	double *lam_lb[N+1]; for(ii=0; ii<=N; ii++) d_zeros(lam_lb+ii, nb[ii], 1);
	double *lam_ub[N+1]; for(ii=0; ii<=N; ii++) d_zeros(lam_ub+ii, nb[ii], 1);
	double *lam_lg[N+1]; for(ii=0; ii<=N; ii++) d_zeros(lam_lg+ii, ng[ii], 1);
	double *lam_ug[N+1]; for(ii=0; ii<=N; ii++) d_zeros(lam_ug+ii, ng[ii], 1);
	double *lam_ls[N+1]; for(ii=0; ii<=N; ii++) d_zeros(lam_ls+ii, ns[ii], 1);
	double *lam_us[N+1]; for(ii=0; ii<=N; ii++) d_zeros(lam_us+ii, ns[ii], 1);

	d_cvt_ocp_qp_sol_to_colmaj(&ocp_qp_sol, u, x, ls, us, pi, lam_lb, lam_ub, lam_lg, lam_ug, lam_ls, lam_us);

#if 1
	printf("\nfull space solution\n\n");
	printf("\nu\n");
	for(ii=0; ii<=N; ii++)
		d_print_mat(1, nu[ii], u[ii], 1);
	printf("\nx\n");
	for(ii=0; ii<=N; ii++)
		d_print_mat(1, nx[ii], x[ii], 1);
	printf("\nls\n");
	for(ii=0; ii<=N; ii++)
		d_print_mat(1, ns[ii], ls[ii], 1);
	printf("\nus\n");
	for(ii=0; ii<=N; ii++)
		d_print_mat(1, ns[ii], us[ii], 1);
	printf("\npi\n");
	for(ii=0; ii<N; ii++)
		d_print_mat(1, nx[ii+1], pi[ii], 1);
	printf("\nlam_lb\n");
	for(ii=0; ii<=N; ii++)
		d_print_mat(1, nb[ii], lam_lb[ii], 1);
	printf("\nlam_ub\n");
	for(ii=0; ii<=N; ii++)
		d_print_mat(1, nb[ii], lam_ub[ii], 1);
	printf("\nlam_lg\n");
	for(ii=0; ii<=N; ii++)
		d_print_mat(1, ng[ii], lam_lg[ii], 1);
	printf("\nlam_ug\n");
	for(ii=0; ii<=N; ii++)
		d_print_mat(1, ng[ii], lam_ug[ii], 1);
	printf("\nlam_ls\n");
	for(ii=0; ii<=N; ii++)
		d_print_mat(1, ns[ii], lam_ls[ii], 1);
	printf("\nlam_us\n");
	for(ii=0; ii<=N; ii++)
		d_print_mat(1, ns[ii], lam_us[ii], 1);

	printf("\nt_lb\n");
	for(ii=0; ii<=N; ii++)
		d_print_mat(1, nb[ii], (ocp_qp_sol.t+ii)->pa, 1);
	printf("\nt_ub\n");
	for(ii=0; ii<=N; ii++)
		d_print_mat(1, nb[ii], (ocp_qp_sol.t+ii)->pa+nb[ii]+ng[ii], 1);
	printf("\nt_lg\n");
	for(ii=0; ii<=N; ii++)
		d_print_mat(1, ng[ii], (ocp_qp_sol.t+ii)->pa+nb[ii], 1);
	printf("\nt_ug\n");
	for(ii=0; ii<=N; ii++)
		d_print_mat(1, ng[ii], (ocp_qp_sol.t+ii)->pa+2*nb[ii]+ng[ii], 1);
	printf("\nt_ls\n");
	for(ii=0; ii<=N; ii++)
		d_print_mat(1, ns[ii], (ocp_qp_sol.t+ii)->pa+2*nb[ii]+2*ng[ii], 1);
	printf("\nt_us\n");
	for(ii=0; ii<=N; ii++)
		d_print_mat(1, ns[ii], (ocp_qp_sol.t+ii)->pa+2*nb[ii]+2*ng[ii]+ns[ii], 1);
#endif

/************************************************
* free memory
************************************************/

	d_free(A);
	d_free(B);
	d_free(b);
	d_free(x0);
	d_free(Q);
	d_free(R);
	d_free(S);
	d_free(q);
	d_free(r);
	d_free(r0);
	int_free(idxb0);
	d_free(d_lb0);
	d_free(d_ub0);
	int_free(idxb1);
	d_free(d_lb1);
	d_free(d_ub1);
	int_free(idxbN);
	d_free(d_lbN);
	d_free(d_ubN);
	d_free(C0);
	d_free(D0);
	d_free(d_lg0);
	d_free(d_ug0);
	d_free(C1);
	d_free(D1);
	d_free(d_lg1);
	d_free(d_ug1);
	d_free(CN);
	d_free(DN);
	d_free(d_lgN);
	d_free(d_ugN);

	d_free(Zl0);
	d_free(Zu0);
	d_free(zl0);
	d_free(zu0);
	int_free(idxs0);
	d_free(d_ls0);
	d_free(d_us0);
	d_free(Zl1);
	d_free(Zu1);
	d_free(zl1);
	d_free(zu1);
	int_free(idxs1);
	d_free(d_ls1);
	d_free(d_us1);
	d_free(ZlN);
	d_free(ZuN);
	d_free(zlN);
	d_free(zuN);
	int_free(idxsN);
	d_free(d_lsN);
	d_free(d_usN);

	free(ocp_qp_mem);
	free(ocp_qp_sol_mem);
	free(dense_qp_mem);
	free(dense_qp_sol_mem);
	free(cond_mem);
	free(cond_arg_mem);
	free(dense_ipm_mem);

/************************************************
* return
************************************************/

	return 0;

	}
