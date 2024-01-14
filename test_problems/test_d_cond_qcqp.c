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
#include <math.h>
#include <sys/time.h>

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_v_aux_ext_dep.h>
#include <blasfeo_d_aux_ext_dep.h>
#include <blasfeo_i_aux_ext_dep.h>
#include <blasfeo_d_aux.h>
#include <blasfeo_d_blas.h>

#include <hpipm_d_ocp_qcqp_dim.h>
#include <hpipm_d_ocp_qcqp.h>
#include <hpipm_d_ocp_qcqp_sol.h>
#include <hpipm_d_ocp_qcqp_utils.h>
#include <hpipm_d_ocp_qcqp_red.h>
#include <hpipm_d_dense_qcqp.h>
#include <hpipm_d_dense_qcqp_sol.h>
#include <hpipm_d_dense_qcqp_res.h>
#include <hpipm_d_dense_qcqp_ipm.h>
#include <hpipm_d_dense_qcqp_utils.h>
#include <hpipm_d_cond_qcqp.h>
#include <hpipm_timing.h>

#include "d_tools.h"



// printing
#ifndef PRINT
#define PRINT 1
#endif



// remove initial state x0 from optimization variables
#define REMOVE_X0 1



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

	int ii, jj, kk, ll;

	int hpipm_status; // 0 normal; 1 max iter; 2 min alpha; 3 NaN

	int rep, nrep=1000;

	hpipm_timer timer;

/************************************************
* problem size
************************************************/

	int nx_ = 4; // number of states (it has to be even for the mass-spring system test problem)
	int nu_ = 1; // number of inputs (controllers) (it has to be at least 1 and at most nx/2 for the mass-spring system test problem)
	int N  = 15; // horizon lenght



	// stage-wise variant size

	int nx[N+1];
	for(ii=0; ii<=N; ii++)
		nx[ii] = nx_;

	int nu[N+1];
	for(ii=0; ii<N; ii++)
		nu[ii] = nu_;
	nu[N] = 0;

	int nbu[N+1];
	for (ii=0; ii<=N; ii++)
//		nbu[ii] = 0;
		nbu[ii] = nu[ii];

	int nbx[N+1];
	nbx[0] = nx[0];
	for(ii=1; ii<=N; ii++)
		nbx[ii] = 0;

	// mark initial state bounds as equalities to be removed
	int nbxe[N+1];
#if REMOVE_X0
	nbxe[0] = nbx[0]; // remove all states at stage 0
#else
	nbxe[0] = 0; // do not remove anything
#endif
	for(ii=1; ii<=N; ii++)
		nbxe[ii] = 0;

	int nb[N+1];
	for (ii=0; ii<=N; ii++)
		nb[ii] = nbu[ii]+nbx[ii];

	int ng[N+1];
	for(ii=0; ii<N; ii++)
		ng[ii] = 0;
	ng[N] = 0;

	int nq[N+1];
	for(ii=0; ii<N; ii++)
//		nq[ii] = 1;
		nq[ii] = 0;
	nq[N] = 1;

	int nsbx[N+1];
	for(ii=0; ii<=N; ii++)
		nsbx[ii] = 0;

	int nsbu[N+1];
	for(ii=0; ii<=N; ii++)
		nsbu[ii] = 0;

	int nsg[N+1];
	for(ii=0; ii<=N; ii++)
		nsg[ii] = 0;

	int nsq[N+1];
	for(ii=0; ii<N; ii++)
		nsq[ii] = 0;
	nsq[N] = 1;

	int ns[N+1];
	for(ii=0; ii<=N; ii++)
		ns[ii] = nsbx[ii] + nsbu[ii] + nsg[ii] + nsq[ii];

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

#if PRINT
	d_print_mat(nx_, nx_, A, nx_);
	d_print_mat(nx_, nu_, B, nu_);
	d_print_mat(1, nx_, b, 1);
	d_print_mat(1, nx_, x0, 1);
#endif

/************************************************
* cost function
************************************************/

	double *Q; d_zeros(&Q, nx_, nx_);
	for(ii=0; ii<nx_; ii++) Q[ii*(nx_+1)] = 1.0;

	double *R; d_zeros(&R, nu_, nu_);
	for(ii=0; ii<nu_; ii++) R[ii*(nu_+1)] = 2.0;

	double *S; d_zeros(&S, nu_, nx_);

	double *q; d_zeros(&q, nx_, 1);
	for(ii=0; ii<nx_; ii++) q[ii] = 0.0;

	double *r; d_zeros(&r, nu_, 1);
	for(ii=0; ii<nu_; ii++) r[ii] = 0.0;

#if PRINT
	d_print_mat(nx_, nx_, Q, nx_);
	d_print_mat(nu_, nu_, R, nu_);
	d_print_mat(nu_, nx_, S, nu_);
	d_print_mat(1, nx_, q, 1);
	d_print_mat(1, nu_, r, 1);
#endif

/************************************************
* box & general constraints
************************************************/

	// initial state
	int *idxbx0; int_zeros(&idxbx0, nx[0], 1);
	double *lbx0; d_zeros(&lbx0, nx[0], 1);
	double *ubx0; d_zeros(&ubx0, nx[0], 1);
	int *idxbxe0; int_zeros(&idxbxe0, nx[0], 1);
	for(ii=0; ii<nx[0]; ii++)
		{
		lbx0[ii] = x0[ii]; // xmin
		ubx0[ii] = x0[ii]; // xmax
		idxbx0[ii] = ii;
#if REMOVE_X0
		idxbxe0[ii] = ii; // mark state bounds as equalities to be removed
#endif
		}

	int *idxbx1; int_zeros(&idxbx1, nbx[1], 1);
	double *lbx1; d_zeros(&lbx1, nbx[1], 1);
	double *ubx1; d_zeros(&ubx1, nbx[1], 1);
	int *idxbu1; int_zeros(&idxbu1, nbu[1], 1);
	double *lbu1; d_zeros(&lbu1, nbu[1], 1);
	double *ubu1; d_zeros(&ubu1, nbu[1], 1);
	double *lg1; d_zeros(&lg1, ng[1], 1);
	double *ug1; d_zeros(&ug1, ng[1], 1);
	for(ii=0; ii<nbu[1]; ii++)
		{
		lbu1[ii] = - 0.5; // umin
		ubu1[ii] =   0.5; // umax
		idxbu1[ii] = ii;
		}
	for(ii=0; ii<nbx[1]; ii++)
		{
		lbx1[ii] = - 4.0; // xmin
		ubx1[ii] =   4.0; // xmax
		idxbx1[ii] = ii;
		}
	for(ii=0; ii<ng[1]; ii++)
		{
		if(ii<nu[1]-nb[1]) // input
			{
			lg1[ii] = - 0.5; // umin
			ug1[ii] =   0.5; // umax
			}
		else // state
			{
			lg1[ii] = - 4.0; // xmin
			ug1[ii] =   4.0; // xmax
			}
		}


	int *idxbxN; int_zeros(&idxbxN, nbx[N], 1);
	double *lbxN; d_zeros(&lbxN, nbx[N], 1);
	double *ubxN; d_zeros(&ubxN, nbx[N], 1);
	double *lgN; d_zeros(&lgN, ng[N], 1);
	double *ugN; d_zeros(&ugN, ng[N], 1);
	for(ii=0; ii<nbx[N]; ii++)
		{
		lbxN[ii] = - 4.0; // xmin
		ubxN[ii] =   4.0; // xmax
		idxbxN[ii] = ii;
		}
	for(ii=0; ii<ng[N]; ii++)
		{
		lgN[ii] = - 4.0; // dmin
		ugN[ii] =   4.0; // dmax
		}

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
	int_print_mat(1, nbx[0], idxbx0, 1);
	d_print_mat(1, nbx[0], lbx0, 1);
	d_print_mat(1, nbx[0], ubx0, 1);
	int_print_mat(1, nx[0], idxbxe0, 1);
	int_print_mat(1, nbx[1], idxbx1, 1);
	d_print_mat(1, nbx[1], lbx1, 1);
	d_print_mat(1, nbx[1], ubx1, 1);
	int_print_mat(1, nbu[1], idxbu1, 1);
	d_print_mat(1, nbu[1], lbu1, 1);
	d_print_mat(1, nbu[1], ubu1, 1);
	int_print_mat(1, nbx[N], idxbxN, 1);
	d_print_mat(1, nbx[N], lbxN, 1);
	d_print_mat(1, nbx[N], ubxN, 1);
	// general constraints
	d_print_mat(1, ng[1], lg1, 1);
	d_print_mat(1, ng[1], ug1, 1);
	d_print_mat(ng[1], nu[1], D1, ng[1]);
	d_print_mat(ng[1], nx[1], C1, ng[1]);
	d_print_mat(1, ng[N], lgN, 1);
	d_print_mat(1, ng[N], ugN, 1);
	d_print_mat(ng[N], nu[N], DN, ng[N]);
	d_print_mat(ng[N], nx[N], CN, ng[N]);
#endif

/************************************************
* quadratic constraints
************************************************/
	
	double *Rq1; d_zeros(&Rq1, nu[1], nu[1]*nq[1]);
	if(nq[1]>0)
		for(ii=0; ii<nu[1]; ii++)
			Rq1[ii*(nu[1]+1)] = 2*1.0;

	double *Rq2; d_zeros(&Rq2, nu[2], nu[2]*nq[2]);
	if(nq[2]>0)
		for(ii=0; ii<nu[2]; ii++)
			Rq2[ii*(nu[2]+1)] = 2*1.0/2;

	double *Rq3; d_zeros(&Rq3, nu[3], nu[3]*nq[3]);
	if(nq[3]>0)
		for(ii=0; ii<nu[3]; ii++)
			Rq3[ii*(nu[3]+1)] = 2*1.0/3;

	double *Qq1; d_zeros(&Qq1, nx[1], nx[1]*nq[1]);
//	for(ii=0; ii<nx[1]/2; ii++)
//		Qq1[(nx[1]/2+ii)*(nx[1]+1)] = 0.0;
//	d_print_mat(nx[1], nx[1], Qq1, nx[1]);

	double *qq1; d_zeros(&qq1, nx[1], nq[1]);
//	qq1[0*nx[1]+0] = -1;
//	qq1[1*nx[1]+0] =  1;

	double *uq1; d_zeros(&uq1, nq[1], 1);
	if(nq[1]>0)
		uq1[0] =  0.5*0.5;
//	uq1[0] =  4.0;
//	uq1[1] =  4.0;

	double *uq1_mask; d_zeros(&uq1_mask, nq[1], 1);
//	uq1_mask[0] = 1.0;
//	uq1_mask[1] = 1.0;


//	double *QqNm1; d_zeros(&QqNm1, nx[N-1], nx[N-1]*nq[N-1]);
//	for(ii=0; ii<nx[N-1]; ii++)
//		QqNm1[ii*(nx[N-1]+1)] = 1.0;

//	double *qqNm1; d_zeros(&qqNm1, nx[N-1], nq[N-1]);
//	qqN[0*nx[N]+0] = 0.0;

//	double *uqNm1; d_zeros(&uqNm1, nq[N-1], 1);
//	uqNm1[0] = 1.5;

//	double *uqNm1_mask; d_zeros(&uqNm1_mask, nq[N-1], 1);
//	uqNm1_mask[0] = 1.0;


	double *QqN; d_zeros(&QqN, nx[N], nx[N]*nq[N]);
	if(nq[N]>0)
		for(ii=0; ii<nx[N]; ii++)
			QqN[ii*(nx[N]+1)] = 1.0;

	double *qqN; d_zeros(&qqN, nx[N], nq[N]);
//	qqN[0*nx[N]+0] = 0.0;

	double *uqN; d_zeros(&uqN, nq[N], 1);
	if(nq[N]>0)
		uqN[0] = 1.5;

	double *uqN_mask; d_zeros(&uqN_mask, nq[N], 1);
	if(nq[N]>0)
		uqN_mask[0] = 1.0;

/************************************************
* soft constraints
************************************************/

#if 0
	double *Zl0; d_zeros(&Zl0, ns[0], 1);
	for(ii=0; ii<ns[0]; ii++)
		Zl0[ii] = 0e3;
	double *Zu0; d_zeros(&Zu0, ns[0], 1);
	for(ii=0; ii<ns[0]; ii++)
		Zu0[ii] = 0e3;
	double *zl0; d_zeros(&zl0, ns[0], 1);
	for(ii=0; ii<ns[0]; ii++)
		zl0[ii] = 1e2;
	double *zu0; d_zeros(&zu0, ns[0], 1);
	for(ii=0; ii<ns[0]; ii++)
		zu0[ii] = 1e2;
	int *idxs0; int_zeros(&idxs0, ns[0], 1);
	for(ii=0; ii<ns[0]; ii++)
		idxs0[ii] = nu[0]+ii;
	double *d_ls0; d_zeros(&d_ls0, ns[0], 1);
	for(ii=0; ii<ns[0]; ii++)
		d_ls0[ii] = 0.0; //-1.0;
	double *d_us0; d_zeros(&d_us0, ns[0], 1);
	for(ii=0; ii<ns[0]; ii++)
		d_us0[ii] = 0.0;

	double *Zl1; d_zeros(&Zl1, ns[1], 1);
	for(ii=0; ii<ns[1]; ii++)
		Zl1[ii] = 0e3;
	double *Zu1; d_zeros(&Zu1, ns[1], 1);
	for(ii=0; ii<ns[1]; ii++)
		Zu1[ii] = 0e3;
	double *zl1; d_zeros(&zl1, ns[1], 1);
	for(ii=0; ii<ns[1]; ii++)
		zl1[ii] = 1e2;
	double *zu1; d_zeros(&zu1, ns[1], 1);
	for(ii=0; ii<ns[1]; ii++)
		zu1[ii] = 1e2;
	int *idxs1; int_zeros(&idxs1, ns[1], 1);
	for(ii=0; ii<ns[1]; ii++)
		idxs1[ii] = nu[1]+ii;
	double *d_ls1; d_zeros(&d_ls1, ns[1], 1);
	for(ii=0; ii<ns[1]; ii++)
		d_ls1[ii] = 0.0; //-1.0;
	double *d_us1; d_zeros(&d_us1, ns[1], 1);
	for(ii=0; ii<ns[1]; ii++)
		d_us1[ii] = 0.0;
#endif

	double *ZlN; d_zeros(&ZlN, ns[N], 1);
	for(ii=0; ii<ns[N]; ii++)
		ZlN[ii] = 1e2;
	double *ZuN; d_zeros(&ZuN, ns[N], 1);
	for(ii=0; ii<ns[N]; ii++)
		ZuN[ii] = 1e2;
	double *zlN; d_zeros(&zlN, ns[N], 1);
	for(ii=0; ii<ns[N]; ii++)
		zlN[ii] = 0e2;
	double *zuN; d_zeros(&zuN, ns[N], 1);
	for(ii=0; ii<ns[N]; ii++)
		zuN[ii] = 0e2;
	int *idxsN; int_zeros(&idxsN, ns[N], 1);
	for(ii=0; ii<ns[N]; ii++)
		idxsN[ii] = nb[N]+ng[N]+ii; // softed quadratic constr
	double *llsN; d_zeros(&llsN, ns[N], 1);
	for(ii=0; ii<ns[N]; ii++)
		llsN[ii] = 0.0; //-1.0;
	double *lusN; d_zeros(&lusN, ns[N], 1);
	for(ii=0; ii<ns[N]; ii++)
		lusN[ii] = 0.0;

#if 0
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
	d_print_mat(1, ns[N], llsN, 1);
	d_print_mat(1, ns[N], lusN, 1);
#endif

/************************************************
* ocp qp dim
************************************************/

	hpipm_size_t ocp_dim_size = d_ocp_qcqp_dim_memsize(N);
#if PRINT
	printf("\ndim size = %d\n", (int) ocp_dim_size);
#endif
	void *ocp_dim_mem = malloc(ocp_dim_size);

	struct d_ocp_qcqp_dim ocp_dim;
	d_ocp_qcqp_dim_create(N, &ocp_dim, ocp_dim_mem);
	
	for(ii=0; ii<=N; ii++)
		{
		d_ocp_qcqp_dim_set_nx(ii, nx[ii], &ocp_dim);
		d_ocp_qcqp_dim_set_nu(ii, nu[ii], &ocp_dim);
		d_ocp_qcqp_dim_set_nbx(ii, nbx[ii], &ocp_dim);
		d_ocp_qcqp_dim_set_nbu(ii, nbu[ii], &ocp_dim);
		d_ocp_qcqp_dim_set_ng(ii, ng[ii], &ocp_dim);
		d_ocp_qcqp_dim_set_nq(ii, nq[ii], &ocp_dim);
		d_ocp_qcqp_dim_set_ns(ii, ns[ii], &ocp_dim);
		d_ocp_qcqp_dim_set_nsbx(ii, nsbx[ii], &ocp_dim);
		d_ocp_qcqp_dim_set_nsbu(ii, nsbu[ii], &ocp_dim);
		d_ocp_qcqp_dim_set_nsg(ii, nsg[ii], &ocp_dim);
		d_ocp_qcqp_dim_set_nsq(ii, nsq[ii], &ocp_dim);
		d_ocp_qcqp_dim_set_nbxe(ii, nbxe[ii], &ocp_dim); // state bounds to be removed
		}
	
#if PRINT
	d_ocp_qcqp_dim_print(&ocp_dim);
#endif

/************************************************
* ocp qp dim red eq dof (reduce equation dof, i.e. x0 elimination)
************************************************/

	hpipm_size_t ocp_dim_size2 = d_ocp_qcqp_dim_memsize(N);
#if PRINT
	printf("\ndim size red = %d\n", (int) ocp_dim_size2);
#endif
	void *ocp_dim_mem2 = malloc(ocp_dim_size2);

	struct d_ocp_qcqp_dim ocp_dim2;
	d_ocp_qcqp_dim_create(N, &ocp_dim2, ocp_dim_mem2);

#if REMOVE_X0
	d_ocp_qcqp_dim_reduce_eq_dof(&ocp_dim, &ocp_dim2);
#else
	d_ocp_qcqp_dim_copy_all(&ocp_dim, &ocp_dim2);
#endif

#if PRINT
	d_ocp_qcqp_dim_print(&ocp_dim2);
#endif

/************************************************
* dense qp dim
************************************************/

	hpipm_size_t dense_dim_size = d_dense_qcqp_dim_memsize();
#if PRINT
	printf("\ndense dim size = %d\n", (int) dense_dim_size);
#endif
	void *dense_dim_mem = malloc(dense_dim_size);

	struct d_dense_qcqp_dim dense_dim;
	d_dense_qcqp_dim_create(&dense_dim, dense_dim_mem);

	d_cond_qcqp_compute_dim(&ocp_dim2, &dense_dim);

#if PRINT
	d_dense_qcqp_dim_print(&dense_dim);
#endif

/************************************************
* ocp qp
************************************************/

	hpipm_size_t ocp_qp_size = d_ocp_qcqp_memsize(&ocp_dim);
#if PRINT
	printf("\nqp size = %d\n", (int) ocp_qp_size);
#endif
	void *ocp_qp_mem = malloc(ocp_qp_size);

	struct d_ocp_qcqp ocp_qp;
	d_ocp_qcqp_create(&ocp_dim, &ocp_qp, ocp_qp_mem);

	// dynamics
	for(ii=0; ii<N; ii++)
		{
		d_ocp_qcqp_set_A(ii, A, &ocp_qp);
		d_ocp_qcqp_set_B(ii, B, &ocp_qp);
		d_ocp_qcqp_set_b(ii, b, &ocp_qp);
		}
	
	// cost
	for(ii=0; ii<N; ii++)
		{
		d_ocp_qcqp_set_Q(ii, Q, &ocp_qp);
		d_ocp_qcqp_set_S(ii, S, &ocp_qp);
		d_ocp_qcqp_set_R(ii, R, &ocp_qp);
		d_ocp_qcqp_set_q(ii, q, &ocp_qp);
		d_ocp_qcqp_set_r(ii, r, &ocp_qp);
		}
	ii = N;
	d_ocp_qcqp_set_Q(ii, Q, &ocp_qp);
	d_ocp_qcqp_set_q(ii, q, &ocp_qp);
	
	// constraints
	ii = 0;
	d_ocp_qcqp_set_idxbx(ii, idxbx0, &ocp_qp); // initial state
	d_ocp_qcqp_set_lbx(ii, lbx0, &ocp_qp); // initial state
	d_ocp_qcqp_set_ubx(ii, ubx0, &ocp_qp); // initial state
	d_ocp_qcqp_set_idxbxe(ii, idxbxe0, &ocp_qp); // initial state
	d_ocp_qcqp_set_idxbu(ii, idxbu1, &ocp_qp);
	d_ocp_qcqp_set_lbu(ii, lbu1, &ocp_qp);
	d_ocp_qcqp_set_ubu(ii, ubu1, &ocp_qp);
	d_ocp_qcqp_set_C(ii, C1, &ocp_qp);
	d_ocp_qcqp_set_D(ii, D1, &ocp_qp);
	d_ocp_qcqp_set_lg(ii, lg1, &ocp_qp);
	d_ocp_qcqp_set_ug(ii, ug1, &ocp_qp);
	d_ocp_qcqp_set_Rq(ii, Rq1, &ocp_qp);
	d_ocp_qcqp_set_uq(ii, uq1, &ocp_qp);
	for(ii=1; ii<N; ii++)
		{
		d_ocp_qcqp_set_idxbx(ii, idxbx1, &ocp_qp);
		d_ocp_qcqp_set_lbx(ii, lbx1, &ocp_qp);
		d_ocp_qcqp_set_ubx(ii, ubx1, &ocp_qp);
		d_ocp_qcqp_set_idxbu(ii, idxbu1, &ocp_qp);
		d_ocp_qcqp_set_lbu(ii, lbu1, &ocp_qp);
		d_ocp_qcqp_set_ubu(ii, ubu1, &ocp_qp);
		d_ocp_qcqp_set_C(ii, C1, &ocp_qp);
		d_ocp_qcqp_set_D(ii, D1, &ocp_qp);
		d_ocp_qcqp_set_lg(ii, lg1, &ocp_qp);
		d_ocp_qcqp_set_ug(ii, ug1, &ocp_qp);
		d_ocp_qcqp_set_Rq(ii, Rq1, &ocp_qp);
		d_ocp_qcqp_set_Qq(ii, Qq1, &ocp_qp);
		d_ocp_qcqp_set_qq(ii, qq1, &ocp_qp);
		d_ocp_qcqp_set_uq(ii, uq1, &ocp_qp);
//		d_ocp_qcqp_set_uq_mask(ii, uq1_mask, &ocp_qp);
		}
//if(2<N)
//	d_ocp_qcqp_set_Rq(2, Rq2, &ocp_qp);
//if(3<N)
//	d_ocp_qcqp_set_Rq(3, Rq3, &ocp_qp);
//d_ocp_qcqp_set_Qq(N-1, QqNm1, &ocp_qp);
//d_ocp_qcqp_set_uq(N-1, uqNm1, &ocp_qp);
	ii = N;
	d_ocp_qcqp_set_idxbx(ii, idxbxN, &ocp_qp);
	d_ocp_qcqp_set_lbx(ii, lbxN, &ocp_qp);
	d_ocp_qcqp_set_ubx(ii, ubxN, &ocp_qp);
	d_ocp_qcqp_set_C(ii, CN, &ocp_qp);
	d_ocp_qcqp_set_lg(ii, lgN, &ocp_qp);
	d_ocp_qcqp_set_ug(ii, ugN, &ocp_qp);
	d_ocp_qcqp_set_Qq(ii, QqN, &ocp_qp);
	d_ocp_qcqp_set_qq(ii, qqN, &ocp_qp);
	d_ocp_qcqp_set_uq(ii, uqN, &ocp_qp);
	d_ocp_qcqp_set_uq_mask(ii, uqN_mask, &ocp_qp);
	// soft constraints
	d_ocp_qcqp_set_Zl(ii, ZlN, &ocp_qp);
	d_ocp_qcqp_set_Zu(ii, ZuN, &ocp_qp);
	d_ocp_qcqp_set_zl(ii, zlN, &ocp_qp);
	d_ocp_qcqp_set_zu(ii, zuN, &ocp_qp);
	d_ocp_qcqp_set_lls(ii, llsN, &ocp_qp);
	d_ocp_qcqp_set_lus(ii, lusN, &ocp_qp);
	d_ocp_qcqp_set_idxs(ii, idxsN, &ocp_qp);

	// dynamic constraints removal
	double *lbu_mask; d_zeros(&lbu_mask, nbu[0], 1);
	double *ubu_mask; d_zeros(&ubu_mask, nbu[0], 1);
	double *lbx_mask; d_zeros(&lbx_mask, nbx[0], 1);
	double *ubx_mask; d_zeros(&ubx_mask, nbx[0], 1);
//	d_ocp_qp_set("lbu_mask", 0, lbu_mask, &qp);
//	d_ocp_qp_set("ubu_mask", 0, ubu_mask, &qp);
//	d_ocp_qp_set("lbx_mask", N, lbx_mask, &qp);
//	d_ocp_qp_set("ubx_mask", N, ubx_mask, &qp);

#if PRINT
	d_ocp_qcqp_print(&ocp_dim, &ocp_qp);
#endif

/************************************************
* ocp qp red eq dof
************************************************/

	hpipm_size_t ocp_qp_size2 = d_ocp_qcqp_memsize(&ocp_dim2);
#if PRINT
	printf("\nqp size red = %d\n", (int) ocp_qp_size2);
#endif
	void *ocp_qp_mem2 = malloc(ocp_qp_size2);

	struct d_ocp_qcqp ocp_qp2;
	d_ocp_qcqp_create(&ocp_dim2, &ocp_qp2, ocp_qp_mem2);

/************************************************
* dense qp
************************************************/

	hpipm_size_t dense_qp_size = d_dense_qcqp_memsize(&dense_dim);
#if PRINT
	printf("\nqp size = %d\n", (int) dense_qp_size);
#endif
	void *dense_qp_mem = malloc(dense_qp_size);

	struct d_dense_qcqp dense_qp;
	d_dense_qcqp_create(&dense_dim, &dense_qp, dense_qp_mem);

/************************************************
* ocp qp sol
************************************************/

	hpipm_size_t ocp_qp_sol_size = d_ocp_qcqp_sol_memsize(&ocp_dim);
#if PRINT
	printf("\nocp qp sol size = %d\n", (int) ocp_qp_sol_size);
#endif
	void *ocp_qp_sol_mem = malloc(ocp_qp_sol_size);

	struct d_ocp_qcqp_sol ocp_qp_sol;
	d_ocp_qcqp_sol_create(&ocp_dim, &ocp_qp_sol, ocp_qp_sol_mem);

/************************************************
* ocp qp sol red eq dof
************************************************/

	hpipm_size_t ocp_qp_sol_size2 = d_ocp_qcqp_sol_memsize(&ocp_dim2);
#if PRINT
	printf("\nocp qp sol size = %d\n", (int) ocp_qp_sol_size2);
#endif
	void *ocp_qp_sol_mem2 = malloc(ocp_qp_sol_size2);

	struct d_ocp_qcqp_sol ocp_qp_sol2;
	d_ocp_qcqp_sol_create(&ocp_dim2, &ocp_qp_sol2, ocp_qp_sol_mem2);

/************************************************
* dense qp sol
************************************************/

	hpipm_size_t dense_qp_sol_size = d_dense_qcqp_sol_memsize(&dense_dim);
#if PRINT
	printf("\ndense qp sol size = %d\n", (int) dense_qp_sol_size);
#endif
	void *dense_qp_sol_mem = malloc(dense_qp_sol_size);

	struct d_dense_qcqp_sol dense_qp_sol;
	d_dense_qcqp_sol_create(&dense_dim, &dense_qp_sol, dense_qp_sol_mem);

/************************************************
* red eq dof arg
************************************************/

	hpipm_size_t ocp_qp_red_arg_size = d_ocp_qcqp_reduce_eq_dof_arg_memsize();
	void *ocp_qp_red_arg_mem = malloc(ocp_qp_red_arg_size);

	struct d_ocp_qcqp_reduce_eq_dof_arg ocp_qp_red_arg;
	d_ocp_qcqp_reduce_eq_dof_arg_create(&ocp_qp_red_arg, ocp_qp_red_arg_mem);

	int tmp_i1 = 1;
	d_ocp_qcqp_reduce_eq_dof_arg_set_default(&ocp_qp_red_arg);
	d_ocp_qcqp_reduce_eq_dof_arg_set_alias_unchanged(&tmp_i1, &ocp_qp_red_arg);
	d_ocp_qcqp_reduce_eq_dof_arg_set_comp_dual_sol_eq(&tmp_i1, &ocp_qp_red_arg);
	d_ocp_qcqp_reduce_eq_dof_arg_set_comp_dual_sol_ineq(&tmp_i1, &ocp_qp_red_arg);

/************************************************
* cond arguments
************************************************/

	hpipm_size_t cond_arg_size = d_cond_qcqp_arg_memsize();
#if PRINT
	printf("\ncond_arg size = %d\n", (int) cond_arg_size);
#endif
	void *cond_arg_mem = malloc(cond_arg_size);

	struct d_cond_qcqp_arg cond_arg;
	d_cond_qcqp_arg_create(&cond_arg, cond_arg_mem);
	d_cond_qcqp_arg_set_default(&cond_arg);

/************************************************
* ipm arg
************************************************/

	hpipm_size_t ipm_arg_size = d_dense_qcqp_ipm_arg_memsize(&dense_dim);
#if PRINT
	printf("\nipm arg size = %d\n", (int) ipm_arg_size);
#endif
	void *ipm_arg_mem = malloc(ipm_arg_size);

	struct d_dense_qcqp_ipm_arg ipm_arg;
	d_dense_qcqp_ipm_arg_create(&dense_dim, &ipm_arg, ipm_arg_mem);
//	enum hpipm_mode mode = SPEED_ABS;
	enum hpipm_mode mode = SPEED;
//	enum hpipm_mode mode = BALANCE;
//	enum hpipm_mode mode = ROBUST;
	d_dense_qcqp_ipm_arg_set_default(mode, &ipm_arg);

	double mu0 = 1e2;
	int split_step = 0;

	d_dense_qcqp_ipm_arg_set_mu0(&mu0, &ipm_arg);
	d_dense_qcqp_ipm_arg_set_split_step(&split_step, &ipm_arg);

/************************************************
* red eq dof workspace
************************************************/

	hpipm_size_t ocp_qp_red_work_size = d_ocp_qcqp_reduce_eq_dof_ws_memsize(&ocp_dim);
	void *ocp_qp_red_work_mem = malloc(ocp_qp_red_work_size);

	struct d_ocp_qcqp_reduce_eq_dof_ws ocp_qp_red_work;
	d_ocp_qcqp_reduce_eq_dof_ws_create(&ocp_dim, &ocp_qp_red_work, ocp_qp_red_work_mem);

/************************************************
* cond workspace
************************************************/

	hpipm_size_t cond_size = d_cond_qcqp_ws_memsize(&ocp_dim2, &cond_arg);
#if PRINT
	printf("\ncond size = %d\n", (int) cond_size);
#endif
	void *cond_mem = malloc(cond_size);

	struct d_cond_qcqp_ws cond_ws;
	d_cond_qcqp_ws_create(&ocp_dim2, &cond_arg, &cond_ws, cond_mem);

/************************************************
* ipm workspace
************************************************/

	hpipm_size_t dense_ipm_size = d_dense_qcqp_ipm_ws_memsize(&dense_dim, &ipm_arg);
#if PRINT
	printf("\ndense ipm size = %d\n", (int) dense_ipm_size);
#endif
	void *ipm_work_mem = malloc(dense_ipm_size);

	struct d_dense_qcqp_ipm_ws ipm_work;
	d_dense_qcqp_ipm_ws_create(&dense_dim, &ipm_arg, &ipm_work, ipm_work_mem);

/************************************************
* reduce equation dof (i.e. x0 elimination)
************************************************/
	
	/* red all */

#if REMOVE_X0
	hpipm_tic(&timer);

	for(rep=0; rep<nrep; rep++)
		{
		d_ocp_qcqp_reduce_eq_dof(&ocp_qp, &ocp_qp2, &ocp_qp_red_arg, &ocp_qp_red_work);
		}

	double time_red_eq_dof_all = hpipm_toc(&timer) / nrep;
#else
	d_ocp_qcqp_copy_all(&ocp_qp, &ocp_qp2);

	double time_red_eq_dof_all = 0.0;
#endif

	/* red rhs */

#if REMOVE_X0
	hpipm_tic(&timer);

	for(rep=0; rep<nrep; rep++)
		{
		d_ocp_qcqp_reduce_eq_dof_rhs(&ocp_qp, &ocp_qp2, &ocp_qp_red_arg, &ocp_qp_red_work);
		}

	double time_red_eq_dof_rhs = hpipm_toc(&timer) / nrep;
#else
	d_ocp_qcqp_copy_all(&ocp_qp, &ocp_qp2);

	double time_red_eq_dof_rhs = 0.0;
#endif

	/* red lhs */

#if REMOVE_X0
	hpipm_tic(&timer);

	for(rep=0; rep<nrep; rep++)
		{
		d_ocp_qcqp_reduce_eq_dof_lhs(&ocp_qp, &ocp_qp2, &ocp_qp_red_arg, &ocp_qp_red_work);
		}

	double time_red_eq_dof_lhs = hpipm_toc(&timer) / nrep;
#else
	d_ocp_qcqp_copy_all(&ocp_qp, &ocp_qp2);

	double time_red_eq_dof_lhs = 0.0;
#endif

#if PRINT
	d_ocp_qcqp_print(&ocp_dim2, &ocp_qp2);
#endif

/************************************************
* cond
************************************************/

	/* cond all */

	hpipm_tic(&timer);

	for(rep=0; rep<nrep; rep++)
		{
		d_cond_qcqp_cond(&ocp_qp2, &dense_qp, &cond_arg, &cond_ws);
		}

	double time_cond_all = hpipm_toc(&timer) / nrep;

#if PRINT
	d_dense_qcqp_print(&dense_dim, &dense_qp);
#endif


	/* cond lhs */

	hpipm_tic(&timer);

	for(rep=0; rep<nrep; rep++)
		{
		d_cond_qcqp_cond_lhs(&ocp_qp2, &dense_qp, &cond_arg, &cond_ws);
		}

	double time_cond_lhs = hpipm_toc(&timer) / nrep;

#if PRINT
	d_dense_qcqp_print(&dense_dim, &dense_qp);
#endif

	/* cond rhs */

	hpipm_tic(&timer);

	for(rep=0; rep<nrep; rep++)
		{
		d_cond_qcqp_cond_rhs(&ocp_qp2, &dense_qp, &cond_arg, &cond_ws);
		}

	double time_cond_rhs = hpipm_toc(&timer) / nrep;

#if PRINT
	d_dense_qcqp_print(&dense_dim, &dense_qp);
#endif

/************************************************
* ipm
************************************************/

	hpipm_tic(&timer);

	for(rep=0; rep<nrep; rep++)
		{
		d_dense_qcqp_ipm_solve(&dense_qp, &dense_qp_sol, &ipm_arg, &ipm_work);
		d_dense_qcqp_ipm_get_status(&ipm_work, &hpipm_status);
		}

	double time_dense_ipm = hpipm_toc(&timer) / nrep;

#if PRINT
	d_dense_qcqp_sol_print(&dense_dim, &dense_qp_sol);
#endif

/************************************************
* expand solution
************************************************/

	hpipm_tic(&timer);

	for(rep=0; rep<nrep; rep++)
		{
		d_cond_qcqp_expand_sol(&ocp_qp2, &dense_qp_sol, &ocp_qp_sol2, &cond_arg, &cond_ws);
		}

	double time_expa = hpipm_toc(&timer) / nrep;

#if PRINT
	d_ocp_qcqp_sol_print(&ocp_dim2, &ocp_qp_sol2);
#endif

/************************************************
* restore equation dof
************************************************/

#if REMOVE_X0
	hpipm_tic(&timer);

	for(rep=0; rep<nrep; rep++)
		{
		d_ocp_qcqp_restore_eq_dof(&ocp_qp, &ocp_qp_sol2, &ocp_qp_sol, &ocp_qp_red_arg, &ocp_qp_red_work);
		}

	double time_res_eq_dof = hpipm_toc(&timer) / nrep;
#else
	d_ocp_qcqp_sol_copy_all(&ocp_qp_sol2, &ocp_qp_sol);

	double time_res_eq_dof = 0.0;
#endif

#if PRINT
	d_ocp_qcqp_sol_print(&ocp_dim, &ocp_qp_sol);
#endif

/************************************************
* print ipm statistics
************************************************/

	int iter; d_dense_qcqp_ipm_get_iter(&ipm_work, &iter);

	double max_res_stat; d_dense_qcqp_ipm_get_max_res_stat(&ipm_work, &max_res_stat);
	double max_res_eq  ; d_dense_qcqp_ipm_get_max_res_eq(&ipm_work, &max_res_eq);
	double max_res_ineq; d_dense_qcqp_ipm_get_max_res_ineq(&ipm_work, &max_res_ineq);
	double max_res_comp; d_dense_qcqp_ipm_get_max_res_comp(&ipm_work, &max_res_comp);
	double *stat; d_dense_qcqp_ipm_get_stat(&ipm_work, &stat);
	int stat_m;  d_dense_qcqp_ipm_get_stat_m(&ipm_work, &stat_m);

#if PRINT
	printf("\nipm return = %d\n", hpipm_status);
	printf("\nipm iter = %d\n", iter);
	printf("\nipm max res: stat = %e, eq =  %e, ineq =  %e, comp = %e\n", max_res_stat, max_res_eq, max_res_ineq, max_res_comp);

	printf("\nalpha_aff\tmu_aff\t\tsigma\t\talpha_prim\talpha_dual\tmu\t\tres_stat\tres_eq\t\tres_ineq\tres_comp\tobj\t\tlq fact\t\titref pred\titref corr\tlin res stat\tlin res eq\tlin res ineq\tlin res comp\n");
	d_print_exp_tran_mat(stat_m, iter+1, stat, stat_m);

	printf("\nred eq dof all time = %e [s]\n", time_red_eq_dof_all);
	printf("\nred eq dof lhs time = %e [s]\n", time_red_eq_dof_lhs);
	printf("\nred eq dof rhs time = %e [s]\n", time_red_eq_dof_rhs);
	printf("\ncond all time       = %e [s]\n", time_cond_all);
	printf("\ncond lhs time       = %e [s]\n", time_cond_lhs);
	printf("\ncond rhs time       = %e [s]\n", time_cond_rhs);
	printf("\ndense ipm time      = %e [s]\n", time_dense_ipm);
	printf("\nexpand time         = %e [s]\n", time_expa);
	printf("\nres eq dof time     = %e [s]\n\n", time_res_eq_dof);
	printf("\ntotal solution time = %e [s]\n\n", time_red_eq_dof_all+time_cond_all+time_dense_ipm+time_expa+time_res_eq_dof);
#endif

/************************************************
* free memory
************************************************/

	d_free(A);
	d_free(B);
	d_free(b);
	d_free(x0);
	d_free(Q);
//	d_free(QN);
	d_free(R);
	d_free(S);
	d_free(q);
//	d_free(qN);
	d_free(r);
	int_free(idxbx0);
	d_free(lbx0);
	d_free(ubx0);
	int_free(idxbxe0);
	int_free(idxbx1);
	d_free(lbx1);
	d_free(ubx1);
	int_free(idxbu1);
	d_free(lbu1);
	d_free(ubu1);
	int_free(idxbxN);
	d_free(lbxN);
	d_free(ubxN);
	d_free(C1);
	d_free(D1);
	d_free(lg1);
	d_free(ug1);
	d_free(CN);
	d_free(DN);
	d_free(lgN);
	d_free(ugN);
	d_free(Rq1);
	d_free(Rq2);
	d_free(Rq3);
	d_free(Qq1);
	d_free(QqN);
	d_free(qqN);
	d_free(uq1);
	d_free(uqN);
	d_free(uq1_mask);
	d_free(uqN_mask);
	d_free(lbu_mask);
	d_free(ubu_mask);
	d_free(lbx_mask);
	d_free(ubx_mask);

#if 0
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
#endif
	d_free(ZlN);
	d_free(ZuN);
	d_free(zlN);
	d_free(zuN);
	int_free(idxsN);
	d_free(llsN);
	d_free(lusN);

	free(ocp_dim_mem);
	free(ocp_qp_mem);
	free(ocp_qp_sol_mem);
	free(dense_dim_mem);
	free(dense_qp_mem);
	free(dense_qp_sol_mem);
	free(cond_mem);
	free(cond_arg_mem);
	free(ipm_arg_mem);
	free(ipm_work_mem);

	free(ocp_dim_mem2);
	free(ocp_qp_mem2);
	free(ocp_qp_sol_mem2);
	free(ocp_qp_red_arg_mem);
	free(ocp_qp_red_work_mem);

/************************************************
* return
************************************************/

	return hpipm_status;

	}
