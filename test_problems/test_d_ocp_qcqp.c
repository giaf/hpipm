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

#include <hpipm_d_ocp_qcqp_dim.h>
#include <hpipm_d_ocp_qcqp.h>
#include <hpipm_d_ocp_qcqp_sol.h>
#include <hpipm_d_ocp_qcqp_ipm.h>
#include <hpipm_d_ocp_qcqp_utils.h>
#include <hpipm_d_ocp_qcqp_red.h>
#include <hpipm_timing.h>

#include "d_tools.h"



// printing
#ifndef PRINT
#define PRINT 0
#endif


// remove initial state x0 from optimization variables
#define REMOVE_X0 1



#if ! defined(EXT_DEP)
/* creates a zero matrix */
void d_zeros(double **pA, int row, int col)
	{
	*pA = malloc((row*col)*sizeof(double));
	double *A = *pA;
	int i;
	for(i=0; i<row*col; i++) A[i] = 0.0;
	}
/* frees matrix */
void d_free(double *pA)
	{
	free( pA );
	}
/* prints a matrix in column-major format */
void d_print_mat(int m, int n, double *A, int lda)
	{
	int i, j;
	for(i=0; i<m; i++)
		{
		for(j=0; j<n; j++)
			{
			printf("%9.5f ", A[i+lda*j]);
			}
		printf("\n");
		}
	printf("\n");
	}
/* prints the transposed of a matrix in column-major format */
void d_print_tran_mat(int row, int col, double *A, int lda)
	{
	int i, j;
	for(j=0; j<col; j++)
		{
		for(i=0; i<row; i++)
			{
			printf("%9.5f ", A[i+lda*j]);
			}
		printf("\n");
		}
	printf("\n");
	}
/* prints a matrix in column-major format (exponential notation) */
void d_print_exp_mat(int m, int n, double *A, int lda)
	{
	int i, j;
	for(i=0; i<m; i++)
		{
		for(j=0; j<n; j++)
			{
			printf("%e\t", A[i+lda*j]);
			}
		printf("\n");
		}
	printf("\n");
	}
/* prints the transposed of a matrix in column-major format (exponential notation) */
void d_print_exp_tran_mat(int row, int col, double *A, int lda)
	{
	int i, j;
	for(j=0; j<col; j++)
		{
		for(i=0; i<row; i++)
			{
			printf("%e\t", A[i+lda*j]);
			}
		printf("\n");
		}
	printf("\n");
	}
/* creates a zero matrix aligned */
void int_zeros(int **pA, int row, int col)
	{
	void *temp = malloc((row*col)*sizeof(int));
	*pA = temp;
	int *A = *pA;
	int i;
	for(i=0; i<row*col; i++) A[i] = 0;
	}
/* frees matrix */
void int_free(int *pA)
	{
	free( pA );
	}
/* prints a matrix in column-major format */
void int_print_mat(int row, int col, int *A, int lda)
	{
	int i, j;
	for(i=0; i<row; i++)
		{
		for(j=0; j<col; j++)
			{
			printf("%d ", A[i+lda*j]);
			}
		printf("\n");
		}
	printf("\n");
	}
#endif



/************************************************
Mass-spring system: nx/2 masses connected each other with springs (in a row), and the first and the last one to walls. nu (<=nx) controls act on the first nu masses. The system is sampled with sampling time Ts.
************************************************/
void mass_spring_system(double Ts, int nx, int nu, double *A, double *B, double *b, double *x0)
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
		nbu[ii] = 0;
//		nbu[ii] = nu[ii];

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
		nq[ii] = 1;
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
	mass_spring_system(Ts, nx_, nu_, A, B, b, x0);

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

#if 0
	double *QN; d_zeros(&QN, nx_, nx_);
	for(ii=0; ii<2; ii++) QN[ii*(nx_+1)] = 1e15;
	for(ii=0; ii<nx_; ii++) QN[ii*(nx_+1)] += Q[ii*(nx_+1)];
	double *qN; d_zeros(&qN, nx_, 1);
	qN[0] = - 0.1;
	qN[1] = - 0.1;
	for(ii=0; ii<2; ii++) qN[ii] *= 1e15;
	for(ii=0; ii<nx_; ii++) qN[ii] += q[ii];
#endif

#if PRINT
	d_print_mat(nx_, nx_, Q, nx_);
	d_print_mat(nu_, nu_, R, nu_);
	d_print_mat(nu_, nx_, S, nu_);
	d_print_mat(1, nx_, q, 1);
	d_print_mat(1, nu_, r, 1);
//	d_print_mat(nx_, nx_, QN, nx_);
//	d_print_mat(1, nx_, qN, 1);
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
	int_print_mat(1, nx[0], idxbx0, 1);
	d_print_mat(1, nx[0], lbx0, 1);
	d_print_mat(1, nx[0], ubx0, 1);
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
//	if(nq[1]>0)
//		uq1_mask[0] = 1.0;
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
* ocp qcqp dim
************************************************/

	hpipm_size_t dim_size = d_ocp_qcqp_dim_memsize(N);
#if PRINT
	printf("\ndim size = %d\n", (int) dim_size);
#endif
	void *dim_mem = malloc(dim_size);

	struct d_ocp_qcqp_dim dim;
	d_ocp_qcqp_dim_create(N, &dim, dim_mem);
	
	for(ii=0; ii<=N; ii++)
		{
		d_ocp_qcqp_dim_set_nx(ii, nx[ii], &dim);
		d_ocp_qcqp_dim_set_nu(ii, nu[ii], &dim);
		d_ocp_qcqp_dim_set_nbx(ii, nbx[ii], &dim);
		d_ocp_qcqp_dim_set_nbu(ii, nbu[ii], &dim);
		d_ocp_qcqp_dim_set_ng(ii, ng[ii], &dim);
		d_ocp_qcqp_dim_set_nq(ii, nq[ii], &dim);
		d_ocp_qcqp_dim_set_nsbx(ii, nsbx[ii], &dim);
		d_ocp_qcqp_dim_set_nsq(ii, nsq[ii], &dim);
		d_ocp_qcqp_dim_set_nbxe(ii, nbxe[ii], &dim); // state bounds to be removed
		}
	
#if PRINT
	d_ocp_qcqp_dim_print(&dim);
#endif

//	d_ocp_qcqp_dim_codegen("mass_spring_qcqp_data.c", "w", &dim);

/************************************************
* ocp qp dim red eq dof (reduce equation dof, i.e. x0 elimination)
************************************************/

	hpipm_size_t dim_size2 = d_ocp_qcqp_dim_memsize(N);
#if PRINT
	printf("\ndim size red = %d\n", (int) dim_size2);
#endif
	void *dim_mem2 = malloc(dim_size2);

	struct d_ocp_qcqp_dim dim2;
	d_ocp_qcqp_dim_create(N, &dim2, dim_mem2);

#if REMOVE_X0
	d_ocp_qcqp_dim_reduce_eq_dof(&dim, &dim2);
#else
	d_ocp_qcqp_dim_copy_all(&dim, &dim2);
#endif

#if PRINT
	d_ocp_qcqp_dim_print(&dim2);
#endif

/************************************************
* ocp qp
************************************************/

	hpipm_size_t qcqp_size = d_ocp_qcqp_memsize(&dim);
#if PRINT
	printf("\nqp size = %d\n", (int) qcqp_size);
#endif
	void *qcqp_mem = malloc(qcqp_size);

	struct d_ocp_qcqp qcqp;
	d_ocp_qcqp_create(&dim, &qcqp, qcqp_mem);

	// dynamics
	for(ii=0; ii<N; ii++)
		{
		d_ocp_qcqp_set_A(ii, A, &qcqp);
		d_ocp_qcqp_set_B(ii, B, &qcqp);
		d_ocp_qcqp_set_b(ii, b, &qcqp);
		}
	
	// cost
	for(ii=0; ii<N; ii++)
		{
		d_ocp_qcqp_set_Q(ii, Q, &qcqp);
		d_ocp_qcqp_set_S(ii, S, &qcqp);
		d_ocp_qcqp_set_R(ii, R, &qcqp);
		d_ocp_qcqp_set_q(ii, q, &qcqp);
		d_ocp_qcqp_set_r(ii, r, &qcqp);
		}
	ii = N;
	d_ocp_qcqp_set_Q(ii, Q, &qcqp);
	d_ocp_qcqp_set_q(ii, q, &qcqp);
	
	// constraints
	ii = 0;
	d_ocp_qcqp_set_idxbx(ii, idxbx0, &qcqp); // initial state
	d_ocp_qcqp_set_lbx(ii, lbx0, &qcqp); // initial state
	d_ocp_qcqp_set_ubx(ii, ubx0, &qcqp); // initial state
	d_ocp_qcqp_set_idxbxe(ii, idxbxe0, &qcqp); // initial state
	d_ocp_qcqp_set_idxbu(ii, idxbu1, &qcqp);
	d_ocp_qcqp_set_lbu(ii, lbu1, &qcqp);
	d_ocp_qcqp_set_ubu(ii, ubu1, &qcqp);
	d_ocp_qcqp_set_C(ii, C1, &qcqp);
	d_ocp_qcqp_set_D(ii, D1, &qcqp);
	d_ocp_qcqp_set_lg(ii, lg1, &qcqp);
	d_ocp_qcqp_set_ug(ii, ug1, &qcqp);
	d_ocp_qcqp_set_Rq(ii, Rq1, &qcqp);
	d_ocp_qcqp_set_uq(ii, uq1, &qcqp);
//	d_ocp_qcqp_set_uq_mask(ii, uq1_mask, &qcqp);
	for(ii=1; ii<N; ii++)
		{
		d_ocp_qcqp_set_idxbx(ii, idxbx1, &qcqp);
		d_ocp_qcqp_set_lbx(ii, lbx1, &qcqp);
		d_ocp_qcqp_set_ubx(ii, ubx1, &qcqp);
		d_ocp_qcqp_set_idxbu(ii, idxbu1, &qcqp);
		d_ocp_qcqp_set_lbu(ii, lbu1, &qcqp);
		d_ocp_qcqp_set_ubu(ii, ubu1, &qcqp);
		d_ocp_qcqp_set_C(ii, C1, &qcqp);
		d_ocp_qcqp_set_D(ii, D1, &qcqp);
		d_ocp_qcqp_set_lg(ii, lg1, &qcqp);
		d_ocp_qcqp_set_ug(ii, ug1, &qcqp);
		d_ocp_qcqp_set_Rq(ii, Rq1, &qcqp);
		d_ocp_qcqp_set_Qq(ii, Qq1, &qcqp);
		d_ocp_qcqp_set_qq(ii, qq1, &qcqp);
		d_ocp_qcqp_set_uq(ii, uq1, &qcqp);
//		d_ocp_qcqp_set_uq_mask(ii, uq1_mask, &qcqp);
		}
//if(2<N)
//	d_ocp_qcqp_set_Rq(2, Rq2, &qcqp);
//if(3<N)
//	d_ocp_qcqp_set_Rq(3, Rq3, &qcqp);
//d_ocp_qcqp_set_Qq(N-1, QqNm1, &qcqp);
//d_ocp_qcqp_set_uq(N-1, uqNm1, &qcqp);
	ii = N;
	d_ocp_qcqp_set_idxbx(ii, idxbxN, &qcqp);
	d_ocp_qcqp_set_lbx(ii, lbxN, &qcqp);
	d_ocp_qcqp_set_ubx(ii, ubxN, &qcqp);
	d_ocp_qcqp_set_C(ii, CN, &qcqp);
	d_ocp_qcqp_set_lg(ii, lgN, &qcqp);
	d_ocp_qcqp_set_ug(ii, ugN, &qcqp);
	d_ocp_qcqp_set_Qq(ii, QqN, &qcqp);
	d_ocp_qcqp_set_qq(ii, qqN, &qcqp);
	d_ocp_qcqp_set_uq(ii, uqN, &qcqp);
	d_ocp_qcqp_set_uq_mask(ii, uqN_mask, &qcqp);
	// soft constraints
	d_ocp_qcqp_set_Zl(ii, ZlN, &qcqp);
	d_ocp_qcqp_set_Zu(ii, ZuN, &qcqp);
	d_ocp_qcqp_set_zl(ii, zlN, &qcqp);
	d_ocp_qcqp_set_zu(ii, zuN, &qcqp);
	d_ocp_qcqp_set_lls(ii, llsN, &qcqp);
	d_ocp_qcqp_set_lus(ii, lusN, &qcqp);
	d_ocp_qcqp_set_idxs(ii, idxsN, &qcqp);

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
	d_ocp_qcqp_print(&dim, &qcqp);
#endif

//	d_ocp_qcqp_codegen("mass_spring_qcqp_data.c", "a", &dim, &qcqp);

/************************************************
* ocp qp red eq dof
************************************************/

	hpipm_size_t qcqp_size2 = d_ocp_qcqp_memsize(&dim2);
#if PRINT
	printf("\nqp size red = %d\n", (int) qcqp_size2);
#endif
	void *qcqp_mem2 = malloc(qcqp_size2);

	struct d_ocp_qcqp qcqp2;
	d_ocp_qcqp_create(&dim2, &qcqp2, qcqp_mem2);

/************************************************
* ocp qp sol
************************************************/

	hpipm_size_t qcqp_sol_size = d_ocp_qcqp_sol_memsize(&dim);
#if PRINT
	printf("\nqcqp sol size = %d\n", (int) qcqp_sol_size);
#endif
	void *qcqp_sol_mem = malloc(qcqp_sol_size);

	struct d_ocp_qcqp_sol qcqp_sol;
	d_ocp_qcqp_sol_create(&dim, &qcqp_sol, qcqp_sol_mem);

/************************************************
* ocp qp sol red eq dof
************************************************/

	hpipm_size_t qcqp_sol_size2 = d_ocp_qcqp_sol_memsize(&dim2);
#if PRINT
	printf("\nqcqp sol size red = %d\n", (int) qcqp_sol_size2);
#endif
	void *qcqp_sol_mem2 = malloc(qcqp_sol_size2);

	struct d_ocp_qcqp_sol qcqp_sol2;
	d_ocp_qcqp_sol_create(&dim2, &qcqp_sol2, qcqp_sol_mem2);

/************************************************
* red eq dof arg
************************************************/

	hpipm_size_t qp_red_arg_size = d_ocp_qcqp_reduce_eq_dof_arg_memsize();
	void *qp_red_arg_mem = malloc(qp_red_arg_size);

	struct d_ocp_qcqp_reduce_eq_dof_arg qp_red_arg;
	d_ocp_qcqp_reduce_eq_dof_arg_create(&qp_red_arg, qp_red_arg_mem);

	int tmp_i1 = 1;
	d_ocp_qcqp_reduce_eq_dof_arg_set_default(&qp_red_arg);
	d_ocp_qcqp_reduce_eq_dof_arg_set_alias_unchanged(&tmp_i1, &qp_red_arg);
	d_ocp_qcqp_reduce_eq_dof_arg_set_comp_dual_sol_eq(&tmp_i1, &qp_red_arg);
	d_ocp_qcqp_reduce_eq_dof_arg_set_comp_dual_sol_ineq(&tmp_i1, &qp_red_arg);

/************************************************
* ipm arg
************************************************/

	hpipm_size_t ipm_arg_size = d_ocp_qcqp_ipm_arg_memsize(&dim2);
	void *ipm_arg_mem = malloc(ipm_arg_size);

	struct d_ocp_qcqp_ipm_arg ipm_arg;
	d_ocp_qcqp_ipm_arg_create(&dim2, &ipm_arg, ipm_arg_mem);

//	enum hpipm_mode mode = SPEED_ABS;
	enum hpipm_mode mode = SPEED;
//	enum hpipm_mode mode = BALANCE;
//	enum hpipm_mode mode = ROBUST;
	d_ocp_qcqp_ipm_arg_set_default(mode, &ipm_arg);

	double mu0 = 1e2;
	int iter_max = 30;
	double alpha_min = 1e-8;
	double tol_stat = 1e-8;
	double tol_eq = 1e-8;
	double tol_ineq = 1e-8;
	double tol_comp = 1e-8;
	double reg_prim = 1e-12;
	int warm_start = 0;
	int pred_corr = 1;
	int ric_alg = 1;
	int comp_res_exit = 1;
	int split_step = 0;

	d_ocp_qcqp_ipm_arg_set_mu0(&mu0, &ipm_arg);
	d_ocp_qcqp_ipm_arg_set_iter_max(&iter_max, &ipm_arg);
	d_ocp_qcqp_ipm_arg_set_alpha_min(&alpha_min, &ipm_arg);
	d_ocp_qcqp_ipm_arg_set_tol_stat(&tol_stat, &ipm_arg);
	d_ocp_qcqp_ipm_arg_set_tol_eq(&tol_eq, &ipm_arg);
	d_ocp_qcqp_ipm_arg_set_tol_ineq(&tol_ineq, &ipm_arg);
	d_ocp_qcqp_ipm_arg_set_tol_comp(&tol_comp, &ipm_arg);
	d_ocp_qcqp_ipm_arg_set_reg_prim(&reg_prim, &ipm_arg);
	d_ocp_qcqp_ipm_arg_set_warm_start(&warm_start, &ipm_arg);
	d_ocp_qcqp_ipm_arg_set_pred_corr(&pred_corr, &ipm_arg);
	d_ocp_qcqp_ipm_arg_set_ric_alg(&ric_alg, &ipm_arg);
	d_ocp_qcqp_ipm_arg_set_comp_res_exit(&comp_res_exit, &ipm_arg);
	d_ocp_qcqp_ipm_arg_set_split_step(&split_step, &ipm_arg);

//	d_ocp_qcqp_ipm_arg_codegen("mass_spring_qcqp_data.c", "a", &dim, &ipm_arg);

/************************************************
* red eq dof workspace
************************************************/

	hpipm_size_t qp_red_work_size = d_ocp_qcqp_reduce_eq_dof_ws_memsize(&dim);
	void *qp_red_work_mem = malloc(qp_red_work_size);

	struct d_ocp_qcqp_reduce_eq_dof_ws qp_red_work;
	d_ocp_qcqp_reduce_eq_dof_ws_create(&dim, &qp_red_work, qp_red_work_mem);

/************************************************
* ipm workspace
************************************************/

	hpipm_size_t ipm_size = d_ocp_qcqp_ipm_ws_memsize(&dim2, &ipm_arg);
#if PRINT
	printf("\nipm size = %d\n", (int) ipm_size);
#endif
	void *ipm_mem = malloc(ipm_size);

	struct d_ocp_qcqp_ipm_ws ipm_work;
	d_ocp_qcqp_ipm_ws_create(&dim2, &ipm_arg, &ipm_work, ipm_mem);

/************************************************
* reduce equation dof (i.e. x0 elimination)
************************************************/

#if REMOVE_X0
	hpipm_tic(&timer);

	for(rep=0; rep<nrep; rep++)
		{
		d_ocp_qcqp_reduce_eq_dof(&qcqp, &qcqp2, &qp_red_arg, &qp_red_work);
		}

	double time_red_eq_dof = hpipm_toc(&timer) / nrep;
#else
	d_ocp_qcqp_copy_all(&qcqp, &qcqp2);

	double time_red_eq_dof = 0.0;
#endif

//	d_ocp_qcqp_print(&dim3, &qp3);

/************************************************
* ipm solution
************************************************/

	hpipm_tic(&timer);

	for(rep=0; rep<nrep; rep++)
		{
		d_ocp_qcqp_ipm_solve(&qcqp2, &qcqp_sol2, &ipm_arg, &ipm_work);
		d_ocp_qcqp_ipm_get_status(&ipm_work, &hpipm_status);
		}

	double time_ipm = hpipm_toc(&timer) / nrep;

/************************************************
* restore equation dof
************************************************/

#if REMOVE_X0
	hpipm_tic(&timer);

	for(rep=0; rep<nrep; rep++)
		{
		d_ocp_qcqp_restore_eq_dof(&qcqp, &qcqp_sol2, &qcqp_sol, &qp_red_arg, &qp_red_work);
		}

	double time_res_eq_dof = hpipm_toc(&timer) / nrep;
#else
	d_ocp_qcqp_sol_copy_all(&qcqp_sol2, &qcqp_sol);

	double time_res_eq_dof = 0.0;
#endif

//	d_ocp_qcqp_sol_print(&dim, &qp_sol);

/************************************************
* print solution
************************************************/
	
#if PRINT
	d_ocp_qcqp_sol_print(&dim, &qcqp_sol);
#endif

	double *tmp_nx; d_zeros(&tmp_nx, nx_, 1);
	d_ocp_qcqp_sol_get_x(N, &qcqp_sol, tmp_nx);
//	d_print_mat(1, nx_, tmp_nx, 1);

	// check quadratic constr
	double *tmp_nux; d_zeros(&tmp_nux, nu_+nx_, 1);

	double tmp_qc;

#if 1
	for(ii=N; ii<=N; ii++)
		{
		for(jj=0; jj<nq[ii]; jj++)
			{
			for(kk=0; kk<nu[ii]+nx[ii]; kk++)
				{
				tmp_nux[kk] = 0.0;
				for(ll=0; ll<nu[ii]+nx[ii]; ll++)
					{
					tmp_nux[kk] += QqN[kk+ll*nx[N]]*tmp_nx[ll];
					}
				tmp_nux[kk] *= 0.5;
				tmp_nux[kk] += qqN[kk];
				}
			tmp_qc = 0;
			for(kk=0; kk<nu[ii]+nx[ii]; kk++)
				{
				tmp_qc += tmp_nux[kk]*tmp_nx[kk];
				}
#if PRINT
			printf("\nquadr constr %f\n\n", tmp_qc);
#endif
			}
		}
#endif

	d_free(tmp_nx);
	d_free(tmp_nux);

/************************************************
* print residuals
************************************************/
	
#if PRINT
	d_ocp_qcqp_res_print(&dim, ipm_work.qcqp_res);
#endif

/************************************************
* print ipm statistics
************************************************/

	int iter; d_ocp_qcqp_ipm_get_iter(&ipm_work, &iter);
	double res_stat; d_ocp_qcqp_ipm_get_max_res_stat(&ipm_work, &res_stat);
	double res_eq; d_ocp_qcqp_ipm_get_max_res_eq(&ipm_work, &res_eq);
	double res_ineq; d_ocp_qcqp_ipm_get_max_res_ineq(&ipm_work, &res_ineq);
	double res_comp; d_ocp_qcqp_ipm_get_max_res_comp(&ipm_work, &res_comp);
	double *stat; d_ocp_qcqp_ipm_get_stat(&ipm_work, &stat);
	int stat_m; d_ocp_qcqp_ipm_get_stat_m(&ipm_work, &stat_m);

#if PRINT
	printf("\nipm return = %d\n", hpipm_status);
	printf("\nipm residuals max: res_g = %e, res_b = %e, res_d = %e, res_m = %e\n", res_stat, res_eq, res_ineq, res_comp);

	printf("\nipm iter = %d\n", iter);
	printf("\nalpha_aff\tmu_aff\t\tsigma\t\talpha_prim\talpha_dual\tmu\t\tres_stat\tres_eq\t\tres_ineq\tres_comp\tobj\t\tlq fact\t\titref pred\titref corr\tlin res stat\tlin res eq\tlin res ineq\tlin res comp\n");
	d_print_exp_tran_mat(stat_m, iter+1, stat, stat_m);

	printf("\nred eq for time     = %e [s]\n", time_red_eq_dof);
	printf("\nocp ipm time        = %e [s]\n", time_ipm);
	printf("\nres eq for time     = %e [s]\n\n", time_res_eq_dof);
	printf("\ntotal solution time = %e [s]\n\n", time_red_eq_dof+time_ipm+time_res_eq_dof);
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

	free(dim_mem);
	free(qcqp_mem);
	free(qcqp_sol_mem);
	free(ipm_arg_mem);
	free(ipm_mem);

	free(dim_mem2);
	free(qcqp_mem2);
	free(qcqp_sol_mem2);
	free(qp_red_arg_mem);
	free(qp_red_work_mem);

/************************************************
* return
************************************************/

	return hpipm_status;

	}

