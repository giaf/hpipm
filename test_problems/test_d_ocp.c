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

#include "../include/hpipm_d_ocp_qp.h"
#include "../include/hpipm_d_ipm_hard_ocp_qp.h"

#include "tools.h"



#define KEEP_X0 0

// printing
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
	


	// problem size

	int nx_ = 8; // number of states (it has to be even for the mass-spring system test problem)
	int nu_ = 3; // number of inputs (controllers) (it has to be at least 1 and at most nx/2 for the mass-spring system test problem)
	int N  = 5; // horizon lenght



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
	ng[0] = 0;
	for(ii=1; ii<N; ii++)
		ng[ii] = 0;
	ng[N] = 0;

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
		b[jj] = 0.1;
	
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
	for(ii=0; ii<nx_; ii++) Q[ii*(nx_+1)] = 1.0;

	double *R; d_zeros(&R, nu_, nu_);
	for(ii=0; ii<nu_; ii++) R[ii*(nu_+1)] = 2.0;

	double *S; d_zeros(&S, nu_, nx_);

	double *q; d_zeros(&q, nx_, 1);
	for(ii=0; ii<nx_; ii++) q[ii] = 0.1;

	double *r; d_zeros(&r, nu_, 1);
	for(ii=0; ii<nu_; ii++) r[ii] = 0.2;

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
	double mu0 = 2.0;

/************************************************
* box & general constraints
************************************************/	

	int *idxb0; int_zeros(&idxb0, nb[0], 1);
	double *lb0; d_zeros(&lb0, nb[0], 1);
	double *ub0; d_zeros(&ub0, nb[0], 1);
	double *lg0; d_zeros(&lg0, ng[0], 1);
	double *ug0; d_zeros(&ug0, ng[0], 1);
	for(ii=0; ii<nb[0]; ii++)
		{
		if(ii<nu[0]) // input
			{
			lb0[ii] = - 0.5; // umin
			ub0[ii] =   0.5; // umax
			}
		else // state
			{
			lb0[ii] = - 4.0; // xmin
			ub0[ii] =   4.0; // xmax
			}
		idxb0[ii] = ii;
		}
	for(ii=0; ii<ng[0]; ii++)
		{
		if(ii<nu[0]) // input
			{
			lg0[ii] = - 0.5; // umin
			ug0[ii] =   0.5; // umax
			}
		else // state
			{
			lg0[ii] = - 4.0; // xmin
			ug0[ii] =   4.0; // xmax
			}
		}

	int *idxb1; int_zeros(&idxb1, nb[1], 1);
	double *lb1; d_zeros(&lb1, nb[1], 1);
	double *ub1; d_zeros(&ub1, nb[1], 1);
	double *lg1; d_zeros(&lg1, ng[1], 1);
	double *ug1; d_zeros(&ug1, ng[1], 1);
	for(ii=0; ii<nb[1]; ii++)
		{
		if(ii<nu[1]) // input
			{
			lb1[ii] = - 0.5; // umin
			ub1[ii] =   0.5; // umax
			}
		else // state
			{
			lb1[ii] = - 4.0; // xmin
			ub1[ii] =   4.0; // xmax
			}
		idxb1[ii] = ii;
		}
	for(ii=0; ii<ng[1]; ii++)
		{
		if(ii<nu[1]) // input
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


	int *idxbN; int_zeros(&idxbN, nb[N], 1);
	double *lbN; d_zeros(&lbN, nb[N], 1);
	double *ubN; d_zeros(&ubN, nb[N], 1);
	double *lgN; d_zeros(&lgN, ng[N], 1);
	double *ugN; d_zeros(&ugN, ng[N], 1);
	for(ii=0; ii<nb[N]; ii++)
		{
		lbN[ii] = - 4.0; // xmin
		ubN[ii] =   4.0; // xmax
		idxbN[ii] = ii;
		}
	for(ii=0; ii<ng[N]; ii++)
		{
		lgN[ii] = - 4.0; // dmin
		ugN[ii] =   4.0; // dmax
		}

	double *DC0; d_zeros(&DC0, ng[0], nu[0]+nx[0]);
//	for(ii=0; ii<ng[0]; ii++)
//		DC0[ii*(ng[0]+1)] = 1.0;

	double *DC1; d_zeros(&DC1, ng[1], nu[1]+nx[1]);
//	for(ii=0; ii<ng[1]; ii++)
//		DC1[ii*(ng[1]+1)] = 1.0;

	double *DCN; d_zeros(&DCN, ng[N], nx[N]);
//	for(ii=0; ii<ng[N]; ii++)
//		DCN[ii*(ng[N]+1)] = 1.0;

	double *C;
	double *D;

#if PRINT
	int_print_mat(1, nb[0], idxb0, 1);
	d_print_mat(1, nb[0], lb0, 1);
	d_print_mat(1, nb[0], ub0, 1);
	int_print_mat(1, nb[1], idxb1, 1);
	d_print_mat(1, nb[1], lb1, 1);
	d_print_mat(1, nb[1], ub1, 1);
	int_print_mat(1, nb[N], idxbN, 1);
	d_print_mat(1, nb[N], lbN, 1);
	d_print_mat(1, nb[N], ubN, 1);
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
	double *hlb[N+1];
	double *hub[N+1];
	double *hlg[N+1];
	double *hug[N+1];
	double *hC[N+1];
	double *hD[N+1];
	int *hidxb[N+1];

	hA[0] = A;
	hB[0] = B;
	hb[0] = b0;
	hQ[0] = Q;
	hS[0] = S;
	hR[0] = R;
	hq[0] = q;
	hr[0] = r0;
	hidxb[0] = idxb0;
	hlb[0] = lb0;
	hub[0] = ub0;
	hlg[0] = lg0;
	hug[0] = ug0;
	hC[0] = C;
	hD[0] = D;
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
		hlb[ii] = lb1;
		hub[ii] = ub1;
		hlg[ii] = lg1;
		hug[ii] = ug1;
		hC[ii] = C;
		hD[ii] = D;
		}
	hQ[N] = Q;
	hS[N] = S;
	hR[N] = R;
	hq[N] = q;
	hr[N] = r;
	hidxb[N] = idxbN;
	hlb[N] = lbN;
	hub[N] = ubN;
	hlg[N] = lgN;
	hug[N] = ugN;
	hC[N] = C;
	hD[N] = D;
	
/************************************************
* ocp qp
************************************************/	
	
	int qp_size = d_memsize_ocp_qp(N, nx, nu, nb, ng);
	printf("\nqp size = %d\n", qp_size);
	void *qp_mem = malloc(qp_size);

	struct d_ocp_qp qp;
	d_create_ocp_qp(N, nx, nu, nb, ng, &qp, qp_mem);
	d_cvt_colmaj_to_ocp_qp(hA, hB, hb, hQ, hS, hR, hq, hr, hidxb, hlb, hub, hC, hD, hlg, hug, &qp);
#if 0
	printf("\nN = %d\n", qp.N);
	for(ii=0; ii<N; ii++)
		d_print_strmat(qp.nu[ii]+qp.nx[ii]+1, qp.nx[ii+1], qp.BAbt+ii, 0, 0);
	for(ii=0; ii<N; ii++)
		d_print_tran_strvec(qp.nx[ii+1], qp.b+ii, 0);
	for(ii=0; ii<=N; ii++)
		d_print_strmat(qp.nu[ii]+qp.nx[ii]+1, qp.nu[ii]+qp.nx[ii], qp.RSQrq+ii, 0, 0);
	for(ii=0; ii<=N; ii++)
		d_print_tran_strvec(qp.nu[ii]+qp.nx[ii], qp.rq+ii, 0);
	for(ii=0; ii<=N; ii++)
		int_print_mat(1, nb[ii], qp.idxb[ii], 1);
	for(ii=0; ii<=N; ii++)
		d_print_tran_strvec(qp.nb[ii], qp.d_lb+ii, 0);
	for(ii=0; ii<=N; ii++)
		d_print_tran_strvec(qp.nb[ii], qp.d_ub+ii, 0);
	for(ii=0; ii<=N; ii++)
		d_print_strmat(qp.nu[ii]+qp.nx[ii], qp.ng[ii], qp.DCt+ii, 0, 0);
	for(ii=0; ii<=N; ii++)
		d_print_tran_strvec(qp.ng[ii], qp.d_lg+ii, 0);
	for(ii=0; ii<=N; ii++)
		d_print_tran_strvec(qp.ng[ii], qp.d_ug+ii, 0);
	return;
#endif

/************************************************
* ipm
************************************************/	

	struct d_ipm_hard_ocp_qp_arg arg;
	arg.alpha_min = 1e-8;
	arg.mu_max = 1e-12;
	arg.iter_max = 20;
	arg.mu0 = 2.0;

	int ipm_size = d_memsize_ipm_hard_ocp_qp(&qp, &arg);
	printf("\nipm size = %d\n", ipm_size);
	void *ipm_mem = malloc(ipm_size);

	struct d_ipm_hard_ocp_qp_workspace workspace;
	d_create_ipm_hard_ocp_qp(&qp, &arg, &workspace, ipm_mem);

	int rep, nrep=1000;

	struct timeval tv0, tv1;

	gettimeofday(&tv0, NULL); // start

	for(rep=0; rep<nrep; rep++)
		{
		d_solve_ipm_hard_ocp_qp(&qp, &workspace);
		}

	gettimeofday(&tv1, NULL); // stop

	double time_ocp_ipm = (tv1.tv_sec-tv0.tv_sec)/(nrep+0.0)+(tv1.tv_usec-tv0.tv_usec)/(nrep*1e6);

#if 1
	printf("\nsolution\n\n");
	printf("\nux\n");
	for(ii=0; ii<=N; ii++)
		d_print_tran_strvec(nu[ii]+nx[ii], workspace.ux+ii, 0);
	printf("\npi\n");
	for(ii=0; ii<N; ii++)
		d_print_tran_strvec(nx[ii+1], workspace.pi+ii, 0);
	printf("\nlam_lb\n");
	for(ii=0; ii<=N; ii++)
		d_print_tran_strvec(nb[ii], workspace.lam_lb+ii, 0);
	printf("\nlam_ub\n");
	for(ii=0; ii<=N; ii++)
		d_print_tran_strvec(nb[ii], workspace.lam_ub+ii, 0);
	printf("\nlam_lg\n");
	for(ii=0; ii<=N; ii++)
		d_print_tran_strvec(ng[ii], workspace.lam_lg+ii, 0);
	printf("\nlam_ug\n");
	for(ii=0; ii<=N; ii++)
		d_print_tran_strvec(ng[ii], workspace.lam_ug+ii, 0);
	printf("\nt_lb\n");
	for(ii=0; ii<=N; ii++)
		d_print_tran_strvec(nb[ii], workspace.t_lb+ii, 0);
	printf("\nt_ub\n");
	for(ii=0; ii<=N; ii++)
		d_print_tran_strvec(nb[ii], workspace.t_ub+ii, 0);
	printf("\nt_lg\n");
	for(ii=0; ii<=N; ii++)
		d_print_tran_strvec(ng[ii], workspace.t_lg+ii, 0);
	printf("\nt_ug\n");
	for(ii=0; ii<=N; ii++)
		d_print_tran_strvec(ng[ii], workspace.t_ug+ii, 0);

	printf("\nresiduals\n\n");
	printf("\nres_g\n");
	for(ii=0; ii<=N; ii++)
		d_print_e_tran_strvec(nu[ii]+nx[ii], workspace.res_g+ii, 0);
	printf("\nres_b\n");
	for(ii=0; ii<N; ii++)
		d_print_e_tran_strvec(nx[ii+1], workspace.res_b+ii, 0);
	printf("\nres_m_lb\n");
	for(ii=0; ii<=N; ii++)
		d_print_e_tran_strvec(nb[ii], workspace.res_m_lb+ii, 0);
	printf("\nres_m_ub\n");
	for(ii=0; ii<=N; ii++)
		d_print_e_tran_strvec(nb[ii], workspace.res_m_ub+ii, 0);
	printf("\nres_m_lg\n");
	for(ii=0; ii<=N; ii++)
		d_print_e_tran_strvec(ng[ii], workspace.res_m_lg+ii, 0);
	printf("\nres_m_ug\n");
	for(ii=0; ii<=N; ii++)
		d_print_e_tran_strvec(ng[ii], workspace.res_m_ug+ii, 0);
	printf("\nres_d_lb\n");
	for(ii=0; ii<=N; ii++)
		d_print_e_tran_strvec(nb[ii], workspace.res_d_lb+ii, 0);
	printf("\nres_d_ub\n");
	for(ii=0; ii<=N; ii++)
		d_print_e_tran_strvec(nb[ii], workspace.res_d_ub+ii, 0);
	printf("\nres_d_lg\n");
	for(ii=0; ii<=N; ii++)
		d_print_e_tran_strvec(ng[ii], workspace.res_d_lg+ii, 0);
	printf("\nres_d_ug\n");
	for(ii=0; ii<=N; ii++)
		d_print_e_tran_strvec(ng[ii], workspace.res_d_ug+ii, 0);
	printf("\nres_mu\n");
	printf("\n%e\n\n", workspace.res_mu);
#endif

	printf("\nipm iter = %d\n", workspace.iter);
	printf("\nsigma\t\talpha_aff\tmu_aff\t\talpha\t\tmu\n");
	d_print_e_tran_mat(5, workspace.iter, workspace.stat, 5);

	printf("\nocp ipm time = %e [s]\n\n", time_ocp_ipm);

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
	d_free(lb0);
	d_free(ub0);
	int_free(idxb1);
	d_free(lb1);
	d_free(ub1);
	int_free(idxbN);
	d_free(lbN);
	d_free(ubN);
	d_free(DC0);
	d_free(lg0);
	d_free(ug0);
	d_free(DC1);
	d_free(lg1);
	d_free(ug1);
	d_free(DCN);
	d_free(lgN);
	d_free(ugN);

	free(qp_mem);
	free(ipm_mem);

/************************************************
* return
************************************************/	

	return 0;

	}
