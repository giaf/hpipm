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

#include "../include/hpipm_tree.h"
#include "../include/hpipm_scenario_tree.h"
#include "../include/hpipm_d_tree_ocp_qp.h"

#include "d_tools.h"



#define KEEP_X0 0

// printing
#define PRINT 0



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
void d_print_e_mat(int m, int n, double *A, int lda)
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
void d_print_e_tran_mat(int row, int col, double *A, int lda)
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

	int ii, jj;
	int stage;

	int nx_ = 8;
	int nu_ = 3;

	int md = 3;
	int Nr = 2;
	int Nh = 3;

	// stage-wise size
	int nx[Nh+1];
	int nu[Nh+1];
	int nb[Nh+1];
	int ng[Nh+1];

	nx[0] = 0;
	nu[0] = nu_;
	nb[0] = nu[0]+nx[0];
	ng[0] = 0;
	for(ii=1; ii<Nh; ii++)
		{
		nx[ii] = nx_;
		nu[ii] = nu_;
		nb[ii] = nu[ii]+nx[ii];
		ng[ii] = 0;
		}
	nx[Nh] = nx_;
	nu[Nh] = 0;
	nb[Nh] = nu[Nh]+nx[Nh];
	ng[Nh] = 0;

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
		b[jj] = 0.1;
	
	for(jj=0; jj<nx_; jj++)
		x0[jj] = 0;
	x0[0] = 2.5;
	x0[1] = 2.5;

	double *b0; d_zeros(&b0, nx_, 1);
	dgemv_n_3l(nx_, nx_, A, nx_, x0, b0);
	daxpy_3l(nx_, 1.0, b, b0);

#if 0
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

#if 0
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
			d_lb1[ii] = - 4.0; // xmin
			d_ub1[ii] =   4.0; // xmax
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


	int *idxbN; int_zeros(&idxbN, nb[Nh], 1);
	double *d_lbN; d_zeros(&d_lbN, nb[Nh], 1);
	double *d_ubN; d_zeros(&d_ubN, nb[Nh], 1);
	double *d_lgN; d_zeros(&d_lgN, ng[Nh], 1);
	double *d_ugN; d_zeros(&d_ugN, ng[Nh], 1);
	for(ii=0; ii<nb[Nh]; ii++)
		{
		d_lbN[ii] = - 4.0; // xmin
		d_ubN[ii] =   4.0; // xmax
		idxbN[ii] = ii;
		}
	for(ii=0; ii<ng[Nh]; ii++)
		{
		d_lgN[ii] =   0.1; // dmin
		d_ugN[ii] =   0.1; // dmax
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

	double *CN; d_zeros(&CN, ng[Nh], nx[Nh]);
	double *DN; d_zeros(&DN, ng[Nh], nu[Nh]);
	for(ii=0; ii<nu[Nh]-nb[Nh] & ii<ng[Nh]; ii++)
		DN[ii+(nb[Nh]+ii)*ng[Nh]] = 1.0;
	for(; ii<ng[Nh]; ii++)
		CN[ii+(nb[Nh]+ii-nu[Nh])*ng[Nh]] = 1.0;

#if 0
	// box constraints
	int_print_mat(1, nb[0], idxb0, 1);
	d_print_mat(1, nb[0], d_lb0, 1);
	d_print_mat(1, nb[0], d_ub0, 1);
	int_print_mat(1, nb[1], idxb1, 1);
	d_print_mat(1, nb[1], d_lb1, 1);
	d_print_mat(1, nb[1], d_ub1, 1);
	int_print_mat(1, nb[Nh], idxbN, 1);
	d_print_mat(1, nb[Nh], d_lbN, 1);
	d_print_mat(1, nb[Nh], d_ubN, 1);
	// general constraints
	d_print_mat(1, ng[0], d_lg0, 1);
	d_print_mat(1, ng[0], d_ug0, 1);
	d_print_mat(ng[0], nu[0], D0, ng[0]);
	d_print_mat(ng[0], nx[0], C0, ng[0]);
	d_print_mat(1, ng[1], d_lg1, 1);
	d_print_mat(1, ng[1], d_ug1, 1);
	d_print_mat(ng[1], nu[1], D1, ng[1]);
	d_print_mat(ng[1], nx[1], C1, ng[1]);
	d_print_mat(1, ng[Nh], d_lgN, 1);
	d_print_mat(1, ng[Nh], d_ugN, 1);
	d_print_mat(ng[Nh], nu[Nh], DN, ng[Nh]);
	d_print_mat(ng[Nh], nx[Nh], CN, ng[Nh]);
#endif

/************************************************
* create scenario tree
************************************************/	

	int tree_memory_size = memsize_sctree(md, Nr, Nh);
	printf("\ntree memsize = %d\n", tree_memory_size);
	void *tree_memory = malloc(tree_memory_size);

	struct sctree st;
	create_sctree(md, Nr, Nh, &st, tree_memory);

	int Nn = st.Nn;

#if 0
	int Nn = st.Nn;
	printf("\nscenario tree\n");
	for(ii=0; ii<Nn; ii++)
		{
		printf("\n");
		printf("idx = %d\n", (st.root+ii)->idx);
		printf("stage = %d\n", (st.root+ii)->stage);
		printf("real = %d\n", (st.root+ii)->real);
		printf("idxkid = %d\n", (st.root+ii)->idxkid);
		printf("dad = %d\n", (st.root+ii)->dad);
		printf("nkids = %d\n", (st.root+ii)->nkids);
		printf("kids =");
		for(jj=0; jj<(st.root+ii)->nkids; jj++)
			printf(" %d", (st.root+ii)->kids[jj]);
		printf("\n\n");
		}
#endif

/************************************************
* cast scenario tree into tree
************************************************/	

	struct tree tt;
	cast_sctree2tree(&st, &tt);

#if 0
	Nn = tt.Nn;
	printf("\ntree\n");
	for(ii=0; ii<Nn; ii++)
		{
		printf("\n");
		printf("idx = %d\n", (tt.root+ii)->idx);
		printf("stage = %d\n", (tt.root+ii)->stage);
		printf("real = %d\n", (tt.root+ii)->real);
		printf("idxkid = %d\n", (tt.root+ii)->idxkid);
		printf("dad = %d\n", (tt.root+ii)->dad);
		printf("nkids = %d\n", (tt.root+ii)->nkids);
		printf("kids =");
		for(jj=0; jj<(tt.root+ii)->nkids; jj++)
			printf(" %d", (tt.root+ii)->kids[jj]);
		printf("\n\n");
		}
#endif

/************************************************
* tree ocp problem size
************************************************/	

	// node-wise size
	int nxt[Nn];
	int nut[Nn];
	int nbt[Nn];
	int ngt[Nn];

	for(ii=0; ii<Nn; ii++)
		{
		stage = (tt.root+ii)->stage;
		nxt[ii] = nx[stage];
		nut[ii] = nu[stage];
		nbt[ii] = nb[stage];
		ngt[ii] = ng[stage];
		}
	
#if 0
	for(ii=0; ii<Nn; ii++)
		{
		printf("\n%d %d %d %d\n", nxt[ii], nut[ii], nbt[ii], ngt[ii]);
		}
#endif

/************************************************
* tree ocp data
************************************************/	

	// stage-wise data

	double *hA[Nh];
	double *hB[Nh];
	double *hb[Nh];
	double *hQ[Nh+1];
	double *hS[Nh+1];
	double *hR[Nh+1];
	double *hq[Nh+1];
	double *hr[Nh+1];
	double *hd_lb[Nh+1];
	double *hd_ub[Nh+1];
	double *hd_lg[Nh+1];
	double *hd_ug[Nh+1];
	double *hC[Nh+1];
	double *hD[Nh+1];
	int *hidxb[Nh+1];

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
	hd_lg[0] = d_lg0;
	hd_ug[0] = d_ug0;
	hC[0] = C0;
	hD[0] = D0;
	for(ii=1; ii<Nh; ii++)
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
		}
	hQ[Nh] = Q;
	hS[Nh] = S;
	hR[Nh] = R;
	hq[Nh] = q;
	hr[Nh] = r;
	hidxb[Nh] = idxbN;
	hd_lb[Nh] = d_lbN;
	hd_ub[Nh] = d_ubN;
	hd_lg[Nh] = d_lgN;
	hd_ug[Nh] = d_ugN;
	hC[Nh] = CN;
	hD[Nh] = DN;
	
	// node-wise data

	double *hAt[Nn-1];
	double *hBt[Nn-1];
	double *hbt[Nn-1];
	double *hQt[Nn];
	double *hSt[Nn];
	double *hRt[Nn];
	double *hqt[Nn];
	double *hrt[Nn];
	double *hd_lbt[Nn];
	double *hd_ubt[Nn];
	double *hd_lgt[Nn];
	double *hd_ugt[Nn];
	double *hCt[Nn];
	double *hDt[Nn];
	int *hidxbt[Nn];

	for(ii=0; ii<Nn-1; ii++)
		{
		stage = (tt.root+ii+1)->stage-1;
		hAt[ii] = hA[stage];
		hBt[ii] = hB[stage];
		hbt[ii] = hb[stage];
		}

	for(ii=0; ii<Nn; ii++)
		{
		stage = (tt.root+ii)->stage;
		hQt[ii] = hQ[stage];
		hRt[ii] = hR[stage];
		hSt[ii] = hS[stage];
		hqt[ii] = hq[stage];
		hrt[ii] = hr[stage];
		hd_lbt[ii] = hd_lb[stage];
		hd_ubt[ii] = hd_ub[stage];
		hd_lgt[ii] = hd_lg[stage];
		hd_ugt[ii] = hd_ug[stage];
		hidxbt[ii] = hidxb[stage];
		}

/************************************************
* create tree ocp qp
************************************************/	

	int tree_ocp_qp_memory_size = d_memsize_tree_ocp_qp(&tt, nxt, nut, nbt, ngt);
	printf("\ntree ocp qp memsize = %d\n", tree_ocp_qp_memory_size);
	void *tree_ocp_qp_memory = malloc(tree_ocp_qp_memory_size);

	struct d_tree_ocp_qp qp;
	d_create_tree_ocp_qp(&tt, nxt, nut, nbt, ngt, &qp, tree_ocp_qp_memory);
	d_cvt_colmaj_to_tree_ocp_qp(hAt, hBt, hbt, hQt, hSt, hRt, hqt, hrt, hidxbt, hd_lbt, hd_ubt, hCt, hDt, hd_lgt, hd_ugt, &qp);

#if 0
	struct d_strmat *tmat;
	struct d_strvec *tvec;
	for(ii=0; ii<Nn-1; ii++)
		{
		tmat = qp.BAbt+ii;
		d_print_strmat(tmat->m, tmat->n, tmat, 0, 0);
		}
	for(ii=0; ii<Nn-1; ii++)
		{
		tvec = qp.b+ii;
		d_print_tran_strvec(tvec->m, tvec, 0);
		}
	for(ii=0; ii<Nn; ii++)
		{
		tmat = qp.RSQrq+ii;
		d_print_strmat(tmat->m, tmat->n, tmat, 0, 0);
		}
	for(ii=0; ii<Nn; ii++)
		{
		tvec = qp.rq+ii;
		d_print_tran_strvec(tvec->m, tvec, 0);
		}
	for(ii=0; ii<Nn; ii++)
		{
		tvec = qp.d_lb+ii;
		d_print_tran_strvec(tvec->m, tvec, 0);
		}
	for(ii=0; ii<Nn; ii++)
		{
		tvec = qp.d_ub+ii;
		d_print_tran_strvec(tvec->m, tvec, 0);
		}
	for(ii=0; ii<Nn; ii++)
		{
		tvec = qp.d_lg+ii;
		d_print_tran_strvec(tvec->m, tvec, 0);
		}
	for(ii=0; ii<Nn; ii++)
		{
		tvec = qp.d_ug+ii;
		d_print_tran_strvec(tvec->m, tvec, 0);
		}
	for(ii=0; ii<Nn; ii++)
		{
		int_print_mat(1, qp.nb[ii], qp.idxb[ii], 1);
		}
#endif

/************************************************
* ocp qp sol
************************************************/	
	
/************************************************
* ipm
************************************************/	

/************************************************
* extract and print solution
************************************************/	

/************************************************
* free memory
************************************************/	

	free(A);
	free(B);
	free(b);
	free(x0);
	free(b0);
	free(Q);
	free(R);
	free(S);
	free(q);
	free(r);
	free(r0);
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

	free(tree_memory);
	free(tree_ocp_qp_memory);

	return 0;

	}
