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



void PRINT_OCP_QP_DIM(struct OCP_QP_DIM *qp_dim)
	{
	int ii;

	int N   = qp_dim->N;
	int *nx = qp_dim->nx;
	int *nu = qp_dim->nu;
	int *nbx = qp_dim->nbx;
	int *nbu = qp_dim->nbu;
	int *ng = qp_dim->ng;
	int *ns = qp_dim->ns;

	printf("N = %d\n\n", N);

	printf("nx =\n");
	for (ii = 0; ii <= N; ii++)
		printf("\t%d", nx[ii]);
	printf("\n\n");

	printf("nu =\n");
	for (ii = 0; ii <= N; ii++)
		printf("\t%d", nu[ii]);
	printf("\n\n");

	printf("nbx =\n");
	for (ii = 0; ii <= N; ii++)
		printf("\t%d", nbx[ii]);
	printf("\n\n");

	printf("nbu =\n");
	for (ii = 0; ii <= N; ii++)
		printf("\t%d", nbu[ii]);
	printf("\n\n");

	printf("ng =\n");
	for (ii = 0; ii <= N; ii++)
		printf("\t%d", ng[ii]);
	printf("\n\n");

	printf("ns =\n");
	for (ii = 0; ii <= N; ii++)
		printf("\t%d", ns[ii]);
	printf("\n\n");

	return;
	}



void PRINT_OCP_QP(struct OCP_QP *qp)
	{
	int ii;

	struct OCP_QP_DIM *dim = qp->dim;

	int N   = dim->N;
	int *nx = dim->nx;
	int *nu = dim->nu;
	int *nb = dim->nb;
	int *ng = dim->ng;
	int *ns = dim->ns;

	printf("BAt =\n");
	for (ii = 0; ii < N; ii++)
		BLASFEO_PRINT_MAT(nu[ii]+nx[ii], nx[ii+1], qp->BAbt+ii, 0, 0);

	printf("b =\n");
	for (ii = 0; ii < N; ii++)
		BLASFEO_PRINT_TRAN_VEC(nx[ii+1], qp->b+ii, 0);

	printf("RSQ =\n");
	for (ii = 0; ii <= N; ii++)
		BLASFEO_PRINT_MAT(nu[ii]+nx[ii], nu[ii]+nx[ii], qp->RSQrq+ii, 0, 0);

	printf("Z =\n");
	for (ii = 0; ii <= N; ii++)
		BLASFEO_PRINT_TRAN_VEC(2*ns[ii], qp->Z+ii, 0);

	printf("rqz =\n");
	for (ii = 0; ii <= N; ii++)
		BLASFEO_PRINT_TRAN_VEC(nu[ii]+nx[ii]+ns[ii], qp->rqz+ii, 0);

	printf("idxb = \n");
	for (ii = 0; ii <= N; ii++)
		int_print_mat(1, nb[ii], qp->idxb[ii], 1);

	printf("d =\n");
	for (ii = 0; ii <= N; ii++)
		BLASFEO_PRINT_TRAN_VEC(2*nb[ii]+2*ng[ii]+2*ns[ii], qp->d+ii, 0);

	printf("DCt =\n");
	for (ii = 0; ii <= N; ii++)
		BLASFEO_PRINT_MAT(nu[ii]+nx[ii], ng[ii], qp->DCt+ii, 0, 0);

	printf("m =\n");
	for (ii = 0; ii <= N; ii++)
		BLASFEO_PRINT_TRAN_VEC(2*nb[ii]+2*ng[ii]+2*ns[ii], qp->m+ii, 0);

	return;
	}



void PRINT_OCP_QP_SOL(struct OCP_QP_SOL *qp_sol, struct OCP_QP_DIM *qp_dim)
	{
	int ii;

	int N   = qp_dim->N;
	int *nx = qp_dim->nx;
	int *nu = qp_dim->nu;
	int *nb = qp_dim->nb;
	int *ng = qp_dim->ng;
	int *ns = qp_dim->ns;

	printf("ux =\n");
	for (ii = 0; ii <= N; ii++)
		BLASFEO_PRINT_TRAN_VEC(nu[ii] + nx[ii] + 2 * ns[ii], &qp_sol->ux[ii], 0);

	printf("pi =\n");
	for (ii = 0; ii < N; ii++)
		BLASFEO_PRINT_TRAN_VEC(nx[ii + 1], &qp_sol->pi[ii], 0);

	printf("lam =\n");
	for (ii = 0; ii <= N; ii++)
		BLASFEO_PRINT_TRAN_VEC(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], &qp_sol->lam[ii], 0);

	printf("t =\n");
	for (ii = 0; ii <= N; ii++)
		BLASFEO_PRINT_TRAN_VEC(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], &qp_sol->t[ii], 0);

	return;
	}
