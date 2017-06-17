/**************************************************************************************************
*                                                                                                 *
* This file is part of HPMPC.                                                                     *
*                                                                                                 *
* HPMPC -- Library for High-Performance implementation of solvers for MPC.                        *
* Copyright (C) 2014-2015 by Technical University of Denmark. All rights reserved.                *
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
* Author: Gianluca Frison, giaf (at) dtu.dk                                                       *
*                                                                                                 *
**************************************************************************************************/

void sgemv_n_3l(int m, int n, float *A, int lda , float *x, float *z);
void sgemm_nn_3l(int m, int n, int k, float *A, int lda , float *B, int ldb, float *C, int ldc);
void saxpy_3l(int n, float da, float *dx, float *dy);
void sscal_3l(int n, float da, float *dx);

/* copies a matrix into another matrix */
void smcopy(int row, int col, float *ptrA, int lda, float *ptrB, int ldb);

/* solution of a system of linear equations */
void sgesv_3l(int n, int nrhs, float *A, int lda, int *ipiv, float *B, int ldb, int *info);

/* matrix exponential */
void expm(int row, float *A);
