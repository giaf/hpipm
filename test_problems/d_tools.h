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

void dgemv_n_3l(int m, int n, double *A, int lda , double *x, double *z);
void dgemm_nn_3l(int m, int n, int k, double *A, int lda , double *B, int ldb, double *C, int ldc);
void daxpy_3l(int n, double da, double *dx, double *dy);
void dscal_3l(int n, double da, double *dx);

/* copies a matrix into another matrix */
void dmcopy(int row, int col, double *ptrA, int lda, double *ptrB, int ldb);

/* solution of a system of linear equations */
void dgesv_3l(int n, int nrhs, double *A, int lda, int *ipiv, double *B, int ldb, int *info);

/* matrix exponential */
void expm(int row, double *A);
