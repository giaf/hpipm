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
