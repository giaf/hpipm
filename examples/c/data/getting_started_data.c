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

#include <stdlib.h>



// QP size

// horizon lenght
int N = 5;
// number of input
static int nnu[6] = {1, 1, 1, 1, 1, 0};
// number of states
static int nnx[6] = {2, 2, 2, 2, 2, 2};
// number of input box constraints
static int nnbu[6] = {0, 0, 0, 0, 0, 0};
// number of states box constraints
static int nnbx[6] = {2, 0, 0, 0, 0, 0};
// number of general constraints
static int nng[6] = {0, 0, 0, 0, 0, 0};
// number of softed constraints on state box constraints
static int nnsbx[6] = {0, 0, 0, 0, 0, 0};
// number of softed constraints on input box constraints
static int nnsbu[6] = {0, 0, 0, 0, 0, 0};
// number of softed constraints on general constraints
static int nnsg[6] = {0, 0, 0, 0, 0, 0};


// QP data

//
static double A[] = {1, 0, 1, 1};
//
static double B[] = {0, 1};
//
static double b[] = {0, 0};

//
static double Q[] = {1, 0, 0, 1};
//
static double R[] = {1};
//
static double S[] = {0, 0};
//
static double q[] = {1, 1};
//
static double r[] = {0};

//
static double lbx0[] = {1, 1};
//
static double ubx0[] = {1, 1};
//
static int idxbx0[] = {0, 1};

//
static double u_guess[] = {0};
//
static double x_guess[] = {0, 0};
//
static double sl_guess[] = {};
//
static double su_guess[] = {};

// array of pointers

//
static double *AA[5] = {A, A, A, A, A};
//
static double *BB[5] = {B, B, B, B, B};
//
static double *bb[5] = {b, b, b, b, b};
//
static double *QQ[6] = {Q, Q, Q, Q, Q, Q};
//
static double *RR[6] = {R, R, R, R, R, R};
//
static double *SS[6] = {S, S, S, S, S, S};
//
static double *qq[6] = {q, q, q, q, q, q};
//
static double *rr[6] = {r, r, r, r, r, r};
//
static int *iidxbx[6] = {idxbx0, NULL, NULL, NULL, NULL, NULL};
//
static double *llbx[6] = {lbx0, NULL, NULL, NULL, NULL, NULL};
//
static double *uubx[6] = {ubx0, NULL, NULL, NULL, NULL, NULL};
//
static int *iidxbu[6] = {};
//
static double *llbu[6] = {};
//
static double *uubu[6] = {};
//
static double *CC[6] = {};
//
static double *DD[6] = {};
//
static double *llg[6] = {};
//
static double *uug[6] = {};
//
static double *ZZl[6] = {};
//
static double *ZZu[6] = {};
//
static double *zzl[6] = {};
//
static double *zzu[6] = {};
//
static int *iidxs[6] = {};
//
static double *llls[6] = {};
//
static double *llus[6] = {};

//
static double *uu_guess[6] = {u_guess, u_guess, u_guess, u_guess, u_guess, u_guess};
//
static double *xx_guess[6] = {x_guess, x_guess, x_guess, x_guess, x_guess, x_guess};
//
static double *ssl_guess[6] = {sl_guess, sl_guess, sl_guess, sl_guess, sl_guess, sl_guess};
//
static double *ssu_guess[6] = {su_guess, su_guess, su_guess, su_guess, su_guess, su_guess};



// export as global data

int *nu = nnu;
int *nx = nnx;
int *nbu = nnbu;
int *nbx = nnbx;
int *ng = nng;
int *nsbx = nnsbx;
int *nsbu = nnsbu;
int *nsg = nnsg;

double **hA = AA;
double **hB = BB;
double **hb = bb;
double **hQ = QQ;
double **hR = RR;
double **hS = SS;
double **hq = qq;
double **hr = rr;
int **hidxbx = iidxbx;
double **hlbx = llbx;
double **hubx = uubx;
int **hidxbu = iidxbu;
double **hlbu = llbu;
double **hubu = uubu;
double **hC = CC;
double **hD = DD;
double **hlg = llg;
double **hug = uug;
double **hZl = ZZl;
double **hZu = ZZu;
double **hzl = zzl;
double **hzu = zzu;
int **hidxs = iidxs;
double **hlls = llls;
double **hlus = llus;

double **hu_guess = uu_guess;
double **hx_guess = xx_guess;
double **hsl_guess = ssl_guess;
double **hsu_guess = ssu_guess;

// arg
int mode = 1;
int iter_max = 30;
double alpha_min = 1e-8;
double mu0 = 1e4;
double tol_stat = 1e-4;
double tol_eq = 1e-5;
double tol_ineq = 1e-5;
double tol_comp = 1e-5;
double reg_prim = 1e-12;
int warm_start = 0;
int pred_corr = 1;
int ric_alg = 0;

