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

#ifndef DENSE_QP_IN_COMMON_H_
#define DENSE_QP_IN_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif

#define _SVID_SOURCE
// standard
#include <dirent.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <assert.h>
// blasfeo
#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_v_aux_ext_dep.h>
#include <blasfeo_d_aux_ext_dep.h>
#include <blasfeo_i_aux_ext_dep.h>
#include <blasfeo_d_aux.h>
#include <blasfeo_d_blas.h>
// hpipm
#include "../include/hpipm_d_dense_qp.h"
#include "../include/hpipm_d_dense_qp_sol.h"
#include "../include/hpipm_d_dense_qp_ipm.h"

#define MAX_STRING_LENGTH 160
#define EPSILON 1e-5 // regularization parameter

/* define global */
#ifdef __NO_SNPRINTF__
  #if (!defined(_MSC_VER)) || defined(__DSPACE__) || defined(__XPCTARGET__)
    /* If snprintf is not available, provide an empty implementation... */
    int snprintf( char* s, size_t n, const char* format, ... );
  #else
	/* ... or substitute snprintf by _snprintf for Microsoft compilers. */
    #define snprintf _snprintf
  #endif
#endif /* __NO_SNPRINTF__ */

struct benchmark_qp
  {
    int nv;
    int ne; // equality C
    int nc; // ineqality C
    double *H;
    double *g;
    double *lbx;
    double *ubx;
    double *C;
    double *lbC;
    double *ubC;
   };

struct benchmark_to_hpipm
  {
    int *idxb;
    double *C_eq;
    double *b;
    double *C_ieq;
    double *d_lg0;
    double *d_ug0;
  };

int d_memsize_benchmark_qp(int nv, int ne, int nc);

void d_create_benchmark_qp(int nv, int ne, int nc, struct benchmark_qp *qp, void *mem);

int d_memsize_benchmark_to_hpipm(int nv, int ne, int nc);

void d_create_benchmark_to_hpipm(int nv, int ne, int nc, struct benchmark_to_hpipm *tran_space, void *mem);

int d_benchmark_to_hpipm(struct benchmark_qp *qp_bench,
                         struct d_dense_qp *qpd,
                         struct benchmark_to_hpipm *tran_space);

int readOQPdimensions(const char* path,	/**< Full path of the data files (without trailing slash!). */
      								int* nQP,			    /**< Output: Number of QPs. */
      								int* nV,			    /**< Output: Number of variables. */
      								int* nC,			    /**< Output: Number of constraints. */
      								int* nEC			    /**< Output: Number of equality constraints. */
      								);

int readOQPdata(const char* path,	/**< Full path of the data files (without trailing slash!). */
  							int* nQP,    			/**< Output: Number of QPs. */
  							int* nV,		     	/**< Output: Number of variables. */
  							int* nC,			    /**< Output: Number of constraints. */
  							int* nEC,			    /**< Output: Number of equality constraints. */
  							double* H,		 	  /**< Output: Hessian matrix. */
  							double* g,		 	  /**< Output: Sequence of gradient vectors. */
  							double* A,		 	  /**< Output: Constraint matrix. */
  							double* lb,			  /**< Output: Sequence of lower bound vectors (on variables). */
  							double* ub,			  /**< Output: Sequence of upper bound vectors (on variables). */
  							double* lbA,		  /**< Output: Sequence of lower constraints' bound vectors. */
  							double* ubA,		  /**< Output: Sequence of upper constraints' bound vectors. */
  							double* xOpt,		  /**< Output: Sequence of primal solution vectors
  												 *           (not read if a null pointer is passed). */
  							double* yOpt,		  /**< Output: Sequence of dual solution vectors
  												 *           (not read if a null pointer is passed). */
  							double* objOpt		/**< Output: Sequence of optimal objective function values
  												 *           (not read if a null pointer is passed). */
  							);



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif  // DENSE_QP_IN_COMMON_H_
