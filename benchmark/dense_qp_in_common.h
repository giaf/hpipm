/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    acados is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with acados; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#ifndef DENSE_QP_IN_COMMON_H_
#define DENSE_QP_IN_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif
// blasfeo
#include <blasfeo_target.h>
#include <blasfeo_common.h>
// hpipm
#include "../include/hpipm_d_dense_qp.h"

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
#define QPOASES_MAX_STRING_LENGTH 160
#define ALIGNMENT 64

/*================================================================*/
// struct of qpoases_to_hpipm
typedef struct {
    int *idxb;
    double *C_eq;
    double *b;
    double *C_ieq;
    double *d_lg0;
    double *d_ug0;
} qpoases_to_hpipm;

int qpoases_to_hpipm_calculate_size(int nv,
                                    int ne,
                                    int nc);

char *assign_qpoases_to_hpipm(int nv,
                              int ne,
                              int nc,
                              qpoases_to_hpipm **tran_space,
                              void *ptr);

qpoases_to_hpipm *create_qpoases_to_hpipm(int nv,
                                          int ne,
                                          int nc);

/*================================================================*/
// struct of dense_qp_in_qpoases
typedef struct{
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
} dense_qp_in_qpoases;


int dense_qp_in_qpoases_calculate_size(int nv,
                                       int ne,
                                       int nc);


char *assign_dense_qp_in_qpoases(int nv,
                                 int ne,
                                 int nc,
                                 dense_qp_in_qpoases **qpd_qpoases,
                                 void *ptr);


dense_qp_in_qpoases *create_dense_qp_in_qpoases(int nv,
                                                int ne,
                                                int nc);


/*================================================================*/
int qpd_qpoases_to_hpipm(const dense_qp_in_qpoases *qpd_qpoases,
                         struct d_dense_qp *qpd,
                         void *memory_);



/*================================================================*/
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
