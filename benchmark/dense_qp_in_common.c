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

#include "dense_qp_in_common.h"

#include <dirent.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

// blasfeo
#include <blasfeo_target.h>
#include <blasfeo_common.h>
// hpipm
#include "../include/hpipm_d_dense_qp.h"

/*================================================================*/
int qpoases_to_hpipm_calculate_size(int nv,
                                    int ne,
                                    int nc) {

    int bytes = sizeof(qpoases_to_hpipm);
    // size of pointer
    bytes += 1 * sizeof(int *); // idxb
    bytes += 5 * sizeof(double *); // C_eq, b, C_ieq, d_lg0, d_ug0

    // size of double, int
    bytes += nv * sizeof(int); // idxb
    bytes += ne * nv * sizeof(double);  // C_eq
    bytes += ne * sizeof(double); // b
    bytes += nc * nv * sizeof(double); // C_ieq
    bytes += 2 * nc * sizeof(double); // d_lg0, d_ug0

    bytes = (bytes+ALIGNMENT-1)/ALIGNMENT*ALIGNMENT;
    bytes += ALIGNMENT;

    return bytes;
}


char *assign_qpoases_to_hpipm(int nv,
                              int ne,
                              int nc,
                              qpoases_to_hpipm **tran_space,
                              void *ptr) {

    // pointer to initialize QP data to zero
//    char *c_ptr_QPdata;

    // char pointer
    char *c_ptr = (char *) ptr;

    *tran_space = (qpoases_to_hpipm *) c_ptr;
    c_ptr += sizeof(qpoases_to_hpipm);

    // assign pointers to ints
    (*tran_space)->idxb = (int *) c_ptr;
    c_ptr += nv * sizeof(int);

    // align data
    size_t l_ptr = (size_t) c_ptr;
    l_ptr = (l_ptr+ALIGNMENT-1)/ALIGNMENT*ALIGNMENT;
    c_ptr = (char *) l_ptr;

    // assign pointers to doubles
//    c_ptr_QPdata = c_ptr;

    (*tran_space)->C_eq = (double *) c_ptr;
    c_ptr += ne * nv * sizeof(double);

    (*tran_space)->b = (double *) c_ptr;
    c_ptr += ne * sizeof(double);

    (*tran_space)->C_ieq = (double *) c_ptr;
    c_ptr += nc * nv * sizeof(double);

    (*tran_space)->d_lg0 = (double *) c_ptr;
    c_ptr += nc * sizeof(double);

    (*tran_space)->d_ug0 = (double *) c_ptr;
    c_ptr += nc * sizeof(double);

    // set QP data to zero (mainly for valgrind)
//    for (char *idx = c_ptr_QPdata; idx < c_ptr; idx++)
//        *idx = 0;

    return c_ptr;
}


qpoases_to_hpipm *create_qpoases_to_hpipm(int nv,
                                          int ne,
                                          int nc) {

    qpoases_to_hpipm *tran_space;

    int bytes = qpoases_to_hpipm_calculate_size(nv, ne, nc);

    void *ptr = calloc(bytes, 1);

    // // set a value for debugging
    // char *c_ptr = (char *) ptr;
    // for (int i = 0; i < bytes; i++) c_ptr[i] = 13;

    char *ptr_end = assign_qpoases_to_hpipm(nv, ne, nc, &tran_space, ptr);
    assert((char*)ptr + bytes >= ptr_end); (void) ptr_end;

    // for (int i = 0; i < bytes; i++) printf("%d - ", c_ptr[i]);
    // exit(1);

    return tran_space;
}


/*================================================================*/
int dense_qp_in_qpoases_calculate_size(int nv,
                                       int ne,
                                       int nc) {

    int bytes = sizeof(dense_qp_in_qpoases);
    // size of pointer
    bytes += 1 * sizeof(double *);  // H
    bytes += 3 * sizeof(double *);  // g, lbx, ubx
    bytes += 1 * sizeof(double *);  // C
    bytes += 2 * sizeof(double *);  // lbC, ubC

    // size of double, int
    bytes += nv * nv * sizeof(double);         // H
    bytes += (nv + 2 * nv) * sizeof(double);  // g, lbx, ubx
    bytes += (ne + nc) * nv * sizeof(double); // C
    bytes += 2 * (ne + nc) * sizeof(double);  // lbC, ubC

    bytes = (bytes+ALIGNMENT-1)/ALIGNMENT*ALIGNMENT;
    bytes += ALIGNMENT;

    return bytes;
}



char *assign_dense_qp_in_qpoases(int nv,
                         int ne,
                         int nc,
                         dense_qp_in_qpoases **qpd_qpoases,
                         void *ptr) {

    // pointer to initialize QP data to zero
//    char *c_ptr_QPdata;

    // char pointer
    char *c_ptr = (char *) ptr;

    *qpd_qpoases = (dense_qp_in_qpoases *) c_ptr;
    c_ptr += sizeof(dense_qp_in_qpoases);

    // copy dimensions to workspace
    (*qpd_qpoases)->nv = nv;
    (*qpd_qpoases)->ne = ne;
    (*qpd_qpoases)->nc = nc;

    // align data
    size_t l_ptr = (size_t) c_ptr;
    l_ptr = (l_ptr+ALIGNMENT-1)/ALIGNMENT*ALIGNMENT;
    c_ptr = (char *) l_ptr;

    // assign pointers to doubles
//    c_ptr_QPdata = c_ptr;

    (*qpd_qpoases)->H = (double *) c_ptr;
    c_ptr += nv*nv*sizeof(double);

    (*qpd_qpoases)->g = (double *) c_ptr;
    c_ptr += nv*sizeof(double);

    (*qpd_qpoases)->lbx = (double *) c_ptr;
    c_ptr += nv*sizeof(double);

    (*qpd_qpoases)->ubx = (double *) c_ptr;
    c_ptr += nv*sizeof(double);

    (*qpd_qpoases)->C = (double *) c_ptr;
    c_ptr += (ne+nc)*nv*sizeof(double);

    (*qpd_qpoases)->lbC = (double *) c_ptr;
    c_ptr += (ne+nc)*sizeof(double);

    (*qpd_qpoases)->ubC = (double *) c_ptr;
    c_ptr += (ne+nc)*sizeof(double);

    // set QP data to zero (mainly for valgrind)
//    for (char *idx = c_ptr_QPdata; idx < c_ptr; idx++)
//        *idx = 0;

    return c_ptr;
}



dense_qp_in_qpoases *create_dense_qp_in_qpoases(int nv,
                                                int ne,
                                                int nc) {

    dense_qp_in_qpoases *qpd_qpoases;

    int bytes = dense_qp_in_qpoases_calculate_size(nv, ne, nc);

    void *ptr = calloc(bytes,1);

    // // set a value for debugging
    // char *c_ptr = (char *) ptr;
    // for (int i = 0; i < bytes; i++) c_ptr[i] = 13;

    char *ptr_end = assign_dense_qp_in_qpoases(nv, ne, nc, &qpd_qpoases, ptr);
    assert((char*)ptr + bytes >= ptr_end); (void) ptr_end;

    // for (int i = 0; i < bytes; i++) printf("%d - ", c_ptr[i]);
    // exit(1);

    return qpd_qpoases;
}


/*================================================================*/
int qpd_qpoases_to_hpipm(const dense_qp_in_qpoases *qpd_qpoases,
                         struct d_dense_qp *qpd,
                         void *memory_) {
     qpoases_to_hpipm *memory =
        (qpoases_to_hpipm *)memory_;
     /* extract qpoases dense qp in */
     int nvd = qpd_qpoases->nv;
     int ned = qpd_qpoases->ne;
     int ncd = qpd_qpoases->nc;
     double *H = qpd_qpoases->H;
     double *g = qpd_qpoases->g;
     double *C = qpd_qpoases->C;
     double *d_lb = qpd_qpoases->lbx;
     double *d_ub = qpd_qpoases->ubx;
     double *d_lg = qpd_qpoases->lbC;
     double *d_ug = qpd_qpoases->ubC;

     /* construct transfer workspace */
     int *idxb = memory->idxb;
     double *C_eq = memory->C_eq;
     double *b = memory->b;
     double *C_ieq = memory->C_ieq;
     double *d_lg0 = memory->d_lg0;
     double *d_ug0 = memory->d_ug0;

     int ii,jje,kk,jji;
     jje = 0;
     jji = 0;
     for (ii = 0; ii < nvd; ii++) {
       /* full box constraint */
       idxb[ii] = ii;
     }

     for (ii = 0; ii < ned+ncd; ii++) {
       /* split C into C_ieq and C_eq */
       if ( d_lg[ii] == d_ug[ii]) {
            for (kk = 0; kk < nvd; kk++) {
              C_eq[jje*nvd + kk] = C[ii*nvd + kk];
            }
            b[jje] = d_lg[ii];
            jje += 1;
        }
       else {
            for (kk = 0; kk < nvd; kk++) {
              C_ieq[jji*nvd + kk] = C[ii*nvd + kk];
            }
            d_lg0[jji] = d_lg[ii];
            d_ug0[jji] = d_ug[ii];
            jji += 1;
       }
     }
/*
     printf("A(%d,%d) = \n",jje,nvd);
     int k,j;
     for (j = 0; j < ned; j++) {
       for (k = 0; k < nvd; k++) {
       printf("%f ", C_eq[j*nvd + k]);
       }
       printf("\n");
     }
*/

     d_cvt_rowmaj_to_dense_qp(H, g, C_eq, b, idxb, d_lb, d_ub, C_ieq, d_lg0, d_ug0,
                              NULL, NULL, NULL, NULL, NULL, qpd);

     return 0;
}


/*================================================================*/
/*========== read data from benchmark fold =======================*/

int qpOASES_readFromFileI(int* data, int n,
									        const char* datafilename) {

	int i;
	FILE* datafile;
	char errstr[QPOASES_MAX_STRING_LENGTH];

	/* 1) Open file. */
	if ( ( datafile = fopen( datafilename, "r" ) ) == 0 )
	{
		snprintf( errstr,QPOASES_MAX_STRING_LENGTH,"(%s)",datafilename );
	}

	/* 2) Read data from file. */
	for( i=0; i<n; ++i )
	{
		if ( fscanf( datafile, "%d\n", &(data[i]) ) == 0 )
		{
			fclose( datafile );
			snprintf( errstr,QPOASES_MAX_STRING_LENGTH,"(%s)",datafilename );
		}
	}

	/* 3) Close file. */
	fclose( datafile );

	return 0;

}

int readOQPdimensions(const char* path,
	     	              int* nQP, int* nV,
                      int* nC, int* nEC){
	int dims[4];

	/* 1) Setup file name where dimensions are stored. */
	char filename[QPOASES_MAX_STRING_LENGTH];
	snprintf( filename, QPOASES_MAX_STRING_LENGTH, "%sdims.oqp", path );

	/* 2) Load dimensions from file. */
	qpOASES_readFromFileI( dims, 4, filename );

	*nQP = dims[0];
	*nV  = dims[1];
	*nC  = dims[2];
	*nEC = dims[3];

//  printf( "nQP = %d,  nV = %d,  nC = %d,  nEC = %d\n", *nQP, *nV, *nC, *nEC );

	return 0;
}

int qpOASES_readFromFileM(double* data, int nrow, int ncol,
									        const char* datafilename) {

	int i, j;
	double float_data;
	FILE* datafile;
	char errstr[QPOASES_MAX_STRING_LENGTH];

	/* 1) Open file. */
	if ( ( datafile = fopen( datafilename, "r" ) ) == 0 )
	{
		snprintf( errstr,QPOASES_MAX_STRING_LENGTH,"(%s)",datafilename );
	}

	/* 2) Read data from file. */
	for( i=0; i<nrow; ++i )
	{
		for( j=0; j<ncol; ++j )
		{
			#ifdef __USE_SINGLE_PRECISION__
			if ( fscanf( datafile, "%f ", &float_data ) == 0 )
			#else
			if ( fscanf( datafile, "%lf ", &float_data ) == 0 )
			#endif /* __USE_SINGLE_PRECISION__ */
			{
				fclose( datafile );
				snprintf( errstr,QPOASES_MAX_STRING_LENGTH,"(%s)",datafilename );
			}
			data[i*ncol + j] = ( (double) float_data );
		}
	}

	/* 3) Close file. */
	fclose( datafile );

	return 0;

}

int readOQPdata(const char* path,
					      int* nQP, int* nV,
                int* nC, int* nEC,
						    double* H, double* g,
                double* A,
                double* lb, double* ub,
                double* lbA,
                double* ubA,
							  double* xOpt, double* yOpt,
                double* objOpt) {
	char filename[QPOASES_MAX_STRING_LENGTH];

	/* 1) Obtain OQP dimensions. */
	readOQPdimensions( path, nQP, nV, nC,nEC );

	/* 2) Allocate memory and load OQP data: */
	/* Hessian matrix */
	snprintf( filename,QPOASES_MAX_STRING_LENGTH,"%sH.oqp",path );
  qpOASES_readFromFileM( H,(*nV),(*nV),filename );

	/* gradient vector sequence */
	snprintf( filename,QPOASES_MAX_STRING_LENGTH,"%sg.oqp",path );
	qpOASES_readFromFileM( g,(*nQP),(*nV),filename );

	/* lower bound vector sequence */
	snprintf( filename,QPOASES_MAX_STRING_LENGTH,"%slb.oqp",path );
  qpOASES_readFromFileM( lb,(*nQP),(*nV),filename );

	/* upper bound vector sequence */
	snprintf( filename,QPOASES_MAX_STRING_LENGTH,"%sub.oqp",path );
  qpOASES_readFromFileM( ub,(*nQP),(*nV),filename );

	if ( (*nC) > 0 )
	{
		/* Constraint matrix */
		snprintf( filename,QPOASES_MAX_STRING_LENGTH,"%sA.oqp",path );
		qpOASES_readFromFileM( A,(*nC),(*nV),filename );

		/* lower constraints' bound vector sequence */
		snprintf( filename,QPOASES_MAX_STRING_LENGTH,"%slbA.oqp",path );
		qpOASES_readFromFileM( lbA,(*nQP),(*nC),filename );

		/* upper constraints' bound vector sequence */
		snprintf( filename,QPOASES_MAX_STRING_LENGTH,"%subA.oqp",path );
		qpOASES_readFromFileM( ubA,(*nQP),(*nC),filename );
	}

	if ( xOpt != 0 )
	{
		/* primal solution vector sequence */
		snprintf( filename,QPOASES_MAX_STRING_LENGTH,"%sx_opt.oqp",path );
		qpOASES_readFromFileM( xOpt,(*nQP),(*nV),filename );
	}

	if ( yOpt != 0 )
	{
		/* dual solution vector sequence */
		snprintf( filename,QPOASES_MAX_STRING_LENGTH,"%sy_opt.oqp",path );
		qpOASES_readFromFileM( yOpt,(*nQP),(*nV)+(*nC),filename );
	}

	if ( objOpt != 0 )
	{
		/* dual solution vector sequence */
		snprintf( filename,QPOASES_MAX_STRING_LENGTH,"%sobj_opt.oqp",path );
		qpOASES_readFromFileM( objOpt,(*nQP),1,filename );
	}

	return 0;
}
