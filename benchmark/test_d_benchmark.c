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


#include "d_benchmark.h"

int d_memsize_benchmark_qp(int nv, int ne, int nc)
{
    int size = 0;

    // size of double, int
    size += nv * nv * sizeof(double);         // H
    size += (nv + 2 * nv) * sizeof(double);  // g, lbx, ubx
    size += (ne + nc) * nv * sizeof(double); // C
    size += 2 * (ne + nc) * sizeof(double);  // lbC, ubC

    size += 8;

    return size;
}


void d_create_benchmark_qp(int nv, int ne, int nc, struct benchmark_qp *qp, void *mem)
{

    // problem size
    qp->nv = nv;
    qp->ne = ne;
    qp->nc = nc;

    // char pointer
    double *c_ptr = (double *) mem;

    qp->H = c_ptr;
    c_ptr += nv * nv;

    qp->g = c_ptr;
    c_ptr += nv;

    qp->lbx = c_ptr;
    c_ptr += nv;

    qp->ubx = c_ptr;
    c_ptr += nv;

    qp->C = c_ptr;
    c_ptr += (ne+nc)*nv;

    qp->lbC = c_ptr;
    c_ptr += ne + nc;

    qp->ubC = c_ptr;
    c_ptr += ne + nc;

}

int d_memsize_benchmark_to_hpipm(int nv, int ne, int nc)
{

    int size = 0;

    // size of double, int
    size += nv * sizeof(int); // idxb
    size += ne * nv * sizeof(double);  // C_eq
    size += ne * sizeof(double); // b
    size += nc * nv * sizeof(double); // C_ieq
    size += 2 * nc * sizeof(double); // d_lg0, d_ug0

    size = (size+8-1)/8*8;
    size += 8;

    return size;
}


void d_create_benchmark_to_hpipm(int nv, int ne, int nc, struct benchmark_to_hpipm *tran_space, void *mem)
{

    // int pointer
    int *i_ptr = (int *) mem;

    // assign pointers to ints
    tran_space->idxb = i_ptr;
    i_ptr += nv;

    // align data
    size_t s_ptr = (size_t) i_ptr;
  	s_ptr = (s_ptr+7)/8*8;

    // char pointer
  	double *c_ptr = (double *) s_ptr;

    tran_space->C_eq = c_ptr;
    c_ptr += ne * nv;

    tran_space->b = c_ptr;
    c_ptr += ne;

    tran_space->C_ieq = c_ptr;
    c_ptr += nc * nv;

    tran_space->d_lg0 = c_ptr;
    c_ptr += nc;

    tran_space->d_ug0 = c_ptr;
    c_ptr += nc;

}

int d_benchmark_to_hpipm(struct benchmark_qp *qp_bench,
                         struct d_dense_qp *qpd,
                         struct benchmark_to_hpipm *tran_space)
{

     /* extract benchmark qp */
     int nvd = qp_bench->nv;
     int ned = qp_bench->ne;
     int ncd = qp_bench->nc;
     double *H = qp_bench->H;
     double *g = qp_bench->g;
     double *C = qp_bench->C;
     double *d_lb = qp_bench->lbx;
     double *d_ub = qp_bench->ubx;
     double *d_lg = qp_bench->lbC;
     double *d_ug = qp_bench->ubC;

     /* construct transfer workspace */
     int *idxb = tran_space->idxb;
     double *C_eq = tran_space->C_eq;
     double *b = tran_space->b;
     double *C_ieq = tran_space->C_ieq;
     double *d_lg0 = tran_space->d_lg0;
     double *d_ug0 = tran_space->d_ug0;

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
/*---------- read data from benchmark folder ---------------------*/

int readFromFileI(int* data, int n, const char* datafilename)
{

	int i;
	FILE* datafile;
	char errstr[MAX_STRING_LENGTH];

	/* 1) Open file. */
	if ( ( datafile = fopen( datafilename, "r" ) ) == 0 )
	{
		snprintf( errstr,MAX_STRING_LENGTH,"(%s)",datafilename );
	}

	/* 2) Read data from file. */
	for( i=0; i<n; ++i )
	{
		if ( fscanf( datafile, "%d\n", &(data[i]) ) == 0 )
		{
			fclose( datafile );
			snprintf( errstr,MAX_STRING_LENGTH,"(%s)",datafilename );
		}
	}

	/* 3) Close file. */
	fclose( datafile );

	return 0;

}

int readOQPdimensions(const char* path,
	     	              int* nQP, int* nV,
                      int* nC, int* nEC)
{

	int dims[4];

	/* 1) Setup file name where dimensions are stored. */
	char filename[MAX_STRING_LENGTH];
	snprintf( filename, MAX_STRING_LENGTH, "%sdims.oqp", path );

	/* 2) Load dimensions from file. */
	readFromFileI( dims, 4, filename );

	*nQP = dims[0];
	*nV  = dims[1];
	*nC  = dims[2];
	*nEC = dims[3];

//  printf( "nQP = %d,  nV = %d,  nC = %d,  nEC = %d\n", *nQP, *nV, *nC, *nEC );

	return 0;

}

int readFromFileM(double* data, int nrow, int ncol, const char* datafilename)
{
	int i, j;
	double float_data;
	FILE* datafile;
	char errstr[MAX_STRING_LENGTH];

	/* 1) Open file. */
	if ( ( datafile = fopen( datafilename, "r" ) ) == 0 )
	{
		snprintf( errstr,MAX_STRING_LENGTH,"(%s)",datafilename );
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
				snprintf( errstr,MAX_STRING_LENGTH,"(%s)",datafilename );
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
                double* objOpt)
{
	char filename[MAX_STRING_LENGTH];

	/* 1) Obtain OQP dimensions. */
	readOQPdimensions( path, nQP, nV, nC,nEC );

	/* 2) Allocate memory and load OQP data: */
	/* Hessian matrix */
	snprintf( filename,MAX_STRING_LENGTH,"%sH.oqp",path );
  readFromFileM( H,(*nV),(*nV),filename );

	/* gradient vector sequence */
	snprintf( filename,MAX_STRING_LENGTH,"%sg.oqp",path );
	readFromFileM( g,(*nQP),(*nV),filename );

	/* lower bound vector sequence */
	snprintf( filename,MAX_STRING_LENGTH,"%slb.oqp",path );
  readFromFileM( lb,(*nQP),(*nV),filename );

	/* upper bound vector sequence */
	snprintf( filename,MAX_STRING_LENGTH,"%sub.oqp",path );
  readFromFileM( ub,(*nQP),(*nV),filename );

	if ( (*nC) > 0 )
	{
		/* Constraint matrix */
		snprintf( filename,MAX_STRING_LENGTH,"%sA.oqp",path );
		readFromFileM( A,(*nC),(*nV),filename );

		/* lower constraints' bound vector sequence */
		snprintf( filename,MAX_STRING_LENGTH,"%slbA.oqp",path );
		readFromFileM( lbA,(*nQP),(*nC),filename );

		/* upper constraints' bound vector sequence */
		snprintf( filename,MAX_STRING_LENGTH,"%subA.oqp",path );
		readFromFileM( ubA,(*nQP),(*nC),filename );
	}

	if ( xOpt != 0 )
	{
		/* primal solution vector sequence */
		snprintf( filename,MAX_STRING_LENGTH,"%sx_opt.oqp",path );
		readFromFileM( xOpt,(*nQP),(*nV),filename );
	}

	if ( yOpt != 0 )
	{
		/* dual solution vector sequence */
		snprintf( filename,MAX_STRING_LENGTH,"%sy_opt.oqp",path );
		readFromFileM( yOpt,(*nQP),(*nV)+(*nC),filename );
	}

	if ( objOpt != 0 )
	{
		/* dual solution vector sequence */
		snprintf( filename,MAX_STRING_LENGTH,"%sobj_opt.oqp",path );
		readFromFileM( objOpt,(*nQP),1,filename );
	}

	return 0;
}

/************************************************
        test benchmark problem
************************************************/
int main()
{
    printf("\n");
    printf("\n");
    printf("\n");
    printf(
        " HPIPM -- High-Performance Interior Point Method.\n");
    printf("\n");
    printf("\n");
    printf("\n");

    int nQP = 0 , nvc = 0, nec = 0, ngc = 0;
    int nproblems, i, j;
    struct dirent **namelist;
    char resstr[200], OQPproblem[200];
    char *problem;
    nproblems = scandir("./problems", &namelist, NULL, alphasort);
    /*
    int ii,jj;
    char filename[1024];
    FILE * pFile;
    */
    /************** benchmark loop *********************/
    for (i = 0; i < nproblems; i++) {

        /************************************************
        * bechmark data setting
        ************************************************/

        /* skip special directories and zip file cuter.*bz2 */
        if (namelist[i]->d_name[0] == '.' || namelist[i]->d_name[0] == 'c') {
            free(namelist[i]);
            continue;
        }
        problem = namelist[i]->d_name;
        snprintf(OQPproblem, 199, "./problems/%s/", problem);

        /* read dimensions */
        readOQPdimensions( OQPproblem, &nQP, &nvc, &ngc, &nec );

        /************************************************
        * dense qp benchmark
        ************************************************/

        int nc = ngc-nec; // inequality constraint
        int benchmark_size = d_memsize_benchmark_qp(nvc, nec, nc);
        void *benchmark_mem = calloc(benchmark_size,1);

        struct benchmark_qp qp_bench;
        d_create_benchmark_qp(nvc, nec, nc, &qp_bench, benchmark_mem);

        /* read data */
        readOQPdata(OQPproblem, &nQP, &nvc, &ngc, &nec, qp_bench.H, qp_bench.g, qp_bench.C, qp_bench.lbx, qp_bench.ubx, qp_bench.lbC, qp_bench.ubC, NULL, NULL, NULL);

        // print data to text files
        /*
        snprintf(filename, sizeof(filename), "matrixH%d.txt", i);
        pFile = fopen(filename,"w");
        for (ii = 0; ii < nv; ii++){
           for (jj = 0; jj < nv; jj++)
           {
               fprintf(pFile, "%e ", H[ii*nv+jj]);
           }
           fputc('\n', pFile);
        }
        fclose(pFile);
        */

        /************************************************
        * benchmark to hpipm workspace
        ************************************************/

        int tran_size = d_memsize_benchmark_to_hpipm(nvc, nec, nc);
        void *tran_mem = calloc(tran_size,1);

        struct benchmark_to_hpipm tran_space;
        d_create_benchmark_to_hpipm(nvc, nec, nc, &tran_space, tran_mem);

        /************************************************
        * dense qp dim
        ************************************************/

        int nsc = 0;

		int qp_dim_size = d_memsize_dense_qp_dim();
		void *qp_dim_mem = calloc(qp_dim_size, 1);

		struct d_dense_qp_dim dim;
		d_create_dense_qp_dim(&dim, qp_dim_mem);

		d_cvt_int_to_dense_qp_dim(nvc, nec, nvc, nc, nsc, &dim);

        /************************************************
        * dense qp
        ************************************************/

        int qp_size = d_memsize_dense_qp(&dim);
        void *qp_mem = calloc(qp_size,1);

        struct d_dense_qp qpd_hpipm;
        d_create_dense_qp(&dim, &qpd_hpipm, qp_mem);

        /* qp_benchmark -> qpd_hpipm */
        d_benchmark_to_hpipm(&qp_bench, &qpd_hpipm, &tran_space);

        /************************************************
        * dense sol
        ************************************************/

        int qp_sol_size = d_memsize_dense_qp_sol(&dim);
        void *qp_sol_mem = calloc(qp_sol_size,1);

        struct d_dense_qp_sol qpd_sol;
        d_create_dense_qp_sol(&dim, &qpd_sol, qp_sol_mem);

        /************************************************
        * ipm arg
        ************************************************/

        int ipm_arg_size = d_memsize_dense_qp_ipm_arg(&dim);
        void *ipm_arg_mem = calloc(ipm_arg_size,1);

        struct d_dense_qp_ipm_arg argd;
        d_create_dense_qp_ipm_arg(&dim, &argd, ipm_arg_mem);
        d_set_default_dense_qp_ipm_arg(&argd);
        /* consistent with setting in acore */
        argd.res_g_max = 1e-6;
        argd.res_b_max = 1e-8;
        argd.res_d_max = 1e-8;
        argd.res_m_max = 1e-8;
        argd.iter_max = 50;
        argd.stat_max = 50;
        argd.alpha_min = 1e-8;
        argd.mu0 = 1;

        /************************************************
        * dense ipm
        ************************************************/
        int ipm_size = d_memsize_dense_qp_ipm(&dim, &argd);
        void *ipm_mem = calloc(ipm_size,1);

        struct d_dense_qp_ipm_workspace workspace;
        d_create_dense_qp_ipm(&dim, &argd, &workspace, ipm_mem);

        int hpipm_return; // 0 normal; 1 max iter

        hpipm_return = d_solve_dense_qp_ipm(&qpd_hpipm, &qpd_sol, &argd, &workspace);

        /************************************************
        * print ipm statistics
        ************************************************/
        printf("Problem %d\n", i-1);

        if (hpipm_return == 1) {
           /* print original H*/
        //            printf("\nH_org =\n");
        //            d_print_strmat(nvc, nvc, qpd_hpipm.Hg, 0, 0);
            for (j = 0; j < nvc; j++) {
                 qp_bench.H[j*nvc+j] = qp_bench.H[j*nvc+j] + EPSILON;
            }
            /* print modefied H*/
        //            printf("\nH_reg =\n");
        //            d_print_strmat(nvc, nvc, qpd_hpipm.Hg, 0, 0);

            d_benchmark_to_hpipm(&qp_bench, &qpd_hpipm, &tran_space);
            hpipm_return = d_solve_dense_qp_ipm(&qpd_hpipm, &qpd_sol, &argd, &workspace);

            /* print primal solution */
            printf("\n\nipm return = %d\n", hpipm_return);
        //            printf("\nnew_primal_sol = \n");
        //            d_print_strvec(nvc, qpd_sol.v, 0);
         }

        printf("\nipm iter = %d\n", workspace.iter);
        printf(" inf norm res: %e, %e, %e, %e, %e\n", workspace.qp_res[0],
                                 workspace.qp_res[1], workspace.qp_res[2],
                                 workspace.qp_res[3], workspace.res_mu);
        printf("\n\n\n\n");

        /************************************************
        * free memory
        ************************************************/

        free(benchmark_mem);
        free(tran_mem);
        free(qp_mem);
      	free(qp_sol_mem);
      	free(ipm_mem);
        free(ipm_arg_mem);

    }

    /************************************************
    * return
    ************************************************/

  	return 0;

}
