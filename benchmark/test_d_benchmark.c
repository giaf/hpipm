/**************************************************************************************************
*                                                                                                 *
* This file is part of HPIPM.                                                                     *
*                                                                                                 *
* HPIPM -- High-Performance Interior Point Method.                                                *
* Copyright (C) 2019 by Gianluca Frison.                                                          *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
*                                                                                                 *
* The 2-Clause BSD License                                                                        *
*                                                                                                 *
* Redistribution and use in source and binary forms, with or without                              *
* modification, are permitted provided that the following conditions are met:                     *
*                                                                                                 *
* 1. Redistributions of source code must retain the above copyright notice, this                  *
*    list of conditions and the following disclaimer.                                             *
* 2. Redistributions in binary form must reproduce the above copyright notice,                    *
*    this list of conditions and the following disclaimer in the documentation                    *
*    and/or other materials provided with the distribution.                                       *
*                                                                                                 *
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND                 *
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED                   *
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE                          *
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR                 *
* ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES                  *
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;                    *
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND                     *
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT                      *
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS                   *
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                                    *
*                                                                                                 *
* Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             *
*                                                                                                 *
**************************************************************************************************/

#include <sys/time.h>



#include "d_benchmark.h"



#define REG 1e-8



//#define SINGLE
#define DOUBLE



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



int cvt_benchmark_to_hpipm(struct benchmark_qp *qp_bench, struct d_dense_qp *qpd, struct benchmark_to_hpipm *tran_space)
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

	int ii, jje, kk, jji;

	// box constraint on all variables
	for (ii = 0; ii < nvd; ii++)
		{
		idxb[ii] = ii;
		}

	// split C into C_ieq and C_eq
	jje = 0;
	jji = 0;
	for (ii = 0; ii < ned+ncd; ii++) {
	if ( d_lg[ii] == d_ug[ii])
		{
		for (kk = 0; kk < nvd; kk++)
			{
			C_eq[jje*nvd + kk] = C[ii*nvd + kk];
			}
		b[jje] = d_lg[ii];
		jje += 1;
		}
	else
		{
		for (kk = 0; kk < nvd; kk++)
			{
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
						  NULL, NULL, NULL, NULL, NULL, NULL, NULL, qpd);

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
    printf(" HPIPM -- High-Performance Interior Point Method.\n");
    printf("\n");
    printf("\n");
    printf("\n");

    int nQP, nv, ne, ng;
    int nproblems, i, j;
    struct dirent **namelist;
    char resstr[200], OQPproblem[200];
    char *problem;
    nproblems = scandir("./projects.coin-or.org/svn/qpOASES/misc/testingdata/cpp/problems", &namelist, NULL, alphasort);
    /*
    int ii,jj;
    char filename[1024];
    FILE * pFile;
    */

	printf("probl\tnv\tne\tnc\tdp\treturn\titer\tres_g\t\tres_b\t\tres_d\t\tres_m\t\tmu\t\ttime\t\ttime[ms]\n");

	int npass = 0;
	int nfail = 0;
	int s_npass = 0;
	int s_nfail = 0;
	int dp;

	int ii;

    /************** benchmark loop *********************/

//    for (i = 0; i < nproblems; i++)
    for (i = 0; i < nproblems-1; i++)
//    for (i = 44; i < 45; i++)
		{

        /************************************************
        * bechmark data setting
        ************************************************/

        /* skip special directories and zip file cuter.*bz2 */
        if (namelist[i]->d_name[0] == '.' || namelist[i]->d_name[0] == 'c')
			{
            free(namelist[i]);
            continue;
			}
        problem = namelist[i]->d_name;
        snprintf(OQPproblem, 199, "./projects.coin-or.org/svn/qpOASES/misc/testingdata/cpp/problems/%s/", problem);

        /* read dimensions */
        readOQPdimensions( OQPproblem, &nQP, &nv, &ng, &ne );

        /************************************************
        * dense qp benchmark
        ************************************************/

        int nc = ng-ne; // inequality constraint
        int benchmark_size = d_memsize_benchmark_qp(nv, ne, nc);
        void *benchmark_mem = calloc(benchmark_size,1);

        struct benchmark_qp qp_bench;
        d_create_benchmark_qp(nv, ne, nc, &qp_bench, benchmark_mem);

        /* read data */
        readOQPdata(OQPproblem, &nQP, &nv, &ng, &ne, qp_bench.H, qp_bench.g, qp_bench.C, qp_bench.lbx, qp_bench.ubx, qp_bench.lbC, qp_bench.ubC, NULL, NULL, NULL);

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

        int tran_size = d_memsize_benchmark_to_hpipm(nv, ne, nc);
        void *tran_mem = calloc(tran_size,1);

        struct benchmark_to_hpipm tran_space;
        d_create_benchmark_to_hpipm(nv, ne, nc, &tran_space, tran_mem);

        /************************************************
        * dense qp dim
        ************************************************/

		// double

		int qp_dim_size = d_memsize_dense_qp_dim();
		void *qp_dim_mem = calloc(qp_dim_size, 1);

		struct d_dense_qp_dim dim;
		d_create_dense_qp_dim(&dim, qp_dim_mem);

		d_cvt_int_to_dense_qp_dim(nv, ne, nv, nc, 0, 0, &dim);

		// single

		int s_qp_dim_size = s_memsize_dense_qp_dim();
		void *s_qp_dim_mem = calloc(s_qp_dim_size, 1);

		struct s_dense_qp_dim s_dim;
		s_create_dense_qp_dim(&s_dim, s_qp_dim_mem);

		cvt_d2s_dense_qp_dim(&dim, &s_dim);

        /************************************************
        * dense qp
        ************************************************/

		// double

        int qp_size = d_memsize_dense_qp(&dim);
        void *qp_mem = calloc(qp_size, 1);

        struct d_dense_qp qpd_hpipm;
        d_create_dense_qp(&dim, &qpd_hpipm, qp_mem);

        cvt_benchmark_to_hpipm(&qp_bench, &qpd_hpipm, &tran_space);

		double reg = 0e-8;
		blasfeo_ddiare(nv, reg, qpd_hpipm.Hv, 0, 0);

		int H_fact_size = blasfeo_memsize_dmat(nv, nv);
		void *H_fact_mem; v_zeros_align(&H_fact_mem, H_fact_size);

		struct blasfeo_dmat H_fact;
		blasfeo_create_dmat(nv, nv, &H_fact, H_fact_mem);

		blasfeo_dpotrf_l(nv, qpd_hpipm.Hv, 0, 0, &H_fact, 0, 0);

		dp = 1;
		for(ii=0; ii<nv; ii++)
			if(H_fact.dA[ii]<=0.0)
				dp = 0;

		// single

        int s_qp_size = s_memsize_dense_qp(&s_dim);
        void *s_qp_mem = calloc(s_qp_size, 1);

        struct s_dense_qp s_qpd_hpipm;
        s_create_dense_qp(&s_dim, &s_qpd_hpipm, s_qp_mem);

		cvt_d2s_dense_qp(&qpd_hpipm, &s_qpd_hpipm);

        /************************************************
        * dense sol
        ************************************************/

		// double

        int qp_sol_size = d_memsize_dense_qp_sol(&dim);
        void *qp_sol_mem = calloc(qp_sol_size,1);

        struct d_dense_qp_sol qpd_sol;
        d_create_dense_qp_sol(&dim, &qpd_sol, qp_sol_mem);

		// single

        int s_qp_sol_size = s_memsize_dense_qp_sol(&s_dim);
        void *s_qp_sol_mem = calloc(s_qp_sol_size,1);

        struct s_dense_qp_sol s_qpd_sol;
        s_create_dense_qp_sol(&s_dim, &s_qpd_sol, s_qp_sol_mem);

        /************************************************
        * ipm arg
        ************************************************/

//		enum dense_qp_ipm_mode mode = SPEED_ABS;
//		enum dense_qp_ipm_mode mode = SPEED;
		enum dense_qp_ipm_mode mode = BALANCE;
//		enum dense_qp_ipm_mode mode = ROBUST;

		// double

        int ipm_arg_size = d_memsize_dense_qp_ipm_arg(&dim);
        void *ipm_arg_mem = calloc(ipm_arg_size,1);

        struct d_dense_qp_ipm_arg argd;
        d_create_dense_qp_ipm_arg(&dim, &argd, ipm_arg_mem);
        d_set_default_dense_qp_ipm_arg(mode, &argd);
		// overwirte default
        argd.res_g_max = 1e-6;
        argd.res_b_max = 1e-6;
        argd.res_d_max = 1e-6;
        argd.res_m_max = 1e-6;
        argd.iter_max = 200;
        argd.stat_max = 200;
//        argd.alpha_min = 1e-12;
        argd.mu0 = 1e1;
//		argd.pred_corr = 1;
//		argd.cond_pred_corr = 1;
//		argd.scale = 1;
//		argd.itref_pred_max = 0;
//		argd.itref_corr_max = 4;
//		argd.reg_prim = 1e-15;
//		argd.reg_dual = 1e-15;
//		argd.lq_fact = 1;

		// single

        int s_ipm_arg_size = s_memsize_dense_qp_ipm_arg(&s_dim);
        void *s_ipm_arg_mem = calloc(s_ipm_arg_size,1);

        struct s_dense_qp_ipm_arg s_argd;
        s_create_dense_qp_ipm_arg(&s_dim, &s_argd, s_ipm_arg_mem);
        s_set_default_dense_qp_ipm_arg(mode, &s_argd);
		// overwirte default
        s_argd.res_g_max = 1e-3;
        s_argd.res_b_max = 1e-3;
        s_argd.res_d_max = 1e-3;
        s_argd.res_m_max = 1e-3;
        s_argd.iter_max = 200;
        s_argd.stat_max = 200;
//		s_argd.alpha_min = 1e-12;
        s_argd.mu0 = 1e1;
//		s_argd.pred_corr = 1;
//		s_argd.cond_pred_corr = 1;
//		s_argd.scale = 1;
//		s_argd.itref_pred_max = 0;
		s_argd.itref_corr_max = 4;
		s_argd.reg_prim = 1e-7;
		s_argd.reg_dual = 1e-7;
//		s_argd.lq_fact = 1;

        /************************************************
        * dense ipm
        ************************************************/

		struct timeval tv0, tv1;

		int rep, nrep=1;

		// double

        int ipm_size = d_memsize_dense_qp_ipm(&dim, &argd);

        void *ipm_mem = calloc(ipm_size,1);
        struct d_dense_qp_ipm_workspace workspace;
        d_create_dense_qp_ipm(&dim, &argd, &workspace, ipm_mem);

        int hpipm_return; // 0 normal; 1 max iter; 2 alpha_minl; 3 NaN

		double sol_time;

		gettimeofday(&tv0, NULL); // start

		for(rep=0; rep<nrep; rep++)
			{
#ifdef DOUBLE
			hpipm_return = d_solve_dense_qp_ipm(&qpd_hpipm, &qpd_sol, &argd, &workspace);
#endif
			}

		gettimeofday(&tv1, NULL); // stop

		sol_time = (tv1.tv_sec-tv0.tv_sec)/(nrep+0.0)+(tv1.tv_usec-tv0.tv_usec)/(nrep*1e6);

		// single

        int s_ipm_size = s_memsize_dense_qp_ipm(&s_dim, &s_argd);

        void *s_ipm_mem = calloc(s_ipm_size,1);
        struct s_dense_qp_ipm_workspace s_workspace;
        s_create_dense_qp_ipm(&s_dim, &s_argd, &s_workspace, s_ipm_mem);

        int s_hpipm_return; // 0 normal; 1 max iter; 2 alpha_minl; 3 NaN

		double s_sol_time;

		gettimeofday(&tv0, NULL); // start

		for(rep=0; rep<nrep; rep++)
			{
#ifdef SINGLE
			s_hpipm_return = s_solve_dense_qp_ipm(&s_qpd_hpipm, &s_qpd_sol, &s_argd, &s_workspace);
#endif
			}

		gettimeofday(&tv1, NULL); // stop

		s_sol_time = (tv1.tv_sec-tv0.tv_sec)/(nrep+0.0)+(tv1.tv_usec-tv0.tv_usec)/(nrep*1e6);

        /************************************************
        * print ipm statistics
        ************************************************/
#if 0
//		if(i==17)
		if(1)
			{
			printf("\nipm iter = %d\n", workspace.iter);
			printf("\nalpha_aff\tmu_aff\t\tsigma\t\talpha\t\tmu\n");
			d_print_exp_tran_mat(5, workspace.iter, workspace.stat, 5);
			printf("\n\n\n\n");
			}
#endif

		// number of passed and failed problems
		if(hpipm_return==0)
			npass++;
		else
			nfail++;

#ifdef DOUBLE
		printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%e\t%e\t%e\t%e\t%e\t%e%12.4f\n", i-1, nv, ne, nc, dp, hpipm_return, workspace.iter, workspace.qp_res[0], workspace.qp_res[1], workspace.qp_res[2], workspace.qp_res[3], workspace.res->res_mu, sol_time, sol_time*1000);
#endif

		// number of passed and failed problems
		if(s_hpipm_return==0)
			s_npass++;
		else
			s_nfail++;

#ifdef SINGLE
		printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%e\t%e\t%e\t%e\t%e\t%e%12.4f\n", i-1, nv, ne, nc, dp, s_hpipm_return, s_workspace.iter, s_workspace.qp_res[0], s_workspace.qp_res[1], s_workspace.qp_res[2], s_workspace.qp_res[3], s_workspace.res->res_mu, s_sol_time, s_sol_time*1000);
#endif

        /************************************************
        * free memory
        ************************************************/

		// double
        free(benchmark_mem);
        free(tran_mem);
        free(qp_mem);
		free(H_fact_mem);
      	free(qp_sol_mem);
      	free(ipm_mem);
        free(ipm_arg_mem);
		// single
        free(s_qp_mem);
      	free(s_qp_sol_mem);
      	free(s_ipm_mem);
        free(s_ipm_arg_mem);

		}

    /************************************************
    * overall statistics
    ************************************************/

#ifdef DOUBLE
	printf("\n\nTestbench results (double precision):\n");
	printf("=====================\n");
//	printf("\nTotal:  %d\n", nproblems-2);
	printf("Pass:\t%3d\n", npass);
	printf("Fail:\t%3d\n", nfail);
	printf("Ratio:\t%5.1f\n", 100.0 * (double)npass / (double) (npass+nfail));
	printf("\n\n");
#endif

#ifdef SINGLE
	printf("\n\nTestbench results (single precision):\n");
	printf("=====================\n");
//	printf("\nTotal:  %d\n", nproblems-2);
	printf("Pass:\t%3d\n", s_npass);
	printf("Fail:\t%3d\n", s_nfail);
	printf("Ratio:\t%5.1f\n", 100.0 * (double)s_npass / (double) (s_npass+s_nfail));
	printf("\n\n");
#endif

    /************************************************
    * return
    ************************************************/

  	return 0;

	}
