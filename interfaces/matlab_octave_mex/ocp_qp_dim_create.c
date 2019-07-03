// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// hpipm
#include "hpipm_d_ocp_qp_dim.h"
// mex
#include "mex.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{

//	printf("\nin ocp_qp_dim_create\n");

	mxArray *tmp_mat;
	long long *l_ptr;
	char *c_ptr;

	/* RHS */

	int N = mxGetScalar( prhs[0] );

	/* body */

	int dim_size = sizeof(struct d_ocp_qp_dim) + d_ocp_qp_dim_memsize(N);
	void *dim_mem = malloc(dim_size);

	c_ptr = dim_mem;

	struct d_ocp_qp_dim *dim = (struct d_ocp_qp_dim *) c_ptr;
	c_ptr += sizeof(struct d_ocp_qp_dim);

	d_ocp_qp_dim_create(N, dim, c_ptr);
	c_ptr += dim->memsize;

	/* LHS */

	tmp_mat = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
	l_ptr = mxGetData(tmp_mat);
	l_ptr[0] = (long long) dim_mem;
	plhs[0] = tmp_mat;

	return;

	}


