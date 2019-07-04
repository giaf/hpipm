// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// hpipm
#include "hpipm_d_ocp_qp_dim.h"
#include "hpipm_d_ocp_qp.h"
// mex
#include "mex.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{

//	printf("\nin ocp_qp_create\n");

	mxArray *tmp_mat;
	long long *l_ptr;
	char *c_ptr;

	/* RHS */

	// dim
	l_ptr = mxGetData( prhs[0] );
	struct d_ocp_qp_dim *dim = (struct d_ocp_qp_dim *) *l_ptr;

	/* body */

	int qp_size = sizeof(struct d_ocp_qp) + d_ocp_qp_memsize(dim);
	void *qp_mem = malloc(qp_size);

	c_ptr = qp_mem;

	struct d_ocp_qp *qp = (struct d_ocp_qp *) c_ptr;
	c_ptr += sizeof(struct d_ocp_qp);

	d_ocp_qp_create(dim, qp, c_ptr);
	c_ptr += d_ocp_qp_memsize(dim);

	/* LHS */

	tmp_mat = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
	l_ptr = mxGetData(tmp_mat);
	l_ptr[0] = (long long) qp_mem;
	plhs[0] = tmp_mat;

	return;

	}



