// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// hpipm
#include "hpipm_d_ocp_qp_dim.h"
#include "hpipm_d_ocp_qp_ipm.h"
// mex
#include "mex.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{

//	printf("\nin ocp_solver_create\n");

	mxArray *tmp_mat;
	long long *l_ptr;
	char *c_ptr;

	/* RHS */

	// dim
	l_ptr = mxGetData( prhs[0] );
	struct d_ocp_qp_dim *dim = (struct d_ocp_qp_dim *) *l_ptr;

	// arg
	l_ptr = mxGetData( prhs[1] );
	struct d_ocp_qp_ipm_arg *arg = (struct d_ocp_qp_ipm_arg *) *l_ptr;

	/* body */

	int ws_size = sizeof(struct d_ocp_qp_ipm_ws) + d_ocp_qp_ipm_ws_memsize(dim, arg);
	void *ws_mem = malloc(ws_size);

	c_ptr = ws_mem;

	struct d_ocp_qp_ipm_ws *ws = (struct d_ocp_qp_ipm_ws *) c_ptr;
	c_ptr += sizeof(struct d_ocp_qp_ipm_ws);

	d_ocp_qp_ipm_ws_create(dim, arg, ws, c_ptr);
	c_ptr += d_ocp_qp_ipm_ws_memsize(dim, arg);

	/* LHS */

	tmp_mat = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
	l_ptr = mxGetData(tmp_mat);
	l_ptr[0] = (long long) ws_mem;
	plhs[0] = tmp_mat;

	return;

	}





