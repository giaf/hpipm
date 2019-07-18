// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// hpipm
#include "hpipm_timing.h"
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

	// qp
	l_ptr = mxGetData( prhs[0] );
	struct d_ocp_qp *qp = (struct d_ocp_qp *) *l_ptr;

	// qp_sol
	l_ptr = mxGetData( prhs[1] );
	struct d_ocp_qp_sol *qp_sol = (struct d_ocp_qp_sol *) *l_ptr;

	// arg
	l_ptr = mxGetData( prhs[2] );
	struct d_ocp_qp_ipm_arg *arg = (struct d_ocp_qp_ipm_arg *) *l_ptr;

	// ws
	l_ptr = mxGetData( prhs[3] );
	struct d_ocp_qp_ipm_ws *ws = (struct d_ocp_qp_ipm_ws *) *l_ptr;

	/* RHS */

	// time_ext
	plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
	double *mat_ptr = mxGetPr( plhs[0] );

	/* body */

	hpipm_timer timer;
	hpipm_tic(&timer);

	d_ocp_qp_ipm_solve(qp, qp_sol, arg, ws);

	double time_ext = hpipm_toc(&timer);
	*mat_ptr = time_ext;

	return;

	}






