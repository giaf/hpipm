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

//	mexPrintf("\nin ocp_qp_set\n");

	long long *l_ptr;

	int ii;

	/* RHS */

	// arg
	l_ptr = mxGetData( prhs[0] );
	struct d_ocp_qp_ipm_arg *arg = (struct d_ocp_qp_ipm_arg *) *l_ptr;

	// field
	char *field = mxArrayToString( prhs[1] );

	// manually set integers
	if(!strcmp(field, "iter_max") | !strcmp(field, "warm_start") | !strcmp(field, "pred_corr") | !strcmp(field, "ric_alg"))
		{
		int value = mxGetScalar( prhs[2] );
		d_ocp_qp_ipm_arg_set(field, &value, arg);
		}
	else // real
		{
		// value
		double *value = mxGetData( prhs[2] );
		d_ocp_qp_ipm_arg_set(field, value, arg);
		}

	return;

	}




