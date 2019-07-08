// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// hpipm
#include "hpipm_d_ocp_qp_ipm.h"
// mex
#include "mex.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{

//	mexPrintf("\nin ocp_solver_arg_destroy\n");

	long long *ptr;

	/* RHS */

	// arg_mem
	ptr = (long long *) mxGetData( prhs[0] );
	free( (void *) ptr[0] );

	return;

	}



