// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// hpipm
#include "hpipm_d_ocp_qp.h"
// mex
#include "mex.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{

//	mexPrintf("\nin ocp_qp_destroy\n");

	long long *ptr;

	/* RHS */

	// qp_mem
	ptr = (long long *) mxGetData( prhs[0] );
	free( (void *) ptr[0] );

	return;

	}


