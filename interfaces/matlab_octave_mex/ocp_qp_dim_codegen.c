// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// hpipm
#include "hpipm_d_ocp_qp_dim.h"
#include "hpipm_d_ocp_qp_utils.h"
// mex
#include "mex.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{

//	mexPrintf("\nin ocp_qp_dim_print\n");

	long long *l_ptr;

	/* RHS */

	// dim
	l_ptr = mxGetData( prhs[0] );
	struct d_ocp_qp_dim *dim = (struct d_ocp_qp_dim *) *l_ptr;

	// field
	char *file_name = mxArrayToString( prhs[1] );

	// mode
	char *mode = mxArrayToString( prhs[2] );

	/* body */
	d_codegen_ocp_qp_dim(file_name, mode, dim);

	return;

	}




