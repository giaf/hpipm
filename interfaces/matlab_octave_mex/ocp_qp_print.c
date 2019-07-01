// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// hpipm
#include "hpipm_d_ocp_qp.h"
#include "hpipm_d_ocp_qp_utils.h"
// mex
#include "mex.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{

//	mexPrintf("\nin ocp_qp_dim_print\n");

	long long *l_ptr;

	/* RHS */

	// qp
	l_ptr = mxGetData( prhs[0] );
	struct d_ocp_qp *qp = (struct d_ocp_qp *) *l_ptr;

	d_print_ocp_qp(qp);

	return;

	}




