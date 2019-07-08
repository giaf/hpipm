// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// hpipm
#include "hpipm_d_ocp_qp_dim.h"
#include "hpipm_d_ocp_qp_sol.h"
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

	// sol
	l_ptr = mxGetData( prhs[1] );
	struct d_ocp_qp_sol *sol = (struct d_ocp_qp_sol *) *l_ptr;

	d_ocp_qp_sol_print(dim, sol);

	return;

	}





