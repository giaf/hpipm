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

//	mexPrintf("\nin ocp_qp_dim_set\n");

	long long *l_ptr;

	int ii;

	/* RHS */

	// dim
	l_ptr = mxGetData( prhs[0] );
	struct d_ocp_qp_dim *dim = (struct d_ocp_qp_dim *) *l_ptr;

	// field
	char *field = mxArrayToString( prhs[1] );

	// value
	int value = mxGetScalar( prhs[2] );

	// stage0
	int stage0 = mxGetScalar( prhs[3] );

	// stage1
	int stage1;
	if(nrhs==5)
		{
		stage1 = mxGetScalar( prhs[4] );
		}
	
	if(nrhs==5)
		{
		for(ii=stage0; ii<=stage1; ii++)
			{
			d_set_ocp_qp_dim(field, ii, value, dim);
			}
		}
	else
		{
		d_set_ocp_qp_dim(field, stage0, value, dim);
		}

	return;

	}


