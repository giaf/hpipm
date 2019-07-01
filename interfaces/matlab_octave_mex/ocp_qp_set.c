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

//	mexPrintf("\nin ocp_qp_set\n");

	long long *l_ptr;

	int ii;

	/* RHS */

	// qp
	l_ptr = mxGetData( prhs[0] );
	struct d_ocp_qp *qp = (struct d_ocp_qp *) *l_ptr;

	// field
	char *field = mxArrayToString( prhs[1] );

	// value
	double *value = mxGetData( prhs[2] );

	// stage0
	int stage0 = mxGetScalar( prhs[3] );

	// stage1
	int stage1;
	if(nrhs==5)
		{
		stage1 = mxGetScalar( prhs[4] );
		}
	
	/* body */

	if(nrhs==5)
		{
		for(ii=stage0; ii<=stage1; ii++)
			{
			d_cvt_colmaj_gen_to_ocp_qp(field, ii, value, qp);
			}
		}
	else
		{
		d_cvt_colmaj_gen_to_ocp_qp(field, stage0, value, qp);
		}

	return;

	}



