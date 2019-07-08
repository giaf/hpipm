// system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
// hpipm
#include "hpipm_d_ocp_qp_dim.h"
#include "hpipm_d_ocp_qp_sol.h"
// mex
#include "mex.h"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
	{

//	mexPrintf("\nin ocp_qp_sol_get\n");

	long long *l_ptr;

	int ii;

	/* RHS */

	// sol
	l_ptr = mxGetData( prhs[0] );
	struct d_ocp_qp_sol *sol = (struct d_ocp_qp_sol *) *l_ptr;

	// dim
	struct d_ocp_qp_dim *dim = sol->dim;
	int N = dim->N;
	int *nx = dim->nx;
	int *nu = dim->nu;

	// field
	char *field = mxArrayToString( prhs[1] );

	// stage0
	int stage0 = mxGetScalar( prhs[2] );

	// stage1
	int stage1;
	if(nrhs==4)
		{
		stage1 = mxGetScalar( prhs[3] );
		}
	
	/* body */

	if(!strcmp(field, "x"))
		{
		if(nrhs==4)
			{
			int nx_sum = 0;
			for(ii=stage0; ii<=stage1; ii++)
				{
				nx_sum += nx[ii];
				}
			plhs[0] = mxCreateNumericMatrix(nx_sum, 1, mxDOUBLE_CLASS, mxREAL);
			double *x = mxGetPr( plhs[0] );
			nx_sum = 0;
			for(ii=stage0; ii<=stage1; ii++)
				{
				d_ocp_qp_sol_get(field, ii, sol, x+nx_sum);
				nx_sum += nx[ii];
				}
			}
		else
			{
			plhs[0] = mxCreateNumericMatrix(nx[stage0], 1, mxDOUBLE_CLASS, mxREAL);
			double *x = mxGetPr( plhs[0] );
			d_ocp_qp_sol_get(field, stage0, sol, x);
			}
		}
	else if(!strcmp(field, "u"))
		{
		if(nrhs==4)
			{
			int nu_sum = 0;
			for(ii=stage0; ii<=stage1; ii++)
				{
				nu_sum += nu[ii];
				}
			plhs[0] = mxCreateNumericMatrix(nu_sum, 1, mxDOUBLE_CLASS, mxREAL);
			double *u = mxGetPr( plhs[0] );
			nu_sum = 0;
			for(ii=stage0; ii<=stage1; ii++)
				{
				d_ocp_qp_sol_get(field, ii, sol, u+nu_sum);
				nu_sum += nu[ii];
				}
			}
		else
			{
			plhs[0] = mxCreateNumericMatrix(nu[stage0], 1, mxDOUBLE_CLASS, mxREAL);
			double *u = mxGetPr( plhs[0] );
			d_ocp_qp_sol_get(field, stage0, sol, u);
			}
		}
	else
		{
		mexPrintf("\nocp_qp_sol_get: field not supported: %s\n", field);
		return;
		}

	return;

	}




