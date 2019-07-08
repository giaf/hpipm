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

//	mexPrintf("\nin ocp_qp_sol_get\n");

	long long *l_ptr;

	int ii, jj;

	/* RHS */

	// ws
	l_ptr = mxGetData( prhs[0] );
	struct d_ocp_qp_ipm_ws *ws = (struct d_ocp_qp_ipm_ws *) *l_ptr;

	// field
	char *field = mxArrayToString( prhs[1] );

	if(!strcmp(field, "status") | !strcmp(field, "iter"))
		{
		plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
		double *mat_ptr = mxGetPr( plhs[0] );
		int tmp_int;
		d_ocp_qp_ipm_get(field, ws, &tmp_int);
		*mat_ptr = (double) tmp_int;
		}
	else if(!strcmp(field, "res_stat") | !strcmp(field, "res_eq") | !strcmp(field, "res_ineq") | !strcmp(field, "res_comp"))
		{
		plhs[0] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
		double *mat_ptr = mxGetPr( plhs[0] );
		d_ocp_qp_ipm_get(field, ws, mat_ptr);
		}
	else if(!strcmp(field, "stat"))
		{
		int iter;
		int stat_m;
		double *stat;
		d_ocp_qp_ipm_get("iter", ws, &iter);
		d_ocp_qp_ipm_get("stat_m", ws, &stat_m);
		d_ocp_qp_ipm_get("stat", ws, &stat);
		plhs[0] = mxCreateNumericMatrix(iter+1, stat_m+1, mxDOUBLE_CLASS, mxREAL);
		double *mat_ptr = mxGetPr( plhs[0] );
		for(ii=0; ii<iter+1; ii++)
			{
			mat_ptr[ii+0] = ii;
			for(jj=0; jj<stat_m; jj++)
				{
				mat_ptr[ii+(jj+1)*(iter+1)] = stat[jj+ii*stat_m];
				}
			}
		}
	else
		{
		mexPrintf("\nocp_qp_solver_get: field not supported: %s\n", field);
		return;
		}

	return;

	}


