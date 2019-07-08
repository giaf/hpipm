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

//	printf("\nin ocp_solver_arg_create\n");

	mxArray *tmp_mat;
	long long *l_ptr;
	char *c_ptr;

	/* RHS */

	// dim
	l_ptr = mxGetData( prhs[0] );
	struct d_ocp_qp_dim *dim = (struct d_ocp_qp_dim *) *l_ptr;

	// mode
	char *str_mode = mxArrayToString( prhs[1] );

	int mode;
	if(!strcmp(str_mode, "speed_abs"))
		{
		mode = SPEED_ABS;
		}
	else if(!strcmp(str_mode, "speed"))
		{
		mode = SPEED;
		}
	else if(!strcmp(str_mode, "balance"))
		{
		mode = BALANCE;
		}
	else if(!strcmp(str_mode, "robust"))
		{
		mode = ROBUST;
		}
	else
		{
		mode = SPEED;
		mexPrintf("\nocp_qp_solver_arg_create: mode not supported: %s; speed mode used instead\n", str_mode);
		}

	/* body */

	int arg_size = sizeof(struct d_ocp_qp_ipm_arg) + d_ocp_qp_ipm_arg_memsize(dim);
	void *arg_mem = malloc(arg_size);

	c_ptr = arg_mem;

	struct d_ocp_qp_ipm_arg *arg = (struct d_ocp_qp_ipm_arg *) c_ptr;
	c_ptr += sizeof(struct d_ocp_qp_ipm_arg);

	d_ocp_qp_ipm_arg_create(dim, arg, c_ptr);
	c_ptr += d_ocp_qp_ipm_arg_memsize(dim);

	d_ocp_qp_ipm_arg_set_default(mode, arg);

	/* LHS */

	tmp_mat = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
	l_ptr = mxGetData(tmp_mat);
	l_ptr[0] = (long long) arg_mem;
	plhs[0] = tmp_mat;

	return;

	}




