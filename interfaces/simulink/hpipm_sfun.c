#define S_FUNCTION_NAME hpipm_solver_sfunction
#define S_FUNCTION_LEVEL  2

#define MDL_START

// HPIPM
#include <blasfeo_d_aux_ext_dep.h>

#include "hpipm_d_ocp_qp_ipm.h"
#include "hpipm_d_ocp_qp_dim.h"
#include "hpipm_d_ocp_qp.h"
#include "hpipm_d_ocp_qp_sol.h"

#define SAMPLINGTIME -1

// memory as global data
void *dim_mem;
void *qp_mem;
void *qp_sol_mem;
void *ipm_arg_mem;
void *ipm_mem; 

// qp data as global data
extern int N;
extern int *nx;
extern int *nu;
extern int *nbu;
extern int *nbx;
extern int *ng;
extern int *ns;
extern double **hA;
extern double **hB;
extern double **hb;
extern double **hQ;
extern double **hR;
extern double **hS;
extern double **hq;
extern double **hr;
extern int **hidxb;
extern double **hlb;
extern double **hub;
extern double **hC;
extern double **hD;
extern double **hlg;
extern double **hug;
extern double **hZl;
extern double **hZu;
extern double **hzl;
extern double **hzu;
extern int **hidxs;
extern double **hlls;
extern double **hlus;
extern double **hu_guess;
extern double **hx_guess;
extern double **hsl_guess;
extern double **hsu_guess;

static void mdlInitializeSizes (SimStruct *S)
{
    // specify the number of continuous and discrete states 
    ssSetNumContStates(S, 0);
    ssSetNumDiscStates(S, 0);

    // specify the number of input ports 
    if ( !ssSetNumInputPorts(S, 1) )
        return;

    // specify the number of output ports 
    if ( !ssSetNumOutputPorts(S, 5) )
        return;

    // specify dimension information for the input ports 
    ssSetInputPortVectorDimension(S, 0, NX0);

    // specify dimension information for the output ports 
    ssSetOutputPortVectorDimension(S, 0, NU);  // optimal input
    ssSetOutputPortVectorDimension(S, 1, 1 );  // solver status
    ssSetOutputPortVectorDimension(S, 2, 1 );  // KKT residuals
    ssSetOutputPortVectorDimension(S, 3, NX0); // first state
    ssSetOutputPortVectorDimension(S, 4, 1);   // computation times

    // specify the direct feedthrough status 
    ssSetInputPortDirectFeedThrough(S, 0, 1); // current state x0

    // one sample time 
    ssSetNumSampleTimes(S, 1);
    }


#if defined(MATLAB_MEX_FILE)

#define MDL_SET_INPUT_PORT_DIMENSION_INFO
#define MDL_SET_OUTPUT_PORT_DIMENSION_INFO

static void mdlSetInputPortDimensionInfo(SimStruct *S, int_T port, const DimsInfo_T *dimsInfo)
{
    if ( !ssSetInputPortDimensionInfo(S, port, dimsInfo) )
         return;
}

static void mdlSetOutputPortDimensionInfo(SimStruct *S, int_T port, const DimsInfo_T *dimsInfo)
{
    if ( !ssSetOutputPortDimensionInfo(S, port, dimsInfo) )
         return;
}

    #endif /* MATLAB_MEX_FILE */


static void mdlInitializeSampleTimes(SimStruct *S)
{
    ssSetSampleTime(S, 0, SAMPLINGTIME);
    ssSetOffsetTime(S, 0, 0.0);
}


static void mdlStart(SimStruct *S)
{
	int ii;

	int hpipm_return;

	int rep, nrep=10;

	struct timeval tv0, tv1;

	/************************************************
	* ocp qp dim
	************************************************/

	int dim_size = d_memsize_ocp_qp_dim(N);
	void *dim_mem = malloc(dim_size);

	struct d_ocp_qp_dim dim;
	d_create_ocp_qp_dim(N, &dim, dim_mem);

	d_cvt_int_to_ocp_qp_dim(N, nx, nu, nbx, nbu, ng, ns, &dim);

	/************************************************
	* ocp qp
	************************************************/

	int qp_size = d_memsize_ocp_qp(&dim);
	void *qp_mem = malloc(qp_size);

	struct d_ocp_qp qp;
	d_create_ocp_qp(&dim, &qp, qp_mem);

	d_cvt_colmaj_to_ocp_qp(hA, hB, hb, hQ, hS, hR, hq, hr, hidxb, hlb, hub, hC, hD, hlg, hug, hZl, hZu, hzl, hzu, hidxs, hlls, hlus, &qp);

	/************************************************
	* ocp qp sol
	************************************************/

	int qp_sol_size = d_memsize_ocp_qp_sol(&dim);
	void *qp_sol_mem = malloc(qp_sol_size);

	struct d_ocp_qp_sol qp_sol;
	d_create_ocp_qp_sol(&dim, &qp_sol, qp_sol_mem);

	/************************************************
	* ipm arg
	************************************************/

	int ipm_arg_size = d_memsize_ocp_qp_ipm_arg(&dim);
	void *ipm_arg_mem = malloc(ipm_arg_size);

	struct d_ocp_qp_ipm_arg arg;
	d_create_ocp_qp_ipm_arg(&dim, &arg, ipm_arg_mem);

	// enum hpipm_mode mode = SPEED_ABS;
	enum hpipm_mode mode = SPEED;
	// enum hpipm_mode mode = BALANCE;
	// enum hpipm_mode mode = ROBUST;
	d_set_default_ocp_qp_ipm_arg(mode, &arg);

	d_set_ocp_qp_ipm_arg_mu0(1e4, &arg);
	d_set_ocp_qp_ipm_arg_iter_max(30, &arg);
	d_set_ocp_qp_ipm_arg_tol_stat(1e-4, &arg);
	d_set_ocp_qp_ipm_arg_tol_eq(1e-5, &arg);
	d_set_ocp_qp_ipm_arg_tol_ineq(1e-5, &arg);
	d_set_ocp_qp_ipm_arg_tol_comp(1e-5, &arg);
	d_set_ocp_qp_ipm_arg_reg_prim(1e-12, &arg);
	d_set_ocp_qp_ipm_arg_warm_start(0, &arg);

	/************************************************
	* ipm workspace
	************************************************/

	int ipm_size = d_memsize_ocp_qp_ipm(&dim, &arg);
	void *ipm_mem = malloc(ipm_size);

	struct d_ocp_qp_ipm_workspace workspace;
	d_create_ocp_qp_ipm(&dim, &arg, &workspace, ipm_mem);

}

static void mdlOutputs(SimStruct *S, int_T tid)
{
    // get input signals
    InputRealPtrsType in_x0_sign;
    InputRealPtrsType in_y_ref_sign;
    InputRealPtrsType in_y_ref_e_sign;
    {% if ocp.dims.np > 0 %}
    InputRealPtrsType in_p_sign;
    {% endif %}
    
    // local buffers
    real_t in_x0[{{ ocp.dims.nx }}];
    real_t in_y_ref[{{ ocp.dims.ny }}];
    real_t in_y_ref_e[{{ ocp.dims.ny_e }}];
    {% if ocp.dims.np > 0 %}
    real_t in_p[{{ ocp.dims.np }}];
    {% endif %}

    in_x0_sign = ssGetInputPortRealSignalPtrs(S, 0);
    in_y_ref_sign = ssGetInputPortRealSignalPtrs(S, 1);
    in_y_ref_e_sign = ssGetInputPortRealSignalPtrs(S, 2);
    {% if ocp.dims.np > 0 %}
    in_p_sign = ssGetInputPortRealSignalPtrs(S, 3);
    {% endif %}

    // copy signals into local buffers
    for (int i = 0; i < {{ ocp.dims.nx }}; i++) in_x0[i] = (double)(*in_x0_sign[i]);
    for (int i = 0; i < {{ ocp.dims.ny }}; i++) in_y_ref[i] = (double)(*in_y_ref_sign[i]);
    for (int i = 0; i < {{ ocp.dims.ny_e }}; i++) in_y_ref_e[i] = (double)(*in_y_ref_e_sign[i]);
    {% if ocp.dims.np > 0 %}
    for (int i = 0; i < {{ ocp.dims.np }}; i++) in_p[i] = (double)(*in_p_sign[i]);
    {% endif %}

    // for (int i = 0; i < 4; i++) ssPrintf("x0[%d] = %f\n", i, in_x0[i]);
    // ssPrintf("\n");

    // set initial condition
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "lbx", in_x0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "ubx", in_x0);

    // update reference
    for (int ii = 0; ii < {{ocp.dims.N}}; ii++)
        ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, ii, "yref", (void *) in_y_ref);

    ocp_nlp_cost_model_set(nlp_config, nlp_dims, nlp_in, {{ocp.dims.N}}, "yref", (void *) in_y_ref_e);

    // update value of parameters
    {% if ocp.dims.np > 0%}
    {% if ocp.solver_config.integrator_type == 'IRK' %}
    for (int ii = 0; ii < {{ocp.dims.N}}; ii++) {
    impl_dae_fun[ii].set_param(impl_dae_fun+ii, in_p);
    impl_dae_fun_jac_x_xdot_z[ii].set_param(impl_dae_fun_jac_x_xdot_z+ii, in_p);
    impl_dae_jac_x_xdot_u_z[ii].set_param(impl_dae_jac_x_xdot_u_z+ii, in_p);
    }
    {% else %}
    for (int ii = 0; ii < {{ocp.dims.N}}; ii++) {
    expl_vde_for[ii].set_param(expl_vde_for+ii, in_p);
    }
    {% endif %}
    {% endif %}
    
    // assign pointers to output signals 
    real_t *out_u0, *out_status, *out_KKT_res, *out_x1, *out_cpu_time;

    out_u0          = ssGetOutputPortRealSignal(S, 0);
    out_status      = ssGetOutputPortRealSignal(S, 1);
    out_KKT_res     = ssGetOutputPortRealSignal(S, 2);
    out_x1          = ssGetOutputPortRealSignal(S, 3);
    out_cpu_time    = ssGetOutputPortRealSignal(S, 4);
    
    // call acados_solve()
    int acados_status = acados_solve();

    *out_status = (real_t) acados_status;
    *out_KKT_res = (real_t) nlp_out->inf_norm_res;
    *out_cpu_time = (real_t) nlp_out->total_time;
    
    // get solution
    ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, 0, "u", (void *) out_u0);

    // get next state
    ocp_nlp_out_get(nlp_config, nlp_dims, nlp_out, 1, "x", (void *) out_x1);

	/************************************************
	* extract and print solution
	************************************************/

	// u

	int nu_max = nu[0];
	for(ii=1; ii<=N; ii++)
		if(nu[ii]>nu_max)
			nu_max = nu[ii];

	double *u = malloc(nu_max*sizeof(double));

	printf("\nu = \n");
	for(ii=0; ii<=N; ii++)
		{
		d_cvt_ocp_qp_sol_to_colmaj_u(ii, &qp_sol, u);
		d_print_mat(1, nu[ii], u, 1);
		}

	// x

	int nx_max = nx[0];
	for(ii=1; ii<=N; ii++)
		if(nx[ii]>nx_max)
			nx_max = nx[ii];

	double *x = malloc(nx_max*sizeof(double));

	printf("\nx = \n");
	for(ii=0; ii<=N; ii++)
		{
		d_cvt_ocp_qp_sol_to_colmaj_x(ii, &qp_sol, x);
		d_print_mat(1, nx[ii], x, 1);
		}

}

static void mdlTerminate(SimStruct *S)
{
    acados_free();
}


#ifdef  MATLAB_MEX_FILE
#include "simulink.c"
#else
#include "cg_sfun.h"
#endif
