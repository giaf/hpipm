/**************************************************************************************************
*                                                                                                 *
* This file is part of HPIPM.                                                                     *
*                                                                                                 *
* HPIPM -- High-Performance Interior Point Method.                                                *
* Copyright (C) 2019 by Gianluca Frison.                                                          *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
*                                                                                                 *
* The 2-Clause BSD License                                                                        *
*                                                                                                 *
* Redistribution and use in source and binary forms, with or without                              *
* modification, are permitted provided that the following conditions are met:                     *
*                                                                                                 *
* 1. Redistributions of source code must retain the above copyright notice, this                  *
*    list of conditions and the following disclaimer.                                             *
* 2. Redistributions in binary form must reproduce the above copyright notice,                    *
*    this list of conditions and the following disclaimer in the documentation                    *
*    and/or other materials provided with the distribution.                                       *
*                                                                                                 *
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND                 *
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED                   *
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE                          *
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR                 *
* ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES                  *
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;                    *
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND                     *
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT                      *
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS                   *
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                                    *
*                                                                                                 *
* Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             *
*                                                                                                 *
**************************************************************************************************/






hpipm_size_t OCP_QP_SOLVER_ARG_STRSIZE()
	{
	hpipm_size_t size = 0;
	size += sizeof(struct d_ocp_qp_solver_arg);
	return size;
	}



hpipm_size_t OCP_QP_SOLVER_ARG_MEMSIZE(struct d_ocp_qp_dim *ocp_dim)
	{
	hpipm_size_t size = 0;

	size += sizeof(struct OCP_QP_IPM_ARG);
	size += 1*OCP_QP_IPM_ARG_MEMSIZE(ocp_dim);

	size += 1*64; // align once to typical cache line size
	size = (size+63)/64*64; // make multiple of typical cache line size

	return size;
	}



void OCP_QP_SOLVER_ARG_CREATE(struct OCP_QP_DIM *ocp_dim, struct OCP_QP_SOLVER_ARG *arg, void *mem)
	{
	// zero memory (to avoid corrupted memory like e.g. NaN)
	hpipm_size_t memsize = OCP_QP_SOLVER_ARG_MEMSIZE(ocp_dim);
	hpipm_zero_memset(memsize, mem);

	// ipm arg struct
	struct OCP_QP_IPM_ARG *ipm_arg_ptr = mem;
	arg->ipm_arg = ipm_arg_ptr;
	ipm_arg_ptr += 1;

	// align to typical cache line size
	hpipm_size_t s_ptr = (hpipm_size_t) ipm_arg_ptr;
	s_ptr = (s_ptr+63)/64*64;

	// void stuff
	char *c_ptr = (char *) s_ptr;

	OCP_QP_IPM_ARG_CREATE(ocp_dim, arg->ipm_arg, c_ptr);
	c_ptr += arg->ipm_arg->memsize;

	arg->memsize = memsize;

	return;
	}



void OCP_QP_SOLVER_ARG_SET_DEFAULT(enum hpipm_mode mode, struct OCP_QP_DIM *ocp_dim, struct OCP_QP_SOLVER_ARG *arg)
	{
	OCP_QP_IPM_ARG_SET_DEFAULT(mode, arg->ipm_arg);
	return;
	}



void OCP_QP_SOLVER_ARG_SET(char *field, void *value, struct OCP_QP_SOLVER_ARG *arg)
	{
	OCP_QP_IPM_ARG_SET(field, value, arg->ipm_arg);
	return;
	}



void OCP_QP_SOLVER_ARG_DEEPCOPY(struct OCP_QP_SOLVER_ARG *arg_s, struct OCP_QP_SOLVER_ARG *arg_d)
	{
	OCP_QP_IPM_ARG_DEEPCOPY(arg_s->ipm_arg, arg_d->ipm_arg);
	return;
	}



hpipm_size_t OCP_QP_SOLVER_WS_STRSIZE()
	{
	hpipm_size_t size = 0;
	size += sizeof(struct OCP_QP_SOLVER_WS);
	return size;
	}



hpipm_size_t OCP_QP_SOLVER_WS_MEMSIZE(struct OCP_QP_DIM *ocp_dim, struct OCP_QP_SOLVER_ARG *arg)
	{
	hpipm_size_t size = 0;

	size += sizeof(struct OCP_QP_IPM_WS);
	size += 1*OCP_QP_IPM_WS_MEMSIZE(ocp_dim, arg->ipm_arg);

	size += sizeof(struct OCP_QP_SOLVER_ARG);
	size += 1*OCP_QP_SOLVER_ARG_MEMSIZE(ocp_dim);

	size += 1*64; // align once to typical cache line size
	size = (size+63)/64*64; // make multiple of typical cache line size

	return size;
	}



void OCP_QP_SOLVER_WS_CREATE(struct OCP_QP_DIM *ocp_dim, struct OCP_QP_SOLVER_ARG *arg, struct OCP_QP_SOLVER_WS *ws, void *mem)
	{
	// zero memory (to avoid corrupted memory like e.g. NaN)
	hpipm_size_t memsize = OCP_QP_SOLVER_WS_MEMSIZE(ocp_dim, arg);
	hpipm_zero_memset(memsize, mem);

	// XXX in nested create routines, multiple zeroing out of same memory !!!

	// arg struct
	struct OCP_QP_SOLVER_ARG *arg_ptr = mem;
	ws->arg = arg_ptr;
	arg_ptr += 1;

	// ipm ws struct
	struct OCP_QP_IPM_WS *ipm_ws_ptr = (struct OCP_QP_IPM_WS *) arg_ptr;
	ws->ipm_ws = ipm_ws_ptr;
	ipm_ws_ptr += 1;

	// align to typical cache line size
	hpipm_size_t s_ptr = (hpipm_size_t) ipm_ws_ptr;
	s_ptr = (s_ptr+63)/64*64;

	// void stuff
	char *c_ptr = (char *) s_ptr;

	OCP_QP_IPM_WS_CREATE(ocp_dim, arg->ipm_arg, ws->ipm_ws, c_ptr);
	c_ptr += ws->ipm_ws->memsize;

	OCP_QP_SOLVER_ARG_CREATE(ocp_dim, ws->arg, c_ptr);
	c_ptr += ws->arg->memsize;
	// deep copy of arg
	OCP_QP_SOLVER_ARG_DEEPCOPY(arg, ws->arg);

	arg->memsize = memsize;

	return;
	}



void OCP_QP_SOLVER_GET(char *field, struct OCP_QP_SOLVER_WS *ws, void *value)
	{
	if(hpipm_strcmp(field, "status"))
		{
		OCP_QP_SOLVER_GET_STATUS(ws, value);
		}
	else
		{
		OCP_QP_IPM_GET(field, ws->ipm_ws, value);
		}
	return;
	}


void OCP_QP_SOLVER_GET_STATUS(struct OCP_QP_SOLVER_WS *ws, int *value)
	{
	OCP_QP_IPM_GET_STATUS(ws->ipm_ws, value);
	return;
	}


// XXX set selected args that can safely change
void OCP_QP_SOLVER_SET(char *field, void *value, struct OCP_QP_SOLVER_WS *ws)
	{
	// TODO
	if(hpipm_strcmp(field, "iter_max"))
		{
		OCP_QP_SOLVER_SET_ITER_MAX(value, ws);
		}
	else if(hpipm_strcmp(field, "alpha_min"))
		{
		OCP_QP_SOLVER_SET_ALPHA_MIN(value, ws);
		}
	else if(hpipm_strcmp(field, "mu0"))
		{
		OCP_QP_SOLVER_SET_MU0(value, ws);
		}
	else if(hpipm_strcmp(field, "tol_stat"))
		{
		OCP_QP_SOLVER_SET_TOL_STAT(value, ws);
		}
	else if(hpipm_strcmp(field, "tol_eq"))
		{
		OCP_QP_SOLVER_SET_TOL_EQ(value, ws);
		}
	else if(hpipm_strcmp(field, "tol_ineq"))
		{
		OCP_QP_SOLVER_SET_TOL_INEQ(value, ws);
		}
	else if(hpipm_strcmp(field, "tol_comp"))
		{
		OCP_QP_SOLVER_SET_TOL_COMP(value, ws);
		}
	else
		{
		printf("error: OCP_QP_SOLVER_ARG_SET: wrong field %s\n", field);
		exit(1);	
		}
	return;
	}



void OCP_QP_SOLVER_SET_ITER_MAX(int *value, struct OCP_QP_SOLVER_WS *ws)
	{
	OCP_QP_IPM_ARG_SET_ITER_MAX(value, ws->arg->ipm_arg);
	return;
	}



void OCP_QP_SOLVER_SET_ALPHA_MIN(REAL *value, struct OCP_QP_SOLVER_WS *ws)
	{
	OCP_QP_IPM_ARG_SET_ALPHA_MIN(value, ws->arg->ipm_arg);
	return;
	}



void OCP_QP_SOLVER_SET_MU0(REAL *value, struct OCP_QP_SOLVER_WS *ws)
	{
	OCP_QP_IPM_ARG_SET_MU0(value, ws->arg->ipm_arg);
	return;
	}



void OCP_QP_SOLVER_SET_TOL_STAT(REAL *value, struct OCP_QP_SOLVER_WS *ws)
	{
	OCP_QP_IPM_ARG_SET_TOL_STAT(value, ws->arg->ipm_arg);
	return;
	}



void OCP_QP_SOLVER_SET_TOL_EQ(REAL *value, struct OCP_QP_SOLVER_WS *ws)
	{
	OCP_QP_IPM_ARG_SET_TOL_EQ(value, ws->arg->ipm_arg);
	return;
	}



void OCP_QP_SOLVER_SET_TOL_INEQ(REAL *value, struct OCP_QP_SOLVER_WS *ws)
	{
	OCP_QP_IPM_ARG_SET_TOL_INEQ(value, ws->arg->ipm_arg);
	return;
	}



void OCP_QP_SOLVER_SET_TOL_COMP(REAL *value, struct OCP_QP_SOLVER_WS *ws)
	{
	OCP_QP_IPM_ARG_SET_TOL_COMP(value, ws->arg->ipm_arg);
	return;
	}



// XXX no arg
void OCP_QP_SOLVER_SOLVE(struct OCP_QP *qp, struct OCP_QP_SOL *qp_sol, struct OCP_QP_SOLVER_WS *ws)
	{
	OCP_QP_IPM_SOLVE(qp, qp_sol, ws->arg->ipm_arg, ws->ipm_ws);
	return;
	}


