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
	size += sizeof(struct OCP_QP_SOLVER_ARG);
	return size;
	}



hpipm_size_t OCP_QP_SOLVER_ARG_MEMSIZE(struct d_ocp_qp_dim *ocp_dim)
	{
	hpipm_size_t size = 0;

	size += sizeof(struct OCP_QP_IPM_ARG);
	size += 1*OCP_QP_IPM_ARG_MEMSIZE(NULL); // XXX dim is not used in there

	size += sizeof(struct OCP_QP_REDUCE_EQ_DOF_ARG);
	size += 1*OCP_QP_REDUCE_EQ_DOF_ARG_MEMSIZE();

	size += sizeof(struct OCP_QP_DIM);
	size += 1*OCP_QP_DIM_MEMSIZE(ocp_dim->N);

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
	// red arg struct
	struct OCP_QP_REDUCE_EQ_DOF_ARG *red_arg_ptr = (struct OCP_QP_REDUCE_EQ_DOF_ARG *) ipm_arg_ptr;
	arg->red_arg = red_arg_ptr;
	red_arg_ptr += 1;
	// dim struct
	struct OCP_QP_DIM *red_dim_ptr = (struct OCP_QP_DIM *) red_arg_ptr;
	arg->red_dim = red_dim_ptr;
	red_dim_ptr += 1;

	// align to typical cache line size
	hpipm_size_t s_ptr = (hpipm_size_t) red_dim_ptr;
	s_ptr = (s_ptr+63)/64*64;

	// void stuff
	char *c_ptr = (char *) s_ptr;

	//
	OCP_QP_IPM_ARG_CREATE(NULL, arg->ipm_arg, c_ptr); // XXX dim is not used in there
	c_ptr += arg->ipm_arg->memsize;
	//
	OCP_QP_REDUCE_EQ_DOF_ARG_CREATE(arg->red_arg, c_ptr);
	c_ptr += arg->red_arg->memsize;
	//
	OCP_QP_DIM_CREATE(ocp_dim->N, arg->red_dim, c_ptr);
	c_ptr += arg->red_dim->memsize;

	arg->valid_red_dim = 0;

	arg->memsize = memsize;

	return;
	}



void OCP_QP_SOLVER_ARG_SET_DEFAULT(enum hpipm_mode mode, struct OCP_QP_DIM *ocp_dim, struct OCP_QP_SOLVER_ARG *arg)
	{
	OCP_QP_IPM_ARG_SET_DEFAULT(mode, arg->ipm_arg);
	arg->reduce_eq_dof = 0;
	OCP_QP_REDUCE_EQ_DOF_ARG_SET_DEFAULT(arg->red_arg);
	int alias_unchanged = 1;
	OCP_QP_REDUCE_EQ_DOF_ARG_SET_ALIAS_UNCHANGED(&alias_unchanged, arg->red_arg);
	return;
	}



void OCP_QP_SOLVER_ARG_SET(char *field, void *value, struct OCP_QP_SOLVER_ARG *arg)
	{
	if(hpipm_strcmp(field, "reduce_eq_dof"))
		{
		OCP_QP_SOLVER_ARG_SET_REDUCE_EQ_DOF(value, arg);
		}
	else
		{
		OCP_QP_IPM_ARG_SET(field, value, arg->ipm_arg);
		}
	// TODO common setters for common stuff like dual solution computation
	return;
	}



void OCP_QP_SOLVER_ARG_SET_REDUCE_EQ_DOF(int *value, struct OCP_QP_SOLVER_ARG *arg)
	{
	arg->reduce_eq_dof = *value;
	return;
	}



void OCP_QP_SOLVER_ARG_DEEPCOPY(struct OCP_QP_SOLVER_ARG *arg_s, struct OCP_QP_SOLVER_ARG *arg_d)
	{
	OCP_QP_IPM_ARG_DEEPCOPY(arg_s->ipm_arg, arg_d->ipm_arg);
	OCP_QP_REDUCE_EQ_DOF_ARG_DEEPCOPY(arg_s->red_arg, arg_d->red_arg);
	if(arg_s->valid_red_dim)
		{
		OCP_QP_DIM_DEEPCOPY(arg_s->red_dim, arg_d->red_dim);
		}
	arg_d->reduce_eq_dof = arg_s->reduce_eq_dof;
	return;
	}



void OCP_QP_SOLVER_ARG_GET_REDUCE_EQ_DOF(struct OCP_QP_SOLVER_ARG *arg, int *value)
	{
	*value = arg->reduce_eq_dof;
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

	// compute red_dim
	if(arg->reduce_eq_dof & arg->valid_red_dim==0)
		{
		OCP_QP_DIM_REDUCE_EQ_DOF(ocp_dim, arg->red_dim);
		arg->valid_red_dim = 1;
		}

	size += sizeof(struct OCP_QP_IPM_WS);
	if(arg->reduce_eq_dof)
		{
		size += 1*OCP_QP_IPM_WS_MEMSIZE(arg->red_dim, arg->ipm_arg);
		}
	else
		{
		size += 1*OCP_QP_IPM_WS_MEMSIZE(ocp_dim, arg->ipm_arg);
		}

	size += sizeof(struct OCP_QP_SOLVER_ARG);
	size += 1*OCP_QP_SOLVER_ARG_MEMSIZE(ocp_dim);

	if(arg->reduce_eq_dof)
		{
		size += sizeof(struct OCP_QP_REDUCE_EQ_DOF_WS);
		size += 1*OCP_QP_REDUCE_EQ_DOF_WS_MEMSIZE(ocp_dim);

		size += sizeof(struct OCP_QP);
		size += 1*OCP_QP_MEMSIZE(arg->red_dim);

		size += 2*sizeof(struct OCP_QP_SOL);
		size += 2*OCP_QP_SOL_MEMSIZE(arg->red_dim);

		size += sizeof(struct OCP_QP_SEED);
		size += 1*OCP_QP_SEED_MEMSIZE(arg->red_dim);
		}

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

	hpipm_size_t s_ptr;

	// ipm ws struct
	struct OCP_QP_IPM_WS *ipm_ws_ptr = mem;
	ws->ipm_ws = ipm_ws_ptr;
	ipm_ws_ptr += 1;

	// arg struct
	struct OCP_QP_SOLVER_ARG *arg_ptr = (struct OCP_QP_SOLVER_ARG *) ipm_ws_ptr;
	ws->arg = arg_ptr;
	arg_ptr += 1;

	if(arg->reduce_eq_dof)
		{
		// red ws struct
		struct OCP_QP_REDUCE_EQ_DOF_WS *red_ws_ptr = (struct OCP_QP_REDUCE_EQ_DOF_WS *) arg_ptr;
		ws->red_ws = red_ws_ptr;
		red_ws_ptr += 1;

		// qp struct
		struct OCP_QP *qp_ptr = (struct OCP_QP *) red_ws_ptr;
		ws->red_qp = qp_ptr;
		qp_ptr += 1;

		// sol struct
		struct OCP_QP_SOL *sol_ptr = (struct OCP_QP_SOL *) qp_ptr;
		ws->red_sol = sol_ptr;
		sol_ptr += 1;
		ws->red_sens = sol_ptr;
		sol_ptr += 1;

		// seed struct
		struct OCP_QP_SEED *seed_ptr = (struct OCP_QP_SEED *) sol_ptr;
		ws->red_seed = seed_ptr;
		seed_ptr += 1;

		s_ptr = (hpipm_size_t) sol_ptr;
		}
	else
		{
		s_ptr = (hpipm_size_t) arg_ptr;
		}

	// align to typical cache line size
	//hpipm_size_t s_ptr = (hpipm_size_t) sol_ptr;
	s_ptr = (s_ptr+63)/64*64;

	// void stuff
	char *c_ptr = (char *) s_ptr;

	if(arg->reduce_eq_dof)
		{
		OCP_QP_IPM_WS_CREATE(arg->red_dim, arg->ipm_arg, ws->ipm_ws, c_ptr);
		}
	else
		{
		OCP_QP_IPM_WS_CREATE(ocp_dim, arg->ipm_arg, ws->ipm_ws, c_ptr);
		}
	c_ptr += ws->ipm_ws->memsize;

	OCP_QP_SOLVER_ARG_CREATE(ocp_dim, ws->arg, c_ptr);
	c_ptr += ws->arg->memsize;
	// deep copy of arg
	OCP_QP_SOLVER_ARG_DEEPCOPY(arg, ws->arg);

	if(arg->reduce_eq_dof)
		{
		OCP_QP_REDUCE_EQ_DOF_WS_CREATE(ocp_dim, ws->red_ws, c_ptr);
		c_ptr += ws->red_ws->memsize;

		OCP_QP_CREATE(arg->red_dim, ws->red_qp, c_ptr);
		c_ptr += ws->red_qp->memsize;

		OCP_QP_SOL_CREATE(arg->red_dim, ws->red_sol, c_ptr);
		c_ptr += ws->red_sol->memsize;

		OCP_QP_SOL_CREATE(arg->red_dim, ws->red_sens, c_ptr);
		c_ptr += ws->red_sens->memsize;

		OCP_QP_SEED_CREATE(arg->red_dim, ws->red_seed, c_ptr);
		c_ptr += ws->red_seed->memsize;
		}

	ws->memsize = memsize;

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
	else if(hpipm_strcmp(field, "pred_corr"))
		{
		OCP_QP_SOLVER_SET_PRED_CORR(value, ws);
		}
	else if(hpipm_strcmp(field, "split_step"))
		{
		OCP_QP_SOLVER_SET_SPLIT_STEP(value, ws);
		}
	else if(hpipm_strcmp(field, "reg_prim"))
		{
		OCP_QP_SOLVER_SET_REG_PRIM(value, ws);
		}
	else
		{
#ifdef EXT_DEP
		printf("error: OCP_QP_SOLVER_ARG_SET: wrong field %s\n", field);
#endif
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



void OCP_QP_SOLVER_SET_PRED_CORR(int *value, struct OCP_QP_SOLVER_WS *ws)
	{
	OCP_QP_IPM_ARG_SET_PRED_CORR(value, ws->arg->ipm_arg);
	return;
	}



void OCP_QP_SOLVER_SET_SPLIT_STEP(int *value, struct OCP_QP_SOLVER_WS *ws)
	{
	OCP_QP_IPM_ARG_SET_SPLIT_STEP(value, ws->arg->ipm_arg);
	return;
	}



void OCP_QP_SOLVER_SET_REG_PRIM(REAL *value, struct OCP_QP_SOLVER_WS *ws)
	{
	OCP_QP_IPM_ARG_SET_REG_PRIM(value, ws->arg->ipm_arg);
	return;
	}



void OCP_QP_SOLVER_GET(char *field, struct OCP_QP_SOLVER_WS *ws, void *value)
	{
	if(hpipm_strcmp(field, "status"))
		{
		OCP_QP_SOLVER_GET_STATUS(ws, value);
		}
	else if(hpipm_strcmp(field, "iter"))
		{
		OCP_QP_SOLVER_GET_ITER(ws, value);
		}
	else if(hpipm_strcmp(field, "stat_m"))
		{
		OCP_QP_SOLVER_GET_STAT_M(ws, value);
		}
	else if(hpipm_strcmp(field, "stat"))
		{
		OCP_QP_SOLVER_GET_STAT(ws, value);
		}
	else if(hpipm_strcmp(field, "reduce_eq_dof"))
		{
		OCP_QP_SOLVER_GET_REDUCE_EQ_DOF(ws, value);
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


void OCP_QP_SOLVER_GET_ITER(struct OCP_QP_SOLVER_WS *ws, int *value)
	{
	OCP_QP_IPM_GET_ITER(ws->ipm_ws, value);
	return;
	}


void OCP_QP_SOLVER_GET_STAT_M(struct OCP_QP_SOLVER_WS *ws, int *value)
	{
	OCP_QP_IPM_GET_STAT_M(ws->ipm_ws, value);
	return;
	}


void OCP_QP_SOLVER_GET_STAT(struct OCP_QP_SOLVER_WS *ws, REAL **value)
	{
	OCP_QP_IPM_GET_STAT(ws->ipm_ws, value);
	return;
	}


void OCP_QP_SOLVER_GET_REDUCE_EQ_DOF(struct OCP_QP_SOLVER_WS *ws, int *value)
	{
	OCP_QP_SOLVER_ARG_GET_REDUCE_EQ_DOF(ws->arg, value);
	return;
	}


// XXX no arg
void OCP_QP_SOLVER_SOLVE(struct OCP_QP *qp, struct OCP_QP_SOL *qp_sol, struct OCP_QP_SOLVER_WS *ws)
	{
	if(ws->arg->reduce_eq_dof)
		{
		OCP_QP_REDUCE_EQ_DOF(qp, ws->red_qp, ws->arg->red_arg, ws->red_ws);
		OCP_QP_IPM_SOLVE(ws->red_qp, ws->red_sol, ws->arg->ipm_arg, ws->ipm_ws);
		OCP_QP_RESTORE_EQ_DOF(qp, ws->red_sol, qp_sol, ws->arg->red_arg, ws->red_ws);
		}
	else
		{
		OCP_QP_IPM_SOLVE(qp, qp_sol, ws->arg->ipm_arg, ws->ipm_ws);
		}
	return;
	}


// XXX no arg
void OCP_QP_SOLVER_SENS_FRW(struct OCP_QP *qp, struct OCP_QP_SEED *qp_seed, struct OCP_QP_SOL *qp_sens, struct OCP_QP_SOLVER_WS *ws)
	{
	if(ws->arg->reduce_eq_dof)
		{
		OCP_QP_REDUCE_EQ_DOF_SEED(qp, qp_seed, ws->red_seed, ws->arg->red_arg, ws->red_ws); // red_seed
		OCP_QP_IPM_SENS_FRW(ws->red_qp, ws->red_seed, ws->red_sens, ws->arg->ipm_arg, ws->ipm_ws); // red_sens
		OCP_QP_RESTORE_EQ_DOF_SEED(qp, qp_seed, ws->red_sens, qp_sens, ws->arg->red_arg, ws->red_ws);
		}
	else
		{
		OCP_QP_IPM_SENS_FRW(qp, qp_seed, qp_sens, ws->arg->ipm_arg, ws->ipm_ws);
		}
	return;
	}


