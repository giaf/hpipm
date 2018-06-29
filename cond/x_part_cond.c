/**************************************************************************************************
*                                                                                                 *
* This file is part of HPIPM.                                                                     *
*                                                                                                 *
* HPIPM -- High Performance Interior Point Method.                                                *
* Copyright (C) 2017 by Gianluca Frison.                                                          *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
*                                                                                                 *
* HPIPM is free software; you can redistribute it and/or                                          *
* modify it under the terms of the GNU Lesser General Public                                      *
* License as published by the Free Software Foundation; either                                    *
* version 2.1 of the License, or (at your option) any later version.                              *
*                                                                                                 *
* HPIPM is distributed in the hope that it will be useful,                                        *
* but WITHOUT ANY WARRANTY; without even the implied warranty of                                  *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                                            *
* See the GNU Lesser General Public License for more details.                                     *
*                                                                                                 *
* You should have received a copy of the GNU Lesser General Public                                *
* License along with HPIPM; if not, write to the Free Software                                    *
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA                  *
*                                                                                                 *
* Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             *
*                                                                                                 *
**************************************************************************************************/

void COMPUTE_BLOCK_SIZE_COND_QP_OCP2OCP(int N, int N2, int *block_size)
	{

	int ii;

	int bs0 = N/N2; // (floor) size of small blocks

	// the first blocks have size bs0+1
	for(ii=0; ii<N-N2*bs0; ii++)
		block_size[ii] = bs0+1;
	// the following blocks have size bs0
	for(; ii<N2; ii++)
		block_size[ii] = bs0;
	// the last block has size 0
	block_size[N2] = 0;

	return;

	}



void COMPUTE_QP_DIM_OCP2OCP(struct OCP_QP_DIM *ocp_dim, int *block_size, struct OCP_QP_DIM *part_dense_dim)
	{

	// TODO run time check on sum(block_size) = N

	int N = ocp_dim->N;
	int *nx = ocp_dim->nx;
	int *nu = ocp_dim->nu;
	int *nb = ocp_dim->nb;
	int *nbx = ocp_dim->nbx;
	int *nbu = ocp_dim->nbu;
	int *ng = ocp_dim->ng;
	int *ns = ocp_dim->ns;
	int *nsbx = ocp_dim->nsbx;
	int *nsbu = ocp_dim->nsbu;
	int *nsg = ocp_dim->nsg;

	int N2 = part_dense_dim->N;
	int *nx2 = part_dense_dim->nx;
	int *nu2 = part_dense_dim->nu;
	int *nb2 = part_dense_dim->nb;
	int *nbx2 = part_dense_dim->nbx;
	int *nbu2 = part_dense_dim->nbu;
	int *ng2 = part_dense_dim->ng;
	int *ns2 = part_dense_dim->ns;
	int *nsbx2 = part_dense_dim->nsbx;
	int *nsbu2 = part_dense_dim->nsbu;
	int *nsg2 = part_dense_dim->nsg;

	int ii, jj;

	int nbb; // box constr that remain box constr
	int nbg; // box constr that becomes general constr
	int N_tmp = 0; // temporary sum of block size
	// first stages
	for(ii=0; ii<N2; ii++)
		{
		nx2[ii] = nx[N_tmp+0];
		nu2[ii] = nu[N_tmp+0];
		nbx2[ii] = nbx[N_tmp+0];
		nbu2[ii] = nbu[N_tmp+0];
		nb2[ii] = nb[N_tmp+0];
		ng2[ii] = ng[N_tmp+0];
		ns2[ii] = ns[N_tmp+0];
		nsbx2[ii] = nsbx[N_tmp+0];
		nsbu2[ii] = nsbu[N_tmp+0];
		nsg2[ii] = nsg[N_tmp+0];
		for(jj=1; jj<block_size[ii]; jj++)
			{
			nx2[ii] += 0;
			nu2[ii] += nu[N_tmp+jj];
			nbx2[ii] += 0;
			nbu2[ii] += nbu[N_tmp+jj];
			nb2[ii] += nbu[N_tmp+jj];
			ng2[ii] += ng[N_tmp+jj] + nbx[N_tmp+jj];
			ns2[ii] += ns[N_tmp+jj];
			nsbx2[ii] += 0;
			nsbu2[ii] += nsbu[N_tmp+jj];
			nsg2[ii] += nsg[N_tmp+jj] + nsbx[N_tmp+jj];
			}
		N_tmp += block_size[ii];
		}
	// last stage: condense also following stage
	ii = N2;
	nx2[ii] = nx[N_tmp+0];
	nu2[ii] = nu[N_tmp+0];
	nbx2[ii] = nbx[N_tmp+0];
	nbu2[ii] = nbu[N_tmp+0];
	nb2[ii] = nb[N_tmp+0];
	ng2[ii] = ng[N_tmp+0];
	ns2[ii] = ns[N_tmp+0];
	nsbx2[ii] = nsbx[N_tmp+0];
	nsbu2[ii] = nsbu[N_tmp+0];
	nsg2[ii] = nsg[N_tmp+0];
	for(jj=1; jj<block_size[ii]+1; jj++)
		{
		nx2[ii] += 0;
		nu2[ii] += nu[N_tmp+jj];
		nbx2[ii] += 0;
		nbu2[ii] += nbu[N_tmp+jj];
		nb2[ii] += nbu[N_tmp+jj];
		ng2[ii] += ng[N_tmp+jj] + nbx[N_tmp+jj];
		ns2[ii] += ns[N_tmp+jj];

		nsbx2[ii] = nsbx[N_tmp+0];
		nsbu2[ii] = nsbu[N_tmp+0];
		nsg2[ii] += nsg[N_tmp+jj] + nsbx[N_tmp+jj];
		}

	return;

	}



int MEMSIZE_COND_QP_OCP2OCP_ARG(int N2)
	{

	int ii;

	int size = 0;

	size += (N2+1)*sizeof(struct COND_QP_OCP2DENSE_ARG);

	for(ii=0; ii<=N2; ii++)
		{

		size += MEMSIZE_COND_QP_OCP2DENSE_ARG();

		}

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

	}



void CREATE_COND_QP_OCP2OCP_ARG(int N2, struct COND_QP_OCP2OCP_ARG *part_cond_arg, void *mem)
	{

	int ii;

	// cond workspace struct
	struct COND_QP_OCP2DENSE_ARG *cws_ptr = mem;
	part_cond_arg->cond_arg = cws_ptr;
	cws_ptr += N2+1;

	// align to typicl cache line size
	size_t s_ptr = (size_t) cws_ptr;
	s_ptr = (s_ptr+63)/64*64;

	char *c_ptr = (char *) s_ptr;

	for(ii=0; ii<=N2; ii++)
		{

		CREATE_COND_QP_OCP2DENSE_ARG(part_cond_arg->cond_arg+ii, c_ptr);
		c_ptr += (part_cond_arg->cond_arg+ii)->memsize;

		}

	part_cond_arg->memsize = MEMSIZE_COND_QP_OCP2OCP_ARG(N2);

#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + part_cond_arg->memsize)
		{
		printf("\nCreate_cond_qp_ocp2ocp_arg: outside memory bounds!\n\n");
		exit(1);
		}
#endif

return;

	}



void SET_DEFAULT_COND_QP_OCP2OCP_ARG(int N2, struct COND_QP_OCP2OCP_ARG *part_cond_arg)
	{

	int ii;

	for(ii=0; ii<=N2; ii++)
		{

		SET_DEFAULT_COND_QP_OCP2DENSE_ARG(part_cond_arg->cond_arg+ii);
		(part_cond_arg->cond_arg+ii)->cond_last_stage = 0;

		}
	// cond_last_stage at last stage
	part_cond_arg->cond_arg[N2].cond_last_stage = 1;

return;

	}



int MEMSIZE_COND_QP_OCP2OCP(struct OCP_QP_DIM *ocp_dim, int *block_size, struct OCP_QP_DIM *part_dense_dim, struct COND_QP_OCP2OCP_ARG *part_cond_arg)
	{

	struct OCP_QP_DIM tmp_ocp_dim;

	int ii;

	int N = ocp_dim->N;
	int N2 = part_dense_dim->N;

	int size = 0;

	size += (N2+1)*sizeof(struct COND_QP_OCP2DENSE_WORKSPACE);

	int N_tmp = 0; // temporary sum of horizons
	for(ii=0; ii<=N2; ii++)
		{

		// alias ocp_dim
		tmp_ocp_dim.N = block_size[ii];
		tmp_ocp_dim.nx = ocp_dim->nx+N_tmp;
		tmp_ocp_dim.nu = ocp_dim->nu+N_tmp;
		tmp_ocp_dim.nbx = ocp_dim->nbx+N_tmp;
		tmp_ocp_dim.nbu = ocp_dim->nbu+N_tmp;
		tmp_ocp_dim.nb = ocp_dim->nb+N_tmp;
		tmp_ocp_dim.ng = ocp_dim->ng+N_tmp;

		size += MEMSIZE_COND_QP_OCP2DENSE(&tmp_ocp_dim, part_cond_arg->cond_arg+ii);

		N_tmp += block_size[ii];

		}

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

	}



void CREATE_COND_QP_OCP2OCP(struct OCP_QP_DIM *ocp_dim, int *block_size, struct OCP_QP_DIM *part_dense_dim, struct COND_QP_OCP2OCP_ARG *part_cond_arg, struct COND_QP_OCP2OCP_WORKSPACE *part_cond_ws, void *mem)
	{

	struct OCP_QP_DIM tmp_ocp_dim;

	int ii;

	int N = ocp_dim->N;
	int N2 = part_dense_dim->N;

	// cond workspace struct
	struct COND_QP_OCP2DENSE_WORKSPACE *cws_ptr = mem;
	part_cond_ws->cond_workspace = cws_ptr;
	cws_ptr += N2+1;

	// align to typicl cache line size
	size_t s_ptr = (size_t) cws_ptr;
	s_ptr = (s_ptr+63)/64*64;

	char *c_ptr = (char *) s_ptr;

	int N_tmp = 0; // temporary sum of horizons
	for(ii=0; ii<=N2; ii++)
		{

		// alias ocp_dim
		tmp_ocp_dim.N = block_size[ii];
		tmp_ocp_dim.nx = ocp_dim->nx+N_tmp;
		tmp_ocp_dim.nu = ocp_dim->nu+N_tmp;
		tmp_ocp_dim.nbx = ocp_dim->nbx+N_tmp;
		tmp_ocp_dim.nbu = ocp_dim->nbu+N_tmp;
		tmp_ocp_dim.nb = ocp_dim->nb+N_tmp;
		tmp_ocp_dim.ng = ocp_dim->ng+N_tmp;

		CREATE_COND_QP_OCP2DENSE(&tmp_ocp_dim, part_cond_arg->cond_arg+ii, part_cond_ws->cond_workspace+ii, c_ptr);
		c_ptr += (part_cond_ws->cond_workspace+ii)->memsize;

		N_tmp += block_size[ii];

		}

	part_cond_ws->memsize = MEMSIZE_COND_QP_OCP2OCP(ocp_dim, block_size, part_dense_dim, part_cond_arg);

#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + part_cond_ws->memsize)
		{
		printf("\nCreate_cond_qp_ocp2ocp: outside memory bounds!\n\n");
		exit(1);
		}
#endif

return;

	}



void COND_QP_OCP2OCP(struct OCP_QP *ocp_qp, struct OCP_QP *part_dense_qp, struct COND_QP_OCP2OCP_ARG *part_cond_arg, struct COND_QP_OCP2OCP_WORKSPACE *part_cond_ws)
	{

	struct OCP_QP_DIM tmp_ocp_dim;
	struct OCP_QP tmp_ocp_qp;

	int ii;

	int N = ocp_qp->dim->N;
	int N2 = part_dense_qp->dim->N;
	int bs; // horizon of current block

	int N_tmp = 0; // temporary sum of horizons
	for(ii=0; ii<=N2; ii++)
		{

		bs = part_cond_ws->cond_workspace[ii].bs;

		// alias ocp_dim
		tmp_ocp_dim.N = bs;
		tmp_ocp_dim.nx = ocp_qp->dim->nx+N_tmp;
		tmp_ocp_dim.nu = ocp_qp->dim->nu+N_tmp;
		tmp_ocp_dim.nbx = ocp_qp->dim->nbx+N_tmp;
		tmp_ocp_dim.nbu = ocp_qp->dim->nbu+N_tmp;
		tmp_ocp_dim.nb = ocp_qp->dim->nb+N_tmp;
		tmp_ocp_dim.ng = ocp_qp->dim->ng+N_tmp;
		tmp_ocp_dim.ns = ocp_qp->dim->ns+N_tmp;

		// alias ocp_qp
		tmp_ocp_qp.dim = &tmp_ocp_dim;
		tmp_ocp_qp.idxb = ocp_qp->idxb+N_tmp;
		tmp_ocp_qp.BAbt = ocp_qp->BAbt+N_tmp;
		tmp_ocp_qp.b = ocp_qp->b+N_tmp;
		tmp_ocp_qp.RSQrq = ocp_qp->RSQrq+N_tmp;
		tmp_ocp_qp.rqz = ocp_qp->rqz+N_tmp;
		tmp_ocp_qp.DCt = ocp_qp->DCt+N_tmp;
		tmp_ocp_qp.d = ocp_qp->d+N_tmp;
		tmp_ocp_qp.Z = ocp_qp->Z+N_tmp;
		tmp_ocp_qp.idxs = ocp_qp->idxs+N_tmp;

		COND_BABT(&tmp_ocp_qp, part_dense_qp->BAbt+ii, part_dense_qp->b+ii, part_cond_arg->cond_arg+ii, part_cond_ws->cond_workspace+ii);

		COND_RSQRQ_N2NX3(&tmp_ocp_qp, part_dense_qp->RSQrq+ii, part_dense_qp->rqz+ii, part_cond_arg->cond_arg+ii, part_cond_ws->cond_workspace+ii);

		COND_DCTD(&tmp_ocp_qp, part_dense_qp->idxb[ii], part_dense_qp->DCt+ii, part_dense_qp->d+ii, part_dense_qp->idxs[ii], part_dense_qp->Z+ii, part_dense_qp->rqz+ii, part_cond_arg->cond_arg+ii, part_cond_ws->cond_workspace+ii);

		N_tmp += bs;

		}

#if 0
	// copy last stage
	int *nx = ocp_qp->dim->nx;
	int *nu = ocp_qp->dim->nu;
	int *nb = ocp_qp->dim->nb;
	int *ng = ocp_qp->dim->ng;
	int *ns = ocp_qp->dim->ns;

	GECP_LIBSTR(nu[N]+nx[N]+1, nu[N]+nx[N], ocp_qp->RSQrq+N, 0, 0, part_dense_qp->RSQrq+N2, 0, 0);
	VECCP_LIBSTR(nu[N]+nx[N], ocp_qp->rq+N, 0, part_dense_qp->rq+N2, 0);
	GECP_LIBSTR(nu[N]+nx[N], ng[N], ocp_qp->DCt+N, 0, 0, part_dense_qp->DCt+N2, 0, 0);
	VECCP_LIBSTR(2*nb[N]+2*ng[N], ocp_qp->d+N, 0, part_dense_qp->d+N2, 0);
	for(ii=0; ii<nb[N]; ii++) part_dense_qp->idxb[N2][ii] = ocp_qp->idxb[N][ii];
	VECCP_LIBSTR(2*ns[N], ocp_qp->Z+N, 0, part_dense_qp->Z+N2, 0);
	VECCP_LIBSTR(2*ns[N], ocp_qp->z+N, 0, part_dense_qp->z+N2, 0);
	for(ii=0; ii<ns[N]; ii++) part_dense_qp->idxs[N2][ii] = ocp_qp->idxs[N][ii];
#endif

	return;

	}



void COND_RHS_QP_OCP2OCP(struct OCP_QP *ocp_qp, struct OCP_QP *part_dense_qp, struct COND_QP_OCP2OCP_ARG *part_cond_arg, struct COND_QP_OCP2OCP_WORKSPACE *part_cond_ws)
	{

	struct OCP_QP_DIM tmp_ocp_dim;
	struct OCP_QP tmp_ocp_qp;

	int ii;

	int N = ocp_qp->dim->N;
	int N2 = part_dense_qp->dim->N;
	int bs; // horizon of current block

	int N_tmp = 0; // temporary sum of horizons
	for(ii=0; ii<=N2; ii++)
		{

		bs = part_cond_ws->cond_workspace[ii].bs;

		// alias ocp_dim
		tmp_ocp_dim.N = bs;
		tmp_ocp_dim.nx = ocp_qp->dim->nx+N_tmp;
		tmp_ocp_dim.nu = ocp_qp->dim->nu+N_tmp;
		tmp_ocp_dim.nbx = ocp_qp->dim->nbx+N_tmp;
		tmp_ocp_dim.nbu = ocp_qp->dim->nbu+N_tmp;
		tmp_ocp_dim.nb = ocp_qp->dim->nb+N_tmp;
		tmp_ocp_dim.ng = ocp_qp->dim->ng+N_tmp;
		tmp_ocp_dim.ns = ocp_qp->dim->ns+N_tmp;

		// alias ocp_qp
		tmp_ocp_qp.dim = &tmp_ocp_dim;
		tmp_ocp_qp.idxb = ocp_qp->idxb+N_tmp;
		tmp_ocp_qp.BAbt = ocp_qp->BAbt+N_tmp;
		tmp_ocp_qp.b = ocp_qp->b+N_tmp;
		tmp_ocp_qp.RSQrq = ocp_qp->RSQrq+N_tmp;
		tmp_ocp_qp.rqz = ocp_qp->rqz+N_tmp;
		tmp_ocp_qp.DCt = ocp_qp->DCt+N_tmp;
		tmp_ocp_qp.d = ocp_qp->d+N_tmp;
		tmp_ocp_qp.Z = ocp_qp->Z+N_tmp;
		tmp_ocp_qp.idxs = ocp_qp->idxs+N_tmp;

		COND_B(&tmp_ocp_qp, part_dense_qp->b+ii, part_cond_arg->cond_arg+ii, part_cond_ws->cond_workspace+ii);

		COND_RQ_N2NX3(&tmp_ocp_qp, part_dense_qp->rqz+ii, part_cond_arg->cond_arg+ii, part_cond_ws->cond_workspace+ii);

		COND_D(&tmp_ocp_qp, part_dense_qp->d+ii, part_dense_qp->rqz+ii, part_cond_arg->cond_arg+ii, part_cond_ws->cond_workspace+ii);

		N_tmp += bs;

		}

#if 0
	// copy last stage
	int *nx = ocp_qp->dim->nx;
	int *nu = ocp_qp->dim->nu;
	int *nb = ocp_qp->dim->nb;
	int *ng = ocp_qp->dim->ng;
	int *ns = ocp_qp->dim->ns;

	VECCP_LIBSTR(nu[N]+nx[N], ocp_qp->rq+N, 0, part_dense_qp->rq+N2, 0);
	VECCP_LIBSTR(2*nb[N]+2*ng[N], ocp_qp->d+N, 0, part_dense_qp->d+N2, 0);
	VECCP_LIBSTR(2*ns[N], ocp_qp->z+N, 0, part_dense_qp->z+N2, 0);
#endif

	return;

	}



void EXPAND_SOL_OCP2OCP(struct OCP_QP *ocp_qp, struct OCP_QP *part_dense_qp, struct OCP_QP_SOL *part_dense_qp_sol, struct OCP_QP_SOL *ocp_qp_sol, struct COND_QP_OCP2OCP_ARG *part_cond_arg, struct COND_QP_OCP2OCP_WORKSPACE *part_cond_ws)
	{

	struct OCP_QP_DIM tmp_ocp_dim;
	struct OCP_QP tmp_ocp_qp;
	struct OCP_QP_SOL tmp_ocp_qp_sol;
	struct DENSE_QP_SOL dense_qp_sol;

	int *nx = ocp_qp->dim->nx;
	int *nu = ocp_qp->dim->nu;
	int *nb = ocp_qp->dim->nb;
	int *ng = ocp_qp->dim->ng;
	int *ns = ocp_qp->dim->ns;

	int ii;

	int N = ocp_qp->dim->N;
	int N2 = part_dense_qp->dim->N;
	int bs; // horizon of current block

	int N_tmp = 0; // temporary sum of horizons
	for(ii=0; ii<=N2; ii++)
		{

		bs = part_cond_ws->cond_workspace[ii].bs;

		// alias ocp_dim
		tmp_ocp_dim.N = bs;
		tmp_ocp_dim.nx = ocp_qp->dim->nx+N_tmp;
		tmp_ocp_dim.nu = ocp_qp->dim->nu+N_tmp;
		tmp_ocp_dim.nbx = ocp_qp->dim->nbx+N_tmp;
		tmp_ocp_dim.nbu = ocp_qp->dim->nbu+N_tmp;
		tmp_ocp_dim.nb = ocp_qp->dim->nb+N_tmp;
		tmp_ocp_dim.ng = ocp_qp->dim->ng+N_tmp;
		tmp_ocp_dim.ns = ocp_qp->dim->ns+N_tmp;

		// alias ocp_qp
		tmp_ocp_qp.dim = &tmp_ocp_dim;
		tmp_ocp_qp.idxb = ocp_qp->idxb+N_tmp;
		tmp_ocp_qp.BAbt = ocp_qp->BAbt+N_tmp;
		tmp_ocp_qp.b = ocp_qp->b+N_tmp;
		tmp_ocp_qp.RSQrq = ocp_qp->RSQrq+N_tmp;
		tmp_ocp_qp.rqz = ocp_qp->rqz+N_tmp;
		tmp_ocp_qp.DCt = ocp_qp->DCt+N_tmp;
		tmp_ocp_qp.d = ocp_qp->d+N_tmp;
		tmp_ocp_qp.Z = ocp_qp->Z+N_tmp;
		tmp_ocp_qp.idxs = ocp_qp->idxs+N_tmp;

		// alias ocp qp sol
		tmp_ocp_qp_sol.ux = ocp_qp_sol->ux+N_tmp;
		tmp_ocp_qp_sol.pi = ocp_qp_sol->pi+N_tmp;
		tmp_ocp_qp_sol.lam = ocp_qp_sol->lam+N_tmp;
		tmp_ocp_qp_sol.t = ocp_qp_sol->t+N_tmp;

		// alias ocp qp sol
		dense_qp_sol.v = part_dense_qp_sol->ux+ii;
		dense_qp_sol.pi = part_dense_qp_sol->pi+ii;
		dense_qp_sol.lam = part_dense_qp_sol->lam+ii;
		dense_qp_sol.t = part_dense_qp_sol->t+ii;

		EXPAND_SOL(&tmp_ocp_qp, &dense_qp_sol, &tmp_ocp_qp_sol, part_cond_arg->cond_arg+ii, part_cond_ws->cond_workspace+ii);

		N_tmp += bs;

		}

#if 0
	// copy last stage
	VECCP_LIBSTR(nu[N]+nx[N]+2*ns[N], part_dense_qp_sol->ux+N2, 0, ocp_qp_sol->ux+N, 0);
	VECCP_LIBSTR(2*nb[N]+2*ng[N]+2*ns[N], part_dense_qp_sol->lam+N2, 0, ocp_qp_sol->lam+N, 0);
	VECCP_LIBSTR(2*nb[N]+2*ng[N]+2*ns[N], part_dense_qp_sol->t+N2, 0, ocp_qp_sol->t+N, 0);
#endif

	return;

	}

/************************************************
* update cond
************************************************/

void UPDATE_COND_QP_OCP2OCP(int *idxc, struct OCP_QP *ocp_qp, struct OCP_QP *part_dense_qp, struct COND_QP_OCP2OCP_ARG *part_cond_arg, struct COND_QP_OCP2OCP_WORKSPACE *part_cond_ws)
	{

	struct OCP_QP_DIM tmp_ocp_dim;
	struct OCP_QP tmp_ocp_qp;

	int ii;

	int N = ocp_qp->dim->N;
	int N2 = part_dense_qp->dim->N;
	int bs; // horizon of current block

	int N_tmp = 0; // temporary sum of horizons
	for(ii=0; ii<=N2; ii++)
		{

		bs = part_cond_ws->cond_workspace[ii].bs;

		// alias ocp_dim
		tmp_ocp_dim.N = bs;
		tmp_ocp_dim.nx = ocp_qp->dim->nx+N_tmp;
		tmp_ocp_dim.nu = ocp_qp->dim->nu+N_tmp;
		tmp_ocp_dim.nbx = ocp_qp->dim->nbx+N_tmp;
		tmp_ocp_dim.nbu = ocp_qp->dim->nbu+N_tmp;
		tmp_ocp_dim.nb = ocp_qp->dim->nb+N_tmp;
		tmp_ocp_dim.ng = ocp_qp->dim->ng+N_tmp;
		tmp_ocp_dim.ns = ocp_qp->dim->ns+N_tmp;

		// alias ocp_qp
		tmp_ocp_qp.dim = &tmp_ocp_dim;
		tmp_ocp_qp.idxb = ocp_qp->idxb+N_tmp;
		tmp_ocp_qp.BAbt = ocp_qp->BAbt+N_tmp;
		tmp_ocp_qp.b = ocp_qp->b+N_tmp;
		tmp_ocp_qp.RSQrq = ocp_qp->RSQrq+N_tmp;
		tmp_ocp_qp.rqz = ocp_qp->rqz+N_tmp;
		tmp_ocp_qp.DCt = ocp_qp->DCt+N_tmp;
		tmp_ocp_qp.d = ocp_qp->d+N_tmp;
		tmp_ocp_qp.Z = ocp_qp->Z+N_tmp;
		tmp_ocp_qp.idxs = ocp_qp->idxs+N_tmp;

		UPDATE_COND_BABT(idxc+N_tmp, &tmp_ocp_qp, part_dense_qp->BAbt+ii, part_dense_qp->b+ii, part_cond_arg->cond_arg+ii, part_cond_ws->cond_workspace+ii);

		UPDATE_COND_RSQRQ_N2NX3(idxc+N_tmp, &tmp_ocp_qp, part_dense_qp->RSQrq+ii, part_dense_qp->rqz+ii, part_cond_arg->cond_arg+ii, part_cond_ws->cond_workspace+ii);

		UPDATE_COND_DCTD(idxc+N_tmp, &tmp_ocp_qp, part_dense_qp->idxb[ii], part_dense_qp->DCt+ii, part_dense_qp->d+ii, part_dense_qp->idxs[ii], part_dense_qp->Z+ii, part_dense_qp->rqz+ii, part_cond_arg->cond_arg+ii, part_cond_ws->cond_workspace+ii);

		N_tmp += bs;

		}

#if 0
	// copy last stage
	int *nx = ocp_qp->dim->nx;
	int *nu = ocp_qp->dim->nu;
	int *nb = ocp_qp->dim->nb;
	int *ng = ocp_qp->dim->ng;
	int *ns = ocp_qp->dim->ns;

	GECP_LIBSTR(nu[N]+nx[N]+1, nu[N]+nx[N], ocp_qp->RSQrq+N, 0, 0, part_dense_qp->RSQrq+N2, 0, 0);
	VECCP_LIBSTR(nu[N]+nx[N], ocp_qp->rq+N, 0, part_dense_qp->rq+N2, 0);
	GECP_LIBSTR(nu[N]+nx[N], ng[N], ocp_qp->DCt+N, 0, 0, part_dense_qp->DCt+N2, 0, 0);
	VECCP_LIBSTR(2*nb[N]+2*ng[N], ocp_qp->d+N, 0, part_dense_qp->d+N2, 0);
	for(ii=0; ii<nb[N]; ii++) part_dense_qp->idxb[N2][ii] = ocp_qp->idxb[N][ii];
	VECCP_LIBSTR(2*ns[N], ocp_qp->Z+N, 0, part_dense_qp->Z+N2, 0);
	VECCP_LIBSTR(2*ns[N], ocp_qp->z+N, 0, part_dense_qp->z+N2, 0);
	for(ii=0; ii<ns[N]; ii++) part_dense_qp->idxs[N2][ii] = ocp_qp->idxs[N][ii];
#endif

	return;

	}



