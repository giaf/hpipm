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

void PART_COND_QP_COMPUTE_BLOCK_SIZE(int N, int N2, int *block_size)
	{

	int ii;

	if (N2 == 0 && N == 0)
		{
		block_size[0] = 0;
		return;
		}

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



void PART_COND_QP_COMPUTE_DIM(struct OCP_QP_DIM *ocp_dim, int *block_size, struct OCP_QP_DIM *part_dense_dim)
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

	// TODO equality constraints !!!!!!!!!

//	int nbb; // box constr that remain box constr
//	int nbg; // box constr that becomes general constr
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
		nsbx2[ii] += 0;
		nsbu2[ii] += nsbu[N_tmp+jj];
//		nsbx2[ii] = nsbx[N_tmp+0];
//		nsbu2[ii] = nsbu[N_tmp+0];
		nsg2[ii] += nsg[N_tmp+jj] + nsbx[N_tmp+jj];
		}

	return;

	}



hpipm_size_t PART_COND_QP_ARG_MEMSIZE(int N2)
	{

	int ii;

	hpipm_size_t size = 0;

	size += (N2+1)*sizeof(struct COND_QP_ARG);

	for(ii=0; ii<=N2; ii++)
		{

		size += COND_QP_ARG_MEMSIZE();

		}

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

	}



void PART_COND_QP_ARG_CREATE(int N2, struct PART_COND_QP_ARG *part_cond_arg, void *mem)
	{

	int ii;

	// cond workspace struct
	struct COND_QP_ARG *cws_ptr = mem;
	part_cond_arg->cond_arg = cws_ptr;
	cws_ptr += N2+1;

	// align to typical cache line size
	hpipm_size_t s_ptr = (hpipm_size_t) cws_ptr;
	s_ptr = (s_ptr+63)/64*64;

	char *c_ptr = (char *) s_ptr;

	for(ii=0; ii<=N2; ii++)
		{

		COND_QP_ARG_CREATE(part_cond_arg->cond_arg+ii, c_ptr);
		c_ptr += (part_cond_arg->cond_arg+ii)->memsize;

		}
	
	part_cond_arg->N2 = N2;

	part_cond_arg->memsize = PART_COND_QP_ARG_MEMSIZE(N2);

#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + part_cond_arg->memsize)
		{
#ifdef EXT_DEP
		printf("\nCreate_cond_qp_ocp2ocp_arg: outside memory bounds!\n\n");
#endif
		exit(1);
		}
#endif

	return;

	}



void PART_COND_QP_ARG_SET_DEFAULT(struct PART_COND_QP_ARG *part_cond_arg)
	{

	int ii;

	int N2 = part_cond_arg->N2;

	for(ii=0; ii<=N2; ii++)
		{

		COND_QP_ARG_SET_DEFAULT(part_cond_arg->cond_arg+ii);
		COND_QP_ARG_SET_COND_LAST_STAGE(0, part_cond_arg->cond_arg+ii);

		}
	// cond_last_stage at last stage
	COND_QP_ARG_SET_COND_LAST_STAGE(1, part_cond_arg->cond_arg+N2);

	return;

	}



void PART_COND_QP_ARG_SET_RIC_ALG(int ric_alg, struct PART_COND_QP_ARG *part_cond_arg)
	{

	int ii;

	int N2 = part_cond_arg->N2;

	for(ii=0; ii<=N2; ii++)
		{
		COND_QP_ARG_SET_RIC_ALG(ric_alg, part_cond_arg->cond_arg+ii);
		}

	return;

	}



void PART_COND_QP_ARG_SET_COMP_PRIM_SOL(int value, struct PART_COND_QP_ARG *part_cond_arg)
	{

	int ii;

	int N2 = part_cond_arg->N2;

	for(ii=0; ii<=N2; ii++)
		{
		COND_QP_ARG_SET_COMP_PRIM_SOL(value, part_cond_arg->cond_arg+ii);
		}

	return;

	}



void PART_COND_QP_ARG_SET_COMP_DUAL_SOL_EQ(int value, struct PART_COND_QP_ARG *part_cond_arg)
	{

	int ii;

	int N2 = part_cond_arg->N2;

	for(ii=0; ii<=N2; ii++)
		{
		COND_QP_ARG_SET_COMP_DUAL_SOL_EQ(value, part_cond_arg->cond_arg+ii);
		}

	return;

	}



void PART_COND_QP_ARG_SET_COMP_DUAL_SOL_INEQ(int value, struct PART_COND_QP_ARG *part_cond_arg)
	{

	int ii;

	int N2 = part_cond_arg->N2;

	for(ii=0; ii<=N2; ii++)
		{
		COND_QP_ARG_SET_COMP_DUAL_SOL_INEQ(value, part_cond_arg->cond_arg+ii);
		}

	return;

	}



hpipm_size_t PART_COND_QP_WS_MEMSIZE(struct OCP_QP_DIM *ocp_dim, int *block_size, struct OCP_QP_DIM *part_dense_dim, struct PART_COND_QP_ARG *part_cond_arg)
	{

	struct OCP_QP_DIM tmp_ocp_dim;

	int ii;

	int N = ocp_dim->N;
	int N2 = part_dense_dim->N;

	hpipm_size_t size = 0;

	size += (N2+1)*sizeof(struct COND_QP_ARG_WS);

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
		tmp_ocp_dim.nsbx = ocp_dim->nsbx+N_tmp;
		tmp_ocp_dim.nsbu = ocp_dim->nsbu+N_tmp;
		tmp_ocp_dim.nsg = ocp_dim->nsg+N_tmp;
		tmp_ocp_dim.ns = ocp_dim->ns+N_tmp;

		size += COND_QP_WS_MEMSIZE(&tmp_ocp_dim, part_cond_arg->cond_arg+ii);

		N_tmp += block_size[ii];

		}

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

	}



void PART_COND_QP_WS_CREATE(struct OCP_QP_DIM *ocp_dim, int *block_size, struct OCP_QP_DIM *part_dense_dim, struct PART_COND_QP_ARG *part_cond_arg, struct PART_COND_QP_WS *part_cond_ws, void *mem)
	{

	struct OCP_QP_DIM tmp_ocp_dim;

	int ii;

	int N = ocp_dim->N;
	int N2 = part_dense_dim->N;

	// cond workspace struct
	struct COND_QP_ARG_WS *cws_ptr = mem;
	part_cond_ws->cond_workspace = cws_ptr;
	cws_ptr += N2+1;

	// align to typical cache line size
	hpipm_size_t s_ptr = (hpipm_size_t) cws_ptr;
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
		tmp_ocp_dim.nsbx = ocp_dim->nsbx+N_tmp;
		tmp_ocp_dim.nsbu = ocp_dim->nsbu+N_tmp;
		tmp_ocp_dim.nsg = ocp_dim->nsg+N_tmp;
		tmp_ocp_dim.ns = ocp_dim->ns+N_tmp;
		// TODO equality constraints !!!!!!!!!!!!!!!!!!!!!!!!!!!

		COND_QP_WS_CREATE(&tmp_ocp_dim, part_cond_arg->cond_arg+ii, part_cond_ws->cond_workspace+ii, c_ptr);
		c_ptr += (part_cond_ws->cond_workspace+ii)->memsize;

		N_tmp += block_size[ii];

		}

	part_cond_ws->memsize = PART_COND_QP_WS_MEMSIZE(ocp_dim, block_size, part_dense_dim, part_cond_arg);

#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + part_cond_ws->memsize)
		{
#ifdef EXT_DEP
		printf("\nCreate_cond_qp_ocp2ocp: outside memory bounds!\n\n");
#endif
		exit(1);
		}
#endif

return;

	}



void PART_COND_QP_COND_LHS(struct OCP_QP *ocp_qp, struct OCP_QP *part_dense_qp, struct PART_COND_QP_ARG *part_cond_arg, struct PART_COND_QP_WS *part_cond_ws)
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
		tmp_ocp_dim.nsbx = ocp_qp->dim->nsbx+N_tmp;
		tmp_ocp_dim.nsbu = ocp_qp->dim->nsbu+N_tmp;
		tmp_ocp_dim.nsg = ocp_qp->dim->nsg+N_tmp;
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
		tmp_ocp_qp.d_mask = ocp_qp->d_mask+N_tmp;
		tmp_ocp_qp.Z = ocp_qp->Z+N_tmp;
		tmp_ocp_qp.idxs_rev = ocp_qp->idxs_rev+N_tmp;
		tmp_ocp_qp.diag_H_flag = ocp_qp->diag_H_flag+N_tmp;

		COND_BAT(&tmp_ocp_qp, part_dense_qp->BAbt+ii, part_cond_arg->cond_arg+ii, part_cond_ws->cond_workspace+ii);

		COND_RSQ(&tmp_ocp_qp, part_dense_qp->RSQrq+ii, part_cond_arg->cond_arg+ii, part_cond_ws->cond_workspace+ii);

		COND_DCT(&tmp_ocp_qp, part_dense_qp->idxb[ii], part_dense_qp->DCt+ii, part_dense_qp->idxs_rev[ii], part_dense_qp->Z+ii, part_cond_arg->cond_arg+ii, part_cond_ws->cond_workspace+ii);

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
	for(ii=0; ii<ns[N]; ii++) part_dense_qp->idxs_rev[N2][ii] = ocp_qp->idxs_rev[N][ii];
#endif

	return;

	}



void PART_COND_QP_COND(struct OCP_QP *ocp_qp, struct OCP_QP *part_dense_qp, struct PART_COND_QP_ARG *part_cond_arg, struct PART_COND_QP_WS *part_cond_ws)
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
		tmp_ocp_dim.nsbx = ocp_qp->dim->nsbx+N_tmp;
		tmp_ocp_dim.nsbu = ocp_qp->dim->nsbu+N_tmp;
		tmp_ocp_dim.nsg = ocp_qp->dim->nsg+N_tmp;
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
		tmp_ocp_qp.d_mask = ocp_qp->d_mask+N_tmp;
		tmp_ocp_qp.Z = ocp_qp->Z+N_tmp;
		tmp_ocp_qp.idxs_rev = ocp_qp->idxs_rev+N_tmp;
		tmp_ocp_qp.diag_H_flag = ocp_qp->diag_H_flag+N_tmp;

		COND_BABT(&tmp_ocp_qp, part_dense_qp->BAbt+ii, part_dense_qp->b+ii, part_cond_arg->cond_arg+ii, part_cond_ws->cond_workspace+ii);

		COND_RSQRQ(&tmp_ocp_qp, part_dense_qp->RSQrq+ii, part_dense_qp->rqz+ii, part_cond_arg->cond_arg+ii, part_cond_ws->cond_workspace+ii);

		COND_DCTD(&tmp_ocp_qp, part_dense_qp->idxb[ii], part_dense_qp->DCt+ii, part_dense_qp->d+ii, part_dense_qp->d_mask+ii, part_dense_qp->idxs_rev[ii], part_dense_qp->Z+ii, part_dense_qp->rqz+ii, part_cond_arg->cond_arg+ii, part_cond_ws->cond_workspace+ii);

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
	for(ii=0; ii<ns[N]; ii++) part_dense_qp->idxs_rev[N2][ii] = ocp_qp->idxs_rev[N][ii];
#endif

	return;

	}



void PART_COND_QP_COND_RHS(struct OCP_QP *ocp_qp, struct OCP_QP *part_dense_qp, struct PART_COND_QP_ARG *part_cond_arg, struct PART_COND_QP_WS *part_cond_ws)
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
		tmp_ocp_dim.nsbx = ocp_qp->dim->nsbx+N_tmp;
		tmp_ocp_dim.nsbu = ocp_qp->dim->nsbu+N_tmp;
		tmp_ocp_dim.nsg = ocp_qp->dim->nsg+N_tmp;
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
		tmp_ocp_qp.d_mask = ocp_qp->d_mask+N_tmp;
		tmp_ocp_qp.Z = ocp_qp->Z+N_tmp;
		tmp_ocp_qp.idxs_rev = ocp_qp->idxs_rev+N_tmp;
		tmp_ocp_qp.diag_H_flag = ocp_qp->diag_H_flag+N_tmp;

		COND_B(&tmp_ocp_qp, part_dense_qp->b+ii, part_cond_arg->cond_arg+ii, part_cond_ws->cond_workspace+ii);

		COND_RQ(&tmp_ocp_qp, part_dense_qp->rqz+ii, part_cond_arg->cond_arg+ii, part_cond_ws->cond_workspace+ii);

		// TODO cond m !!!!!!!!!!!
		COND_D(&tmp_ocp_qp, part_dense_qp->d+ii, part_dense_qp->d_mask+ii, part_dense_qp->rqz+ii, part_cond_arg->cond_arg+ii, part_cond_ws->cond_workspace+ii);

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



void PART_COND_QP_COND_SEED(struct OCP_QP *ocp_qp, struct OCP_QP_SEED *ocp_qp_seed, struct OCP_QP_SEED *part_dense_qp_seed, struct PART_COND_QP_ARG *part_cond_arg, struct PART_COND_QP_WS *part_cond_ws)
	{

	struct OCP_QP_DIM tmp_ocp_dim;
	struct OCP_QP tmp_ocp_qp;

	int ii;

	int N = ocp_qp->dim->N;
	int *nb = ocp_qp->dim->nb;
	int *ng = ocp_qp->dim->ng;

	int N2 = part_dense_qp_seed->dim->N;
	int *nb2 = part_dense_qp_seed->dim->nb;
	int *ng2 = part_dense_qp_seed->dim->ng;

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
		tmp_ocp_dim.nsbx = ocp_qp->dim->nsbx+N_tmp;
		tmp_ocp_dim.nsbu = ocp_qp->dim->nsbu+N_tmp;
		tmp_ocp_dim.nsg = ocp_qp->dim->nsg+N_tmp;
		tmp_ocp_dim.ns = ocp_qp->dim->ns+N_tmp;

		// alias ocp_qp
		tmp_ocp_qp.dim = &tmp_ocp_dim;
		tmp_ocp_qp.idxb = ocp_qp->idxb+N_tmp;
		tmp_ocp_qp.BAbt = ocp_qp->BAbt+N_tmp;
		tmp_ocp_qp.b = ocp_qp_seed->seed_b+N_tmp; // XXX
		tmp_ocp_qp.RSQrq = ocp_qp->RSQrq+N_tmp;
		tmp_ocp_qp.rqz = ocp_qp_seed->seed_g+N_tmp; //XXX
		tmp_ocp_qp.DCt = ocp_qp->DCt+N_tmp;
		tmp_ocp_qp.d = ocp_qp_seed->seed_d+N_tmp; // XXX
		tmp_ocp_qp.d_mask = ocp_qp->d_mask+N_tmp;
		tmp_ocp_qp.Z = ocp_qp->Z+N_tmp;
		tmp_ocp_qp.idxs_rev = ocp_qp->idxs_rev+N_tmp;
		tmp_ocp_qp.diag_H_flag = ocp_qp->diag_H_flag+N_tmp;
		// TODO cond m !!!!!!!!!!!

		COND_B(&tmp_ocp_qp, part_dense_qp_seed->seed_b+ii, part_cond_arg->cond_arg+ii, part_cond_ws->cond_workspace+ii);

		COND_RQ(&tmp_ocp_qp, part_dense_qp_seed->seed_g+ii, part_cond_arg->cond_arg+ii, part_cond_ws->cond_workspace+ii);

		COND_D(&tmp_ocp_qp, part_dense_qp_seed->seed_d+ii, NULL, part_dense_qp_seed->seed_g+ii, part_cond_arg->cond_arg+ii, part_cond_ws->cond_workspace+ii);

		N_tmp += bs;

		}

	return;

	}



void PART_COND_QP_COND_SOL(struct OCP_QP *ocp_qp, struct OCP_QP *part_dense_qp, struct OCP_QP_SOL *ocp_qp_sol, struct OCP_QP_SOL *part_dense_qp_sol, struct PART_COND_QP_ARG *part_cond_arg, struct PART_COND_QP_WS *part_cond_ws)
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
		tmp_ocp_dim.nsbx = ocp_qp->dim->nsbx+N_tmp;
		tmp_ocp_dim.nsbu = ocp_qp->dim->nsbu+N_tmp;
		tmp_ocp_dim.nsg = ocp_qp->dim->nsg+N_tmp;
		tmp_ocp_dim.ns = ocp_qp->dim->ns+N_tmp;

		// alias ocp_qp
		tmp_ocp_qp.dim = &tmp_ocp_dim;
		tmp_ocp_qp.idxb = ocp_qp->idxb+N_tmp;
		//tmp_ocp_qp.BAbt = ocp_qp->BAbt+N_tmp;
		//tmp_ocp_qp.b = ocp_qp->b+N_tmp;
		//tmp_ocp_qp.RSQrq = ocp_qp->RSQrq+N_tmp;
		//tmp_ocp_qp.rqz = ocp_qp->rqz+N_tmp;
		//tmp_ocp_qp.DCt = ocp_qp->DCt+N_tmp;
		//tmp_ocp_qp.d = ocp_qp->d+N_tmp;
		//tmp_ocp_qp.d_mask = ocp_qp->d_mask+N_tmp;
		//tmp_ocp_qp.Z = ocp_qp->Z+N_tmp;
		tmp_ocp_qp.idxs_rev = ocp_qp->idxs_rev+N_tmp;
		//tmp_ocp_qp.diag_H_flag = ocp_qp->diag_H_flag+N_tmp;

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

		COND_SOL(&tmp_ocp_qp, &tmp_ocp_qp_sol, &dense_qp_sol, part_cond_arg->cond_arg+ii, part_cond_ws->cond_workspace+ii);

		N_tmp += bs;

		}

#if 0
	// copy last stage
	VECCP_LIBSTR(nu[N]+nx[N]+2*ns[N], ocp_qp_sol->ux+N, 0, part_dense_qp_sol->ux+N2, 0);
	VECCP_LIBSTR(2*nb[N]+2*ng[N]+2*ns[N], 0, ocp_qp_sol->lam+N, 0, part_dense_qp_sol->lam+N2);
	VECCP_LIBSTR(2*nb[N]+2*ng[N]+2*ns[N], 0, ocp_qp_sol->t+N, 0, part_dense_qp_sol->t+N2);
#endif

	return;

	}



void PART_COND_QP_EXPAND_SOL(struct OCP_QP *ocp_qp, struct OCP_QP_SOL *part_dense_qp_sol, struct OCP_QP_SOL *ocp_qp_sol, struct PART_COND_QP_ARG *part_cond_arg, struct PART_COND_QP_WS *part_cond_ws)
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
	int N2 = part_dense_qp_sol->dim->N;
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
		tmp_ocp_dim.nsbx = ocp_qp->dim->nsbx+N_tmp;
		tmp_ocp_dim.nsbu = ocp_qp->dim->nsbu+N_tmp;
		tmp_ocp_dim.nsg = ocp_qp->dim->nsg+N_tmp;
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
		tmp_ocp_qp.d_mask = ocp_qp->d_mask+N_tmp;
		tmp_ocp_qp.Z = ocp_qp->Z+N_tmp;
		tmp_ocp_qp.idxs_rev = ocp_qp->idxs_rev+N_tmp;
		tmp_ocp_qp.diag_H_flag = ocp_qp->diag_H_flag+N_tmp;

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



void PART_COND_QP_EXPAND_SOL_SEED(struct OCP_QP *ocp_qp, struct OCP_QP_SEED *ocp_qp_seed, struct OCP_QP_SOL *part_dense_qp_sol, struct OCP_QP_SOL *ocp_qp_sol, struct PART_COND_QP_ARG *part_cond_arg, struct PART_COND_QP_WS *part_cond_ws)
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
	int N2 = part_dense_qp_sol->dim->N;
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
		tmp_ocp_dim.nsbx = ocp_qp->dim->nsbx+N_tmp;
		tmp_ocp_dim.nsbu = ocp_qp->dim->nsbu+N_tmp;
		tmp_ocp_dim.nsg = ocp_qp->dim->nsg+N_tmp;
		tmp_ocp_dim.ns = ocp_qp->dim->ns+N_tmp;

		// alias ocp_qp
		tmp_ocp_qp.dim = &tmp_ocp_dim;
		tmp_ocp_qp.idxb = ocp_qp->idxb+N_tmp;
		tmp_ocp_qp.BAbt = ocp_qp->BAbt+N_tmp;
		tmp_ocp_qp.b = ocp_qp_seed->seed_b+N_tmp; // XXX
		tmp_ocp_qp.RSQrq = ocp_qp->RSQrq+N_tmp;
		tmp_ocp_qp.rqz = ocp_qp_seed->seed_g+N_tmp; // XXX
		tmp_ocp_qp.DCt = ocp_qp->DCt+N_tmp;
		tmp_ocp_qp.d = ocp_qp_seed->seed_d+N_tmp; // XXX
		tmp_ocp_qp.d_mask = ocp_qp->d_mask+N_tmp;
		tmp_ocp_qp.Z = ocp_qp->Z+N_tmp;
		tmp_ocp_qp.idxs_rev = ocp_qp->idxs_rev+N_tmp;
		tmp_ocp_qp.diag_H_flag = ocp_qp->diag_H_flag+N_tmp;

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

void PART_COND_QP_UPDATE(int *idxc, struct OCP_QP *ocp_qp, struct OCP_QP *part_dense_qp, struct PART_COND_QP_ARG *part_cond_arg, struct PART_COND_QP_WS *part_cond_ws)
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
		tmp_ocp_dim.nsbx = ocp_qp->dim->nsbx+N_tmp;
		tmp_ocp_dim.nsbu = ocp_qp->dim->nsbu+N_tmp;
		tmp_ocp_dim.nsg = ocp_qp->dim->nsg+N_tmp;
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
		tmp_ocp_qp.d_mask = ocp_qp->d_mask+N_tmp;
		tmp_ocp_qp.Z = ocp_qp->Z+N_tmp;
		tmp_ocp_qp.idxs_rev = ocp_qp->idxs_rev+N_tmp;
		tmp_ocp_qp.diag_H_flag = ocp_qp->diag_H_flag+N_tmp;

		UPDATE_COND_BABT(idxc+N_tmp, &tmp_ocp_qp, part_dense_qp->BAbt+ii, part_dense_qp->b+ii, part_cond_arg->cond_arg+ii, part_cond_ws->cond_workspace+ii);

		UPDATE_COND_RSQRQ_N2NX3(idxc+N_tmp, &tmp_ocp_qp, part_dense_qp->RSQrq+ii, part_dense_qp->rqz+ii, part_cond_arg->cond_arg+ii, part_cond_ws->cond_workspace+ii);

		UPDATE_COND_DCTD(idxc+N_tmp, &tmp_ocp_qp, part_dense_qp->idxb[ii], part_dense_qp->DCt+ii, part_dense_qp->d+ii, part_dense_qp->idxs_rev[ii], part_dense_qp->Z+ii, part_dense_qp->rqz+ii, part_cond_arg->cond_arg+ii, part_cond_ws->cond_workspace+ii);

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
	for(ii=0; ii<ns[N]; ii++) part_dense_qp->idxs_rev[N2][ii] = ocp_qp->idxs_rev[N][ii];
#endif

	return;

	}



