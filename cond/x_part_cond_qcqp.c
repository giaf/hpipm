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

void PART_COND_QCQP_COMPUTE_BLOCK_SIZE(int N, int N2, int *block_size)
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



void PART_COND_QCQP_COMPUTE_DIM(struct OCP_QCQP_DIM *ocp_dim, int *block_size, struct OCP_QCQP_DIM *part_dense_dim)
	{

	// TODO run time check on sum(block_size) = N

	int N = ocp_dim->N;
	int *nx = ocp_dim->nx;
	int *nu = ocp_dim->nu;
	int *nb = ocp_dim->nb;
	int *nbx = ocp_dim->nbx;
	int *nbu = ocp_dim->nbu;
	int *ng = ocp_dim->ng;
	int *nq = ocp_dim->nq;
	int *ns = ocp_dim->ns;
	int *nsbx = ocp_dim->nsbx;
	int *nsbu = ocp_dim->nsbu;
	int *nsg = ocp_dim->nsg;
	int *nsq = ocp_dim->nsq;

	int N2 = part_dense_dim->N;
//	int *nx2 = part_dense_dim->nx;
//	int *nu2 = part_dense_dim->nu;
//	int *nb2 = part_dense_dim->nb;
//	int *nbx2 = part_dense_dim->nbx;
//	int *nbu2 = part_dense_dim->nbu;
//	int *ng2 = part_dense_dim->ng;
//	int *nq2 = part_dense_dim->nq;
//	int *ns2 = part_dense_dim->ns;
//	int *nsbx2 = part_dense_dim->nsbx;
//	int *nsbu2 = part_dense_dim->nsbu;
//	int *nsg2 = part_dense_dim->nsg;
//	int *nsq2 = part_dense_dim->nsq;

	int nx2, nu2, nb2, nbx2, nbu2, ng2, nq2, ns2, nsbu2, nsbx2, nsg2, nsq2;

	int ii, jj;

	int nbb; // box constr that remain box constr
	int nbg; // box constr that becomes general constr
	int N_tmp = 0; // temporary sum of block size
	// first stages
	for(ii=0; ii<N2; ii++)
		{
		nx2 = nx[N_tmp+0];
		nu2 = nu[N_tmp+0];
		nbx2 = nbx[N_tmp+0];
		nbu2 = nbu[N_tmp+0];
		nb2 = nb[N_tmp+0];
		ng2 = ng[N_tmp+0];
		nq2 = nq[N_tmp+0];
		ns2 = ns[N_tmp+0];
		nsbx2 = nsbx[N_tmp+0];
		nsbu2 = nsbu[N_tmp+0];
		nsg2 = nsg[N_tmp+0];
		nsq2 = nsq[N_tmp+0];
		for(jj=1; jj<block_size[ii]; jj++)
			{
			nx2 += 0;
			nu2 += nu[N_tmp+jj];
			nbx2 += 0;
			nbu2 += nbu[N_tmp+jj];
			nb2 += nbu[N_tmp+jj];
			ng2 += ng[N_tmp+jj] + nbx[N_tmp+jj];
			nq2 += nq[N_tmp+jj];
			ns2 += ns[N_tmp+jj];
			nsbx2 += 0;
			nsbu2 += nsbu[N_tmp+jj];
			nsg2 += nsg[N_tmp+jj] + nsbx[N_tmp+jj];
			nsq2 += nsq[N_tmp+jj];
			}
		N_tmp += block_size[ii];
		// XXX must use setters to correctly set qp ones too !
		OCP_QCQP_DIM_SET_NX(nx2, ii, part_dense_dim);
		OCP_QCQP_DIM_SET_NU(nu2, ii, part_dense_dim);
		OCP_QCQP_DIM_SET_NBX(nbx2, ii, part_dense_dim);
		OCP_QCQP_DIM_SET_NBU(nbu2, ii, part_dense_dim);
		OCP_QCQP_DIM_SET_NB(nb2, ii, part_dense_dim);
		OCP_QCQP_DIM_SET_NG(ng2, ii, part_dense_dim);
		OCP_QCQP_DIM_SET_NQ(nq2, ii, part_dense_dim);
		OCP_QCQP_DIM_SET_NS(ns2, ii, part_dense_dim);
		OCP_QCQP_DIM_SET_NSBX(nsbx2, ii, part_dense_dim);
		OCP_QCQP_DIM_SET_NSBU(nsbu2, ii, part_dense_dim);
		OCP_QCQP_DIM_SET_NSB(nsb2, ii, part_dense_dim);
		OCP_QCQP_DIM_SET_NSG(nsg2, ii, part_dense_dim);
		OCP_QCQP_DIM_SET_NSQ(nsq2, ii, denspart_e_dim);
		}
	// last stage: condense also following stage
	ii = N2;
	nx2 = nx[N_tmp+0];
	nu2 = nu[N_tmp+0];
	nbx2 = nbx[N_tmp+0];
	nbu2 = nbu[N_tmp+0];
	nb2 = nb[N_tmp+0];
	ng2 = ng[N_tmp+0];
	nq2 = nq[N_tmp+0];
	ns2 = ns[N_tmp+0];
	nsbx2 = nsbx[N_tmp+0];
	nsbu2 = nsbu[N_tmp+0];
	nsg2 = nsg[N_tmp+0];
	nsq2 = nsq[N_tmp+0];
	for(jj=1; jj<block_size[ii]+1; jj++)
		{
		nx2 += 0;
		nu2 += nu[N_tmp+jj];
		nbx2 += 0;
		nbu2 += nbu[N_tmp+jj];
		nb2 += nbu[N_tmp+jj];
		ng2 += ng[N_tmp+jj] + nbx[N_tmp+jj];
		nq2 += nq[N_tmp+jj];
		ns2 += ns[N_tmp+jj];
		nsbx2 += 0;
		nsbu2 += nsbu[N_tmp+jj];
//		nsbx2 = nsbx[N_tmp+0];
//		nsbu2 = nsbu[N_tmp+0];
		nsg2 += nsg[N_tmp+jj] + nsbx[N_tmp+jj];
		nsq2 += nsq[N_tmp+jj];
		}
	// XXX must use setters to correctly set qp ones too !
	OCP_QCQP_DIM_SET_NX(nx2, ii, part_dense_dim);
	OCP_QCQP_DIM_SET_NU(nu2, ii, part_dense_dim);
	OCP_QCQP_DIM_SET_NBX(nbx2, ii, part_dense_dim);
	OCP_QCQP_DIM_SET_NBU(nbu2, ii, part_dense_dim);
	OCP_QCQP_DIM_SET_NB(nb2, ii, part_dense_dim);
	OCP_QCQP_DIM_SET_NG(ng2, ii, part_dense_dim);
	OCP_QCQP_DIM_SET_NQ(nq2, ii, part_dense_dim);
	OCP_QCQP_DIM_SET_NS(ns2, ii, part_dense_dim);
	OCP_QCQP_DIM_SET_NSBX(nsbx2, ii, part_dense_dim);
	OCP_QCQP_DIM_SET_NSBU(nsbu2, ii, part_dense_dim);
	OCP_QCQP_DIM_SET_NSB(nsb2, ii, part_dense_dim);
	OCP_QCQP_DIM_SET_NSG(nsg2, ii, part_dense_dim);
	OCP_QCQP_DIM_SET_NSQ(nsq2, ii, denspart_e_dim);

	return;

	}




