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



void OCP_QCQP_DIM_REDUCE_EQ_DOF(struct OCP_QCQP_DIM *dim, struct OCP_QCQP_DIM *dim_red)
	{
	int ii;

	// XXX must use setters to correctly set qp ones too !
	// first stage: DOF: inputs and states
	OCP_QCQP_DIM_SET_NU(0, dim->nu[0] - dim->nbue[0], dim_red);
	OCP_QCQP_DIM_SET_NX(0, dim->nx[0] - dim->nbxe[0], dim_red);
	OCP_QCQP_DIM_SET_NBU(0, dim->nbu[0] - dim->nbue[0], dim_red);
	OCP_QCQP_DIM_SET_NBX(0, dim->nbx[0] - dim->nbxe[0], dim_red);
	OCP_QCQP_DIM_SET_NG(0, dim->ng[0], dim_red);
	OCP_QCQP_DIM_SET_NQ(0, dim->nq[0], dim_red);
	OCP_QCQP_DIM_SET_NS(0, dim->ns[0], dim_red);
	OCP_QCQP_DIM_SET_NSBU(0, dim->nsbu[0], dim_red);
	OCP_QCQP_DIM_SET_NSBX(0, dim->nsbx[0], dim_red);
	OCP_QCQP_DIM_SET_NSG(0, dim->nsg[0], dim_red);
	OCP_QCQP_DIM_SET_NSQ(0, dim->nsq[0], dim_red);
	OCP_QCQP_DIM_SET_NBUE(0, 0, dim_red);
	OCP_QCQP_DIM_SET_NBXE(0, 0, dim_red);
	OCP_QCQP_DIM_SET_NGE(0, dim->nge[0], dim_red);
	OCP_QCQP_DIM_SET_NQE(0, dim->nqe[0], dim_red);
	// other stages: DOF: inputs
	for(ii=1; ii<=dim->N; ii++)
		{
		OCP_QCQP_DIM_SET_NU(ii, dim->nu[ii] - dim->nbue[ii], dim_red);
		OCP_QCQP_DIM_SET_NX(ii, dim->nx[ii], dim_red);
		OCP_QCQP_DIM_SET_NBU(ii, dim->nbu[ii] - dim->nbue[ii], dim_red);
		OCP_QCQP_DIM_SET_NBX(ii, dim->nbx[ii], dim_red);
		OCP_QCQP_DIM_SET_NG(ii, dim->ng[ii], dim_red);
		OCP_QCQP_DIM_SET_NQ(ii, dim->nq[ii], dim_red);
		OCP_QCQP_DIM_SET_NS(ii, dim->ns[ii], dim_red);
		OCP_QCQP_DIM_SET_NSBU(ii, dim->nsbu[ii], dim_red);
		OCP_QCQP_DIM_SET_NSBX(ii, dim->nsbx[ii], dim_red);
		OCP_QCQP_DIM_SET_NSG(ii, dim->nsg[ii], dim_red);
		OCP_QCQP_DIM_SET_NSQ(ii, dim->nsq[ii], dim_red);
		OCP_QCQP_DIM_SET_NBUE(ii, 0, dim_red);
		OCP_QCQP_DIM_SET_NBXE(ii, dim->nbxe[ii], dim_red);
		OCP_QCQP_DIM_SET_NGE(ii, dim->nge[ii], dim_red);
		OCP_QCQP_DIM_SET_NQE(ii, dim->nqe[ii], dim_red);
		}
	
	return;
	}



hpipm_size_t OCP_QCQP_REDUCE_EQ_DOF_ARG_MEMSIZE()
	{

	return 0;

	}



void OCP_QCQP_REDUCE_EQ_DOF_ARG_CREATE(struct OCP_QCQP_REDUCE_EQ_DOF_ARG *arg, void *mem)
	{

	arg->memsize = OCP_QCQP_REDUCE_EQ_DOF_ARG_MEMSIZE();

	return;

	}



void OCP_QCQP_REDUCE_EQ_DOF_ARG_SET_DEFAULT(struct OCP_QCQP_REDUCE_EQ_DOF_ARG *arg)
	{

	arg->lam_min = 1e-16;
	arg->t_min = 1e-16;
	arg->alias_unchanged = 0;
	arg->comp_prim_sol = 1;
	arg->comp_dual_sol_eq = 1;
	arg->comp_dual_sol_ineq = 1;

	return;

	}



void OCP_QCQP_REDUCE_EQ_DOF_ARG_SET_ALIAS_UNCHANGED(int *value, struct OCP_QCQP_REDUCE_EQ_DOF_ARG *arg)
	{

	arg->alias_unchanged = *value;

	return;

	}



void OCP_QCQP_REDUCE_EQ_DOF_ARG_SET_COMP_PRIM_SOL(int *value, struct OCP_QCQP_REDUCE_EQ_DOF_ARG *arg)
	{

	arg->comp_prim_sol = *value;

	return;

	}



void OCP_QCQP_REDUCE_EQ_DOF_ARG_SET_COMP_DUAL_SOL_EQ(int *value, struct OCP_QCQP_REDUCE_EQ_DOF_ARG *arg)
	{

	arg->comp_dual_sol_eq = *value;

	return;

	}



void OCP_QCQP_REDUCE_EQ_DOF_ARG_SET_COMP_DUAL_SOL_INEQ(int *value, struct OCP_QCQP_REDUCE_EQ_DOF_ARG *arg)
	{

	arg->comp_dual_sol_ineq = *value;

	return;

	}



hpipm_size_t OCP_QCQP_REDUCE_EQ_DOF_WS_MEMSIZE(struct OCP_QCQP_DIM *dim)
	{

	int ii;

	// extract dim
	int N = dim->N;
	int *nu = dim->nu;
	int *nx = dim->nx;
	int *nb = dim->nb;
	int *ng = dim->ng;
	int *nq = dim->nq;

	int nuxM = nu[0]+nx[0];
	int nbgqM = nb[0]+ng[0]+nq[0];
	for(ii=1; ii<=N; ii++)
		{
		nuxM = nu[ii]+nx[ii]>nuxM ? nu[ii]+nx[ii] : nuxM;
		nbgqM = nb[ii]+ng[ii]+nq[ii]>nbgqM ? nb[ii]+ng[ii]+nq[ii] : nbgqM;
		}

	hpipm_size_t size = 0;

	size += nuxM*sizeof(int); // e_imask_ux
	size += nbgqM*sizeof(int); // e_imask_d

	size += 3*sizeof(struct STRVEC); // tmp_nuxM(0,1) tmp_nbgqM

	size += 2*SIZE_STRVEC(nuxM); // tmp_nuxM(0,1)
	size += 1*SIZE_STRVEC(nbgqM); // tmp_nbgqM

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 64; // align to typical cache line size

	return size;

	}



void OCP_QCQP_REDUCE_EQ_DOF_WS_CREATE(struct OCP_QCQP_DIM *dim, struct OCP_QCQP_REDUCE_EQ_DOF_WS *work, void *mem)
	{

	int ii;

	// zero memory (to avoid corrupted memory like e.g. NaN)
	hpipm_size_t memsize = OCP_QCQP_REDUCE_EQ_DOF_WS_MEMSIZE(dim);
	hpipm_zero_memset(memsize, mem);

	// extract dim
	int N = dim->N;
	int *nu = dim->nu;
	int *nx = dim->nx;
	int *nb = dim->nb;
	int *ng = dim->ng;
	int *nq = dim->nq;

	int nuxM = nu[0]+nx[0];
	int nbgqM = nb[0]+ng[0]+nq[0];
	for(ii=1; ii<=N; ii++)
		{
		nuxM = nu[ii]+nx[ii]>nuxM ? nu[ii]+nx[ii] : nuxM;
		nbgqM = nb[ii]+ng[ii]+nq[ii]>nbgqM ? nb[ii]+ng[ii]+nq[ii] : nbgqM;
		}

	// vector struct stuff
	struct STRVEC *sv_ptr = (struct STRVEC *) mem;

	// tmp_nuxM
	work->tmp_nuxM = sv_ptr;
	sv_ptr += 2;
	// tmp_nbgqM
	work->tmp_nbgqM = sv_ptr;
	sv_ptr += 1;

	// integer stuff
	int *i_ptr;
	i_ptr = (int *) sv_ptr;

	// e_imask_ux
	work->e_imask_ux = i_ptr;
	i_ptr += nuxM;
	// e_imask_d
	work->e_imask_d = i_ptr;
	i_ptr += nbgqM;

	// align to typical cache line size
	hpipm_size_t l_ptr = (hpipm_size_t) i_ptr;
	l_ptr = (l_ptr+63)/64*64;

	// floating point stuff
	char *c_ptr;
	c_ptr = (char *) l_ptr;

	// tmp_nuxM
	CREATE_STRVEC(nuxM, work->tmp_nuxM+0, c_ptr);
	c_ptr += (work->tmp_nuxM+0)->memsize;
	CREATE_STRVEC(nuxM, work->tmp_nuxM+1, c_ptr);
	c_ptr += (work->tmp_nuxM+1)->memsize;
	// tmp_nbgqM
	CREATE_STRVEC(nbgqM, work->tmp_nbgqM, c_ptr);
	c_ptr += (work->tmp_nbgqM)->memsize;

	work->memsize = memsize;

	return;

	}



void OCP_QCQP_REDUCE_EQ_DOF(struct OCP_QCQP *qp, struct OCP_QCQP *qp_red, struct OCP_QCQP_REDUCE_EQ_DOF_ARG *arg, struct OCP_QCQP_REDUCE_EQ_DOF_WS *work)
	{

	int ii, jj, kk, ll, idx0, idx1;

	struct OCP_QCQP_DIM *dim = qp->dim;
	int N = dim->N;
	int *nx = dim->nx;
	int *nu = dim->nu;
	int *nb = dim->nb;
	int *nbx = dim->nbx;
	int *nbu = dim->nbu;
	int *ng = dim->ng;
	int *nq = dim->nq;
	int *ns = dim->ns;
	int *nbue = dim->nbue;
	int *nbxe = dim->nbxe;
	int *nge = dim->nge;
	int *nqe = dim->nqe;

	struct OCP_QCQP_DIM *dim_red = qp_red->dim;
	int *nx_red = dim_red->nx;
	int *nu_red = dim_red->nu;
	int *nb_red = dim_red->nb;
	int *ng_red = dim_red->ng;
	int *nq_red = dim_red->nq;
	int *ns_red = dim_red->ns;

	// TODO handle case of softed equalities !!!!!!!!!!!!!!!!

	int ne_thr;

	for(ii=0; ii<=N; ii++)
		{
		if(ii==0)
			ne_thr = nbue[ii]+nbxe[ii];
		else
			ne_thr = nbue[ii];
		if(ne_thr>0) // reduce inputs and/or states
			{
			VECSE(nu[ii]+nx[ii], 0.0, work->tmp_nuxM+0, 0);
			for(jj=0; jj<nu[ii]+nx[ii]; jj++)
				work->e_imask_ux[jj] = 0;
			for(jj=0; jj<nbu[ii]+nbx[ii]; jj++)
				work->e_imask_d[jj] = 0;
			for(jj=0; jj<ne_thr; jj++) // set 1s for both inputs and states
				{
				VECEL(work->tmp_nuxM+0, qp->idxb[ii][qp->idxe[ii][jj]]) = VECEL(qp->d+ii, qp->idxe[ii][jj]);
				work->e_imask_ux[qp->idxb[ii][qp->idxe[ii][jj]]] = 1;
				work->e_imask_d[qp->idxe[ii][jj]] = 1;
				}
			// TODO check first and last non-zero in e_mask and only multiply between them
			if(ii<N)
				{
				// BAt
				idx0 = 0;
				for(jj=0; jj<nu[ii]+nx[ii]; jj++)
					{
					if(work->e_imask_ux[jj]==0)
						{
						GECP(1, nx[ii+1], qp->BAbt+ii, jj, 0, qp_red->BAbt+ii, idx0, 0);
						idx0++;
						}
					}
				// b
				GEMV_T(nu[ii]+nx[ii], nx[ii+1], 1.0, qp->BAbt+ii, 0, 0, work->tmp_nuxM+0, 0, 1.0, qp->b+ii, 0, qp_red->b+ii, 0);
				}
			// RSQ & Hq
			idx0 = 0;
			for(jj=0; jj<nu[ii]+nx[ii]; jj++)
				{
				if(work->e_imask_ux[jj]==0)
					{
//					idx1 = 0;
					idx1 = idx0;
//					for(kk=0; kk<nu[ii]+nx[ii]; kk++)
					for(kk=jj; kk<nu[ii]+nx[ii]; kk++)
						{
						if(work->e_imask_ux[kk]==0)
							{
							MATEL(qp_red->RSQrq+ii, idx1, idx0) = MATEL(qp->RSQrq+ii, kk, jj);
							for(ll=0; ll<nq[ii]; ll++)
								{
								MATEL(qp_red->Hq[ii]+ll, idx1, idx0) = MATEL(qp->Hq[ii]+ll, kk, jj);
								}
//#endif
							idx1++;
							}
						}
					idx0++;
					}
				}
			// (at least) a sub-matrix of a zero matrix is a zero matrix
			// TODO do better Hq_nzero ???
			for(jj=0; jj<nq[ii]; jj++)
				{
				qp_red->Hq_nzero[ii][jj] = qp->Hq_nzero[ii][jj];
				}
			// rq
			SYMV_L(nu[ii]+nx[ii], 1.0, qp->RSQrq+ii, 0, 0, work->tmp_nuxM+0, 0, 1.0, qp->rqz+ii, 0, work->tmp_nuxM+1, 0);
			idx0 = 0;
			for(jj=0; jj<nu[ii]+nx[ii]; jj++)
				{
				if(work->e_imask_ux[jj]==0)
					{
					VECEL(qp_red->rqz+ii, idx0) = VECEL(work->tmp_nuxM+1, jj);
					idx0++;
					}
				}
			// gq
			for(ll=0; ll<nq[ii]; ll++)
				{
				COLEX(nu[ii]+nx[ii], qp->DCt+ii, 0, ng[ii]+ll, work->tmp_nuxM+1, 0);
				SYMV_L(nu[ii]+nx[ii], 1.0, qp->Hq[ii]+ll, 0, 0, work->tmp_nuxM+0, 0, 1.0, work->tmp_nuxM+1, 0, work->tmp_nuxM+1, 0);
				VECEL(work->tmp_nbgqM, nb[ii]+ng[ii]+ll) = DOT(nu[ii]+nx[ii], work->tmp_nuxM+0, 0, work->tmp_nuxM+1, 0);
				idx0 = 0;
				for(jj=0; jj<nu[ii]+nx[ii]; jj++)
					{
					if(work->e_imask_ux[jj]==0)
						{
						MATEL(qp_red->DCt+ii, idx0, ng[ii]+ll) = VECEL(work->tmp_nuxM+1, jj);
						idx0++;
						}
					}
				}
			VECCP(nq[ii], qp->d+ii, nb[ii]+ng[ii], qp_red->d+ii, nb_red[ii]+ng_red[ii]);
			VECCP(nq[ii], qp->d+ii, 2*nb[ii]+2*ng[ii]+nq[ii], qp_red->d+ii, 2*nb_red[ii]+2*ng_red[ii]+nq_red[ii]);
			VECCP(nq[ii], qp->d_mask+ii, nb[ii]+ng[ii], qp_red->d_mask+ii, nb_red[ii]+ng_red[ii]);
			VECCP(nq[ii], qp->d_mask+ii, 2*nb[ii]+2*ng[ii], qp_red->d_mask+ii, 2*nb_red[ii]+2*ng_red[ii]+nq_red[ii]);
			VECCP(nq[ii], qp->m+ii, nb[ii]+ng[ii], qp_red->m+ii, nb_red[ii]+ng_red[ii]);
			VECCP(nq[ii], qp->m+ii, 2*nb[ii]+2*ng[ii], qp_red->m+ii, 2*nb_red[ii]+2*ng_red[ii]+nq_red[ii]);
//			AXPY(nq[ii], -1.0, work->tmp_nbgqM, nb[ii]+ng[ii], qp->d+ii, nb[ii]+ng[ii], qp_red->d+ii, nb_red[ii]+ng_red[ii]);
			AXPY(nq[ii], 1.0, work->tmp_nbgqM, 2*nb[ii]+2*ng[ii]+nq[ii], qp->d+ii, 2*nb[ii]+2*ng[ii]+nq[ii], qp_red->d+ii, 2*nb_red[ii]+2*ng_red[ii]+nq_red[ii]);
			for(jj=0; jj<nq[ii]; jj++)
				qp_red->idxs_rev[ii][nb_red[ii]+ng_red[ii]+jj] = qp->idxs_rev[ii][nb[ii]+ng[ii]+jj]; // keep softed inequality constr with same slack
			// d d_mask m idxb idxs_rev
			idx0 = 0;
			for(jj=0; jj<nb[ii]; jj++)
				{
				if(work->e_imask_d[jj]==0)
					{
					VECEL(qp_red->d+ii, idx0) = VECEL(qp->d+ii, jj);
					VECEL(qp_red->d+ii, nb_red[ii]+ng_red[ii]+nq_red[ii]+idx0) = VECEL(qp->d+ii, nb[ii]+ng[ii]+nq[ii]+jj);
					VECEL(qp_red->d_mask+ii, idx0) = VECEL(qp->d_mask+ii, jj);
					VECEL(qp_red->d_mask+ii, nb_red[ii]+ng_red[ii]+nq_red[ii]+idx0) = VECEL(qp->d_mask+ii, nb[ii]+ng[ii]+nq[ii]+jj);
					VECEL(qp_red->m+ii, idx0) = VECEL(qp->m+ii, jj);
					VECEL(qp_red->m+ii, nb_red[ii]+ng_red[ii]+nq_red[ii]+idx0) = VECEL(qp->m+ii, nb[ii]+ng[ii]+nq[ii]+jj);
					qp_red->idxb[ii][idx0] = qp->idxb[ii][jj];
					qp_red->idxs_rev[ii][idx0] = qp->idxs_rev[ii][jj]; // keep softed inequality constr with same slack
					idx0++;
					}
				}
			VECCP(ng[ii], qp->d+ii, nb[ii], qp_red->d+ii, nb_red[ii]);
			VECCP(ng[ii], qp->d+ii, 2*nb[ii]+ng[ii]+nq[ii], qp_red->d+ii, 2*nb_red[ii]+ng_red[ii]+nq_red[ii]);
			VECCP(ng[ii], qp->d_mask+ii, nb[ii], qp_red->d_mask+ii, nb_red[ii]);
			VECCP(ng[ii], qp->d_mask+ii, 2*nb[ii]+ng[ii]+nq[ii], qp_red->d_mask+ii, 2*nb_red[ii]+ng_red[ii]+nq_red[ii]);
			VECCP(ng[ii], qp->m+ii, nb[ii], qp_red->m+ii, nb_red[ii]);
			VECCP(ng[ii], qp->m+ii, 2*nb[ii]+ng[ii]+nq[ii], qp_red->m+ii, 2*nb_red[ii]+ng_red[ii]+nq_red[ii]);
			GEMV_T(nu[ii]+nx[ii], ng[ii], 1.0, qp->DCt+ii, 0, 0, work->tmp_nuxM+0, 0, 0.0, work->tmp_nbgqM, nb[ii], work->tmp_nbgqM, nb[ii]);
			AXPY(ng[ii], -1.0, work->tmp_nbgqM, nb[ii], qp->d+ii, nb[ii], qp_red->d+ii, nb_red[ii]);
			AXPY(ng[ii], 1.0, work->tmp_nbgqM, nb[ii], qp->d+ii, 2*nb[ii]+ng[ii]+nq[ii], qp_red->d+ii, 2*nb_red[ii]+ng_red[ii]+nq_red[ii]);
			// DCt
			idx0 = 0;
			for(jj=0; jj<nu[ii]+nx[ii]; jj++)
				{
				if(work->e_imask_ux[jj]==0)
					{
					GECP(1, ng[ii], qp->DCt+ii, jj, 0, qp_red->DCt+ii, idx0, 0);
					idx0++;
					}
				}
			for(jj=0; jj<ng[ii]; jj++)
				qp_red->idxs_rev[ii][nb_red[ii]+jj] = qp->idxs_rev[ii][nb[ii]+jj]; // keep softed inequality constr with same slack
			// soft constraints
			VECCP(2*ns[ii], qp->d+ii, 2*nb[ii]+2*ng[ii]+2*nq[ii], qp_red->d+ii, 2*nb_red[ii]+2*ng_red[ii]+2*nq_red[ii]);
			VECCP(2*ns[ii], qp->d_mask+ii, 2*nb[ii]+2*ng[ii]+2*nq[ii], qp_red->d_mask+ii, 2*nb_red[ii]+2*ng_red[ii]+2*nq_red[ii]);
			VECCP(2*ns[ii], qp->m+ii, 2*nb[ii]+2*ng[ii]+2*nq[ii], qp_red->m+ii, 2*nb_red[ii]+2*ng_red[ii]+2*nq_red[ii]);
			VECCP(2*ns[ii], qp->Z+ii, 0, qp_red->Z+ii, 0);
			VECCP(2*ns[ii], qp->rqz+ii, nu[ii]+nx[ii], qp_red->rqz+ii, nu_red[ii]+nx_red[ii]);
//			qp_red->diag_H_flag[ii] = qp->diag_H_flag[ii];
			// TODO idxe !!!!!!!!!!!!!!!
			}
		else // copy everything
			{
			// copy vectors which are contiguous in the QP (e.g. to alias to res)
			if(ii<N)
				{
				VECCP(nx[ii+1], qp->b+ii, 0, qp_red->b+ii, 0);
				}
			VECCP(nu[ii]+nx[ii]+2*ns[ii], qp->rqz+ii, 0, qp_red->rqz+ii, 0);
			VECCP(2*nb[ii]+2*ng[ii]+2*nq[ii]+2*ns[ii], qp->d+ii, 0, qp_red->d+ii, 0);
			VECCP(2*nb[ii]+2*ng[ii]+2*nq[ii]+2*ns[ii], qp->d_mask+ii, 0, qp_red->d_mask+ii, 0);
			VECCP(2*nb[ii]+2*ng[ii]+2*nq[ii]+2*ns[ii], qp->m+ii, 0, qp_red->m+ii, 0);
//			qp_red->diag_H_flag[ii] = qp->diag_H_flag[ii];
			for(jj=0; jj<nq[ii]; jj++)
				{
				qp_red->Hq_nzero[ii][jj] = qp->Hq_nzero[ii][jj];
				}
			if(arg->alias_unchanged)
				{
				if(ii<N)
					{
					qp_red->BAbt[ii] = qp->BAbt[ii];
					}
				qp_red->RSQrq[ii] = qp->RSQrq[ii];
				for(ll=0; ll<nq[ii]; ll++)
					{
					qp_red->Hq[ii][ll] = qp->Hq[ii][ll];
					}
				qp_red->Z[ii] = qp->Z[ii];
				qp_red->idxb[ii] = qp->idxb[ii];
				qp_red->DCt[ii] = qp->DCt[ii];
				qp_red->idxs_rev[ii] = qp->idxs_rev[ii];
				qp_red->idxe[ii] = qp->idxe[ii];
				}
			else
				{
				if(ii<N)
					{
					GECP(nu[ii]+nx[ii]+1, nx[ii+1], qp->BAbt+ii, 0, 0, qp_red->BAbt+ii, 0, 0);
					}
				GECP(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], qp->RSQrq+ii, 0, 0, qp_red->RSQrq+ii, 0, 0);
				for(ll=0; ll<nq[ii]; ll++)
					{
					GECP(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], qp->Hq[ii]+ll, 0, 0, qp_red->Hq[ii]+ll, 0, 0);
					}
				VECCP(2*ns[ii], qp->Z+ii, 0, qp_red->Z+ii, 0);
				for(jj=0; jj<nb[ii]; jj++)
					qp_red->idxb[ii][jj] = qp->idxb[ii][jj];
				GECP(nu[ii]+nx[ii], ng[ii]+nq[ii], qp->DCt+ii, 0, 0, qp_red->DCt+ii, 0, 0);
				for(jj=0; jj<nb[ii]+ng[ii]+nq[ii]; jj++)
					qp_red->idxs_rev[ii][jj] = qp->idxs_rev[ii][jj];
				for(jj=0; jj<nbue[ii]+nbxe[ii]+nge[ii]+nqe[ii]; jj++)
					qp_red->idxe[ii][jj] = qp->idxe[ii][jj];
				}
			}
		}

	return;

	}



void OCP_QCQP_REDUCE_EQ_DOF_LHS(struct OCP_QCQP *qp, struct OCP_QCQP *qp_red, struct OCP_QCQP_REDUCE_EQ_DOF_ARG *arg, struct OCP_QCQP_REDUCE_EQ_DOF_WS *work)
	{

	int ii, jj, kk, ll, idx0, idx1;

	struct OCP_QCQP_DIM *dim = qp->dim;
	int N = dim->N;
	int *nx = dim->nx;
	int *nu = dim->nu;
	int *nb = dim->nb;
	int *nbx = dim->nbx;
	int *nbu = dim->nbu;
	int *ng = dim->ng;
	int *nq = dim->nq;
	int *ns = dim->ns;
	int *nbue = dim->nbue;
	int *nbxe = dim->nbxe;
	int *nge = dim->nge;
	int *nqe = dim->nqe;

	struct OCP_QCQP_DIM *dim_red = qp_red->dim;
	int *nx_red = dim_red->nx;
	int *nu_red = dim_red->nu;
	int *nb_red = dim_red->nb;
	int *ng_red = dim_red->ng;
	int *nq_red = dim_red->nq;
	int *ns_red = dim_red->ns;

	// TODO handle case of softed equalities !!!!!!!!!!!!!!!!

	int ne_thr;

	for(ii=0; ii<=N; ii++)
		{
		if(ii==0)
			ne_thr = nbue[ii]+nbxe[ii];
		else
			ne_thr = nbue[ii];
		if(ne_thr>0) // reduce inputs and/or states
			{
			VECSE(nu[ii]+nx[ii], 0.0, work->tmp_nuxM+0, 0);
			for(jj=0; jj<nu[ii]+nx[ii]; jj++)
				work->e_imask_ux[jj] = 0;
			for(jj=0; jj<nbu[ii]+nbx[ii]; jj++)
				work->e_imask_d[jj] = 0;
			for(jj=0; jj<ne_thr; jj++) // set 1s for both inputs and states
				{
				VECEL(work->tmp_nuxM+0, qp->idxb[ii][qp->idxe[ii][jj]]) = VECEL(qp->d+ii, qp->idxe[ii][jj]);
				work->e_imask_ux[qp->idxb[ii][qp->idxe[ii][jj]]] = 1;
				work->e_imask_d[qp->idxe[ii][jj]] = 1;
				}
			// TODO check first and last non-zero in e_mask and only multiply between them
			if(ii<N)
				{
				// BAt
				idx0 = 0;
				for(jj=0; jj<nu[ii]+nx[ii]; jj++)
					{
					if(work->e_imask_ux[jj]==0)
						{
						GECP(1, nx[ii+1], qp->BAbt+ii, jj, 0, qp_red->BAbt+ii, idx0, 0);
						idx0++;
						}
					}
				}
			// RSQ & Hq
			idx0 = 0;
			for(jj=0; jj<nu[ii]+nx[ii]; jj++)
				{
				if(work->e_imask_ux[jj]==0)
					{
//					idx1 = 0;
					idx1 = idx0;
//					for(kk=0; kk<nu[ii]+nx[ii]; kk++)
					for(kk=jj; kk<nu[ii]+nx[ii]; kk++)
						{
						if(work->e_imask_ux[kk]==0)
							{
							MATEL(qp_red->RSQrq+ii, idx1, idx0) = MATEL(qp->RSQrq+ii, kk, jj);
							for(ll=0; ll<nq[ii]; ll++)
								{
								MATEL(qp_red->Hq[ii]+ll, idx1, idx0) = MATEL(qp->Hq[ii]+ll, kk, jj);
								}
//#endif
							idx1++;
							}
						}
					idx0++;
					}
				}
			// (at least) a sub-matrix of a zero matrix is a zero matrix
			// TODO do better Hq_nzero ???
			for(jj=0; jj<nq[ii]; jj++)
				{
				qp_red->Hq_nzero[ii][jj] = qp->Hq_nzero[ii][jj];
				}
			// gq
			for(ll=0; ll<nq[ii]; ll++)
				{
				COLEX(nu[ii]+nx[ii], qp->DCt+ii, 0, ng[ii]+ll, work->tmp_nuxM+1, 0);
				SYMV_L(nu[ii]+nx[ii], 1.0, qp->Hq[ii]+ll, 0, 0, work->tmp_nuxM+0, 0, 1.0, work->tmp_nuxM+1, 0, work->tmp_nuxM+1, 0);
//				VECEL(work->tmp_nbgqM, nb[ii]+ng[ii]+ll) = DOT(nu[ii]+nx[ii], work->tmp_nuxM+0, 0, work->tmp_nuxM+1, 0);
				idx0 = 0;
				for(jj=0; jj<nu[ii]+nx[ii]; jj++)
					{
					if(work->e_imask_ux[jj]==0)
						{
						MATEL(qp_red->DCt+ii, idx0, ng[ii]+ll) = VECEL(work->tmp_nuxM+1, jj);
						idx0++;
						}
					}
				}
			for(jj=0; jj<nq[ii]; jj++)
				qp_red->idxs_rev[ii][nb_red[ii]+ng_red[ii]+jj] = qp->idxs_rev[ii][nb[ii]+ng[ii]+jj]; // keep softed inequality constr with same slack
			// idxb idxs_rev
			idx0 = 0;
			for(jj=0; jj<nb[ii]; jj++)
				{
				if(work->e_imask_d[jj]==0)
					{
					qp_red->idxb[ii][idx0] = qp->idxb[ii][jj];
					qp_red->idxs_rev[ii][idx0] = qp->idxs_rev[ii][jj]; // keep softed inequality constr with same slack
					idx0++;
					}
				}
			// DCt
			idx0 = 0;
			for(jj=0; jj<nu[ii]+nx[ii]; jj++)
				{
				if(work->e_imask_ux[jj]==0)
					{
					GECP(1, ng[ii], qp->DCt+ii, jj, 0, qp_red->DCt+ii, idx0, 0);
					idx0++;
					}
				}
			// soft constraints
			for(jj=0; jj<ng[ii]; jj++)
				qp_red->idxs_rev[ii][nb_red[ii]+jj] = qp->idxs_rev[ii][nb[ii]+jj]; // keep softed inequality constr with same slack
			VECCP(2*ns[ii], qp->Z+ii, 0, qp_red->Z+ii, 0);
//			qp_red->diag_H_flag[ii] = qp->diag_H_flag[ii];
			// TODO idxe !!!!!!!!!!!!!!!
			}
		else // copy everything
			{
//			qp_red->diag_H_flag[ii] = qp->diag_H_flag[ii];
			for(jj=0; jj<nq[ii]; jj++)
				{
				qp_red->Hq_nzero[ii][jj] = qp->Hq_nzero[ii][jj];
				}
			if(arg->alias_unchanged)
				{
				if(ii<N)
					{
					qp_red->BAbt[ii] = qp->BAbt[ii];
					}
				qp_red->RSQrq[ii] = qp->RSQrq[ii];
				for(ll=0; ll<nq[ii]; ll++)
					{
					qp_red->Hq[ii][ll] = qp->Hq[ii][ll];
					}
				qp_red->Z[ii] = qp->Z[ii];
				qp_red->idxb[ii] = qp->idxb[ii];
				qp_red->DCt[ii] = qp->DCt[ii];
				qp_red->idxs_rev[ii] = qp->idxs_rev[ii];
				qp_red->idxe[ii] = qp->idxe[ii];
				}
			else
				{
				if(ii<N)
					{
					GECP(nu[ii]+nx[ii]+1, nx[ii+1], qp->BAbt+ii, 0, 0, qp_red->BAbt+ii, 0, 0);
					}
				GECP(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], qp->RSQrq+ii, 0, 0, qp_red->RSQrq+ii, 0, 0);
				for(ll=0; ll<nq[ii]; ll++)
					{
					GECP(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], qp->Hq[ii]+ll, 0, 0, qp_red->Hq[ii]+ll, 0, 0);
					}
				VECCP(2*ns[ii], qp->Z+ii, 0, qp_red->Z+ii, 0);
				for(jj=0; jj<nb[ii]; jj++)
					qp_red->idxb[ii][jj] = qp->idxb[ii][jj];
				GECP(nu[ii]+nx[ii], ng[ii]+nq[ii], qp->DCt+ii, 0, 0, qp_red->DCt+ii, 0, 0);
				for(jj=0; jj<nb[ii]+ng[ii]; jj++)
					qp_red->idxs_rev[ii][jj] = qp->idxs_rev[ii][jj];
				for(jj=0; jj<nbue[ii]+nbxe[ii]+nge[ii]+nqe[ii]; jj++)
					qp_red->idxe[ii][jj] = qp->idxe[ii][jj];
				}
			}
		}

	return;

	}



#if 1
void OCP_QCQP_REDUCE_EQ_DOF_RHS(struct OCP_QCQP *qp, struct OCP_QCQP *qp_red, struct OCP_QCQP_REDUCE_EQ_DOF_ARG *arg, struct OCP_QCQP_REDUCE_EQ_DOF_WS *work)
	{

	int ii, jj, kk, ll, idx0, idx1;

	struct OCP_QCQP_DIM *dim = qp->dim;
	int N = dim->N;
	int *nx = dim->nx;
	int *nu = dim->nu;
	int *nb = dim->nb;
	int *nbx = dim->nbx;
	int *nbu = dim->nbu;
	int *ng = dim->ng;
	int *nq = dim->nq;
	int *ns = dim->ns;
	int *nbue = dim->nbue;
	int *nbxe = dim->nbxe;
	int *nge = dim->nge;
	int *nqe = dim->nqe;

	struct OCP_QCQP_DIM *dim_red = qp_red->dim;
	int *nx_red = dim_red->nx;
	int *nu_red = dim_red->nu;
	int *nb_red = dim_red->nb;
	int *ng_red = dim_red->ng;
	int *nq_red = dim_red->nq;
	int *ns_red = dim_red->ns;

	// TODO handle case of softed equalities !!!!!!!!!!!!!!!!

	int ne_thr;

	for(ii=0; ii<=N; ii++)
		{
		if(ii==0)
			ne_thr = nbue[ii]+nbxe[ii];
		else
			ne_thr = nbue[ii];
		if(ne_thr>0) // reduce inputs and/or states
			{
			VECSE(nu[ii]+nx[ii], 0.0, work->tmp_nuxM+0, 0);
			for(jj=0; jj<nu[ii]+nx[ii]; jj++)
				work->e_imask_ux[jj] = 0;
			for(jj=0; jj<nbu[ii]+nbx[ii]; jj++)
				work->e_imask_d[jj] = 0;
			for(jj=0; jj<ne_thr; jj++) // set 1s for both inputs and states
				{
				VECEL(work->tmp_nuxM+0, qp->idxb[ii][qp->idxe[ii][jj]]) = VECEL(qp->d+ii, qp->idxe[ii][jj]);
				work->e_imask_ux[qp->idxb[ii][qp->idxe[ii][jj]]] = 1;
				work->e_imask_d[qp->idxe[ii][jj]] = 1;
				}
			// TODO check first and last non-zero in e_mask and only multiply between them
			if(ii<N)
				{
				// b
				GEMV_T(nu[ii]+nx[ii], nx[ii+1], 1.0, qp->BAbt+ii, 0, 0, work->tmp_nuxM+0, 0, 1.0, qp->b+ii, 0, qp_red->b+ii, 0);
				}
			// rq
			SYMV_L(nu[ii]+nx[ii], 1.0, qp->RSQrq+ii, 0, 0, work->tmp_nuxM+0, 0, 1.0, qp->rqz+ii, 0, work->tmp_nuxM+1, 0);
			idx0 = 0;
			for(jj=0; jj<nu[ii]+nx[ii]; jj++)
				{
				if(work->e_imask_ux[jj]==0)
					{
					VECEL(qp_red->rqz+ii, idx0) = VECEL(work->tmp_nuxM+1, jj);
					idx0++;
					}
				}
			// gq
			for(ll=0; ll<nq[ii]; ll++)
				{
				COLEX(nu[ii]+nx[ii], qp->DCt+ii, 0, ng[ii]+ll, work->tmp_nuxM+1, 0);
				SYMV_L(nu[ii]+nx[ii], 1.0, qp->Hq[ii]+ll, 0, 0, work->tmp_nuxM+0, 0, 1.0, work->tmp_nuxM+1, 0, work->tmp_nuxM+1, 0);
				VECEL(work->tmp_nbgqM, nb[ii]+ng[ii]+ll) = DOT(nu[ii]+nx[ii], work->tmp_nuxM+0, 0, work->tmp_nuxM+1, 0);
//				idx0 = 0;
//				for(jj=0; jj<nu[ii]+nx[ii]; jj++)
//					{
//					if(work->e_imask_ux[jj]==0)
//						{
//						MATEL(qp_red->DCt+ii, idx0, ng[ii]+ll) = VECEL(work->tmp_nuxM+1, jj);
//						idx0++;
//						}
//					}
				}
			VECCP(nq[ii], qp->d+ii, nb[ii]+ng[ii], qp_red->d+ii, nb_red[ii]+ng_red[ii]);
			VECCP(nq[ii], qp->d+ii, 2*nb[ii]+2*ng[ii]+nq[ii], qp_red->d+ii, 2*nb_red[ii]+2*ng_red[ii]+nq_red[ii]);
			VECCP(nq[ii], qp->d_mask+ii, nb[ii]+ng[ii], qp_red->d_mask+ii, nb_red[ii]+ng_red[ii]);
			VECCP(nq[ii], qp->d_mask+ii, 2*nb[ii]+2*ng[ii], qp_red->d_mask+ii, 2*nb_red[ii]+2*ng_red[ii]+nq_red[ii]);
			VECCP(nq[ii], qp->m+ii, nb[ii]+ng[ii], qp_red->m+ii, nb_red[ii]+ng_red[ii]);
			VECCP(nq[ii], qp->m+ii, 2*nb[ii]+2*ng[ii], qp_red->m+ii, 2*nb_red[ii]+2*ng_red[ii]+nq_red[ii]);
//			AXPY(nq[ii], -1.0, work->tmp_nbgqM, nb[ii]+ng[ii], qp->d+ii, nb[ii]+ng[ii], qp_red->d+ii, nb_red[ii]+ng_red[ii]);
			AXPY(nq[ii], 1.0, work->tmp_nbgqM, 2*nb[ii]+2*ng[ii]+nq[ii], qp->d+ii, 2*nb[ii]+2*ng[ii]+nq[ii], qp_red->d+ii, 2*nb_red[ii]+2*ng_red[ii]+nq_red[ii]);
			for(jj=0; jj<nq[ii]; jj++)
				qp_red->idxs_rev[ii][nb_red[ii]+ng_red[ii]+jj] = qp->idxs_rev[ii][nb[ii]+ng[ii]+jj]; // keep softed inequality constr with same slack
			// d d_mask m idxb idxs_rev
			idx0 = 0;
			for(jj=0; jj<nb[ii]; jj++)
				{
				if(work->e_imask_d[jj]==0)
					{
					VECEL(qp_red->d+ii, idx0) = VECEL(qp->d+ii, jj);
					VECEL(qp_red->d+ii, nb_red[ii]+ng_red[ii]+nq_red[ii]+idx0) = VECEL(qp->d+ii, nb[ii]+ng[ii]+nq[ii]+jj);
					VECEL(qp_red->d_mask+ii, idx0) = VECEL(qp->d_mask+ii, jj);
					VECEL(qp_red->d_mask+ii, nb_red[ii]+ng_red[ii]+nq_red[ii]+idx0) = VECEL(qp->d_mask+ii, nb[ii]+ng[ii]+nq[ii]+jj);
					VECEL(qp_red->m+ii, idx0) = VECEL(qp->m+ii, jj);
					VECEL(qp_red->m+ii, nb_red[ii]+ng_red[ii]+nq_red[ii]+idx0) = VECEL(qp->m+ii, nb[ii]+ng[ii]+nq[ii]+jj);
					idx0++;
					}
				}
			VECCP(ng[ii], qp->d+ii, nb[ii], qp_red->d+ii, nb_red[ii]);
			VECCP(ng[ii]+2*ns[ii], qp->d+ii, 2*nb[ii]+ng[ii]+nq[ii], qp_red->d+ii, 2*nb_red[ii]+ng_red[ii]+nq_red[ii]);
			VECCP(ng[ii], qp->d_mask+ii, nb[ii], qp_red->d_mask+ii, nb_red[ii]);
			VECCP(ng[ii]+2*ns[ii], qp->d_mask+ii, 2*nb[ii]+ng[ii]+nq[ii], qp_red->d_mask+ii, 2*nb_red[ii]+ng_red[ii]+nq_red[ii]);
			VECCP(ng[ii], qp->m+ii, nb[ii], qp_red->m+ii, nb_red[ii]);
			VECCP(ng[ii]+2*ns[ii], qp->m+ii, 2*nb[ii]+ng[ii]+nq[ii], qp_red->m+ii, 2*nb_red[ii]+ng_red[ii]+nq_red[ii]);
			GEMV_T(nu[ii]+nx[ii], ng[ii], 1.0, qp->DCt+ii, 0, 0, work->tmp_nuxM+0, 0, 0.0, work->tmp_nbgqM, 0, work->tmp_nbgqM, 0);
			AXPY(ng[ii], -1.0, work->tmp_nbgqM, 0, qp->d+ii, nb[ii], qp_red->d+ii, nb_red[ii]);
			AXPY(ng[ii], 1.0, work->tmp_nbgqM, 0, qp->d+ii, 2*nb[ii]+ng[ii]+nq[ii], qp_red->d+ii, 2*nb_red[ii]+ng_red[ii]+nq_red[ii]);
			// soft constraints
			VECCP(2*ns[ii], qp->rqz+ii, nu[ii]+nx[ii], qp_red->rqz+ii, nu_red[ii]+nx_red[ii]);
			// TODO idxe !!!!!!!!!!!!!!!
			}
		else // copy everything
			{
			// copy vectors which are contiguous in the QP (e.g. to alias to res)
			if(ii<N)
				{
				VECCP(nx[ii+1], qp->b+ii, 0, qp_red->b+ii, 0);
				}
			VECCP(nu[ii]+nx[ii]+2*ns[ii], qp->rqz+ii, 0, qp_red->rqz+ii, 0);
			VECCP(2*nb[ii]+2*ng[ii]+2*nq[ii]+2*ns[ii], qp->d+ii, 0, qp_red->d+ii, 0);
			VECCP(2*nb[ii]+2*ng[ii]+2*nq[ii]+2*ns[ii], qp->d_mask+ii, 0, qp_red->d_mask+ii, 0);
			VECCP(2*nb[ii]+2*ng[ii]+2*nq[ii]+2*ns[ii], qp->m+ii, 0, qp_red->m+ii, 0);
			}
		}

	return;

	}
#endif



void OCP_QCQP_RESTORE_EQ_DOF(struct OCP_QCQP *qp, struct OCP_QCQP_SOL *qp_sol_red, struct OCP_QCQP_SOL *qp_sol, struct OCP_QCQP_REDUCE_EQ_DOF_ARG *arg, struct OCP_QCQP_REDUCE_EQ_DOF_WS *work)
	{

	int ii, jj, idx0;

	struct OCP_QCQP_DIM *dim = qp->dim;
	int N = dim->N;
	int *nx = dim->nx;
	int *nu = dim->nu;
	int *nb = dim->nb;
	int *nbx = dim->nbx;
	int *nbu = dim->nbu;
	int *ng = dim->ng;
	int *nq = dim->nq;
	int *ns = dim->ns;
	int *nbue = dim->nbue;
	int *nbxe = dim->nbxe;
	int *nge = dim->nge;
	int *nqe = dim->nqe;

	struct OCP_QCQP_DIM *dim_red = qp_sol_red->dim;
	int *nx_red = dim_red->nx;
	int *nu_red = dim_red->nu;
	int *nb_red = dim_red->nb;
	int *ng_red = dim_red->ng;
	int *nq_red = dim_red->nq;

	int ne_thr;

	REAL tmp;

	for(ii=0; ii<=N; ii++)
		{
		if(ii==0)
			ne_thr = nbue[ii]+nbxe[ii];
		else
			ne_thr = nbue[ii];
		if(ne_thr>0) // restore inputs and/or states
			{
			VECSE(nu[ii]+nx[ii], 0.0, work->tmp_nuxM+0, 0);
			for(jj=0; jj<nu[ii]+nx[ii]; jj++)
				work->e_imask_ux[jj] = 0;
			for(jj=0; jj<nb[ii]; jj++)
				work->e_imask_d[jj] = 0;
			for(jj=0; jj<ne_thr; jj++) // set 1s for both inputs and states
				{
				work->e_imask_ux[qp->idxb[ii][qp->idxe[ii][jj]]] = 1;
				work->e_imask_d[qp->idxe[ii][jj]] = 1;
				}
			// ux
			if(arg->comp_prim_sol)
				{
				idx0 = 0;
				for(jj=0; jj<nu[ii]+nx[ii]; jj++)
					{
					if(work->e_imask_ux[jj]==0)
						{
						VECEL(qp_sol->ux+ii, jj) = VECEL(qp_sol_red->ux+ii, idx0);
						idx0++;
						}
					}
				for(jj=0; jj<ne_thr; jj++)
					{
					VECEL(qp_sol->ux+ii, qp->idxb[ii][qp->idxe[ii][jj]]) = VECEL(qp->d+ii, qp->idxe[ii][jj]);
					}
				// TODO update based on choices on reduce !!!!!!!!!!!!!
				VECCP(2*ns[ii], qp_sol_red->ux+ii, nu_red[ii]+nx_red[ii], qp_sol->ux+ii, nu[ii]+nx[ii]);
				}
			if(arg->comp_dual_sol_eq)
				{
				// pi
				if(ii<N)
					VECCP(nx[ii+1], qp_sol_red->pi+ii, 0, qp_sol->pi+ii, 0);
				}
			if(arg->comp_dual_sol_ineq)
				{
				// lam t
				idx0 = 0;
				for(jj=0; jj<nb[ii]; jj++)
					{
					if(work->e_imask_d[jj]==0)
						{
						VECEL(qp_sol->lam+ii, jj) = VECEL(qp_sol_red->lam+ii, idx0);
						VECEL(qp_sol->lam+ii, nb[ii]+ng[ii]+nq[ii]+jj) = VECEL(qp_sol_red->lam+ii, nb_red[ii]+ng_red[ii]+nq_red[ii]+idx0);
						VECEL(qp_sol->t+ii, jj) = VECEL(qp_sol_red->t+ii, idx0);
						VECEL(qp_sol->t+ii, nb[ii]+ng[ii]+nq[ii]+jj) = VECEL(qp_sol_red->t+ii, nb_red[ii]+ng_red[ii]+nq_red[ii]+idx0);
						idx0++;
						}
					else
						{
						// lam
						// t
						VECEL(qp_sol->lam+ii, jj) = arg->lam_min;
						VECEL(qp_sol->lam+ii, nb[ii]+ng[ii]+nq[ii]+jj) = arg->lam_min;
						VECEL(qp_sol->t+ii, jj) = arg->t_min;
						VECEL(qp_sol->t+ii, nb[ii]+ng[ii]+nq[ii]+jj) = arg->t_min;
						}
					}
				// TODO update based on choices on reduce !!!!!!!!!!!!!
				VECCP(ng[ii]+nq[ii], qp_sol_red->lam+ii, nb_red[ii], qp_sol->lam+ii, nb[ii]);
				VECCP(ng[ii]+nq[ii]+2*ns[ii], qp_sol_red->lam+ii, 2*nb_red[ii]+ng_red[ii]+nq_red[ii], qp_sol->lam+ii, 2*nb[ii]+ng[ii]+nq[ii]);
				VECCP(ng[ii]+nq[ii], qp_sol_red->t+ii, nb_red[ii], qp_sol->t+ii, nb[ii]);
				VECCP(ng[ii]+nq[ii]+2*ns[ii], qp_sol_red->t+ii, 2*nb_red[ii]+ng_red[ii]+nq_red[ii], qp_sol->t+ii, 2*nb[ii]+ng[ii]+nq[ii]);
				// update lam_lb for removed eq if > 0, or lam_ub if < 0 //, keep lam_ub to zero
				VECCP(nu[ii]+nx[ii], qp->rqz+ii, 0, work->tmp_nuxM, 0);
				AXPY(nb[ii]+ng[ii]+nq[ii], -1.0, qp_sol->lam+ii, 0, qp_sol->lam+ii, nb[ii]+ng[ii]+nq[ii], work->tmp_nbgqM, 0);
				VECAD_SP(nb[ii], 1.0, work->tmp_nbgqM, 0, qp->idxb[ii], work->tmp_nuxM, 0);
				SYMV_L(nu[ii]+nx[ii], 1.0, qp->RSQrq+ii, 0, 0, qp_sol->ux+ii, 0, 1.0, work->tmp_nuxM, 0, work->tmp_nuxM, 0);
				if(ii<N)
					GEMV_N(nu[ii]+nx[ii], nx[ii+1], 1.0, qp->BAbt+ii, 0, 0, qp_sol_red->pi+ii, 0, 1.0, work->tmp_nuxM, 0, work->tmp_nuxM, 0);
				GEMV_N(nu[ii]+nx[ii], ng[ii], 1.0, qp->DCt+ii, 0, 0, work->tmp_nbgqM, nb[ii], 1.0, work->tmp_nuxM, 0, work->tmp_nuxM, 0);
				for(jj=0; jj<nq[ii]; jj++)
					{
					COLEX(nu[ii]+nx[ii], qp->DCt+ii, 0, ng[ii]+jj, work->tmp_nuxM+1, 0);
					SYMV_L(nu[ii]+nx[ii], 1.0, qp->Hq[ii]+jj, 0, 0, qp_sol->ux+ii, 0, 1.0, work->tmp_nuxM+1, 0, work->tmp_nuxM+1, 0);
					tmp = VECEL(work->tmp_nbgqM, nb[ii]+ng[ii]+jj); 
					AXPY(nu[ii]+nx[ii], tmp, work->tmp_nuxM+1, 0, work->tmp_nuxM+0, 0, work->tmp_nuxM+0, 0);
					}
				for(jj=0; jj<nb[ii]; jj++)
					{
					if(work->e_imask_d[jj]!=0)
						{
						tmp = VECEL(work->tmp_nuxM, qp->idxb[ii][jj]);
						if(tmp>=0)
							VECEL(qp_sol->lam+ii, jj) = tmp;
						else
							VECEL(qp_sol->lam+ii, nb[ii]+ng[ii]+nq[ii]+jj) = - tmp;
						}
					}
				}
			}
		else // copy
			{
			if(arg->comp_prim_sol)
				{
				// ux
				VECCP(nu[ii]+nx[ii]+2*ns[ii], qp_sol_red->ux+ii, 0, qp_sol->ux+ii, 0);
				}
			if(arg->comp_dual_sol_eq)
				{
				// pi
				if(ii<N)
					VECCP(nx[ii+1], qp_sol_red->pi+ii, 0, qp_sol->pi+ii, 0);
				}
			if(arg->comp_dual_sol_ineq)
				{
				// lam
				VECCP(2*nb[ii]+2*ng[ii]+2*nq[ii]+2*ns[ii], qp_sol_red->lam+ii, 0, qp_sol->lam+ii, 0);
				// t
				VECCP(2*nb[ii]+2*ng[ii]+2*nq[ii]+2*ns[ii], qp_sol_red->t+ii, 0, qp_sol->t+ii, 0);
				}
			}
		}

	return;

	}



void OCP_QCQP_REDUCE_EQ_DOF_SOL(struct OCP_QCQP *qp, struct OCP_QCQP_SOL *qp_sol, struct OCP_QCQP_SOL *qp_sol_red, struct OCP_QCQP_REDUCE_EQ_DOF_ARG *arg, struct OCP_QCQP_REDUCE_EQ_DOF_WS *work)
	{

	int ii, jj, idx0;

	struct OCP_QCQP_DIM *dim = qp->dim;
	int N = dim->N;
	int *nx = dim->nx;
	int *nu = dim->nu;
	int *nb = dim->nb;
	int *nbx = dim->nbx;
	int *nbu = dim->nbu;
	int *ng = dim->ng;
	int *nq = dim->nq;
	int *ns = dim->ns;
	int *nbue = dim->nbue;
	int *nbxe = dim->nbxe;
	int *nge = dim->nge;
	int *nqe = dim->nqe;

	struct OCP_QCQP_DIM *dim_red = qp_sol_red->dim;
	int *nx_red = dim_red->nx;
	int *nu_red = dim_red->nu;
	int *nb_red = dim_red->nb;
	int *ng_red = dim_red->ng;
	int *nq_red = dim_red->nq;

	int ne_thr;

	REAL tmp;

	for(ii=0; ii<=N; ii++)
		{
		if(ii==0)
			ne_thr = nbue[ii]+nbxe[ii];
		else
			ne_thr = nbue[ii];
		if(ne_thr>0) // restore inputs and/or states
			{
			VECSE(nu[ii]+nx[ii], 0.0, work->tmp_nuxM+0, 0);
			for(jj=0; jj<nu[ii]+nx[ii]; jj++)
				work->e_imask_ux[jj] = 0;
			for(jj=0; jj<nb[ii]; jj++)
				work->e_imask_d[jj] = 0;
			for(jj=0; jj<ne_thr; jj++) // set 1s for both inputs and states
				{
				work->e_imask_ux[qp->idxb[ii][qp->idxe[ii][jj]]] = 1;
				work->e_imask_d[qp->idxe[ii][jj]] = 1;
				}
			// ux
			if(arg->comp_prim_sol)
				{
				idx0 = 0;
				for(jj=0; jj<nu[ii]+nx[ii]; jj++)
					{
					if(work->e_imask_ux[jj]==0)
						{
						VECEL(qp_sol_red->ux+ii, idx0) = VECEL(qp_sol->ux+ii, jj);
						idx0++;
						}
					}
				// TODO update based on choices on reduce !!!!!!!!!!!!!
				VECCP(2*ns[ii], qp_sol->ux+ii, nu[ii]+nx[ii], qp_sol_red->ux+ii, nu_red[ii]+nx_red[ii]);
				}
			if(arg->comp_dual_sol_eq)
				{
				// pi
				if(ii<N)
					VECCP(nx[ii+1], qp_sol->pi+ii, 0, qp_sol_red->pi+ii, 0);
				}
			if(arg->comp_dual_sol_ineq)
				{
				// lam t
				idx0 = 0;
				for(jj=0; jj<nb[ii]; jj++)
					{
					if(work->e_imask_d[jj]==0)
						{
						VECEL(qp_sol_red->lam+ii, idx0) = VECEL(qp_sol->lam+ii, jj);
						VECEL(qp_sol_red->lam+ii, nb_red[ii]+ng_red[ii]+nq_red[ii]+idx0) = VECEL(qp_sol->lam+ii, nb[ii]+ng[ii]+nq[ii]+jj);
						VECEL(qp_sol_red->t+ii, idx0) = VECEL(qp_sol->t+ii, jj);
						VECEL(qp_sol_red->t+ii, nb_red[ii]+ng_red[ii]+nq_red[ii]+idx0) = VECEL(qp_sol->t+ii, nb[ii]+ng[ii]+nq[ii]+jj);
						idx0++;
						}
					}
				// TODO update based on choices on reduce !!!!!!!!!!!!!
				VECCP(ng[ii]+nq[ii], qp_sol->lam+ii, nb[ii], qp_sol_red->lam+ii, nb_red[ii]);
				VECCP(ng[ii]+nq[ii]+2*ns[ii], qp_sol->lam+ii, 2*nb[ii]+ng[ii]+nq[ii], qp_sol_red->lam+ii, 2*nb_red[ii]+ng_red[ii]+nq_red[ii]);
				VECCP(ng[ii]+nq[ii], qp_sol->t+ii, nb[ii], qp_sol_red->t+ii, nb_red[ii]);
				VECCP(ng[ii]+nq[ii]+2*ns[ii], qp_sol->t+ii, 2*nb[ii]+ng[ii]+nq[ii], qp_sol_red->t+ii, 2*nb_red[ii]+ng_red[ii]+nq_red[ii]);
				}
			}
		else // copy
			{
			if(arg->comp_prim_sol)
				{
				// ux
				VECCP(nu[ii]+nx[ii]+2*ns[ii], qp_sol->ux+ii, 0, qp_sol_red->ux+ii, 0);
				}
			if(arg->comp_dual_sol_eq)
				{
				// pi
				if(ii<N)
					VECCP(nx[ii+1], qp_sol->pi+ii, 0, qp_sol_red->pi+ii, 0);
				}
			if(arg->comp_dual_sol_ineq)
				{
				// lam
				VECCP(2*nb[ii]+2*ng[ii]+2*nq[ii]+2*ns[ii], qp_sol->lam+ii, 0, qp_sol_red->lam+ii, 0);
				// t
				VECCP(2*nb[ii]+2*ng[ii]+2*nq[ii]+2*ns[ii], qp_sol->t+ii, 0, qp_sol_red->t+ii, 0);
				}
			}
		}

	return;

	}



