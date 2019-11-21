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

void COND_QP_COMPUTE_DIM(struct OCP_QP_DIM *ocp_dim, struct DENSE_QP_DIM *dense_dim)
	{

	int N = ocp_dim->N;
	int *nx = ocp_dim->nx;
	int *nu = ocp_dim->nu;
	int *nbx = ocp_dim->nbx;
	int *nbu = ocp_dim->nbu;
	int *ng = ocp_dim->ng;
	int *ns = ocp_dim->ns;
	int *nsbx = ocp_dim->nsbx;
	int *nsbu = ocp_dim->nsbu;
	int *nsg = ocp_dim->nsg;

	int ii;

	int nvc = 0;
	int nec = 0;
	int nbc = 0;
	int ngc = 0;
	int nsc = 0;
	int nsbc = 0;
	int nsgc = 0;

	// first stage
	nvc += nx[0]+nu[0];
	nbc += nbx[0]+nbu[0];
	ngc += ng[0];
	nsc += ns[0];
	nsbc += nsbx[0]+nsbu[0];
	nsgc += nsg[0];
	// remaining stages
	for(ii=1; ii<=N; ii++)
		{
		nvc += nu[ii];
		nbc += nbu[ii];
		ngc += nbx[ii]+ng[ii];
		nsc += ns[ii];
		nsbc += nsbu[ii];
		nsgc += nsbx[ii]+nsg[ii];
		}

	dense_dim->nv = nvc;
	dense_dim->ne = nec;
	dense_dim->nb = nbc;
	dense_dim->ng = ngc;
	dense_dim->ns = nsc;
	dense_dim->nsb = nsbc;
	dense_dim->nsg = nsgc;

	return;

	}



int COND_QP_ARG_MEMSIZE()
	{

	int size = 0;

	return size;

	}



void COND_QP_ARG_CREATE(struct COND_QP_ARG *cond_arg, void *mem)
	{

	cond_arg->memsize = COND_QP_ARG_MEMSIZE();

	return;

	}



void COND_QP_ARG_SET_DEFAULT(struct COND_QP_ARG *cond_arg)
	{

	cond_arg->cond_last_stage = 1; // condense last stage
	cond_arg->comp_dual_sol = 1; // compute dual solution
	cond_arg->square_root_alg = 1; // square root algorithm (faster but requires RSQ>0)

	return;

	}



void COND_QP_ARG_SET_RIC_ALG(int ric_alg, struct COND_QP_ARG *cond_arg)
	{

	cond_arg->square_root_alg = ric_alg;

	return;

	}



int COND_QP_WS_MEMSIZE(struct OCP_QP_DIM *ocp_dim, struct COND_QP_ARG *cond_arg)
	{

	int ii;

	int N = ocp_dim->N;
	int *nx = ocp_dim->nx;
	int *nu = ocp_dim->nu;
	int *nb = ocp_dim->nb;
	int *ng = ocp_dim->ng;

	// compute core qp size and max size
	int nvt = 0;
	int net = 0;
	int nbt = 0;
	int ngt = 0;
	int nxM = 0;
	int nuM = 0;
	int nbM = 0;
	int ngM = 0;

	for(ii=0; ii<N; ii++)
		{
		nvt += nx[ii]+nu[ii];
		net += nx[ii+1];
		nbt += nb[ii];
		ngt += ng[ii];
		nxM = nx[ii]>nxM ? nx[ii] : nxM;
		nuM = nu[ii]>nuM ? nu[ii] : nuM;
		nbM = nb[ii]>nbM ? nb[ii] : nbM;
		ngM = ng[ii]>ngM ? ng[ii] : ngM;
		}
	ii = N;
	nvt += nx[ii]+nu[ii];
	nbt += nb[ii];
	ngt += ng[ii];
	nxM = nx[ii]>nxM ? nx[ii] : nxM;
	nuM = nu[ii]>nuM ? nu[ii] : nuM;
	nbM = nb[ii]>nbM ? nb[ii] : nbM;
	ngM = ng[ii]>ngM ? ng[ii] : ngM;

	int size = 0;

	size += 2*(N+1)*sizeof(struct STRMAT); // Gamma L
	size += 2*sizeof(struct STRMAT); // Lx AL
	size += 2*(N+1)*sizeof(struct STRVEC); // Gammab l
	size += 2*sizeof(struct STRVEC); // tmp_nbgM tmp_nuxM

	int nu_tmp = 0;
	for(ii=0; ii<N; ii++)
		{
		nu_tmp += nu[ii];
		size += SIZE_STRMAT(nu_tmp+nx[0]+1, nx[ii+1]); // Gamma
		}
	for(ii=0; ii<=N; ii++)
		size += SIZE_STRMAT(nu[ii]+nx[ii]+1, nu[ii]+nx[ii]); // L
	size += SIZE_STRMAT(nxM+1, nxM); // Lx
	size += SIZE_STRMAT(nuM+nxM+1, nxM); // AL
	for(ii=0; ii<N; ii++)
		size += 1*SIZE_STRVEC(nx[ii+1]); // Gammab
	for(ii=0; ii<=N; ii++)
		size += SIZE_STRVEC(nu[ii]+nx[ii]); // l
	size += 1*SIZE_STRVEC(nbM+ngM); // tmp_nbgM
	size += 1*SIZE_STRVEC(nuM+nxM); // tmp_nuxM

	size += 1*(nbM+ngM)*sizeof(int); // idxs_rev

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

	}



void COND_QP_WS_CREATE(struct OCP_QP_DIM *ocp_dim, struct COND_QP_ARG *cond_arg, struct COND_QP_WS *cond_ws, void *mem)
	{

	int ii;

	int N = ocp_dim->N;
	int *nx = ocp_dim->nx;
	int *nu = ocp_dim->nu;
	int *nb = ocp_dim->nb;
	int *ng = ocp_dim->ng;

	// compute core qp dim and max dim
	int nvt = 0;
	int net = 0;
	int nbt = 0;
	int ngt = 0;
	int nxM = 0;
	int nuM = 0;
	int nbM = 0;
	int ngM = 0;

	for(ii=0; ii<N; ii++)
		{
		nvt += nx[ii]+nu[ii];
		net += nx[ii+1];
		nbt += nb[ii];
		ngt += ng[ii];
		nxM = nx[ii]>nxM ? nx[ii] : nxM;
		nuM = nu[ii]>nuM ? nu[ii] : nuM;
		nbM = nb[ii]>nbM ? nb[ii] : nbM;
		ngM = ng[ii]>ngM ? ng[ii] : ngM;
		}
	ii = N;
	nvt += nx[ii]+nu[ii];
	nbt += nb[ii];
	ngt += ng[ii];
	nxM = nx[ii]>nxM ? nx[ii] : nxM;
	nuM = nu[ii]>nuM ? nu[ii] : nuM;
	nbM = nb[ii]>nbM ? nb[ii] : nbM;
	ngM = ng[ii]>ngM ? ng[ii] : ngM;


	// matrix struct
	struct STRMAT *sm_ptr = (struct STRMAT *) mem;

	cond_ws->Gamma = sm_ptr;
	sm_ptr += N+1;
	cond_ws->L = sm_ptr;
	sm_ptr += N+1;
	cond_ws->Lx = sm_ptr;
	sm_ptr += 1;
	cond_ws->AL = sm_ptr;
	sm_ptr += 1;


	// vector struct
	struct STRVEC *sv_ptr = (struct STRVEC *) sm_ptr;

	cond_ws->Gammab = sv_ptr;
	sv_ptr += N+1;
	cond_ws->l = sv_ptr;
	sv_ptr += N+1;
	cond_ws->tmp_nbgM = sv_ptr;
	sv_ptr += 1;
	cond_ws->tmp_nuxM = sv_ptr;
	sv_ptr += 1;


	// int stuff
	int *i_ptr;
	i_ptr = (int *) sv_ptr;

	// idxs_rev
	cond_ws->idxs_rev = i_ptr;
	i_ptr += nbM+ngM;


	// align to typicl cache line size
	size_t s_ptr = (size_t) i_ptr;
	s_ptr = (s_ptr+63)/64*64;


	// void stuf
	char *c_ptr = (char *) s_ptr;
	char *c_tmp;

	int nu_tmp = 0;
	for(ii=0; ii<N; ii++)
		{
		nu_tmp += nu[ii];
		CREATE_STRMAT(nu_tmp+nx[0]+1, nx[ii+1], cond_ws->Gamma+ii, c_ptr);
		c_ptr += (cond_ws->Gamma+ii)->memsize;
		}
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRMAT(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], cond_ws->L+ii, c_ptr);
		c_ptr += (cond_ws->L+ii)->memsize;
		}
	CREATE_STRMAT(nxM+1, nxM, cond_ws->Lx, c_ptr);
	c_ptr += cond_ws->Lx->memsize;
	CREATE_STRMAT(nuM+nxM+1, nxM, cond_ws->AL, c_ptr);
	c_ptr += cond_ws->AL->memsize;
	for(ii=0; ii<N; ii++)
		{
		CREATE_STRVEC(nx[ii+1], cond_ws->Gammab+ii, c_ptr);
		c_ptr += (cond_ws->Gammab+ii)->memsize;
		}
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(nu[ii]+nx[ii], cond_ws->l+ii, c_ptr);
		c_ptr += (cond_ws->l+ii)->memsize;
		}
	CREATE_STRVEC(nbM+ngM, cond_ws->tmp_nbgM, c_ptr);
	c_ptr += cond_ws->tmp_nbgM->memsize;
	c_tmp = c_ptr;
	CREATE_STRVEC(nuM+nxM, cond_ws->tmp_nuxM, c_ptr);
	c_ptr += cond_ws->tmp_nuxM->memsize;

	cond_ws->bs = N;

	cond_ws->memsize = COND_QP_WS_MEMSIZE(ocp_dim, cond_arg);

#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + cond_ws->memsize)
		{
		printf("\nCreate_cond_qp_ocp2dense: outsize memory bounds!\n\n");
		exit(1);
		}
#endif

	return;

	}



void COND_QP_COND(struct OCP_QP *ocp_qp, struct DENSE_QP *dense_qp, struct COND_QP_ARG *cond_arg, struct COND_QP_WS *cond_ws)
	{

	COND_BABT(ocp_qp, NULL, NULL, cond_arg, cond_ws);

	COND_RSQRQ_N2NX3(ocp_qp, dense_qp->Hv, dense_qp->gz, cond_arg, cond_ws);

	COND_DCTD(ocp_qp, dense_qp->idxb, dense_qp->Ct, dense_qp->d, dense_qp->idxs, dense_qp->Z, dense_qp->gz, cond_arg, cond_ws);

	return;

	}



void COND_QP_COND_RHS(struct OCP_QP *ocp_qp, struct DENSE_QP *dense_qp, struct COND_QP_ARG *cond_arg, struct COND_QP_WS *cond_ws)
	{

	COND_B(ocp_qp, NULL, cond_arg, cond_ws);

	COND_RQ_N2NX3(ocp_qp, dense_qp->gz, cond_arg, cond_ws);

	COND_D(ocp_qp, dense_qp->d, dense_qp->gz, cond_arg, cond_ws);

	return;

	}



void COND_QP_EXPAND_SOL(struct OCP_QP *ocp_qp, struct DENSE_QP_SOL *dense_qp_sol, struct OCP_QP_SOL *ocp_qp_sol, struct COND_QP_ARG *cond_arg, struct COND_QP_WS *cond_ws)
	{

	EXPAND_SOL(ocp_qp, dense_qp_sol, ocp_qp_sol, cond_arg, cond_ws);

	return;

	}



// TODO remove
void COND_QP_EXPAND_PRIMAL_SOL(struct OCP_QP *ocp_qp, struct DENSE_QP_SOL *dense_qp_sol, struct OCP_QP_SOL *ocp_qp_sol, struct COND_QP_ARG *cond_arg, struct COND_QP_WS *cond_ws)
	{

	EXPAND_PRIMAL_SOL(ocp_qp, dense_qp_sol, ocp_qp_sol, cond_arg, cond_ws);

	return;

	}



/************************************************
* update cond
************************************************/

void COND_QP_UPDATE(int *idxc, struct OCP_QP *ocp_qp, struct DENSE_QP *dense_qp, struct COND_QP_ARG *cond_arg, struct COND_QP_WS *cond_ws)
	{

	UPDATE_COND_BABT(idxc, ocp_qp, NULL, NULL, cond_arg, cond_ws);

	UPDATE_COND_RSQRQ_N2NX3(idxc, ocp_qp, dense_qp->Hv, dense_qp->gz, cond_arg, cond_ws);

	UPDATE_COND_DCTD(idxc, ocp_qp, dense_qp->idxb, dense_qp->Ct, dense_qp->d, dense_qp->idxs, dense_qp->Z, dense_qp->gz, cond_arg, cond_ws);

	return;

	}




