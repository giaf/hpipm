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

void COMPUTE_QP_SIZE_OCP2DENSE(int N, int *nx, int *nu, int *nb, int **idxb, int *ng, int *ns, int *nvc, int *nec, int *nbc, int *ngc, int *nsc)
	{

	int ii, jj;

	nvc[0] = 0;
	nec[0] = 0;
	nbc[0] = 0;
	ngc[0] = 0;
	nsc[0] = 0;

	// first stage
	nvc[0] += nx[0]+nu[0];
	nbc[0] += nb[0];
	ngc[0] += ng[0];
	nsc[0] += ns[0];
	// remaining stages
	for(ii=1; ii<=N; ii++)
		{
		nvc[0] += nu[ii];
		for(jj=0; jj<nb[ii]; jj++)
			{
			if(idxb[ii][jj]<nu[ii]) // input constraint
				{
				nbc[0]++;
				}
			else // state constraint
				{
				ngc[0]++;
				}
			}
		ngc[0] += ng[ii];
		nsc[0] += ns[ii];
		}

	return;

	}



int MEMSIZE_COND_QP_OCP2DENSE(struct OCP_QP *ocp_qp, struct DENSE_QP *dense_qp) // XXX + args for algorithm type ???
	{

	int ii;

	int N = ocp_qp->N;
	int *nx = ocp_qp->nx;
	int *nu = ocp_qp->nu;
	int *nb = ocp_qp->nb;
	int *ng = ocp_qp->ng;

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
	size += 2*sizeof(struct STRVEC); // tmp_ngM tmp_nuxM

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
	size += SIZE_STRVEC(ngM); // tmp_ngM
	size += 1*SIZE_STRVEC(nuM+nxM); // tmp_nuxM
	size += 1*SIZE_STRVEC(ngM); // tmp_ngM

	size += 1*(nbM+ngM)*sizeof(int); // idxs_rev

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

	}



void CREATE_COND_QP_OCP2DENSE(struct OCP_QP *ocp_qp, struct DENSE_QP *dense_qp, struct COND_QP_OCP2DENSE_WORKSPACE *cond_ws, void *mem)
	{

	int ii;

	int N = ocp_qp->N;
	int *nx = ocp_qp->nx;
	int *nu = ocp_qp->nu;
	int *nb = ocp_qp->nb;
	int *ng = ocp_qp->ng;

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
	cond_ws->tmp_ngM = sv_ptr;
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
		c_ptr += (cond_ws->Gamma+ii)->memory_size;
		}
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRMAT(nu[ii]+nx[ii]+1, nu[ii]+nx[ii], cond_ws->L+ii, c_ptr);
		c_ptr += (cond_ws->L+ii)->memory_size;
		}
	CREATE_STRMAT(nxM+1, nxM, cond_ws->Lx, c_ptr);
	c_ptr += cond_ws->Lx->memory_size;
	CREATE_STRMAT(nuM+nxM+1, nxM, cond_ws->AL, c_ptr);
	c_ptr += cond_ws->AL->memory_size;
	for(ii=0; ii<N; ii++)
		{
		CREATE_STRVEC(nx[ii+1], cond_ws->Gammab+ii, c_ptr);
		c_ptr += (cond_ws->Gammab+ii)->memory_size;
		}
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(nu[ii]+nx[ii], cond_ws->l+ii, c_ptr);
		c_ptr += (cond_ws->l+ii)->memory_size;
		}
	CREATE_STRVEC(ngM, cond_ws->tmp_ngM, c_ptr);
	c_ptr += cond_ws->tmp_ngM->memory_size;
	c_tmp = c_ptr;
	CREATE_STRVEC(nuM+nxM, cond_ws->tmp_nuxM, c_ptr);
	c_ptr += cond_ws->tmp_nuxM->memory_size;

	cond_ws->cond_last_stage = 1; // default: cond last stage

	cond_ws->memsize = MEMSIZE_COND_QP_OCP2DENSE(ocp_qp, dense_qp);

#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + cond_ws->memsize)
		{
		printf("\nCreate_cond_qp_ocp2dense: outsize memory bounds!\n\n");
		exit(1);
		}
#endif

	return;

	}

	

void COND_QP_OCP2DENSE(struct OCP_QP *ocp_qp, struct DENSE_QP *dense_qp, struct COND_QP_OCP2DENSE_WORKSPACE *cond_ws)
	{

	COND_BABT(ocp_qp, NULL, NULL, cond_ws);

	COND_RSQRQ_N2NX3(ocp_qp, dense_qp->Hg, dense_qp->g, cond_ws);

	COND_DCTD(ocp_qp, dense_qp->idxb, dense_qp->Ct, dense_qp->d, dense_qp->idxs, dense_qp->Z, dense_qp->z, cond_ws);

	return;

	}



void COND_RHS_QP_OCP2DENSE(struct OCP_QP *ocp_qp, struct DENSE_QP *dense_qp, struct COND_QP_OCP2DENSE_WORKSPACE *cond_ws)
	{

	COND_B(ocp_qp, NULL, cond_ws);

	COND_RQ_N2NX3(ocp_qp, dense_qp->g, cond_ws);

	COND_D(ocp_qp, dense_qp->d, dense_qp->z, cond_ws);

	return;

	}



void EXPAND_SOL_DENSE2OCP(struct OCP_QP *ocp_qp, struct DENSE_QP_SOL *dense_qp_sol, struct OCP_QP_SOL *ocp_qp_sol, struct COND_QP_OCP2DENSE_WORKSPACE *cond_ws)
	{

	EXPAND_SOL(ocp_qp, dense_qp_sol, ocp_qp_sol, cond_ws);

	return;

	}
