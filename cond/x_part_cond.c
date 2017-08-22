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

void COMPUTE_QP_SIZE_OCP2OCP(int N, int *nx, int *nu, int *nb, int **idxb, int *ng, int *ns, int N2, int *nx2, int *nu2, int *nb2, int *ng2, int *ns2)
	{

	int ii, jj, kk;

	int N1 = N/N2; // (floor) horizon of small blocks
	int R1 = N - N2*N1; // the first R1 blocks have horizon N1+1
	int M1 = R1>0 ? N1+1 : N1; // (ceil) horizon of large blocks
	int T1; // horizon of current block

	int N_tmp = 0; // temporary sum of horizons
	int nbb; // box constr that remain box constr
	int nbg; // box constr that becomes general constr
	for(ii=0; ii<N2; ii++)
		{
		T1 = ii<R1 ? M1 : N1;
		nx2[ii] = nx[N_tmp+0];
		nu2[ii] = nu[N_tmp+0];
		nb2[ii] = nb[N_tmp+0];
		ng2[ii] = ng[N_tmp+0];
		ns2[ii] = ns[N_tmp+0];
		for(jj=1; jj<T1; jj++)
			{
			nbb = 0;
			nbg = 0;
			for(kk=0; kk<nb[N_tmp+jj]; kk++)
				if(idxb[N_tmp+jj][kk]<nu[N_tmp+jj])
					nbb++;
				else
					nbg++;
			nx2[ii] += 0;
			nu2[ii] += nu[N_tmp+jj];
			nb2[ii] += nbb;
			ng2[ii] += ng[N_tmp+jj] + nbg;
			ns2[ii] += ns[N_tmp+jj];
			}
		N_tmp += T1;
		}
	nx2[N2] = nx[N];
	nu2[N2] = nu[N];
	nb2[N2] = nb[N];
	ng2[N2] = ng[N];
	ns2[N2] = ns[N];

	return;

	}



int MEMSIZE_COND_QP_OCP2OCP(struct OCP_QP *ocp_qp, struct OCP_QP *part_dense_qp)
	{

	struct OCP_QP tmp_ocp_qp;
	struct DENSE_QP tmp_dense_qp;

	int ii;

	int N = ocp_qp->N;
	int N2 = part_dense_qp->N;
	int N1 = N/N2; // (floor) horizon of small blocks
	int R1 = N - N2*N1; // the first R1 blocks have horizon N1+1
	int M1 = R1>0 ? N1+1 : N1; // (ceil) horizon of large blocks
	int T1; // horizon of current block

	int size = 0;

	size += N2*sizeof(struct COND_QP_OCP2DENSE_WORKSPACE);

	int N_tmp = 0; // temporary sum of horizons
	for(ii=0; ii<N2; ii++)
		{

		T1 = ii<R1 ? M1 : N1;

		// alias ocp_qp
		tmp_ocp_qp.N = T1;
		tmp_ocp_qp.nx = ocp_qp->nx+N_tmp;
		tmp_ocp_qp.nu = ocp_qp->nu+N_tmp;
		tmp_ocp_qp.nb = ocp_qp->nb+N_tmp;
		tmp_ocp_qp.ng = ocp_qp->ng+N_tmp;
		tmp_ocp_qp.idxb = ocp_qp->idxb+N_tmp;

		size += MEMSIZE_COND_QP_OCP2DENSE(&tmp_ocp_qp, &tmp_dense_qp); // TODO floag to avoid to condense the last stage !!!!!

		N_tmp += T1;

		}

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

	}



void CREATE_COND_QP_OCP2OCP(struct OCP_QP *ocp_qp, struct OCP_QP *part_dense_qp, struct COND_QP_OCP2OCP_WORKSPACE *cond_ws, void *mem)
	{

	struct OCP_QP tmp_ocp_qp;
	struct DENSE_QP tmp_dense_qp;

	int ii;

	int N = ocp_qp->N;
	int N2 = part_dense_qp->N;
	int N1 = N/N2; // (floor) horizon of small blocks
	int R1 = N - N2*N1; // the first R1 blocks have horizon N1+1
	int M1 = R1>0 ? N1+1 : N1; // (ceil) horizon of large blocks
	int T1; // horizon of current block

	// cond workspace struct
	struct COND_QP_OCP2DENSE_WORKSPACE *cws_ptr = mem;
	cond_ws->cond_workspace = cws_ptr;
	cws_ptr += N2;

	// align to typicl cache line size
	size_t s_ptr = (size_t) cws_ptr;
	s_ptr = (s_ptr+63)/64*64;

	char *c_ptr = (char *) s_ptr;

	int N_tmp = 0; // temporary sum of horizons
	for(ii=0; ii<N2; ii++)
		{

		T1 = ii<R1 ? M1 : N1;

		// alias ocp_qp
		tmp_ocp_qp.N = T1;
		tmp_ocp_qp.nx = ocp_qp->nx+N_tmp;
		tmp_ocp_qp.nu = ocp_qp->nu+N_tmp;
		tmp_ocp_qp.nb = ocp_qp->nb+N_tmp;
		tmp_ocp_qp.ng = ocp_qp->ng+N_tmp;
		tmp_ocp_qp.idxb = ocp_qp->idxb+N_tmp;

		CREATE_COND_QP_OCP2DENSE(&tmp_ocp_qp, &tmp_dense_qp, cond_ws->cond_workspace+ii, c_ptr);
		c_ptr += (cond_ws->cond_workspace+ii)->memsize;
		(cond_ws->cond_workspace+ii)->cond_last_stage = 0;

		N_tmp += T1;

		}

	cond_ws->memsize = MEMSIZE_COND_QP_OCP2OCP(ocp_qp, part_dense_qp);

	return;

	}

	

void COND_QP_OCP2OCP(struct OCP_QP *ocp_qp, struct OCP_QP *part_dense_qp, struct COND_QP_OCP2OCP_WORKSPACE *part_cond_ws)
	{

	struct OCP_QP tmp_ocp_qp;

	int ii;

	int N = ocp_qp->N;
	int N2 = part_dense_qp->N;
	int N1 = N/N2; // (floor) horizon of small blocks
	int R1 = N - N2*N1; // the first R1 blocks have horizon N1+1
	int M1 = R1>0 ? N1+1 : N1; // (ceil) horizon of large blocks
	int T1; // horizon of current block

	int N_tmp = 0; // temporary sum of horizons
	for(ii=0; ii<N2; ii++)
		{

		T1 = ii<R1 ? M1 : N1;

		// alias ocp_qp
		tmp_ocp_qp.N = T1;
		tmp_ocp_qp.nx = ocp_qp->nx+N_tmp;
		tmp_ocp_qp.nu = ocp_qp->nu+N_tmp;
		tmp_ocp_qp.nb = ocp_qp->nb+N_tmp;
		tmp_ocp_qp.ng = ocp_qp->ng+N_tmp;
		tmp_ocp_qp.idxb = ocp_qp->idxb+N_tmp;
		tmp_ocp_qp.BAbt = ocp_qp->BAbt+N_tmp;
		tmp_ocp_qp.b = ocp_qp->b+N_tmp;
		tmp_ocp_qp.RSQrq = ocp_qp->RSQrq+N_tmp;
		tmp_ocp_qp.rq = ocp_qp->rq+N_tmp;
		tmp_ocp_qp.DCt = ocp_qp->DCt+N_tmp;
		tmp_ocp_qp.d = ocp_qp->d+N_tmp;

		COND_BABT(&tmp_ocp_qp, part_dense_qp->BAbt+ii, part_dense_qp->b+ii, part_cond_ws->cond_workspace+ii);

		COND_RSQRQ_N2NX3(&tmp_ocp_qp, part_dense_qp->RSQrq+ii, part_dense_qp->rq+ii, part_cond_ws->cond_workspace+ii);

		COND_DCTD(&tmp_ocp_qp, part_dense_qp->idxb[ii], part_dense_qp->DCt+ii, part_dense_qp->d+ii, part_cond_ws->cond_workspace+ii);

		N_tmp += T1;

		}

	// copy last stage
	int *nx = ocp_qp->nx;
	int *nu = ocp_qp->nu;
	int *nb = ocp_qp->nb;
	int *ng = ocp_qp->ng;

	GECP_LIBSTR(nu[N]+nx[N]+1, nu[N]+nx[N], ocp_qp->RSQrq+N, 0, 0, part_dense_qp->RSQrq+N2, 0, 0);
	VECCP_LIBSTR(nu[N]+nx[N], ocp_qp->rq+N, 0, part_dense_qp->rq+N2, 0);
	GECP_LIBSTR(nu[N]+nx[N], ng[N], ocp_qp->DCt+N, 0, 0, part_dense_qp->DCt+N2, 0, 0);
	VECCP_LIBSTR(2*nb[N]+2*ng[N], ocp_qp->d+N, 0, part_dense_qp->d+N2, 0);
	for(ii=0; ii<nb[N]; ii++) part_dense_qp->idxb[N2][ii] = ocp_qp->idxb[N][ii];

	return;

	}



void EXPAND_SOL_OCP2OCP(struct OCP_QP *ocp_qp, struct OCP_QP *part_dense_qp, struct OCP_QP_SOL *part_dense_qp_sol, struct OCP_QP_SOL *ocp_qp_sol, struct COND_QP_OCP2OCP_WORKSPACE *part_cond_ws)
	{

	struct OCP_QP tmp_ocp_qp;
	struct OCP_QP_SOL tmp_ocp_qp_sol;
	struct DENSE_QP_SOL dense_qp_sol;

	int *nx = ocp_qp->nx;
	int *nu = ocp_qp->nu;
	int *nb = ocp_qp->nb;
	int *ng = ocp_qp->ng;

	int ii;

	int tmp_nb, tmp_ng;
	struct STRVEC tmp_lam_lb;
	struct STRVEC tmp_lam_lg;
	struct STRVEC tmp_lam_ub;
	struct STRVEC tmp_lam_ug;
	struct STRVEC tmp_t_lb;
	struct STRVEC tmp_t_lg;
	struct STRVEC tmp_t_ub;
	struct STRVEC tmp_t_ug;

	int N = ocp_qp->N;
	int N2 = part_dense_qp->N;
	int N1 = N/N2; // (floor) horizon of small blocks
	int R1 = N - N2*N1; // the first R1 blocks have horizon N1+1
	int M1 = R1>0 ? N1+1 : N1; // (ceil) horizon of large blocks
	int T1; // horizon of current block

	int N_tmp = 0; // temporary sum of horizons
	for(ii=0; ii<N2; ii++)
		{

		T1 = ii<R1 ? M1 : N1;

		// alias ocp_qp
		tmp_ocp_qp.N = T1;
		tmp_ocp_qp.nx = ocp_qp->nx+N_tmp;
		tmp_ocp_qp.nu = ocp_qp->nu+N_tmp;
		tmp_ocp_qp.nb = ocp_qp->nb+N_tmp;
		tmp_ocp_qp.ng = ocp_qp->ng+N_tmp;
		tmp_ocp_qp.idxb = ocp_qp->idxb+N_tmp;
		tmp_ocp_qp.BAbt = ocp_qp->BAbt+N_tmp;
		tmp_ocp_qp.b = ocp_qp->b+N_tmp;
		tmp_ocp_qp.RSQrq = ocp_qp->RSQrq+N_tmp;
		tmp_ocp_qp.rq = ocp_qp->rq+N_tmp;
		tmp_ocp_qp.DCt = ocp_qp->DCt+N_tmp;
		tmp_ocp_qp.d = ocp_qp->d+N_tmp;

		// alias ocp qp sol
		tmp_ocp_qp_sol.ux = ocp_qp_sol->ux+N_tmp;
		tmp_ocp_qp_sol.pi = ocp_qp_sol->pi+N_tmp;
		tmp_ocp_qp_sol.lam = ocp_qp_sol->lam+N_tmp;
		tmp_ocp_qp_sol.t = ocp_qp_sol->t+N_tmp;

		// alias ocp qp sol
		tmp_nb = part_dense_qp->nb[ii];
		tmp_ng = part_dense_qp->ng[ii];
		CREATE_STRVEC(tmp_nb, &tmp_lam_lb, (part_dense_qp_sol->lam+ii)->pa+0);
		CREATE_STRVEC(tmp_ng, &tmp_lam_lg, (part_dense_qp_sol->lam+ii)->pa+tmp_nb);
		CREATE_STRVEC(tmp_nb, &tmp_lam_ub, (part_dense_qp_sol->lam+ii)->pa+tmp_nb+tmp_ng);
		CREATE_STRVEC(tmp_ng, &tmp_lam_ug, (part_dense_qp_sol->lam+ii)->pa+2*tmp_nb+tmp_ng);
		CREATE_STRVEC(tmp_nb, &tmp_t_lb, (part_dense_qp_sol->t+ii)->pa+0);
		CREATE_STRVEC(tmp_ng, &tmp_t_lg, (part_dense_qp_sol->t+ii)->pa+tmp_nb);
		CREATE_STRVEC(tmp_nb, &tmp_t_ub, (part_dense_qp_sol->t+ii)->pa+tmp_nb+tmp_ng);
		CREATE_STRVEC(tmp_ng, &tmp_t_ug, (part_dense_qp_sol->t+ii)->pa+2*tmp_nb+tmp_ng);

		dense_qp_sol.v = part_dense_qp_sol->ux+ii;
		dense_qp_sol.pi = part_dense_qp_sol->pi+ii;
		dense_qp_sol.lam_lb = &tmp_lam_lb;
		dense_qp_sol.lam_ub = &tmp_lam_ub;
		dense_qp_sol.lam_lg = &tmp_lam_lg;
		dense_qp_sol.lam_ug = &tmp_lam_ug;
		dense_qp_sol.t_lb = &tmp_t_lb;
		dense_qp_sol.t_ub = &tmp_t_ub;
		dense_qp_sol.t_lg = &tmp_t_lg;
		dense_qp_sol.t_ug = &tmp_t_ug;

		EXPAND_SOL(&tmp_ocp_qp, &dense_qp_sol, &tmp_ocp_qp_sol, part_cond_ws->cond_workspace);

		N_tmp += T1;

		}

	// copy last stage
	VECCP_LIBSTR(nu[N]+nx[N], part_dense_qp_sol->ux+N2, 0, ocp_qp_sol->ux+N, 0);
	VECCP_LIBSTR(2*nb[N]+2*ng[N], part_dense_qp_sol->lam+N2, 0, ocp_qp_sol->lam+N, 0);
	VECCP_LIBSTR(2*nb[N]+2*ng[N], part_dense_qp_sol->t+N2, 0, ocp_qp_sol->t+N, 0);

	return;

	}
