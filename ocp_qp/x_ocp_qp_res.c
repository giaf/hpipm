/**************************************************************************************************
*                                                                                                 *
* This file is part of HPIPM.                                                                     *
*                                                                                                 *
* HPIPM -- High-Performance Interior Point Method.                                                *
* Copyright (C) 2017-2018 by Gianluca Frison.                                                     *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
*                                                                                                 *
* This program is free software: you can redistribute it and/or modify                            *
* it under the terms of the GNU General Public License as published by                            *
* the Free Software Foundation, either version 3 of the License, or                               *
* (at your option) any later version                                                              *.
*                                                                                                 *
* This program is distributed in the hope that it will be useful,                                 *
* but WITHOUT ANY WARRANTY; without even the implied warranty of                                  *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                   *
* GNU General Public License for more details.                                                    *
*                                                                                                 *
* You should have received a copy of the GNU General Public License                               *
* along with this program.  If not, see <https://www.gnu.org/licenses/>.                          *
*                                                                                                 *
* The authors designate this particular file as subject to the "Classpath" exception              *
* as provided by the authors in the LICENSE file that accompained this code.                      *
*                                                                                                 *
* Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             *
*                                                                                                 *
**************************************************************************************************/



int OCP_QP_RES_MEMSIZE(struct OCP_QP_DIM *dim)
	{

	// loop index
	int ii;

	// extract ocp qp size
	int N = dim->N;
	int *nx = dim->nx;
	int *nu = dim->nu;
	int *nb = dim->nb;
	int *ng = dim->ng;
	int *ns = dim->ns;

	// compute core qp size
	int nvt = 0;
	int net = 0;
	int nct = 0;
	for(ii=0; ii<N; ii++)
		{
		nvt += nx[ii]+nu[ii]+2*ns[ii];
		net += nx[ii+1];
		nct += 2*nb[ii]+2*ng[ii]+2*ns[ii];
		}
	nvt += nx[ii]+nu[ii]+2*ns[ii];
	nct += 2*nb[ii]+2*ng[ii]+2*ns[ii];

	int size = 0;

	size += 3*(N+1)*sizeof(struct STRVEC); // res_g res_d res_m
	size += 3*N*sizeof(struct STRVEC); // res_b

	size += 1*SIZE_STRVEC(nvt); // res_g
	size += 1*SIZE_STRVEC(net); // res_b
	size += 2*SIZE_STRVEC(nct); // res_d res_m

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

	}



void OCP_QP_RES_CREATE(struct OCP_QP_DIM *dim, struct OCP_QP_RES *res, void *mem)
	{

	// loop index
	int ii;

	// zero memory (to avoid corrupted memory like e.g. NaN)
	int memsize = OCP_QP_RES_MEMSIZE(dim);
	int memsize_m8 = memsize/8; // sizeof(double) is 8
//	int memsize_r8 = memsize - 8*memsize_m8;
	double *double_ptr = mem;
	// XXX exploit that it is multiple of 64 bytes !!!!!
	for(ii=0; ii<memsize_m8-7; ii+=8)
		{
		double_ptr[ii+0] = 0.0;
		double_ptr[ii+1] = 0.0;
		double_ptr[ii+2] = 0.0;
		double_ptr[ii+3] = 0.0;
		double_ptr[ii+4] = 0.0;
		double_ptr[ii+5] = 0.0;
		double_ptr[ii+6] = 0.0;
		double_ptr[ii+7] = 0.0;
		}
//	for(; ii<memsize_m8; ii++)
//		{
//		double_ptr[ii] = 0.0;
//		}
//	char *char_ptr = (char *) (&double_ptr[ii]);
//	for(ii=0; ii<memsize_r8; ii++)
//		{
//		char_ptr[ii] = 0;
//		}

	// extract ocp qp size
	int N = dim->N;
	int *nx = dim->nx;
	int *nu = dim->nu;
	int *nb = dim->nb;
	int *ng = dim->ng;
	int *ns = dim->ns;

	// compute core qp size
	int nvt = 0;
	int net = 0;
	int nct = 0;
	for(ii=0; ii<N; ii++)
		{
		nvt += nx[ii]+nu[ii]+2*ns[ii];
		net += nx[ii+1];
		nct += 2*nb[ii]+2*ng[ii]+2*ns[ii];
		}
	nvt += nx[ii]+nu[ii]+2*ns[ii];
	nct += 2*nb[ii]+2*ng[ii]+2*ns[ii];


	// vector struct
	struct STRVEC *sv_ptr = (struct STRVEC *) mem;

	res->res_g = sv_ptr;
	sv_ptr += N+1;
	res->res_b = sv_ptr;
	sv_ptr += N;
	res->res_d = sv_ptr;
	sv_ptr += N+1;
	res->res_m = sv_ptr;
	sv_ptr += N+1;


	// align to typicl cache line size
	size_t s_ptr = (size_t) sv_ptr;
	s_ptr = (s_ptr+63)/64*64;


	// void stuf
	char *c_ptr = (char *) s_ptr;

	CREATE_STRVEC(nvt, res->res_g, c_ptr);
	c_ptr += SIZE_STRVEC(nvt);

	CREATE_STRVEC(net, res->res_b, c_ptr);
	c_ptr += SIZE_STRVEC(net);

	CREATE_STRVEC(nct, res->res_d, c_ptr);
	c_ptr += SIZE_STRVEC(nct);

	CREATE_STRVEC(nct, res->res_m, c_ptr);
	c_ptr += SIZE_STRVEC(nct);

	// alias
	//
	c_ptr = (char *) res->res_g->pa;
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(nu[ii]+nx[ii]+2*ns[ii], res->res_g+ii, c_ptr);
		c_ptr += nu[ii]*sizeof(REAL);
		c_ptr += nx[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		}
	//
	c_ptr = (char *) res->res_b->pa;
	for(ii=0; ii<N; ii++)
		{
		CREATE_STRVEC(nx[ii+1], res->res_b+ii, c_ptr);
		c_ptr += (nx[ii+1])*sizeof(REAL);
		}
	//
	c_ptr = (char *) res->res_d->pa;
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(2*nb[ii]+2*ng[ii]+2*ns[ii], res->res_d+ii, c_ptr);
		c_ptr += nb[ii]*sizeof(REAL);
		c_ptr += ng[ii]*sizeof(REAL);
		c_ptr += nb[ii]*sizeof(REAL);
		c_ptr += ng[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		}
	//
	c_ptr = (char *) res->res_m->pa;
	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(2*nb[ii]+2*ng[ii]+2*ns[ii], res->res_m+ii, c_ptr);
		c_ptr += nb[ii]*sizeof(REAL);
		c_ptr += ng[ii]*sizeof(REAL);
		c_ptr += nb[ii]*sizeof(REAL);
		c_ptr += ng[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		}



	res->dim = dim;

	res->memsize = OCP_QP_RES_MEMSIZE(dim);


#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + res->memsize)
		{
		printf("\ncreate_ocp_qp_res: outside memory bounds!\n\n");
		exit(1);
		}
#endif


	return;

	}



int OCP_QP_RES_WS_MEMSIZE(struct OCP_QP_DIM *dim)
	{

	// loop index
	int ii;

	// extract ocp qp size
	int N = dim->N;
	int *nx = dim->nx;
	int *nu = dim->nu;
	int *nb = dim->nb;
	int *ng = dim->ng;
	int *ns = dim->ns;

	// compute core qp size and max size
	int nbM = 0;
	int ngM = 0;
	int nsM = 0;
	for(ii=0; ii<N; ii++)
		{
		nbM = nb[ii]>nbM ? nb[ii] : nbM;
		ngM = ng[ii]>ngM ? ng[ii] : ngM;
		nsM = ns[ii]>nsM ? ns[ii] : nsM;
		}
	nbM = nb[ii]>nbM ? nb[ii] : nbM;
	ngM = ng[ii]>ngM ? ng[ii] : ngM;
	nsM = ns[ii]>nsM ? ns[ii] : nsM;

	int size = 0;

	size += 3*sizeof(struct STRVEC); // 2*tmp_nbgM tmp_nsM

	size += 2*SIZE_STRVEC(nbM+ngM); // tmp_nbgM
	size += 1*SIZE_STRVEC(nsM); // tmp_nsM

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

	}



void OCP_QP_RES_WS_CREATE(struct OCP_QP_DIM *dim, struct OCP_QP_RES_WS *ws, void *mem)
	{

	// loop index
	int ii;

	// zero memory (to avoid corrupted memory like e.g. NaN)
	int memsize = OCP_QP_RES_WS_MEMSIZE(dim);
	int memsize_m8 = memsize/8; // sizeof(double) is 8
//	int memsize_r8 = memsize - 8*memsize_m8;
	double *double_ptr = mem;
	// XXX exploit that it is multiple of 64 bytes !!!!!
	for(ii=0; ii<memsize_m8-7; ii+=8)
		{
		double_ptr[ii+0] = 0.0;
		double_ptr[ii+1] = 0.0;
		double_ptr[ii+2] = 0.0;
		double_ptr[ii+3] = 0.0;
		double_ptr[ii+4] = 0.0;
		double_ptr[ii+5] = 0.0;
		double_ptr[ii+6] = 0.0;
		double_ptr[ii+7] = 0.0;
		}
//	for(; ii<memsize_m8; ii++)
//		{
//		double_ptr[ii] = 0.0;
//		}
//	char *char_ptr = (char *) (&double_ptr[ii]);
//	for(ii=0; ii<memsize_r8; ii++)
//		{
//		char_ptr[ii] = 0;
//		}

	// extract ocp qp size
	int N = dim->N;
	int *nx = dim->nx;
	int *nu = dim->nu;
	int *nb = dim->nb;
	int *ng = dim->ng;
	int *ns = dim->ns;


	// compute core qp size and max size
	int nbM = 0;
	int ngM = 0;
	int nsM = 0;
	for(ii=0; ii<N; ii++)
		{
		nbM = nb[ii]>nbM ? nb[ii] : nbM;
		ngM = ng[ii]>ngM ? ng[ii] : ngM;
		nsM = ns[ii]>nsM ? ns[ii] : nsM;
		}
	nbM = nb[ii]>nbM ? nb[ii] : nbM;
	ngM = ng[ii]>ngM ? ng[ii] : ngM;
	nsM = ns[ii]>nsM ? ns[ii] : nsM;


	// vector struct
	struct STRVEC *sv_ptr = (struct STRVEC *) mem;

	ws->tmp_nbgM = sv_ptr;
	sv_ptr += 2;
	ws->tmp_nsM = sv_ptr;
	sv_ptr += 1;


	// align to typicl cache line size
	size_t s_ptr = (size_t) sv_ptr;
	s_ptr = (s_ptr+63)/64*64;


	// void stuf
	char *c_ptr = (char *) s_ptr;


	CREATE_STRVEC(nbM+ngM, ws->tmp_nbgM+0, c_ptr);
	c_ptr += (ws->tmp_nbgM+0)->memsize;

	CREATE_STRVEC(nbM+ngM, ws->tmp_nbgM+1, c_ptr);
	c_ptr += (ws->tmp_nbgM+1)->memsize;

	CREATE_STRVEC(nsM, ws->tmp_nsM+0, c_ptr);
	c_ptr += (ws->tmp_nsM+0)->memsize;

	ws->memsize = OCP_QP_RES_WS_MEMSIZE(dim);


#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + ws->memsize)
		{
		printf("\ncreate_ocp_qp_res_workspace: outside memory bounds!\n\n");
		exit(1);
		}
#endif


	return;

	}



void OCP_QP_RES_GET_ALL(struct OCP_QP_RES *res, REAL **res_r, REAL **res_q, REAL **res_ls, REAL **res_us, REAL **res_b, REAL **res_d_lb, REAL **res_d_ub, REAL **res_d_lg, REAL **res_d_ug, REAL **res_d_ls, REAL **res_d_us, REAL **res_m_lb, REAL **res_m_ub, REAL **res_m_lg, REAL **res_m_ug, REAL **res_m_ls, REAL **res_m_us)
	{

	int N = res->dim->N;
	int *nx = res->dim->nx;
	int *nu = res->dim->nu;
	int *nb = res->dim->nb;
	int *ng = res->dim->ng;
	int *ns = res->dim->ns;

	int ii;

	for(ii=0; ii<N; ii++)
		{
		// cost
		CVT_STRVEC2VEC(nu[ii], res->res_g+ii, 0, res_r[ii]);
		CVT_STRVEC2VEC(nx[ii], res->res_g+ii, nu[ii], res_q[ii]);

		// dynamics
		CVT_STRVEC2VEC(nx[ii+1], res->res_b+ii, 0, res_b[ii]);

		// box constraints
		if(nb[ii]>0)
			{
			CVT_STRVEC2VEC(nb[ii], res->res_d+ii, 0, res_d_lb[ii]);
			CVT_STRVEC2VEC(nb[ii], res->res_d+ii, nb[ii]+ng[ii], res_d_ub[ii]);
			CVT_STRVEC2VEC(nb[ii], res->res_m+ii, 0, res_m_lb[ii]);
			CVT_STRVEC2VEC(nb[ii], res->res_m+ii, nb[ii]+ng[ii], res_m_ub[ii]);
			}

		// general constraints
		if(ng[ii]>0)
			{
			CVT_STRVEC2VEC(ng[ii], res->res_d+ii, nb[ii], res_d_lg[ii]);
			CVT_STRVEC2VEC(ng[ii], res->res_d+ii, 2*nb[ii]+ng[ii], res_d_ug[ii]);
			CVT_STRVEC2VEC(ng[ii], res->res_m+ii, nb[ii], res_m_lg[ii]);
			CVT_STRVEC2VEC(ng[ii], res->res_m+ii, 2*nb[ii]+ng[ii], res_m_ug[ii]);
			}

		// soft constraints
		if(ns[ii]>0)
			{
			CVT_STRVEC2VEC(ns[ii], res->res_g+ii, nu[ii]+nx[ii], res_ls[ii]);
			CVT_STRVEC2VEC(ns[ii], res->res_g+ii, nu[ii]+nx[ii]+ns[ii], res_us[ii]);
			CVT_STRVEC2VEC(ns[ii], res->res_d+ii, 2*nb[ii]+2*ng[ii], res_d_ls[ii]);
			CVT_STRVEC2VEC(ns[ii], res->res_d+ii, 2*nb[ii]+2*ng[ii]+ns[ii], res_d_us[ii]);
			CVT_STRVEC2VEC(ns[ii], res->res_m+ii, 2*nb[ii]+2*ng[ii], res_m_ls[ii]);
			CVT_STRVEC2VEC(ns[ii], res->res_m+ii, 2*nb[ii]+2*ng[ii]+ns[ii], res_m_us[ii]);
			}
		}

	// cost
	CVT_STRVEC2VEC(nu[ii], res->res_g+ii, 0, res_r[ii]);
	CVT_STRVEC2VEC(nx[ii], res->res_g+ii, nu[ii], res_q[ii]);

	// box constraints
	if(nb[ii]>0)
		{
		CVT_STRVEC2VEC(nb[ii], res->res_d+ii, 0, res_d_lb[ii]);
		CVT_STRVEC2VEC(nb[ii], res->res_d+ii, nb[ii]+ng[ii], res_d_ub[ii]);
		CVT_STRVEC2VEC(nb[ii], res->res_m+ii, 0, res_m_lb[ii]);
		CVT_STRVEC2VEC(nb[ii], res->res_m+ii, nb[ii]+ng[ii], res_m_ub[ii]);
		}

	// general constraints
	if(ng[ii]>0)
		{
		CVT_STRVEC2VEC(ng[ii], res->res_d+ii, nb[ii], res_d_lg[ii]);
		CVT_STRVEC2VEC(ng[ii], res->res_d+ii, 2*nb[ii]+ng[ii], res_d_ug[ii]);
		CVT_STRVEC2VEC(ng[ii], res->res_m+ii, nb[ii], res_m_lg[ii]);
		CVT_STRVEC2VEC(ng[ii], res->res_m+ii, 2*nb[ii]+ng[ii], res_m_ug[ii]);
		}

	// soft constraints
	if(ns[ii]>0)
		{
		CVT_STRVEC2VEC(ns[ii], res->res_g+ii, nu[ii]+nx[ii], res_ls[ii]);
		CVT_STRVEC2VEC(ns[ii], res->res_g+ii, nu[ii]+nx[ii]+ns[ii], res_us[ii]);
		CVT_STRVEC2VEC(ns[ii], res->res_d+ii, 2*nb[ii]+2*ng[ii], res_d_ls[ii]);
		CVT_STRVEC2VEC(ns[ii], res->res_d+ii, 2*nb[ii]+2*ng[ii]+ns[ii], res_d_us[ii]);
		CVT_STRVEC2VEC(ns[ii], res->res_m+ii, 2*nb[ii]+2*ng[ii], res_m_ls[ii]);
		CVT_STRVEC2VEC(ns[ii], res->res_m+ii, 2*nb[ii]+2*ng[ii]+ns[ii], res_m_us[ii]);
		}

	return;

	}



