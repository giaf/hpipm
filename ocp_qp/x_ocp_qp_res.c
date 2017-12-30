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



int MEMSIZE_OCP_QP_RES(struct OCP_QP_DIM *dim)
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

	int size = 0;

	size += 3*(N+1)*sizeof(struct STRVEC); // res_g res_d res_m
	size += 3*N*sizeof(struct STRVEC); // res_b

	for(ii=0; ii<=N; ii++) size += 1*SIZE_STRVEC(nu[ii]+nx[ii]+2*ns[ii]); // res_g
	for(ii=0; ii<N; ii++) size += 1*SIZE_STRVEC(nx[ii+1]); // res_b
	for(ii=0; ii<=N; ii++) size += 2*SIZE_STRVEC(2*nb[ii]+2*ng[ii]+2*ns[ii]); // res_d res_m

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

	}



void CREATE_OCP_QP_RES(struct OCP_QP_DIM *dim, struct OCP_QP_RES *res, void *mem)
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

	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(nu[ii]+nx[ii]+2*ns[ii], res->res_g+ii, c_ptr);
		c_ptr += (res->res_g+ii)->memsize;
		}

	for(ii=0; ii<N; ii++)
		{
		CREATE_STRVEC(nx[ii+1], res->res_b+ii, c_ptr);
		c_ptr += (res->res_b+ii)->memsize;
		}

	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(2*nb[ii]+2*ng[ii]+2*ns[ii], res->res_d+ii, c_ptr);
		c_ptr += (res->res_d+ii)->memsize;
		}

	for(ii=0; ii<=N; ii++)
		{
		CREATE_STRVEC(2*nb[ii]+2*ng[ii]+2*ns[ii], res->res_m+ii, c_ptr);
		c_ptr += (res->res_m+ii)->memsize;
		}

	res->dim = dim;

	res->memsize = MEMSIZE_OCP_QP_RES(dim);


#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + res->memsize)
		{
		printf("\ncreate_ocp_qp_res: outsize memory bounds!\n\n");
		exit(1);
		}
#endif


	return;

	}



int MEMSIZE_OCP_QP_RES_WORKSPACE(struct OCP_QP_DIM *dim)
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



void CREATE_OCP_QP_RES_WORKSPACE(struct OCP_QP_DIM *dim, struct OCP_QP_RES_WORKSPACE *ws, void *mem)
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

	ws->memsize = MEMSIZE_OCP_QP_RES(dim);


#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + ws->memsize)
		{
		printf("\ncreate_ocp_qp_res_workspace: outsize memory bounds!\n\n");
		exit(1);
		}
#endif


	return;

	}



void CVT_OCP_QP_RES_TO_COLMAJ(struct OCP_QP_RES *res, REAL **res_r, REAL **res_q, REAL **res_ls, REAL **res_us, REAL **res_b, REAL **res_d_lb, REAL **res_d_ub, REAL **res_d_lg, REAL **res_d_ug, REAL **res_d_ls, REAL **res_d_us, REAL **res_m_lb, REAL **res_m_ub, REAL **res_m_lg, REAL **res_m_ug, REAL **res_m_ls, REAL **res_m_us)
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
		CVT_STRVEC2VEC(nu[ii], res->res_g+ii, 0, res_r[ii]);
		CVT_STRVEC2VEC(nx[ii], res->res_g+ii, nu[ii], res_q[ii]);
		CVT_STRVEC2VEC(ns[ii], res->res_g+ii, nu[ii]+nx[ii], res_ls[ii]);
		CVT_STRVEC2VEC(ns[ii], res->res_g+ii, nu[ii]+nx[ii]+ns[ii], res_us[ii]);

		CVT_STRVEC2VEC(nx[ii+1], res->res_b+ii, 0, res_b[ii]);
		CVT_STRVEC2VEC(nb[ii], res->res_d+ii, 0, res_d_lb[ii]);
		CVT_STRVEC2VEC(ng[ii], res->res_d+ii, nb[ii], res_d_lg[ii]);
		CVT_STRVEC2VEC(nb[ii], res->res_d+ii, nb[ii]+ng[ii], res_d_ub[ii]);
		CVT_STRVEC2VEC(ng[ii], res->res_d+ii, 2*nb[ii]+ng[ii], res_d_ug[ii]);
		CVT_STRVEC2VEC(ns[ii], res->res_d+ii, 2*nb[ii]+2*ng[ii], res_d_ls[ii]);
		CVT_STRVEC2VEC(ns[ii], res->res_d+ii, 2*nb[ii]+2*ng[ii]+ns[ii], res_d_us[ii]);
		CVT_STRVEC2VEC(nb[ii], res->res_m+ii, 0, res_m_lb[ii]);
		CVT_STRVEC2VEC(ng[ii], res->res_m+ii, nb[ii], res_m_lg[ii]);
		CVT_STRVEC2VEC(nb[ii], res->res_m+ii, nb[ii]+ng[ii], res_m_ub[ii]);
		CVT_STRVEC2VEC(ng[ii], res->res_m+ii, 2*nb[ii]+ng[ii], res_m_ug[ii]);
		CVT_STRVEC2VEC(ns[ii], res->res_m+ii, 2*nb[ii]+2*ng[ii], res_m_ls[ii]);
		CVT_STRVEC2VEC(ns[ii], res->res_m+ii, 2*nb[ii]+2*ng[ii]+ns[ii], res_m_us[ii]);
		}
	CVT_STRVEC2VEC(nu[ii], res->res_g+ii, 0, res_r[ii]);
	CVT_STRVEC2VEC(nx[ii], res->res_g+ii, nu[ii], res_q[ii]);
	CVT_STRVEC2VEC(ns[ii], res->res_g+ii, nu[ii]+nx[ii], res_ls[ii]);
	CVT_STRVEC2VEC(ns[ii], res->res_g+ii, nu[ii]+nx[ii]+ns[ii], res_us[ii]);

	CVT_STRVEC2VEC(nb[ii], res->res_d+ii, 0, res_d_lb[ii]);
	CVT_STRVEC2VEC(ng[ii], res->res_d+ii, nb[ii], res_d_lg[ii]);
	CVT_STRVEC2VEC(nb[ii], res->res_d+ii, nb[ii]+ng[ii], res_d_ub[ii]);
	CVT_STRVEC2VEC(ng[ii], res->res_d+ii, 2*nb[ii]+ng[ii], res_d_ug[ii]);
	CVT_STRVEC2VEC(ns[ii], res->res_d+ii, 2*nb[ii]+2*ng[ii], res_d_ls[ii]);
	CVT_STRVEC2VEC(ns[ii], res->res_d+ii, 2*nb[ii]+2*ng[ii]+ns[ii], res_d_us[ii]);
	CVT_STRVEC2VEC(nb[ii], res->res_m+ii, 0, res_m_lb[ii]);
	CVT_STRVEC2VEC(ng[ii], res->res_m+ii, nb[ii], res_m_lg[ii]);
	CVT_STRVEC2VEC(nb[ii], res->res_m+ii, nb[ii]+ng[ii], res_m_ub[ii]);
	CVT_STRVEC2VEC(ng[ii], res->res_m+ii, 2*nb[ii]+ng[ii], res_m_ug[ii]);
	CVT_STRVEC2VEC(ns[ii], res->res_m+ii, 2*nb[ii]+2*ng[ii], res_m_ls[ii]);
	CVT_STRVEC2VEC(ns[ii], res->res_m+ii, 2*nb[ii]+2*ng[ii]+ns[ii], res_m_us[ii]);

	return;

	}



void CVT_OCP_QP_RES_TO_ROWMAJ(struct OCP_QP_RES *res, REAL **res_r, REAL **res_q, REAL **res_ls, REAL **res_us, REAL **res_b, REAL **res_d_lb, REAL **res_d_ub, REAL **res_d_lg, REAL **res_d_ug, REAL **res_d_ls, REAL **res_d_us, REAL **res_m_lb, REAL **res_m_ub, REAL **res_m_lg, REAL **res_m_ug, REAL **res_m_ls, REAL **res_m_us)
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
		CVT_STRVEC2VEC(nu[ii], res->res_g+ii, 0, res_r[ii]);
		CVT_STRVEC2VEC(nx[ii], res->res_g+ii, nu[ii], res_q[ii]);
		CVT_STRVEC2VEC(ns[ii], res->res_g+ii, nu[ii]+nx[ii], res_ls[ii]);
		CVT_STRVEC2VEC(ns[ii], res->res_g+ii, nu[ii]+nx[ii]+ns[ii], res_us[ii]);

		CVT_STRVEC2VEC(nx[ii+1], res->res_b+ii, 0, res_b[ii]);
		CVT_STRVEC2VEC(nb[ii], res->res_d+ii, 0, res_d_lb[ii]);
		CVT_STRVEC2VEC(ng[ii], res->res_d+ii, nb[ii], res_d_lg[ii]);
		CVT_STRVEC2VEC(nb[ii], res->res_d+ii, nb[ii]+ng[ii], res_d_ub[ii]);
		CVT_STRVEC2VEC(ng[ii], res->res_d+ii, 2*nb[ii]+ng[ii], res_d_ug[ii]);
		CVT_STRVEC2VEC(ns[ii], res->res_d+ii, 2*nb[ii]+2*ng[ii], res_d_ls[ii]);
		CVT_STRVEC2VEC(ns[ii], res->res_d+ii, 2*nb[ii]+2*ng[ii]+ns[ii], res_d_us[ii]);
		CVT_STRVEC2VEC(nb[ii], res->res_m+ii, 0, res_m_lb[ii]);
		CVT_STRVEC2VEC(ng[ii], res->res_m+ii, nb[ii], res_m_lg[ii]);
		CVT_STRVEC2VEC(nb[ii], res->res_m+ii, nb[ii]+ng[ii], res_m_ub[ii]);
		CVT_STRVEC2VEC(ng[ii], res->res_m+ii, 2*nb[ii]+ng[ii], res_m_ug[ii]);
		CVT_STRVEC2VEC(ns[ii], res->res_m+ii, 2*nb[ii]+2*ng[ii], res_m_ls[ii]);
		CVT_STRVEC2VEC(ns[ii], res->res_m+ii, 2*nb[ii]+2*ng[ii]+ns[ii], res_m_us[ii]);
		}
	CVT_STRVEC2VEC(nu[ii], res->res_g+ii, 0, res_r[ii]);
	CVT_STRVEC2VEC(nx[ii], res->res_g+ii, nu[ii], res_q[ii]);
	CVT_STRVEC2VEC(ns[ii], res->res_g+ii, nu[ii]+nx[ii], res_ls[ii]);
	CVT_STRVEC2VEC(ns[ii], res->res_g+ii, nu[ii]+nx[ii]+ns[ii], res_us[ii]);

	CVT_STRVEC2VEC(nb[ii], res->res_d+ii, 0, res_d_lb[ii]);
	CVT_STRVEC2VEC(ng[ii], res->res_d+ii, nb[ii], res_d_lg[ii]);
	CVT_STRVEC2VEC(nb[ii], res->res_d+ii, nb[ii]+ng[ii], res_d_ub[ii]);
	CVT_STRVEC2VEC(ng[ii], res->res_d+ii, 2*nb[ii]+ng[ii], res_d_ug[ii]);
	CVT_STRVEC2VEC(ns[ii], res->res_d+ii, 2*nb[ii]+2*ng[ii], res_d_ls[ii]);
	CVT_STRVEC2VEC(ns[ii], res->res_d+ii, 2*nb[ii]+2*ng[ii]+ns[ii], res_d_us[ii]);
	CVT_STRVEC2VEC(nb[ii], res->res_m+ii, 0, res_m_lb[ii]);
	CVT_STRVEC2VEC(ng[ii], res->res_m+ii, nb[ii], res_m_lg[ii]);
	CVT_STRVEC2VEC(nb[ii], res->res_m+ii, nb[ii]+ng[ii], res_m_ub[ii]);
	CVT_STRVEC2VEC(ng[ii], res->res_m+ii, 2*nb[ii]+ng[ii], res_m_ug[ii]);
	CVT_STRVEC2VEC(ns[ii], res->res_m+ii, 2*nb[ii]+2*ng[ii], res_m_ls[ii]);
	CVT_STRVEC2VEC(ns[ii], res->res_m+ii, 2*nb[ii]+2*ng[ii]+ns[ii], res_m_us[ii]);

	return;

	}





