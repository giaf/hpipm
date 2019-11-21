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



int MEMSIZE_TREE_OCP_QP_RES(struct TREE_OCP_QP_DIM *dim)
	{

	// loop index
	int ii, idx;

	// extract ocp qp size
	int Nn = dim->Nn;
	int *nx = dim->nx;
	int *nu = dim->nu;
	int *nb = dim->nb;
	int *ng = dim->ng;
	int *ns = dim->ns;

	// compute core qp size
	int nvt = 0;
	int net = 0;
	int nct = 0;
	for(ii=0; ii<Nn; ii++)
		{
		nvt += nx[ii]+nu[ii]+2*ns[ii];
		nct += 2*nb[ii]+2*ng[ii]+2*ns[ii];
		}
	for(ii=0; ii<Nn-1; ii++)
		{
		idx = ii+1;
		net += nx[idx];
		}

	int size = 0;

	size += 3*Nn*sizeof(struct STRVEC); // res_g res_d res_m
	size += 3*(Nn-1)*sizeof(struct STRVEC); // res_b

	size += 1*SIZE_STRVEC(nvt); // res_g
	size += 1*SIZE_STRVEC(net); // res_b
	size += 2*SIZE_STRVEC(nct); // res_d res_m

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

	}



void CREATE_TREE_OCP_QP_RES(struct TREE_OCP_QP_DIM *dim, struct TREE_OCP_QP_RES *res, void *mem)
	{

	// loop index
	int ii, idx;

	// extract ocp qp size
	int Nn = dim->Nn;
	int *nx = dim->nx;
	int *nu = dim->nu;
	int *nb = dim->nb;
	int *ng = dim->ng;
	int *ns = dim->ns;

	// compute core qp size
	int nvt = 0;
	int net = 0;
	int nct = 0;
	for(ii=0; ii<Nn; ii++)
		{
		nvt += nx[ii]+nu[ii]+2*ns[ii];
		nct += 2*nb[ii]+2*ng[ii]+2*ns[ii];
		}
	for(ii=0; ii<Nn-1; ii++)
		{
		idx = ii+1;
		net += nx[idx];
		}


	// vector struct
	struct STRVEC *sv_ptr = (struct STRVEC *) mem;

	res->res_g = sv_ptr;
	sv_ptr += Nn;
	res->res_b = sv_ptr;
	sv_ptr += Nn-1;
	res->res_d = sv_ptr;
	sv_ptr += Nn;
	res->res_m = sv_ptr;
	sv_ptr += Nn;


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
	for(ii=0; ii<Nn; ii++)
		{
		CREATE_STRVEC(nu[ii]+nx[ii]+2*ns[ii], res->res_g+ii, c_ptr);
		c_ptr += nu[ii]*sizeof(REAL);
		c_ptr += nx[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		c_ptr += ns[ii]*sizeof(REAL);
		}
	//
	c_ptr = (char *) res->res_b->pa;
	for(ii=0; ii<Nn-1; ii++)
		{
		idx = ii+1;
		CREATE_STRVEC(nx[idx], res->res_b+ii, c_ptr);
		c_ptr += (nx[idx])*sizeof(REAL);
		}
	//
	c_ptr = (char *) res->res_d->pa;
	for(ii=0; ii<Nn; ii++)
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
	for(ii=0; ii<Nn; ii++)
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

	res->memsize = MEMSIZE_TREE_OCP_QP_RES(dim);


#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + res->memsize)
		{
		printf("\ncreate_tree_ocp_qp_res: outside memory bounds!\n\n");
		exit(1);
		}
#endif


	return;

	}



int MEMSIZE_TREE_OCP_QP_RES_WORKSPACE(struct TREE_OCP_QP_DIM *dim)
	{

	// loop index
	int ii, idx;

	// extract ocp qp size
	int Nn = dim->Nn;
	int *nx = dim->nx;
	int *nu = dim->nu;
	int *nb = dim->nb;
	int *ng = dim->ng;
	int *ns = dim->ns;

	// compute core qp size and max size
	int nbM = 0;
	int ngM = 0;
	int nsM = 0;
	for(ii=0; ii<Nn; ii++)
		{
		nbM = nb[ii]>nbM ? nb[ii] : nbM;
		ngM = ng[ii]>ngM ? ng[ii] : ngM;
		nsM = ns[ii]>nsM ? ns[ii] : nsM;
		}

	int size = 0;

	size += 3*sizeof(struct STRVEC); // 2*tmp_nbgM tmp_nsM

	size += 2*SIZE_STRVEC(nbM+ngM); // tmp_nbgM
	size += 1*SIZE_STRVEC(nsM); // tmp_nsM

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 1*64; // align once to typical cache line size

	return size;

	}



void CREATE_TREE_OCP_QP_RES_WORKSPACE(struct TREE_OCP_QP_DIM *dim, struct TREE_OCP_QP_RES_WORKSPACE *ws, void *mem)
	{

	// loop index
	int ii, idx;

	// extract ocp qp size
	int Nn = dim->Nn;
	int *nx = dim->nx;
	int *nu = dim->nu;
	int *nb = dim->nb;
	int *ng = dim->ng;
	int *ns = dim->ns;


	// compute core qp size and max size
	int nbM = 0;
	int ngM = 0;
	int nsM = 0;
	for(ii=0; ii<Nn; ii++)
		{
		nbM = nb[ii]>nbM ? nb[ii] : nbM;
		ngM = ng[ii]>ngM ? ng[ii] : ngM;
		nsM = ns[ii]>nsM ? ns[ii] : nsM;
		}


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

	ws->memsize = MEMSIZE_TREE_OCP_QP_RES(dim);


#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + ws->memsize)
		{
		printf("\ncreate_tree_ocp_qp_res_workspace: outside memory bounds!\n\n");
		exit(1);
		}
#endif


	return;

	}



void CVT_TREE_OCP_QP_RES_TO_COLMAJ(struct TREE_OCP_QP_RES *res, REAL **res_r, REAL **res_q, REAL **res_ls, REAL **res_us, REAL **res_b, REAL **res_d_lb, REAL **res_d_ub, REAL **res_d_lg, REAL **res_d_ug, REAL **res_d_ls, REAL **res_d_us, REAL **res_m_lb, REAL **res_m_ub, REAL **res_m_lg, REAL **res_m_ug, REAL **res_m_ls, REAL **res_m_us)
	{

	int Nn = res->dim->Nn;
	int *nx = res->dim->nx;
	int *nu = res->dim->nu;
	int *nb = res->dim->nb;
	int *ng = res->dim->ng;
	int *ns = res->dim->ns;

	int ii, idx;

	for(ii=0; ii<Nn; ii++)
		{
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
		}
	
	for(ii=0; ii<Nn-1; ii++)
		{
		// dynamics
		idx = ii+1;
		CVT_STRVEC2VEC(nx[idx], res->res_b+ii, 0, res_b[ii]);
		}

	return;

	}



void CVT_TREE_OCP_QP_RES_TO_ROWMAJ(struct TREE_OCP_QP_RES *res, REAL **res_r, REAL **res_q, REAL **res_ls, REAL **res_us, REAL **res_b, REAL **res_d_lb, REAL **res_d_ub, REAL **res_d_lg, REAL **res_d_ug, REAL **res_d_ls, REAL **res_d_us, REAL **res_m_lb, REAL **res_m_ub, REAL **res_m_lg, REAL **res_m_ug, REAL **res_m_ls, REAL **res_m_us)
	{

	int Nn = res->dim->Nn;
	int *nx = res->dim->nx;
	int *nu = res->dim->nu;
	int *nb = res->dim->nb;
	int *ng = res->dim->ng;
	int *ns = res->dim->ns;

	int ii, idx;

	for(ii=0; ii<Nn; ii++)
		{
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
		}
	
	for(ii=0; ii<Nn-1; ii++)
		{
		// dynamics
		idx = ii+1;
		CVT_STRVEC2VEC(nx[idx], res->res_b+ii, 0, res_b[ii]);
		}

	return;

	}






