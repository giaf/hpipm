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



int MEMSIZE_TREE_OCP_QP_DIM(int Nn)
	{

	int size = 0;

	size += 10*Nn*sizeof(int);

	size = (size+8-1)/8*8;

	return size;

	}



void CREATE_TREE_OCP_QP_DIM(int Nn, struct TREE_OCP_QP_DIM *dim, void *memory)
	{

	// loop index
	int ii;

	char *c_ptr = memory;

	// nx
	dim->nx = (int *) c_ptr;
	c_ptr += Nn*sizeof(int);
	// nu
	dim->nu = (int *) c_ptr;
	c_ptr += Nn*sizeof(int);
	// nb
	dim->nb = (int *) c_ptr;
	c_ptr += Nn*sizeof(int);
	// nbx
	dim->nbx = (int *) c_ptr;
	c_ptr += Nn*sizeof(int);
	// nbu
	dim->nbu = (int *) c_ptr;
	c_ptr += Nn*sizeof(int);
	// ng
	dim->ng = (int *) c_ptr;
	c_ptr += Nn*sizeof(int);
	// ns
	dim->ns = (int *) c_ptr;
	c_ptr += Nn*sizeof(int);
	// nsbx
	dim->nsbx = (int *) c_ptr;
	c_ptr += Nn*sizeof(int);
	// nsbu
	dim->nsbu = (int *) c_ptr;
	c_ptr += Nn*sizeof(int);
	// nsg
	dim->nsg = (int *) c_ptr;
	c_ptr += Nn*sizeof(int);

	// Nn
	dim->Nn = Nn;

	dim->memsize = MEMSIZE_TREE_OCP_QP_DIM(Nn);

	return;

	}


void CVT_INT_TO_TREE_OCP_QP_DIM(struct tree *ttree, int *nx, int *nu, int *nbx, int *nbu, int *ng, int *nsbx, int *nsbu, int *nsg, struct TREE_OCP_QP_DIM *dim)
	{

	// loop index
	int ii;

	// tree
	dim->ttree = ttree;

	// Nn
	int Nn = ttree->Nn;
	dim->Nn = ttree->Nn;

	// copy qp dim
	for(ii=0; ii<Nn; ii++)
		dim->nx[ii] = nx[ii];
	for(ii=0; ii<Nn; ii++)
		dim->nu[ii] = nu[ii];
	for(ii=0; ii<Nn; ii++)
		dim->nb[ii] = nbx[ii]+nbu[ii];
	for(ii=0; ii<Nn; ii++)
		dim->nbx[ii] = nbx[ii];
	for(ii=0; ii<Nn; ii++)
		dim->nbu[ii] = nbu[ii];
	for(ii=0; ii<Nn; ii++)
		dim->ng[ii] = ng[ii];
	for(ii=0; ii<Nn; ii++)
		dim->ns[ii] = nsbx[ii]+nsbu[ii]+nsg[ii];
	for(ii=0; ii<Nn; ii++)
		dim->nsbx[ii] = nsbx[ii];
	for(ii=0; ii<Nn; ii++)
		dim->nsbu[ii] = nsbu[ii];
	for(ii=0; ii<Nn; ii++)
		dim->nsg[ii] = nsg[ii];

	return;

	}

