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



int MEMSIZE_TREE_OCP_QP_DIM(int Nn)
	{

	int size = 0;

	size += 7*Nn*sizeof(int);

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

	dim->memsize = MEMSIZE_TREE_OCP_QP_DIM(Nn);

	return;

	}


void CVT_INT_TO_TREE_OCP_QP_DIM(struct tree *ttree, int *nx, int *nu, int *nbx, int *nbu, int *ng, int *ns, struct TREE_OCP_QP_DIM *dim)
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
		dim->ns[ii] = ns[ii];

	return;

	}

