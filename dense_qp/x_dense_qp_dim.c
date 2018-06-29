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



int MEMSIZE_DENSE_QP_DIM()
	{

	int size = 0;

	size = (size+8-1)/8*8;

	return size;

	}



void CREATE_DENSE_QP_DIM(struct DENSE_QP_DIM *size, void *memory)
	{

	size->memsize = MEMSIZE_DENSE_QP_DIM();

	return;

	}


void CVT_INT_TO_DENSE_QP_DIM(int nv, int ne, int nb, int ng, int nsb, int nsg, struct DENSE_QP_DIM *size)
	{

	size->nv = nv;
	size->ne = ne;
	size->nb = nb;
	size->ng = ng;
	size->ns = nsb+nsg;
	size->nsb = nsb;
	size->nsg = nsg;

	return;

	}

