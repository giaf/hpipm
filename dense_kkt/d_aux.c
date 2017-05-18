/**************************************************************************************************
*                                                                                                 *
* This file is part of HPIPM.                                                                     *
*                                                                                                 *
* HPIPM -- High Performance Interior Point Method.                                                *
* Copyright (C) 2017 by Gianluca Frison.                                                          *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
*                                                                                                 *
* HPMPC is free software; you can redistribute it and/or                                          *
* modify it under the terms of the GNU Lesser General Public                                      *
* License as published by the Free Software Foundation; either                                    *
* version 2.1 of the License, or (at your option) any later version.                              *
*                                                                                                 *
* HPMPC is distributed in the hope that it will be useful,                                        *
* but WITHOUT ANY WARRANTY; without even the implied warranty of                                  *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                                            *
* See the GNU Lesser General Public License for more details.                                     *
*                                                                                                 *
* You should have received a copy of the GNU Lesser General Public                                *
* License along with HPMPC; if not, write to the Free Software                                    *
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA                  *
*                                                                                                 *
* Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             *
*                                                                                                 *
**************************************************************************************************/



#if defined(RUNTIME_CHECKS)
#include <stdlib.h>
#endif

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_d_aux.h>

#include "../include/hpipm_dense_kkt.h"



int d_size_dense_qp(int nv, int ne, int nb, int ng)
	{

	int ii;

	int size = 0;

	size += nb*sizeof(int);

	size += d_size_strmat(nv, nv); // Q
	size += d_size_strvec(nv); // q
	size += d_size_strmat(ne, nv); // A
	size += d_size_strvec(ne); // b
	size += 2*d_size_strvec(nb); // lb ub
	size += d_size_strmat(nv, ng); // Ct
	size += 2*d_size_strvec(ng); // lg ug

	size = (size+63)/64*64; // make multiple of typical cache line size
	size += 64; // align to typical cache line size
	
	return size;

	}



void d_create_dense_qp(int nv, int ne, int nb, int ng, struct d_dense_qp *str_out, void *memory)
	{

	int ii;


	// int stuff
	int *i_ptr;
	i_ptr = (int *) memory;

	// idxb
	str_out->idxb = i_ptr;
	i_ptr += nb;


	// align to typical cache line size
	long long l_ptr = (long long) i_ptr;
	l_ptr = (l_ptr+63)/64*64;


	// double stuff
	void *v_ptr;
	v_ptr = (void *) l_ptr;

	// Q
	d_create_strmat(nv, nv, &(str_out->sQ), v_ptr);
	v_ptr += str_out->sQ.memory_size;

	// A
	d_create_strmat(ne, nv, &(str_out->sA), v_ptr);
	v_ptr += str_out->sA.memory_size;

	// Ct
	d_create_strmat(nv, ng, &(str_out->sCt), v_ptr);
	v_ptr += str_out->sCt.memory_size;

	// q
	d_create_strvec(nv, &(str_out->sq), v_ptr);
	v_ptr += str_out->sq.memory_size;

	// b
	d_create_strvec(ne, &(str_out->sb), v_ptr);
	v_ptr += str_out->sb.memory_size;

	// lb
	d_create_strvec(nb, &(str_out->slb), v_ptr);
	v_ptr += str_out->slb.memory_size;

	// ub
	d_create_strvec(nb, &(str_out->sub), v_ptr);
	v_ptr += str_out->sub.memory_size;

	// lg
	d_create_strvec(ng, &(str_out->slg), v_ptr);
	v_ptr += str_out->slg.memory_size;

	// ug
	d_create_strvec(ng, &(str_out->sug), v_ptr);
	v_ptr += str_out->sug.memory_size;

	return;

	}



void d_copy_dense_qp(struct d_dense_qp *str_in, struct d_dense_qp *str_out)
	{

	int ii;

#if defined(RUNTIME_CHECKS)
	if(str_out->nv != str_in->nv)
		{
		printf("\nError : d_copy_dense_qp : str_out->nv != str_out->nv : %d != %d\n\n", str_out->nv, str_in->nv);
		exit(1);
		}
	if(str_out->ne != str_in->ne)
		{
		printf("\nError : d_copy_dense_qp : str_out->ne != str_out->ne : %d != %d\n\n", str_out->ne, str_in->ne);
		exit(1);
		}
	if(str_out->nb != str_in->nb)
		{
		printf("\nError : d_copy_dense_qp : str_out->nb != str_out->nb : %d != %d\n\n", str_out->nb, str_in->nb);
		exit(1);
		}
	if(str_out->ng != str_in->ng)
		{
		printf("\nError : d_copy_dense_qp : str_out->ng != str_out->ng : %d != %d\n\n", str_out->ng, str_in->ng);
		exit(1);
		}
#endif

	for(ii=0; ii<str_in->nb; ii++) str_out->idxb[ii] = str_in->idxb[ii];
	dgecp_libstr(str_in->nv, str_in->nv, &(str_out->sQ), 0, 0, &(str_in->sQ), 0, 0);
	dveccp_libstr(str_in->nv, &(str_out->sq), 0, &(str_in->sq), 0);
	dgecp_libstr(str_in->ne, str_in->nv, &(str_out->sA), 0, 0, &(str_in->sA), 0, 0);
	dveccp_libstr(str_in->ne, &(str_out->sb), 0, &(str_in->sb), 0);
	dgecp_libstr(str_in->nv, str_in->ng, &(str_out->sCt), 0, 0, &(str_in->sCt), 0, 0);
	dveccp_libstr(str_in->nb, &(str_out->slb), 0, &(str_in->slb), 0);
	dveccp_libstr(str_in->nb, &(str_out->sub), 0, &(str_in->sub), 0);
	dveccp_libstr(str_in->ng, &(str_out->slg), 0, &(str_in->slg), 0);
	dveccp_libstr(str_in->ng, &(str_out->sug), 0, &(str_in->sug), 0);

	return;

	}

