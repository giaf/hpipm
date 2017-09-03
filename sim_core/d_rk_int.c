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



#if defined(RUNTIME_CHECKS)
#include <stdlib.h>
#include <stdio.h>
#endif

#include "../include/hpipm_d_rk_int.h"



int d_memsize_rk_data(int ns)
	{

	int size = 0;

	size += 1*ns*ns*sizeof(double); // A
	size += 2*ns*sizeof(double); // B C

	return size;

	}



void d_create_rk_data(int ns, struct d_rk_data *rk_data, void *mem)
	{

	rk_data->ns = ns;

	double *d_ptr = mem;

	//
	rk_data->A_rk = d_ptr;
	d_ptr += ns*ns;
	//
	rk_data->B_rk = d_ptr;
	d_ptr += ns;
	//
	rk_data->C_rk = d_ptr;
	d_ptr += ns;


	rk_data->memsize = d_memsize_rk_data(ns);


	char *c_ptr = (char *) d_ptr;


#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + rk_data->memsize)
		{
		printf("\nCreate_irk_int: outsize memory bounds!\n\n");
		exit(1);
		}
#endif


	return;

	}



void d_cvt_colmaj_to_rk_data(double *A, double *B, double *C, struct d_rk_data *rk_data)
	{

	int ii, jj;

	int ns = rk_data->ns;
	double *A_rk = rk_data->A_rk;
	double *B_rk = rk_data->B_rk;
	double *C_rk = rk_data->C_rk;

	for(jj=0; jj<ns; jj++)
		for(ii=0; ii<ns; ii++)
			A_rk[ii+ns*jj] = A[ii+ns*jj];
	for(ii=0; ii<ns; ii++)
		B_rk[ii] = B[ii];
	for(ii=0; ii<ns; ii++)
		C_rk[ii] = C[ii];

	return;

	}



void d_cvt_rowmaj_to_rk_data(double *A, double *B, double *C, struct d_rk_data *rk_data)
	{

	int ii, jj;

	int ns = rk_data->ns;
	double *A_rk = rk_data->A_rk;
	double *B_rk = rk_data->B_rk;
	double *C_rk = rk_data->C_rk;

	for(jj=0; jj<ns; jj++)
		for(ii=0; ii<ns; ii++)
			A_rk[jj+ns*ii] = A[ii+ns*jj];
	for(ii=0; ii<ns; ii++)
		B_rk[ii] = B[ii];
	for(ii=0; ii<ns; ii++)
		C_rk[ii] = C[ii];

	return;

	}




