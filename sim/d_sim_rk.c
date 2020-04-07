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



#if defined(RUNTIME_CHECKS)
#include <stdlib.h>
#include <stdio.h>
#endif

#include <hpipm_d_sim_rk.h>



int d_sim_rk_data_memsize(int ns)
	{

	int size = 0;

	size += 1*ns*ns*sizeof(double); // A
	size += 2*ns*sizeof(double); // B C

	return size;

	}



void d_sim_rk_data_create(int ns, struct d_sim_rk_data *rk_data, void *mem)
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


	rk_data->memsize = d_sim_rk_data_memsize(ns);


	char *c_ptr = (char *) d_ptr;


#if defined(RUNTIME_CHECKS)
	if(c_ptr > ((char *) mem) + rk_data->memsize)
		{
		printf("\nSIM_RK_DATA_CREATE: outsize memory bounds!\n\n");
		exit(1);
		}
#endif


	return;

	}



void d_sim_rk_data_set_all(int expl, double *A, double *B, double *C, struct d_sim_rk_data *rk_data)
	{

	int ii, jj;

	rk_data->expl = expl;

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
