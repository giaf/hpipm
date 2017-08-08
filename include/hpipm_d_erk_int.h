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



struct d_erk_data
	{
	double *A_rk;
	double *B_rk;
	double *C_rk;
	int ns;
	int memsize;
	};



struct d_erk_workspace
	{
	struct d_erk_data *erk_data;
	double *K;
	double *xt;
	int nx;
	int memsize;
	};



struct d_erk_args
	{
	double h;
	int steps;
	};



//
int d_memsize_erk_data(int ns);
//
void d_create_erk_data(int ns, struct d_erk_data *erk_data, void *memory);
//
void d_cvt_colmaj_to_erk_data(double *A_rk, double *B_rk, double *C_rk, struct d_erk_data *erk_data);
//
void d_cvt_rowmaj_to_erk_data(double *A_rk, double *B_rk, double *C_rk, struct d_erk_data *erk_data);
//
int d_memsize_erk_int(struct d_erk_data *erk_data, int nx);
//
void d_create_erk_int(struct d_erk_data *erk_data, int nx, struct d_erk_workspace *workspace, void *memory);
//
void d_erk_int(double *x0, double *p, double *xe, void (*ode)(int t, double *x, double *p, void *ode_args, double *xdot), void *ode_args, struct d_erk_args *erk_args, struct d_erk_workspace *workspace);
