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

#ifdef __cplusplus
extern "C" {
#endif

struct d_rk_data
	{
	double *A_rk; // A in butcher tableau
	double *B_rk; // b in butcher tableau
	double *C_rk; // c in butcher tableau
	int expl; // erk vs irk
	int ns; // number of stages
	int memsize;
	};



//
int d_memsize_rk_data(int ns);
//
void d_create_rk_data(int ns, struct d_rk_data *rk_data, void *memory);
//
void d_cvt_colmaj_to_rk_data(int expl, double *A_rk, double *B_rk, double *C_rk, struct d_rk_data *rk_data);
//
void d_cvt_rowmaj_to_rk_data(int expl, double *A_rk, double *B_rk, double *C_rk, struct d_rk_data *rk_data);

#ifdef __cplusplus
} /* extern "C" */
#endif
