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


void SET_DENSE_QP_DIM(const char *field_name, int value, struct DENSE_QP_DIM *dim)
	{
	if(hpipm_strcmp(field_name, "nv"))
		{ 
		dim->nv = value;
		}
	else if(hpipm_strcmp(field_name, "ne"))
		{ 
		dim->ne = value;
		}
	else if(hpipm_strcmp(field_name, "nb"))
		{
		dim->nb = value;
		}
	else if(hpipm_strcmp(field_name, "ng"))
		{
		dim->ng = value;
		}
	else if(hpipm_strcmp(field_name, "nsb"))
		{
		dim->nsb = value;
		dim->ns = dim->nsb + dim->nsg;
		}
	else if(hpipm_strcmp(field_name, "nsg"))
		{
		dim->nsg = value;
		dim->ns = dim->nsb + dim->nsg;
		}
	else if(hpipm_strcmp(field_name, "ns"))
		{
		dim->ns = value;
		}
	else 
		{
		printf("error [SET_OCP_QP_DIM]: unknown field name '%s'. Exiting.\n", field_name);
		exit(1);
		}
	return;
	}




