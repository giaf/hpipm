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


#ifndef HPIPM_TREE_H_
#define HPIPM_TREE_H_

#ifdef __cplusplus
extern "C" {
#endif



struct node
	{
	int *kids;  // 64 bits
	int idx;    // 32 bits
	int dad;    // 32 bits
	int nkids;  // 32 bits
	int stage;  // 32 bits
	int real;   // 32 bits
	int idxkid; // 32 bits // XXX needed ???
	// total     256 bits
	};



struct tree
	{
	struct node *root; // pointer to root
	int *kids; // pointer to array of kids
	int Nn; // numer of nodes
	int memsize;
	};



#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // HPIPM_TREE_H_
