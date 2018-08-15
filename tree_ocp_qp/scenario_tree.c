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



#include "../include/hpipm_tree.h"
#include "../include/hpipm_scenario_tree.h"



static int ipow(int base, int exp)
	{
	int result = 1;
	while(exp)
		{
		if(exp & 1)
			result *= base;
		exp >>= 1;
		base *= base;
		}
	return result;
	}



static int sctree_node_number(int md, int Nr, int Nh)
	{
	int n_nodes;
	if(md==1) // i.e. standard block-banded structure
		n_nodes = Nh+1;
	else
		n_nodes = (Nh-Nr)*ipow(md,Nr) + (ipow(md,Nr+1)-1)/(md-1);
	return n_nodes;
	}



int memsize_sctree(int md, int Nr, int Nh)
	{

	int Nn = sctree_node_number(md, Nr, Nh);

	int size = 0;

	size += Nn*sizeof(struct node); // root
	size += Nn*sizeof(int); // kids

	return size;

	}



void create_sctree(int md, int Nr, int Nh, struct sctree *st, void *memory)
	{

	st->memsize = memsize_sctree(md, Nr, Nh);

	int Nn = sctree_node_number(md, Nr, Nh);

	st->md = md;
	st->Nr = Nr;
	st->Nh = Nh;
	st->Nn = Nn;

	struct node *n_ptr = (struct node *) memory;
	st->root = n_ptr;
	n_ptr += Nn; // root

	int *i_ptr = (int *) n_ptr;
	st->kids = i_ptr;
	i_ptr += Nn; // kids

	int ii;
	int idx, dad, stage, real, nkids, idxkid;
	int tkids;
	struct node *node0, *node1;

	tkids = 0;
	idxkid = 0;

	// root
	node0 = st->root+0;
	node0->idx = 0;
	node0->stage = 0;
	node0->dad = -1;
	node0->real = -1;
	node0->idxkid = 0;

	// kids
	for(idx=0; idx<Nn; idx++)
		{
		node0 = st->root+idx;
		stage = node0->stage;
		if(stage<Nr)
			nkids = md;
		else if(stage<Nh)
			nkids = 1;
		else 
			nkids = 0;
		node0->nkids = nkids;
		if(nkids>0)
			{
			node0->kids = st->kids+tkids;
			tkids += nkids;
			if(nkids>1)
				{
				for(ii=0; ii<nkids; ii++)
					{
					idxkid++;
					node0->kids[ii] = idxkid;
					node1 = st->root+idxkid;
					node1->idx = idxkid;
					node1->stage = stage+1;
					node1->dad = idx;
					node1->real = ii;
					node1->idxkid = ii;
					}
				}
			else // nkids==1
				{
				idxkid++;
				node0->kids[0] = idxkid;
				node1 = st->root+idxkid;
				node1->idx = idxkid;
				node1->stage = stage+1;
				node1->dad = idx;
				node1->real = node0->real;
				node1->idxkid = 0;
				}
			}
		}

	return;

	}



void cast_sctree2tree(struct sctree *st, struct tree *tt)
	{

	tt->root = st->root;
	tt->kids = st->kids;
	tt->Nn = st->Nn;
	tt->memsize = st->memsize;

	return;

	}
