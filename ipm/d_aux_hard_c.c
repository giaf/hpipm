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



#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_d_aux.h>
#include <blasfeo_d_blas.h>



void d_init_var_hard(int N, int *nx, int *nu, int *nb, int **hidxb, int *ng, struct d_strvec *hsux, struct d_strvec *hspi, struct d_strmat *hsDCt, struct d_strvec *hsdb, struct d_strvec *hst, struct d_strvec *hslam, double mu0, int warm_start)
	{

	int jj, ll, ii;

	double *ptr_ux, *ptr_pi, *ptr_db, *ptr_t, *ptr_lam;

	int nb0, ng0, nt0;
	
	double thr0 = 0.1; // minimum vale of t (minimum distance from a constraint)


	// cold start
	if(warm_start==0)
		{
		for(jj=0; jj<=N; jj++)
			{
			ptr_ux = hsux[jj].pa;
			for(ll=0; ll<nu[jj]+nx[jj]; ll++)
				{
				ptr_ux[ll] = 0.0;
				}
			}
		}


	// check bounds & initialize multipliers
	for(jj=0; jj<=N; jj++)
		{
		nb0 = nb[jj];
		nt0 = nb[jj]+ng[jj];
		ptr_ux = hsux[jj].pa;
		ptr_db = hsdb[jj].pa;
		ptr_lam = hslam[jj].pa;
		ptr_t = hst[jj].pa;
		for(ll=0; ll<nb0; ll++)
			{
			ptr_t[ll]     = - ptr_db[ll]     + ptr_ux[hidxb[jj][ll]];
			ptr_t[nt0+ll] =   ptr_db[nt0+ll] - ptr_ux[hidxb[jj][ll]];
			if(ptr_t[ll] < thr0)
				{
				if(ptr_t[nt0+ll] < thr0)
					{
					ptr_ux[hidxb[jj][ll]] = ( - ptr_db[nt0+ll] + ptr_db[ll])*0.5;
					ptr_t[ll]     = thr0; //- hdb[jj][ll]     + hux[jj][hidxb[jj][ll]];
					ptr_t[nt0+ll] = thr0; //  hdb[jj][nt0+ll] - hux[jj][hidxb[jj][ll]];
					}
				else
					{
					ptr_t[ll] = thr0;
					ptr_ux[hidxb[jj][ll]] = ptr_db[ll] + thr0;
					}
				}
			else if(ptr_t[nt0+ll] < thr0)
				{
				ptr_t[nt0+ll] = thr0;
				ptr_ux[hidxb[jj][ll]] = ptr_db[nt0+ll] - thr0;
				}
			ptr_lam[ll]     = mu0/ptr_t[ll];
			ptr_lam[nt0+ll] = mu0/ptr_t[nt0+ll];
			}
		}


	// initialize pi
	for(jj=1; jj<=N; jj++)
		{
		ptr_pi = hspi[jj].pa;
		for(ll=0; ll<nx[jj]; ll++)
			ptr_pi[ll] = 0.0; // initialize multipliers to zero
		}


	// TODO find a better way to initialize general constraints
	for(jj=0; jj<=N; jj++)
		{
		nb0 = nb[jj];
		ng0 = ng[jj];
		nt0 = nb0 + ng0;
		if(ng0>0)
			{
			ptr_t   = hst[jj].pa;
			ptr_lam = hslam[jj].pa;
			ptr_db  = hsdb[jj].pa;
			dgemv_t_libstr(nu[jj]+nx[jj], ng0, 1.0, &hsDCt[jj], 0, 0, &hsux[jj], 0, 0.0, &hst[jj], nb0, &hst[jj], nb0);
			for(ll=nb0; ll<nb0+ng0; ll++)
				{
				ptr_t[ll+nt0] = - ptr_t[ll];
				ptr_t[ll]     -= ptr_db[ll];
				ptr_t[ll+nt0] += ptr_db[ll+nt0];
				ptr_t[ll]     = fmax( thr0, ptr_t[ll] );
				ptr_t[nt0+ll] = fmax( thr0, ptr_t[nt0+ll] );
				ptr_lam[ll]     = mu0/ptr_t[ll];
				ptr_lam[nt0+ll] = mu0/ptr_t[nt0+ll];
				}
			}
		}

	}




