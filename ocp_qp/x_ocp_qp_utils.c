
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
* Author: Gianluca Frison, Andrea Zanelli: {gianluca.frison, andrea.zanelli} (at) imtek.uni-freiburg.de                             *
*                                                                                                 *
**************************************************************************************************/

void PRINT_OCP_QP_SOL(struct OCP_QP_SOL *qp_sol, struct OCP_QP_DIM *qp_dim)
    {
        int ii;

        int N   = qp_dim->N;
        int *nx = qp_dim->nx;
        int *nu = qp_dim->nu;
        int *nb = qp_dim->nb;
        int *ng = qp_dim->ng;
        int *ns = qp_dim->ns;

        printf("ux =\n");
        for (ii = 0; ii <= N; ii++)
            blasfeo_print_tran_dvec(nu[ii] + nx[ii] + 2 * ns[ii], &qp_sol->ux[ii], 0);

        printf("pi =\n");
        for (ii = 0; ii < N; ii++) blasfeo_print_tran_dvec(nx[ii + 1], &qp_sol->pi[ii], 0);

        printf("lam =\n");
        for (ii = 0; ii <= N; ii++)
            blasfeo_print_tran_dvec(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], &qp_sol->lam[ii], 0);

        printf("t =\n");
        for (ii = 0; ii <= N; ii++)
            blasfeo_print_exp_tran_dvec(2 * nb[ii] + 2 * ng[ii] + 2 * ns[ii], &qp_sol->t[ii], 0);
    }
