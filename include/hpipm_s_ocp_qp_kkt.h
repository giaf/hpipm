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

//
void s_init_var_ocp_qp(struct s_ocp_qp *qp, struct s_ocp_qp_sol *qp_sol, struct s_ocp_qp_ipm_arg *arg, struct s_ocp_qp_ipm_ws *ws);
//
void s_compute_res_ocp_qp(struct s_ocp_qp *qp, struct s_ocp_qp_sol *qp_sol, struct s_ocp_qp_res *res, struct s_ocp_qp_res_ws *ws);
//
void s_compute_lin_res_ocp_qp(struct s_ocp_qp *qp, struct s_ocp_qp_sol *qp_sol, struct s_ocp_qp_sol *qp_step, struct s_ocp_qp_res *res, struct s_ocp_qp_res_ws *ws);
//
void s_fact_solve_kkt_unconstr_ocp_qp(struct s_ocp_qp *qp, struct s_ocp_qp_sol *qp_sol, struct s_ocp_qp_ipm_arg *arg, struct s_ocp_qp_ipm_ws *ws);
//
void s_fact_solve_kkt_step_ocp_qp(struct s_ocp_qp *qp, struct s_ocp_qp_sol *qp_sol, struct s_ocp_qp_ipm_arg *arg, struct s_ocp_qp_ipm_ws *ws);
//
void s_fact_lq_solve_kkt_step_ocp_qp(struct s_ocp_qp *qp, struct s_ocp_qp_sol *qp_sol, struct s_ocp_qp_ipm_arg *arg, struct s_ocp_qp_ipm_ws *ws);
//
void s_solve_kkt_step_ocp_qp(struct s_ocp_qp *qp, struct s_ocp_qp_sol *qp_sol, struct s_ocp_qp_ipm_arg *arg, struct s_ocp_qp_ipm_ws *ws);

#ifdef __cplusplus
} /* extern "C" */
#endif
