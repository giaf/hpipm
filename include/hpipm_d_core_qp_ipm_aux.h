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

#ifndef HPIPM_S_CORE_QP_IPM_AUX_
#define HPIPM_S_CORE_QP_IPM_AUX_

#ifdef __cplusplus
extern "C" {
#endif

//
void d_compute_Gamma_gamma_qp(double *res_d, double *res_m, struct d_core_qp_ipm_workspace *rws);
//
void d_compute_gamma_qp(double *res_d, double *res_m, struct d_core_qp_ipm_workspace *rws);
//
void d_compute_lam_t_qp(double *res_d, double *res_m, double *dlam, double *dt, struct d_core_qp_ipm_workspace *rws);
//
void d_compute_alpha_qp(struct d_core_qp_ipm_workspace *rws);
//
void d_update_var_qp(struct d_core_qp_ipm_workspace *rws);
//
void d_compute_mu_aff_qp(struct d_core_qp_ipm_workspace *rws);
//
void d_backup_res_m(struct d_core_qp_ipm_workspace *rws);
//
void d_compute_centering_correction_qp(struct d_core_qp_ipm_workspace *rws);
//
void d_compute_centering_qp(struct d_core_qp_ipm_workspace *rws);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // HPIPM_S_CORE_QP_IPM_AUX_
