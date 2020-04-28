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



#include <stdlib.h>
#include <stdio.h>

#include <blasfeo_target.h>
#include <blasfeo_common.h>
#include <blasfeo_d_aux.h>
#include <blasfeo_d_blas.h>

#include <hpipm_d_sim_rk.h>
#include <hpipm_d_sim_erk.h>
#include <hpipm_aux_mem.h>



#define BLASFEO_AXPY blasfeo_daxpy
#define BLASFEO_VEC blasfeo_dvec
#define REAL double
#define SIM_ERK_ARG d_sim_erk_arg
#define SIM_RK_DATA d_sim_rk_data
#define SIM_ERK_WS d_sim_erk_ws



#define SIM_ERK_ARG_MEMSIZE d_sim_erk_arg_memsize
#define SIM_ERK_ARG_CREATE d_sim_erk_arg_create
#define SIM_ERK_ARG_SET_ALL d_sim_erk_arg_set_all

#define SIM_ERK_WS_MEMSIZE d_sim_erk_ws_memsize
#define SIM_ERK_WS_CREATE d_sim_erk_ws_create
#define SIM_ERK_WS_SET_ALL d_sim_erk_ws_set_all
#define SIM_ERK_WS_SET_NF d_sim_erk_ws_set_nf
#define SIM_ERK_WS_SET_X d_sim_erk_ws_set_x
#define SIM_ERK_WS_SET_FS d_sim_erk_ws_set_fs
#define SIM_ERK_WS_GET_X d_sim_erk_ws_get_x
#define SIM_ERK_WS_SET_P d_sim_erk_ws_set_p
#define SIM_ERK_WS_SET_ODE d_sim_erk_ws_set_ode
#define SIM_ERK_WS_SET_VDE_FOR d_sim_erk_ws_set_vde_for
#define SIM_ERK_WS_SET_ODE_ARGS d_sim_erk_ws_set_ode_args
#define SIM_ERK_WS_GET_FS d_sim_erk_ws_get_fs
#define SIM_ERK_SOLVE d_sim_erk_solve



#include "x_sim_erk.c"
