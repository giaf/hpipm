/**************************************************************************************************
*                                                                                                 *
* This file is part of HPIPM.                                                                     *
*                                                                                                 *
* HPIPM -- High Performance Interior Point Method.                                                *
* Copyright (C) 2017 by Gianluca Frison.                                                          *
* Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              *
* All rights reserved.                                                                            *
*                                                                                                 *
<<<<<<< HEAD
* HPIPM is free software; you can redistribute it and/or                                          *
=======
* HPMPC is free software; you can redistribute it and/or                                          *
>>>>>>> c012862cce654620c1a17ee48746163031e73a9a
* modify it under the terms of the GNU Lesser General Public                                      *
* License as published by the Free Software Foundation; either                                    *
* version 2.1 of the License, or (at your option) any later version.                              *
*                                                                                                 *
<<<<<<< HEAD
* HPIPM is distributed in the hope that it will be useful,                                        *
=======
* HPMPC is distributed in the hope that it will be useful,                                        *
>>>>>>> c012862cce654620c1a17ee48746163031e73a9a
* but WITHOUT ANY WARRANTY; without even the implied warranty of                                  *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                                            *
* See the GNU Lesser General Public License for more details.                                     *
*                                                                                                 *
* You should have received a copy of the GNU Lesser General Public                                *
<<<<<<< HEAD
* License along with HPIPM; if not, write to the Free Software                                    *
=======
* License along with HPMPC; if not, write to the Free Software                                    *
>>>>>>> c012862cce654620c1a17ee48746163031e73a9a
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA                  *
*                                                                                                 *
* Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             *
*                                                                                                 *
**************************************************************************************************/

#ifndef HPIPM_X_OCP_QP_IPM_H_
#define HPIPM_X_OCP_QP_IPM_H_



#ifdef __cplusplus
extern "C" {
#endif



enum ocp_qp_ipm_mode
	{
	SPEED_ABS, // focus on speed, absolute IPM formulation
	SPEED, // focus on speed, relative IPM formulation
	BALANCE, // balanced mode, relative IPM formulation
	ROBUST, // focus on robustness, relative IPM formulation
	};



#ifdef __cplusplus
} /* extern "C" */
#endif



#endif // HPIPM_X_OCP_QP_IPM_H_
