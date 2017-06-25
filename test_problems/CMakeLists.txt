###################################################################################################
#                                                                                                 #
# This file is part of HPIPM.                                                                     #
#                                                                                                 #
# HPIPM -- High Performance Interior Point Method.                                                #
# Copyright (C) 2017 by Gianluca Frison.                                                          #
# Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              #
# All rights reserved.                                                                            #
#                                                                                                 #
# HPMPC is free software; you can redistribute it and/or                                          #
# modify it under the terms of the GNU Lesser General Public                                      #
# License as published by the Free Software Foundation; either                                    #
# version 2.1 of the License, or (at your option) any later version.                              #
#                                                                                                 #
# HPMPC is distributed in the hope that it will be useful,                                        #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                  #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                                            #
# See the GNU Lesser General Public License for more details.                                     #
#                                                                                                 #
# You should have received a copy of the GNU Lesser General Public                                #
# License along with HPMPC; if not, write to the Free Software                                    #
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA                  #
#                                                                                                 #
# Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             #
#                                                                                                 #
###################################################################################################

add_executable(d_dense_qp test_d_dense.c)
target_link_libraries(d_dense_qp hpipm blasfeo m)

add_executable(s_dense_qp test_s_dense.c)
target_link_libraries(s_dense_qp hpipm blasfeo m)

add_executable(d_ocp_qp test_d_ocp.c d_tools.c)
target_link_libraries(d_ocp_qp hpipm blasfeo m)

add_executable(s_ocp_qp test_s_ocp.c s_tools.c)
target_link_libraries(s_ocp_qp hpipm blasfeo m)

add_executable(m_ocp_qp test_m_ocp.c d_tools.c)
target_link_libraries(m_ocp_qp hpipm blasfeo m)

add_executable(d_cond test_d_cond.c d_tools.c)
target_link_libraries(d_cond hpipm blasfeo m)
