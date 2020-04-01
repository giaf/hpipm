###################################################################################################
#                                                                                                 #
# This file is part of HPIPM.                                                                     #
#                                                                                                 #
# HPIPM -- High-Performance Interior Point Method.                                                #
# Copyright (C) 2019 by Gianluca Frison.                                                          #
# Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              #
# All rights reserved.                                                                            #
#                                                                                                 #
# The 2-Clause BSD License                                                                        #
#                                                                                                 #
# Redistribution and use in source and binary forms, with or without                              #
# modification, are permitted provided that the following conditions are met:                     #
#                                                                                                 #
# 1. Redistributions of source code must retain the above copyright notice, this                  #
#    list of conditions and the following disclaimer.                                             #
# 2. Redistributions in binary form must reproduce the above copyright notice,                    #
#    this list of conditions and the following disclaimer in the documentation                    #
#    and/or other materials provided with the distribution.                                       #
#                                                                                                 #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND                 #
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED                   #
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE                          #
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR                 #
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES                  #
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;                    #
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND                     #
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT                      #
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS                   #
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                                    #
#                                                                                                 #
# Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             #
#                                                                                                 #
###################################################################################################

add_executable(d_cond test_d_cond.c d_tools.c)
target_include_directories(d_cond PRIVATE "${HPIPM_INCLUDE_DIR}" "${BLASFEO_INCLUDE_DIR}")
target_link_libraries(d_cond hpipm blasfeo m)

add_executable(d_part_cond test_d_part_cond.c d_tools.c)
target_include_directories(d_part_cond PRIVATE "${HPIPM_INCLUDE_DIR}" "${BLASFEO_INCLUDE_DIR}")
target_link_libraries(d_part_cond hpipm blasfeo m)

add_executable(d_dense_qp test_d_dense.c)
target_include_directories(d_dense_qp PRIVATE "${HPIPM_INCLUDE_DIR}" "${BLASFEO_INCLUDE_DIR}")
target_link_libraries(d_dense_qp hpipm blasfeo m)

#add_executable(s_dense_qp test_s_dense.c)
#target_link_libraries(s_dense_qp hpipm blasfeo m)

add_executable(d_ocp_qp test_d_ocp.c d_tools.c)
target_include_directories(d_ocp_qp PRIVATE "${HPIPM_INCLUDE_DIR}" "${BLASFEO_INCLUDE_DIR}")
target_link_libraries(d_ocp_qp hpipm blasfeo m)

#add_executable(s_ocp_qp test_s_ocp.c s_tools.c)
#target_link_libraries(s_ocp_qp hpipm blasfeo m)

#add_executable(m_ocp_qp test_m_ocp.c d_tools.c)
#target_link_libraries(m_ocp_qp hpipm blasfeo m)

add_executable(d_tree_ocp_qp test_d_tree_ocp.c d_tools.c)
target_include_directories(d_tree_ocp_qp PRIVATE "${HPIPM_INCLUDE_DIR}" "${BLASFEO_INCLUDE_DIR}")
target_link_libraries(d_tree_ocp_qp hpipm blasfeo m)

#add_executable(s_tree_ocp_qp test_s_tree_ocp.c s_tools.c)
#target_link_libraries(s_tree_ocp_qp hpipm blasfeo m)

