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

# select data set for QP and QCQP examples
set(QP_DATA_SET data/getting_started_data.c)
#set(QP_DATA_SET ../matlab_octave/ocp_qp_data.c)
#set(QP_DATA_SET ../python/ocp_qp_data.c)

set(QCQP_DATA_SET data/pmsm_ocp_qcqp_data.c)
#set(QCQP_DATA_SET data/mass_spring_qcqp_data.c)

add_executable(example_d_ocp_qp_unconstr example_d_ocp_qp_unconstr.c ${QP_DATA_SET})
add_executable(example_d_ocp_qp example_d_ocp_qp.c ${QP_DATA_SET})
add_executable(example_d_ocp_qp_x0emb example_d_ocp_qp_x0emb.c ${QP_DATA_SET})
add_executable(example_d_ocp_qp_x0emb_cond example_d_ocp_qp_x0emb_cond.c ${QP_DATA_SET})
add_executable(example_d_ocp_qp_part_cond example_d_ocp_qp_part_cond.c ${QP_DATA_SET})
add_executable(example_d_ocp_qp_x0emb_part_cond example_d_ocp_qp_x0emb_part_cond.c ${QP_DATA_SET})
add_executable(example_d_ocp_qp_sens example_d_ocp_qp_sens.c ${QP_DATA_SET})
add_executable(example_d_ocp_qp_part_cond_sens example_d_ocp_qp_part_cond_sens.c ${QP_DATA_SET})

add_executable(example_d_ocp_qcqp_part_cond example_d_ocp_qp_part_cond.c ${QCQP_DATA_SET})
add_executable(example_d_ocp_qcqp_x0emb_part_cond example_d_ocp_qp_x0emb_part_cond.c ${QCQP_DATA_SET})

target_include_directories(example_d_ocp_qp_unconstr PRIVATE "${HPIPM_INCLUDE_DIR}" "${BLASFEO_INCLUDE_DIR}")
target_include_directories(example_d_ocp_qp PRIVATE "${HPIPM_INCLUDE_DIR}" "${BLASFEO_INCLUDE_DIR}")
target_include_directories(example_d_ocp_qp_x0emb PRIVATE "${HPIPM_INCLUDE_DIR}" "${BLASFEO_INCLUDE_DIR}")
target_include_directories(example_d_ocp_qp_x0emb_cond PRIVATE "${HPIPM_INCLUDE_DIR}" "${BLASFEO_INCLUDE_DIR}")
target_include_directories(example_d_ocp_qp_part_cond PRIVATE "${HPIPM_INCLUDE_DIR}" "${BLASFEO_INCLUDE_DIR}")
target_include_directories(example_d_ocp_qp_x0emb_part_cond PRIVATE "${HPIPM_INCLUDE_DIR}" "${BLASFEO_INCLUDE_DIR}")
target_include_directories(example_d_ocp_qp_sens PRIVATE "${HPIPM_INCLUDE_DIR}" "${BLASFEO_INCLUDE_DIR}")
target_include_directories(example_d_ocp_qp_part_cond_sens PRIVATE "${HPIPM_INCLUDE_DIR}" "${BLASFEO_INCLUDE_DIR}")

target_include_directories(example_d_ocp_qcqp_part_cond PRIVATE "${HPIPM_INCLUDE_DIR}" "${BLASFEO_INCLUDE_DIR}")
target_include_directories(example_d_ocp_qcqp_x0emb_part_cond PRIVATE "${HPIPM_INCLUDE_DIR}" "${BLASFEO_INCLUDE_DIR}")

if(CMAKE_C_COMPILER_ID MATCHES MSVC) # no explicit math library
	target_link_libraries(example_d_ocp_qp_unconstr hpipm blasfeo)
	target_link_libraries(example_d_ocp_qp hpipm blasfeo)
	target_link_libraries(example_d_ocp_qp_x0emb hpipm blasfeo)
	target_link_libraries(example_d_ocp_qp_x0emb_cond hpipm blasfeo)
	target_link_libraries(example_d_ocp_qp_part_cond hpipm blasfeo)
	target_link_libraries(example_d_ocp_qp_x0emb_part_cond hpipm blasfeo)
	target_link_libraries(example_d_ocp_qp_sens hpipm blasfeo)
	target_link_libraries(example_d_ocp_qp_part_cond_sens hpipm blasfeo)

	target_link_libraries(example_d_ocp_qcqp_part_cond hpipm blasfeo)
	target_link_libraries(example_d_ocp_qcqp_x0emb_part_cond hpipm blasfeo)
else() # add explicit math library
	target_link_libraries(example_d_ocp_qp_unconstr hpipm blasfeo m)
	target_link_libraries(example_d_ocp_qp hpipm blasfeo m)
	target_link_libraries(example_d_ocp_qp_x0emb hpipm blasfeo m)
	target_link_libraries(example_d_ocp_qp_x0emb_cond hpipm blasfeo m)
	target_link_libraries(example_d_ocp_qp_part_cond hpipm blasfeo m)
	target_link_libraries(example_d_ocp_qp_x0emb_part_cond hpipm blasfeo m)
	target_link_libraries(example_d_ocp_qp_sens hpipm blasfeo m)
	target_link_libraries(example_d_ocp_qp_part_cond_sens hpipm blasfeo m)

	target_link_libraries(example_d_ocp_qcqp_part_cond hpipm blasfeo m)
	target_link_libraries(example_d_ocp_qcqp_x0emb_part_cond hpipm blasfeo m)
endif()

