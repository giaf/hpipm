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

include ../Makefile.rule

LIBS =

LIBS += $(BLASFEO_PATH)/lib/libblasfeo.a -lm

OBJS_TEST =
#OBJS_TEST = d_tools.o test_d_cond.o
#OBJS_TEST = d_tools.o test_d_part_cond.o
#OBJS_TEST = test_d_dense.o
#OBJS_TEST = test_s_dense.o
OBJS_TEST += test_d_benchmark.o
#OBJS_TEST = d_tools.o test_d_ocp2.o
#OBJS_TEST = s_tools.o test_s_ocp.o
#OBJS_TEST = d_tools.o test_m_ocp.o
#OBJS_TEST = d_tools.o test_d_tree_ocp.o
#OBJS_TEST = s_tools.o test_s_tree_ocp.o
#OBJS_TEST = d_tools.o test_d_sim.o
#OBJS_TEST = d_tools.o test_d_nlp_sqp.o
#OBJS_TEST = d_tools.o test_d_nlp_ipm.o
#OBJS_TEST = test_van_der_pol_sqp.o
#OBJS_TEST = test_van_der_pol_ipm.o
#OBJS_TEST = test_van_der_pol_hyb.o
#OBJS_TEST = pendulum/vde_forw_model.o test_d_pendulum_sqp.o
#OBJS_TEST = pendulum/vde_forw_model.o test_d_pendulum_ipm.o
#OBJS_TEST = pendulum/ode_model.o pendulum/vde_forw_model.o pendulum/vde_adj_model.o test_d_pendulum_hyb.o

obj: $(OBJS_TEST)
	$(CC) -o test.out $(OBJS_TEST) -L. ../libhpipm.a $(LIBS) #-pg

clean:
	rm -f *.o
	rm -f test.out
	rm -f libhpipm.a
