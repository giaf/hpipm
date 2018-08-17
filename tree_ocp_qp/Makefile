###################################################################################################
#                                                                                                 #
# This file is part of HPIPM.                                                                     #
#                                                                                                 #
# HPIPM -- High-Performance Interior Point Method.                                                #
# Copyright (C) 2017-2018 by Gianluca Frison.                                                     #
# Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              #
# All rights reserved.                                                                            #
#                                                                                                 #
# This program is free software: you can redistribute it and/or modify                            #
# it under the terms of the GNU General Public License as published by                            #
# the Free Software Foundation, either version 3 of the License, or                               #
# (at your option) any later version                                                              #.
#                                                                                                 #
# This program is distributed in the hope that it will be useful,                                 #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                  #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                   #
# GNU General Public License for more details.                                                    #
#                                                                                                 #
# You should have received a copy of the GNU General Public License                               #
# along with this program.  If not, see <https://www.gnu.org/licenses/>.                          #
#                                                                                                 #
# The authors designate this particular file as subject to the "Classpath" exception              #
# as provided by the authors in the LICENSE file that accompained this code.                      #
#                                                                                                 #
# Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             #
#                                                                                                 #
###################################################################################################

include ../Makefile.rule

OBJS =

ifeq ($(TARGET), GENERIC)
OBJS +=
endif

# common
OBJS += scenario_tree.o
# double
OBJS += d_tree_ocp_qp_dim.o
OBJS += d_tree_ocp_qp.o
OBJS += d_tree_ocp_qp_sol.o
OBJS += d_tree_ocp_qp_res.o
OBJS += d_tree_ocp_qp_kkt.o
OBJS += d_tree_ocp_qp_ipm.o
# single
OBJS += s_tree_ocp_qp_dim.o
OBJS += s_tree_ocp_qp.o
OBJS += s_tree_ocp_qp_sol.o
OBJS += s_tree_ocp_qp_res.o
OBJS += s_tree_ocp_qp_kkt.o
OBJS += s_tree_ocp_qp_ipm.o

obj: $(OBJS)

clean:
	rm -f *.o
	rm -f *.s

d_tree_ocp_qp_dim.o: d_tree_ocp_qp_dim.c x_tree_ocp_qp_dim.c
s_tree_ocp_qp_dim.o: s_tree_ocp_qp_dim.c x_tree_ocp_qp_dim.c
d_tree_ocp_qp.o: d_tree_ocp_qp.c x_tree_ocp_qp.c
s_tree_ocp_qp.o: s_tree_ocp_qp.c x_tree_ocp_qp.c
d_tree_ocp_qp_sol.o: d_tree_ocp_qp_sol.c x_tree_ocp_qp_sol.c
s_tree_ocp_qp_sol.o: s_tree_ocp_qp_sol.c x_tree_ocp_qp_sol.c
d_tree_ocp_qp_res.o: d_tree_ocp_qp_res.c x_tree_ocp_qp_res.c
s_tree_ocp_qp_res.o: s_tree_ocp_qp_res.c x_tree_ocp_qp_res.c
d_tree_ocp_qp_kkt.o: d_tree_ocp_qp_kkt.c x_tree_ocp_qp_kkt.c
s_tree_ocp_qp_kkt.o: s_tree_ocp_qp_kkt.c x_tree_ocp_qp_kkt.c
d_tree_ocp_qp_ipm.o: d_tree_ocp_qp_ipm.c x_tree_ocp_qp_ipm.c
s_tree_ocp_qp_ipm.o: s_tree_ocp_qp_ipm.c x_tree_ocp_qp_ipm.c
