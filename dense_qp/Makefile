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

OBJS =

ifeq ($(TARGET), GENERIC)
OBJS +=
endif

OBJS += d_dense_qp.o d_dense_qp_sol.o d_dense_qp_kkt.o d_dense_qp_ipm_hard.o
OBJS += s_dense_qp.o s_dense_qp_sol.o s_dense_qp_kkt.o s_dense_qp_ipm_hard.o

obj: $(OBJS)

clean:
	rm -f *.o
	rm -f *.s

d_dense_qp.o: d_dense_qp.c x_dense_qp.c
s_dense_qp.o: s_dense_qp.c x_dense_qp.c
d_dense_qp_sol.o: d_dense_qp_sol.c x_dense_qp_sol.c
s_dense_qp_sol.o: s_dense_qp_sol.c x_dense_qp_sol.c
d_dense_qp_kkt.o: d_dense_qp_kkt.c x_dense_qp_kkt.c
s_dense_qp_kkt.o: s_dense_qp_kkt.c x_dense_qp_kkt.c
d_dense_qp_ipm_hard.o: d_dense_qp_ipm_hard.c x_dense_qp_ipm_hard.c
s_dense_qp_ipm_hard.o: s_dense_qp_ipm_hard.c x_dense_qp_ipm_hard.c
