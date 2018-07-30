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

# double
OBJS += d_ocp_qp_dim.o
OBJS += d_ocp_qp.o
OBJS += d_ocp_qp_sol.o
OBJS += d_ocp_qp_res.o
OBJS += d_ocp_qp_kkt.o
OBJS += d_ocp_qp_ipm.o
OBJS += d_ocp_qp_utils.o
# single
OBJS += s_ocp_qp_dim.o
OBJS += s_ocp_qp.o
OBJS += s_ocp_qp_sol.o
OBJS += s_ocp_qp_res.o
OBJS += s_ocp_qp_kkt.o
OBJS += s_ocp_qp_ipm.o
OBJS += s_ocp_qp_utils.o
# mixed
#OBJS += m_ocp_qp.o                m_ocp_qp_kkt.o m_ocp_qp_ipm.o

obj: $(OBJS)

clean:
	rm -f *.o
	rm -f *.s

d_ocp_qp_dim.o: d_ocp_qp_dim.c x_ocp_qp_dim.c
s_ocp_qp_dim.o: s_ocp_qp_dim.c x_ocp_qp_dim.c
d_ocp_qp.o: d_ocp_qp.c x_ocp_qp.c
s_ocp_qp.o: s_ocp_qp.c x_ocp_qp.c
d_ocp_qp_sol.o: d_ocp_qp_sol.c x_ocp_qp_sol.c
s_ocp_qp_sol.o: s_ocp_qp_sol.c x_ocp_qp_sol.c
d_ocp_qp_res.o: d_ocp_qp_res.c x_ocp_qp_res.c
s_ocp_qp_res.o: s_ocp_qp_res.c x_ocp_qp_res.c
d_ocp_qp_kkt.o: d_ocp_qp_kkt.c x_ocp_qp_kkt.c
s_ocp_qp_kkt.o: s_ocp_qp_kkt.c x_ocp_qp_kkt.c
d_ocp_qp_ipm.o: d_ocp_qp_ipm.c x_ocp_qp_ipm.c
s_ocp_qp_ipm.o: s_ocp_qp_ipm.c x_ocp_qp_ipm.c
d_ocp_qp_utils.o: d_ocp_qp_utils.c x_ocp_qp_utils.c
s_ocp_qp_utils.o: s_ocp_qp_utils.c x_ocp_qp_utils.c
