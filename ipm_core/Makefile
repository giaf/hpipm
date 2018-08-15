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
OBJS += d_core_qp_ipm_aux.o
OBJS += s_core_qp_ipm_aux.o
endif
ifeq ($(TARGET), AVX)
OBJS += d_core_qp_ipm_aux_avx.o
OBJS += s_core_qp_ipm_aux_avx.o
endif

OBJS += d_core_qp_ipm.o
OBJS += s_core_qp_ipm.o

obj: $(OBJS)

clean:
	rm -f *.o
	rm -f *.s

d_core_qp_ipm.o: d_core_qp_ipm.c x_core_qp_ipm.c
s_core_qp_ipm.o: s_core_qp_ipm.c x_core_qp_ipm.c
d_core_qp_ipm_aux.o: d_core_qp_ipm_aux.c x_core_qp_ipm_aux.c
s_core_qp_ipm_aux.o: s_core_qp_ipm_aux.c x_core_qp_ipm_aux.c
