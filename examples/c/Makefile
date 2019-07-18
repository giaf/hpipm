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

include ../../Makefile.rule

LIBS =

LIBS += $(BLASFEO_PATH)/lib/libblasfeo.a

ifeq ($(REF_BLAS), 0)
LIBS += -lm 
endif
ifeq ($(REF_BLAS), OPENBLAS)
LIBS += /opt/openblas/lib/libopenblas.a -pthread -lgfortran -lm
endif
ifeq ($(REF_BLAS), BLIS)
LIBS += /opt/netlib/liblapack.a /opt/blis/lib/libblis.a -lgfortran -lm -fopenmp
endif
ifeq ($(REF_BLAS), NETLIB)
LIBS += /opt/netlib/liblapack.a /opt/netlib/libblas.a -lgfortran -lm
endif
ifeq ($(REF_BLAS), MKL)
LIBS += -Wl,--start-group /opt/intel/mkl/lib/intel64/libmkl_gf_lp64.a /opt/intel/mkl/lib/intel64/libmkl_core.a /opt/intel/mkl/lib/intel64/libmkl_sequential.a -Wl,--end-group -ldl -lpthread -lm
endif
ifeq ($(REF_BLAS), ATLAS)
LIBS += /opt/atlas/lib/liblapack.a /opt/atlas/lib/libcblas.a /opt/atlas/lib/libf77blas.a /opt/atlas/lib/libatlas.a -lgfortran -lm
endif

OBJS_DATA = data/getting_started_data.o 
#OBJS_DATA = data/ocp_qp_data.o 
#OBJS_DATA = ../matlab_octave_mex/qp_data.o
#OBJS_DATA = ../python/qp_data.o

#OBJS_EXMP = $(OBJS_DATA) example_d_ocp_qp_unconstr.o
OBJS_EXMP = $(OBJS_DATA) example_d_ocp_qp.o
#OBJS_EXMP = $(OBJS_DATA) example_d_ocp_qp_part_cond.o

obj: $(OBJS_EXMP)
	$(CC) -o  example.out $(OBJS_EXMP) -L. libhpipm.a $(LIBS) -I../../include #-pg

clean:
	rm -f *.o
	rm -f data/*.o
	rm -f example.out
	rm -f libhpipm.a


