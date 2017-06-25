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

#OBJS_TEST = d_tools.o test_d_cond.o
#OBJS_TEST = test_d_dense.o
#OBJS_TEST = test_s_dense.o
OBJS_TEST = d_tools.o test_d_ocp.o
#OBJS_TEST = s_tools.o test_s_ocp.o
#OBJS_TEST = d_tools.o test_m_ocp.o

obj: $(OBJS_TEST)
	$(CC) -o test.out $(OBJS_TEST) -L. libhpipm.a $(LIBS) #-pg

clean:
	rm -f *.o
	rm -f test.out
	rm -f libhpipm.a

