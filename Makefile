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

include ./Makefile.rule

OBJS = 

ifeq ($(TARGET), GENERIC)
OBJS +=
endif

# dense qp
OBJS += dense_qp/d_dense_qp.o dense_qp/d_dense_qp_sol.o dense_qp/d_dense_qp_kkt.o dense_qp/d_dense_qp_ipm_hard.o
OBJS += dense_qp/s_dense_qp.o dense_qp/s_dense_qp_sol.o dense_qp/s_dense_qp_kkt.o dense_qp/s_dense_qp_ipm_hard.o
# ocp qp
OBJS += ocp_qp/d_ocp_qp.o ocp_qp/d_ocp_qp_sol.o ocp_qp/d_ocp_qp_kkt.o ocp_qp/d_ocp_qp_ipm_hard.o
OBJS += ocp_qp/s_ocp_qp.o ocp_qp/s_ocp_qp_sol.o ocp_qp/s_ocp_qp_kkt.o ocp_qp/s_ocp_qp_ipm_hard.o
OBJS += ocp_qp/m_ocp_qp.o                       ocp_qp/m_ocp_qp_kkt.o ocp_qp/m_ocp_qp_ipm_hard.o
# core qp
OBJS += core_qp/d_core_qp_ipm_hard_aux.o core_qp/d_core_qp_ipm_hard.o
OBJS += core_qp/s_core_qp_ipm_hard_aux.o core_qp/s_core_qp_ipm_hard.o
# cond
OBJS += cond/d_cond_aux.o cond/d_cond.o cond/d_part_cond.o

all: clean static_library

static_library: target
	( cd cond; $(MAKE) obj)
	( cd core_qp; $(MAKE) obj)
	( cd dense_qp; $(MAKE) obj)
	( cd ocp_qp; $(MAKE) obj)
	ar rcs libhpipm.a $(OBJS) 
	cp libhpipm.a ./lib/
	@echo
	@echo " libhpipm.a static library build complete."
	@echo

shared_library: target
	( cd cond; $(MAKE) obj)
	( cd core_qp; $(MAKE) obj)
	( cd dense_qp; $(MAKE) obj)
	( cd ocp_qp; $(MAKE) obj)
	gcc -shared -o libhpipm.so $(OBJS)
	cp libhpipm.so ./lib/
	@echo
	@echo " libhpipm.so shared library build complete."
	@echo

target:
	touch ./include/hpipm_target.h
ifeq ($(TARGET), GENERIC)
	echo "#ifndef TARGET_GENERIC" > ./include/hpipm_target.h
	echo "#define TARGET_GENERIC" >> ./include/hpipm_target.h
	echo "#endif" >> ./include/hpipm_target.h
endif

install_static:
	mkdir -p $(PREFIX)/hpipm
	mkdir -p $(PREFIX)/hpipm/lib
	cp -f libhpipm.a $(PREFIX)/hpipm/lib/
	mkdir -p $(PREFIX)/hpipm/include
	cp -f ./include/*.h $(PREFIX)/hpipm/include/

install_shared:
	mkdir -p $(PREFIX)/hpipm
	mkdir -p $(PREFIX)/hpipm/lib
	cp -f libhpipm.so $(PREFIX)/hpipm/lib/
	mkdir -p $(PREFIX)/hpipm/include
	cp -f ./include/*.h $(PREFIX)/hpipm/include/

test_problem:
	cp libhpipm.a ./test_problems/libhpipm.a
	make -C test_problems obj
	@echo
	@echo " Test problem build complete."
	@echo

run:
	./test_problems/test.out

clean:
	rm -f libhpipm.a
	rm -f libhpipm.so
	rm -f ./lib/libhpipm.a
	rm -f ./lib/libhpipm.so
	make -C cond clean
	make -C core_qp clean
	make -C dense_qp clean
	make -C ocp_qp clean
	make -C test_problems clean

