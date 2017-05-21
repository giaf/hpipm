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

# dense kkt
OBJS += dense_kkt/d_aux.o
#OBJS += dense_kkt/s_aux.o
# ocp kkt
OBJS += ocp_kkt/d_aux.o ocp_kkt/d_cond.o ocp_kkt/d_part_cond.o
OBJS += ocp_kkt/s_aux.o

all: clean static_library

static_library: target
	( cd auxiliary; $(MAKE) obj)
	( cd dense_kkt; $(MAKE) obj)
	( cd ocp_kkt; $(MAKE) obj)
	ar rcs libhpipm.a $(OBJS) 
	cp libhpipm.a ./lib/
	@echo
	@echo " libhpipm.a static library build complete."
	@echo

shared_library: target
	( cd auxiliary; $(MAKE) obj)
	( cd dense_kkt; $(MAKE) obj)
	( cd ocp_kkt; $(MAKE) obj)
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
	make -C auxiliary clean
	make -C dense_kkt clean
	make -C ocp_kkt clean
	make -C test_problems clean

