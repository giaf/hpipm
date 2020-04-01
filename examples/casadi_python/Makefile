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

PYTHON = python
#PYTHON = python3

PYTHON_VERSION_FULL := $(wordlist 2,4,$(subst ., ,$(shell $(PYTHON) --version 2>&1)))
PYTHON_VERSION_MAJOR := $(word 1,${PYTHON_VERSION_FULL})
PYTHON_VERSION_MINOR := $(word 2,${PYTHON_VERSION_FULL})
#PYTHON_VERSION_PATCH := $(word 3,${PYTHON_VERSION_FULL})

PYTHON_VERSION = $(PYTHON_VERSION_MAJOR)$(PYTHON_VERSION_MINOR)
# TODO fallback to supported versions in casadi: 27, 35, 36, 37
#PYTHON_VERSION = 27
#PYTHON_VERSION = 35
#PYTHON_VERSION = 36
#PYTHON_VERSION = 37

CASADI_VERSION = 3.5.1
CASADI_GITHUB_RELEASES = https://github.com/casadi/casadi/releases/download/$(CASADI_VERSION)
CASADI_SOURCE = casadi-linux-py$(PYTHON_VERSION)-v$(CASADI_VERSION)-64bit.tar.gz
CASADI_URL = $(CASADI_GITHUB_RELEASES)/$(CASADI_SOURCE)

#my_cmd.python.2 := python2 some_script.py2
#my_cmd.python.3 := python3 some_script.py3
#my_cmd := ${my_cmd.python.${PYTHON_VERSION_MAJOR}}

all:
	$(PYTHON) example_PMSM.py

### all individual tests ###

example_PMSM \
: %:
	$(PYTHON) $@.py

get_casadi: | $(CASADI_SOURCE)
	mkdir -p casadi
	tar -xf $(CASADI_SOURCE) -C casadi

$(CASADI_SOURCE):
	wget $(CASADI_URL) -O $(CASADI_SOURCE)

clean:
	rm -rf casadi
	rm -f dummy.txt
	rm -f *.pyc

.PHONY : all


