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
from ctypes import *
import ctypes.util 
import numpy as np
#import faulthandler

#faulthandler.enable()



class hpipm_ocp_qp_solver:
	def __init__(self, qp_dims, arg):

		# load hpipm shared library
		__hpipm   = CDLL('libhpipm.so')
		self.__hpipm = __hpipm

		# set up ipm workspace struct
		sizeof_d_ocp_qp_ipm_workspace = __hpipm.d_ocp_qp_ipm_ws_strsize()
		ipm_ws_struct = cast(create_string_buffer(sizeof_d_ocp_qp_ipm_workspace), c_void_p)
		self.ipm_ws_struct = ipm_ws_struct

		# allocate memory for ipm workspace 
		ipm_size = __hpipm.d_ocp_qp_ipm_ws_memsize(qp_dims.dim_struct, arg.arg_struct)
		ipm_ws_mem = cast(create_string_buffer(ipm_size), c_void_p)
		self.ipm_ws_mem = ipm_ws_mem

		# create C ws
		__hpipm.d_ocp_qp_ipm_ws_create(qp_dims.dim_struct, arg.arg_struct, ipm_ws_struct, ipm_ws_mem)

		self.arg = arg
		self.dim_struct = qp_dims.dim_struct
		

	def solve(self, qp, qp_sol):
		hpipm_return = self.__hpipm.d_ocp_qp_ipm_solve(qp.qp_struct, qp_sol.qp_sol_struct, self.arg.arg_struct, self.ipm_ws_struct)
		return hpipm_return

	# TODO single getter !!!!!!!!!!!!!!!!!!!!!!!!!!

	def get_max_res_stat(self):
		res = np.zeros((1,1));
		tmp = cast(res.ctypes.data, POINTER(c_double))
		self.__hpipm.d_ocp_qp_ipm_get_max_res_stat(self.ipm_ws_struct, tmp)
		return res[0][0]


	def get_max_res_eq(self):
		res = np.zeros((1,1));
		tmp = cast(res.ctypes.data, POINTER(c_double))
		self.__hpipm.d_ocp_qp_ipm_get_max_res_eq(self.ipm_ws_struct, tmp)
		return res[0][0]


	def get_max_res_ineq(self):
		res = np.zeros((1,1));
		tmp = cast(res.ctypes.data, POINTER(c_double))
		self.__hpipm.d_ocp_qp_ipm_get_max_res_ineq(self.ipm_ws_struct, tmp)
		return res[0][0]


	def get_max_res_comp(self):
		res = np.zeros((1,1));
		tmp = cast(res.ctypes.data, POINTER(c_double))
		self.__hpipm.d_ocp_qp_ipm_get_max_res_comp(self.ipm_ws_struct, tmp)
		return res[0][0]


	def get_iter(self):
		res = np.zeros((1,1), dtype=int);
		tmp = cast(res.ctypes.data, POINTER(c_int))
		self.__hpipm.d_ocp_qp_ipm_get_iter(self.ipm_ws_struct, tmp)
		return res[0][0]


	def get_stat(self):
		# get iters
		iters = np.zeros((1,1), dtype=int);
		tmp = cast(iters.ctypes.data, POINTER(c_int))
		self.__hpipm.d_ocp_qp_ipm_get_iter(self.ipm_ws_struct, tmp)
		# get stat_m
		stat_m = np.zeros((1,1), dtype=int);
		tmp = cast(stat_m.ctypes.data, POINTER(c_int))
		self.__hpipm.d_ocp_qp_ipm_get_stat_m(self.ipm_ws_struct, tmp)
		# get stat pointer
		res = np.zeros((iters[0][0]+1, stat_m[0][0]));
		ptr = c_void_p()
		self.__hpipm.d_ocp_qp_ipm_get_stat(self.ipm_ws_struct, byref(ptr))
		tmp = cast(ptr, POINTER(c_double))
		for ii in range(iters[0][0]+1):
			for jj in range(stat_m[0][0]):
				res[ii][jj] = tmp[jj+ii*stat_m[0][0]]
				res[ii][jj] = tmp[jj+ii*stat_m[0][0]]
		return res

