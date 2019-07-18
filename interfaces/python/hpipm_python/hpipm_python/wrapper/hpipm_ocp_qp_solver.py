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

	def get_res_stat(self):
		res = np.zeros((1,1));
		tmp = cast(res.ctypes.data, POINTER(c_double))
		self.__hpipm.d_ocp_qp_ipm_get_res_stat(self.ipm_ws_struct, tmp)
		return res[0][0]


	def get_res_eq(self):
		res = np.zeros((1,1));
		tmp = cast(res.ctypes.data, POINTER(c_double))
		self.__hpipm.d_ocp_qp_ipm_get_res_eq(self.ipm_ws_struct, tmp)
		return res[0][0]


	def get_res_ineq(self):
		res = np.zeros((1,1));
		tmp = cast(res.ctypes.data, POINTER(c_double))
		self.__hpipm.d_ocp_qp_ipm_get_res_ineq(self.ipm_ws_struct, tmp)
		return res[0][0]


	def get_res_comp(self):
		res = np.zeros((1,1));
		tmp = cast(res.ctypes.data, POINTER(c_double))
		self.__hpipm.d_ocp_qp_ipm_get_res_comp(self.ipm_ws_struct, tmp)
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





#	def get_res_eq(self):
#		self.__hpipm.d_get_ocp_qp_ipm_res_eq.restype = c_double
#		res = self.__hpipm.d_get_ocp_qp_ipm_res_eq(self.ipm_ws_struct)
#		return res


#	def get_res_ineq(self):
#		self.__hpipm.d_get_ocp_qp_ipm_res_ineq.restype = c_double
#		res = self.__hpipm.d_get_ocp_qp_ipm_res_ineq(self.ipm_ws_struct)
#		return res


#	def get_res_comp(self):
#		self.__hpipm.d_get_ocp_qp_ipm_res_comp.restype = c_double
#		res = self.__hpipm.d_get_ocp_qp_ipm_res_comp(self.ipm_ws_struct)
#		return res


#	def get_iter(self):
#		iters = self.__hpipm.d_get_ocp_qp_ipm_iter(self.ipm_ws_struct)
#		return iters
	
#	def get_stat(self):
#		iters = self.__hpipm.d_get_ocp_qp_ipm_iter(self.ipm_ws_struct)
#		self.__hpipm.d_get_ocp_qp_ipm_stat.restype = POINTER(c_double)
#		stat = self.__hpipm.d_get_ocp_qp_ipm_stat(self.ipm_ws_struct)
#		stat_py = np.zeros((iters, 5))
#		for ii in range(iters):
#			for jj in range(5):
#				stat_py[ii][jj] = stat[jj+ii*5]
#				stat_py[ii][jj] = stat[jj+ii*5]
#		return stat_py

