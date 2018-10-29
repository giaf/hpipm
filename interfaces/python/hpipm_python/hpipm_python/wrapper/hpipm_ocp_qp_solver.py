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

		# allocate memory for ipm workspace 
		ipm_size = __hpipm.d_memsize_ocp_qp_ipm(qp_dims.dim_struct, arg.arg_struct)
		ipm_mem = cast(create_string_buffer(ipm_size), c_void_p)
		self.ipm_mem = ipm_mem

		# set up ipm workspace
		sizeof_d_ocp_qp_ipm_workspace = __hpipm.d_sizeof_ocp_qp_ipm_workspace()
		workspace = cast(create_string_buffer(sizeof_d_ocp_qp_ipm_workspace), c_void_p)
		self.ocp_qp_ipm_workspace = workspace

		__hpipm.d_create_ocp_qp_ipm(qp_dims.dim_struct, arg.arg_struct, workspace, ipm_mem)

		self.arg = arg
		self.workspace = workspace
		self.dim_struct = qp_dims.dim_struct
		

	def solve(self, qp, qp_sol):
		hpipm_return = self.__hpipm.d_solve_ocp_qp_ipm(qp.qp_struct, qp_sol.qp_sol_struct, self.arg.arg_struct, self.workspace)
		return hpipm_return


	def get_res_stat(self):
		self.__hpipm.d_get_ocp_qp_ipm_res_stat.restype = c_double
		res = self.__hpipm.d_get_ocp_qp_ipm_res_stat(self.workspace)
		return res


	def get_res_eq(self):
		self.__hpipm.d_get_ocp_qp_ipm_res_eq.restype = c_double
		res = self.__hpipm.d_get_ocp_qp_ipm_res_eq(self.workspace)
		return res


	def get_res_ineq(self):
		self.__hpipm.d_get_ocp_qp_ipm_res_ineq.restype = c_double
		res = self.__hpipm.d_get_ocp_qp_ipm_res_ineq(self.workspace)
		return res


	def get_res_comp(self):
		self.__hpipm.d_get_ocp_qp_ipm_res_comp.restype = c_double
		res = self.__hpipm.d_get_ocp_qp_ipm_res_comp(self.workspace)
		return res


	def get_iter(self):
		iters = self.__hpipm.d_get_ocp_qp_ipm_iter(self.workspace)
		return iters
	
	def get_stat(self):
		iters = self.__hpipm.d_get_ocp_qp_ipm_iter(self.workspace)
		self.__hpipm.d_get_ocp_qp_ipm_stat.restype = POINTER(c_double)
		stat = self.__hpipm.d_get_ocp_qp_ipm_stat(self.workspace)
		stat_py = np.zeros((iters, 5))
		for ii in range(iters):
			for jj in range(5):
				stat_py[ii][jj] = stat[jj+ii*5]
				stat_py[ii][jj] = stat[jj+ii*5]
		return stat_py

