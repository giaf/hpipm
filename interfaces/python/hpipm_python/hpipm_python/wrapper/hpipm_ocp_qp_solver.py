from ctypes import *
import ctypes.util 
import numpy as np
#import faulthandler

#faulthandler.enable()



class hpipm_ocp_qp_solver:
	def __init__(self, qp_dims):

		# load hpipm shared library
		__hpipm   = CDLL('libhpipm.so')
		self.__hpipm = __hpipm

		# allocate memory for ipm_arg struct
		ipm_arg_size = __hpipm.d_memsize_ocp_qp_ipm_arg(qp_dims.dim_struct)
		ipm_arg_mem = cast(create_string_buffer(ipm_arg_size), c_void_p)
		self.ipm_arg_mem = ipm_arg_mem

		# set up ipm_arg
		sizeof_d_ocp_qp_ipm_arg = __hpipm.d_sizeof_ocp_qp_ipm_arg()
		arg = cast(create_string_buffer(sizeof_d_ocp_qp_ipm_arg), c_void_p)
		self.ocp_qp_ipm_arg = arg

		__hpipm.d_create_ocp_qp_ipm_arg(qp_dims.dim_struct, arg, ipm_arg_mem)
		__hpipm.d_set_default_ocp_qp_ipm_arg(1, arg)

		# allocate memory for ipm workspace 
		ipm_size = __hpipm.d_memsize_ocp_qp_ipm(qp_dims.dim_struct, arg)
		ipm_mem = cast(create_string_buffer(ipm_size), c_void_p)
		self.ipm_mem = ipm_mem

		# set up ipm workspace
		sizeof_d_ocp_qp_ipm_workspace = __hpipm.d_sizeof_ocp_qp_ipm_workspace()
		workspace = cast(create_string_buffer(sizeof_d_ocp_qp_ipm_workspace), c_void_p)
		self.ocp_qp_ipm_workspace = workspace

		__hpipm.d_create_ocp_qp_ipm(qp_dims.dim_struct, arg, workspace, ipm_mem)

		self.arg = arg
		self.workspace = workspace
		self.dim_struct = qp_dims.dim_struct
		

	def solve(self, qp, qp_sol):
		hpipm_return = self.__hpipm.d_solve_ocp_qp_ipm(qp.qp_struct, qp_sol.qp_sol_struct, self.arg, self.workspace)
		return hpipm_return



