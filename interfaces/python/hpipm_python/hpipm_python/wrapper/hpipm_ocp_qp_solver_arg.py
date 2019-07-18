from ctypes import *
import ctypes.util 
import numpy as np
#import faulthandler

#faulthandler.enable()



class hpipm_ocp_qp_solver_arg:
	def __init__(self, dim):

		# save dim internally
		self.dim = dim

		# load hpipm shared library
		__hpipm   = CDLL('libhpipm.so')
		self.__hpipm = __hpipm

		# C qp struct
		arg_struct_size = __hpipm.d_ocp_qp_ipm_arg_strsize()
		arg_struct = cast(create_string_buffer(arg_struct_size), c_void_p)
		self.arg_struct = arg_struct

		# C qp internal memory
		arg_mem_size = __hpipm.d_ocp_qp_ipm_arg_memsize(dim.dim_struct)
		arg_mem = cast(create_string_buffer(arg_mem_size), c_void_p)
		self.arg_mem = arg_mem

		# create C qp
		__hpipm.d_ocp_qp_ipm_arg_create(dim.dim_struct, arg_struct, arg_mem)

		# initialize default arguments
		__hpipm.d_ocp_qp_ipm_arg_set_default(1, arg_struct) # mode==SPEED
	

	# TODO single setter !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	def set_mu0(self, mu0):
		tmp_in = np.zeros((1,1))
		tmp_in[0][0] = mu0
		tmp = cast(tmp_in.ctypes.data, POINTER(c_double))
		self.__hpipm.d_ocp_qp_ipm_arg_set_mu0.argtypes = [POINTER(c_double), c_void_p]
		self.__hpipm.d_ocp_qp_ipm_arg_set_mu0(tmp, self.arg_struct)
		return


	def set_iter_max(self, iter_max):
		tmp_in = np.zeros((1,1), dtype=int)
		tmp_in[0][0] = iter_max
		tmp = cast(tmp_in.ctypes.data, POINTER(c_int))
		self.__hpipm.d_ocp_qp_ipm_arg_set_iter_max.argtypes = [POINTER(c_int), c_void_p]
		self.__hpipm.d_ocp_qp_ipm_arg_set_iter_max(tmp, self.arg_struct)
		return


	def set_tol_stat(self, tol_stat):
		tmp_in = np.zeros((1,1))
		tmp_in[0][0] = tol_stat
		tmp = cast(tmp_in.ctypes.data, POINTER(c_double))
		self.__hpipm.d_ocp_qp_ipm_arg_set_tol_stat.argtypes = [POINTER(c_double), c_void_p]
		self.__hpipm.d_ocp_qp_ipm_arg_set_tol_stat(tmp, self.arg_struct)
		return


	def set_tol_eq(self, tol_eq):
		tmp_in = np.zeros((1,1))
		tmp_in[0][0] = tol_eq
		tmp = cast(tmp_in.ctypes.data, POINTER(c_double))
		self.__hpipm.d_ocp_qp_ipm_arg_set_tol_eq.argtypes = [POINTER(c_double), c_void_p]
		self.__hpipm.d_ocp_qp_ipm_arg_set_tol_eq(tmp, self.arg_struct)
		return


	def set_tol_ineq(self, tol_ineq):
		tmp_in = np.zeros((1,1))
		tmp_in[0][0] = tol_ineq
		tmp = cast(tmp_in.ctypes.data, POINTER(c_double))
		self.__hpipm.d_ocp_qp_ipm_arg_set_tol_ineq.argtypes = [POINTER(c_double), c_void_p]
		self.__hpipm.d_ocp_qp_ipm_arg_set_tol_ineq(tmp, self.arg_struct)
		return


	def set_tol_comp(self, tol_comp):
		tmp_in = np.zeros((1,1))
		tmp_in[0][0] = tol_comp
		tmp = cast(tmp_in.ctypes.data, POINTER(c_double))
		self.__hpipm.d_ocp_qp_ipm_arg_set_tol_comp.argtypes = [POINTER(c_double), c_void_p]
		self.__hpipm.d_ocp_qp_ipm_arg_set_tol_comp(tmp, self.arg_struct)
		return

	def set_reg_prim(self, reg_prim):
		tmp_in = np.zeros((1,1))
		tmp_in[0][0] = reg_prim
		tmp = cast(tmp_in.ctypes.data, POINTER(c_double))
		self.__hpipm.d_ocp_qp_ipm_arg_set_reg_prim.argtypes = [POINTER(c_double), c_void_p]
		self.__hpipm.d_ocp_qp_ipm_arg_set_reg_prim(tmp, self.arg_struct)
		return



