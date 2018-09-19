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
		arg_struct_size = __hpipm.d_sizeof_ocp_qp_ipm_arg()
		arg_struct = cast(create_string_buffer(arg_struct_size), c_void_p)
		self.arg_struct = arg_struct

		# C qp internal memory
		arg_mem_size = __hpipm.d_memsize_ocp_qp_ipm_arg(dim.dim_struct)
		arg_mem = cast(create_string_buffer(arg_mem_size), c_void_p)
		self.arg_mem = arg_mem

		# create C qp
		__hpipm.d_create_ocp_qp_ipm_arg(dim.dim_struct, arg_struct, arg_mem)

		# initialize default arguments
		__hpipm.d_set_default_ocp_qp_ipm_arg(1, arg_struct) # mode==SPEED
	

	def set_mu0(self, mu0):
		self.__hpipm.d_set_ocp_qp_ipm_arg_mu0.argtypes = [c_double, c_void_p]
		self.__hpipm.d_set_ocp_qp_ipm_arg_mu0(mu0, self.arg_struct)
		return


	def set_iter_max(self, iter_max):
		self.__hpipm.d_set_ocp_qp_ipm_arg_iter_max.argtypes = [c_int, c_void_p]
		self.__hpipm.d_set_ocp_qp_ipm_arg_iter_max(iter_max, self.arg_struct)
		return


	def set_tol_stat(self, tol_stat):
		self.__hpipm.d_set_ocp_qp_ipm_arg_tol_stat.argtypes = [c_double, c_void_p]
		self.__hpipm.d_set_ocp_qp_ipm_arg_tol_stat(tol_stat, self.arg_struct)
		return


	def set_tol_eq(self, tol_eq):
		self.__hpipm.d_set_ocp_qp_ipm_arg_tol_eq.argtypes = [c_double, c_void_p]
		self.__hpipm.d_set_ocp_qp_ipm_arg_tol_eq(tol_eq, self.arg_struct)
		return


	def set_tol_ineq(self, tol_ineq):
		self.__hpipm.d_set_ocp_qp_ipm_arg_tol_ineq.argtypes = [c_double, c_void_p]
		self.__hpipm.d_set_ocp_qp_ipm_arg_tol_ineq(tol_ineq, self.arg_struct)
		return


	def set_tol_comp(self, tol_comp):
		self.__hpipm.d_set_ocp_qp_ipm_arg_tol_comp.argtypes = [c_double, c_void_p]
		self.__hpipm.d_set_ocp_qp_ipm_arg_tol_comp(tol_comp, self.arg_struct)
		return

	def set_reg_prim(self, reg_prim):
		self.__hpipm.d_set_ocp_qp_ipm_arg_reg_prim.argtypes = [c_double, c_void_p]
		self.__hpipm.d_set_ocp_qp_ipm_arg_reg_prim(reg_prim, self.arg_struct)
		return



