from ctypes import *
import ctypes.util 
import numpy as np
#import faulthandler

#faulthandler.enable()

class hpipm_ocp_qp_dim:

	def __init__(self, N):

		# memory for python class
		self.N	= N
		self.nx   = np.zeros(N+1, dtype=np.int32)
		self.nu   = np.zeros(N+1, dtype=np.int32)
		self.nbx  = np.zeros(N+1, dtype=np.int32)
		self.nbu  = np.zeros(N+1, dtype=np.int32)
		self.ng   = np.zeros(N+1, dtype=np.int32)
		self.ns   = np.zeros(N+1, dtype=np.int32)

		# load hpipm shared library
		__hpipm   = CDLL('libhpipm.so')
		self.__hpipm = __hpipm

		# C dim struct
		dim_struct_size = __hpipm.d_sizeof_ocp_qp_dim()
		dim_struct = cast(create_string_buffer(dim_struct_size), c_void_p)
		self.dim_struct = dim_struct

		# C dim internal memory
		dim_mem_size = __hpipm.d_memsize_ocp_qp_dim(N)
		dim_mem = cast(create_string_buffer(dim_mem_size), c_void_p)
		self.dim_mem = dim_mem

		# create C dim
		__hpipm.d_create_ocp_qp_dim(N, self.dim_struct, self.dim_mem)

#		__hpipm.d_cvt_int_to_ocp_qp_dim(N,
#			cast(nx.ctypes.data, POINTER(c_double)),
#			cast(nu.ctypes.data, POINTER(c_double)),
#			cast(nbx.ctypes.data, POINTER(c_double)),
#			cast(nbu.ctypes.data, POINTER(c_double)),
#			cast(ng.ctypes.data, POINTER(c_double)),
#			cast(ns.ctypes.data, POINTER(c_double)),
#			self.dim_struct)


	def set(self, field, value, idx=None):
		self.__hpipm.d_set_ocp_qp_dim.argtypes = [c_char_p, c_int, c_int, c_void_p]
		if idx==None:
			for i in range(value.size):
				field_ = getattr(self, field)
				field_[i] = value[i]
				field_name_b = field.encode('utf-8')
				self.__hpipm.d_set_ocp_qp_dim(c_char_p(field_name_b), i, value[i], self.dim_struct)
		else:
			field_ = getattr(self, field)
			field_[idx] = value
			field_name_b = field.encode('utf-8')
			self.__hpipm.d_set_ocp_qp_dim(c_char_p(field_name_b), idx, value, self.dim_struct)
		return

	def print_C_struct(self):
		self.__hpipm.d_print_ocp_qp_dim(self.dim_struct)
		return 

