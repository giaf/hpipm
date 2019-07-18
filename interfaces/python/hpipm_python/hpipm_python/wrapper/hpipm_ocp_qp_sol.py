from ctypes import *
import ctypes.util 
import numpy as np
#import faulthandler

#faulthandler.enable()



class hpipm_ocp_qp_sol:
	def __init__(self, dim):

		# save dim internally
		self.dim = dim

		# extract dim
		N = dim.N
		nx = dim.nx
		nu = dim.nu
		nbx = dim.nbx
		nbu = dim.nbu
		ng = dim.ng
		ns = dim.ns

		# local qp_sol memory
		self.u = []
		for i in range(N+1):
			self.u.append(np.zeros((nu[i], 1)))

		self.x = []
		for i in range(N+1):
			self.x.append(np.zeros((nx[i], 1)))

		# TODO pi lam t

		# load hpipm shared library
		__hpipm   = CDLL('libhpipm.so')
		self.__hpipm = __hpipm

		# C qp struct
		qp_sol_struct_size = __hpipm.d_ocp_qp_sol_strsize()
		qp_sol_struct = cast(create_string_buffer(qp_sol_struct_size), c_void_p)
		self.qp_sol_struct = qp_sol_struct

		# C qp internal memory
		qp_sol_mem_size = __hpipm.d_ocp_qp_sol_memsize(dim.dim_struct)
		qp_sol_mem = cast(create_string_buffer(qp_sol_mem_size), c_void_p)
		self.qp_sol_mem = qp_sol_mem

		# create C qp
		__hpipm.d_ocp_qp_sol_create(dim.dim_struct, qp_sol_struct, qp_sol_mem)



	# TODO single getter !!!!!!!!!!!!!!!!!!!!!!!!!!!!

	def get_u(self, idx=None):
		# extract dims
		N = self.dim.N
		nu = self.dim.nu
		if idx==None:
			u = []
			for i in range(N+1):
				u.append(np.zeros((nu[i], 1)))
				tmp = cast(u[i].ctypes.data, POINTER(c_double))
				self.__hpipm.d_ocp_qp_sol_get_u(i, self.qp_sol_struct, tmp)
		else:
			u = np.zeros((nu[idx], 1))
			tmp = cast(u.ctypes.data, POINTER(c_double))
			self.__hpipm.d_ocp_qp_sol_get_u(idx, self.qp_sol_struct, tmp)
		return u



	def get_x(self, idx=None):
		# extract dims
		N = self.dim.N
		nx = self.dim.nx
		if idx==None:
			x = []
			for i in range(N+1):
				x.append(np.zeros((nx[i], 1)))
				tmp = cast(x[i].ctypes.data, POINTER(c_double))
				self.__hpipm.d_ocp_qp_sol_get_x(i, self.qp_sol_struct, tmp)
		else:
			x = np.zeros((nx[idx], 1))
			tmp = cast(x.ctypes.data, POINTER(c_double))
			self.__hpipm.d_ocp_qp_sol_get_x(idx, self.qp_sol_struct, tmp)
		return x



	def print_C_struct(self):
		self.__hpipm.d_ocp_qp_sol_print(self.dim.dim_struct, self.qp_sol_struct)
		return 



