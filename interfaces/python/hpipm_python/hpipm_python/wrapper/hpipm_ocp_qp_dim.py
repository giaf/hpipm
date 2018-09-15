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


	def set_nx(self, nx, idx=None):
		self.__hpipm.d_set_ocp_qp_dim_nx.argtypes = [c_int, c_int, c_void_p]
		if idx==None:
			for i in range(nx.size):
				self.nx[i] = nx[i]
				self.__hpipm.d_set_ocp_qp_dim_nx(i, nx[i], self.dim_struct)
		else:
			self.nx[idx] = nx
			self.__hpipm.d_set_ocp_qp_dim_nx(idx, nx, self.dim_struct)
		return


	def set_nu(self, nu, idx=None):
		self.__hpipm.d_set_ocp_qp_dim_nu.argtypes = [c_int, c_int, c_void_p]
		if idx==None:
			for i in range(nu.size):
				self.nu[i] = nu[i]
				self.__hpipm.d_set_ocp_qp_dim_nu(i, nu[i], self.dim_struct)
		else:
			self.nu[idx] = nu
			self.__hpipm.d_set_ocp_qp_dim_nu(idx, nu, self.dim_struct)
		return


	def set_nbx(self, nbx, idx=None):
		self.__hpipm.d_set_ocp_qp_dim_nbx.argtypes = [c_int, c_int, c_void_p]
		if idx==None:
			for i in range(nbx.size):
				self.nbx[i] = nbx[i]
				self.__hpipm.d_set_ocp_qp_dim_nbx(i, nbx[i], self.dim_struct)
		else:
			self.nbx[idx] = nbx
			self.__hpipm.d_set_ocp_qp_dim_nbx(idx, nbx, self.dim_struct)
		return


	def set_nbu(self, nbu, idx=None):
		self.__hpipm.d_set_ocp_qp_dim_nbu.argtypes = [c_int, c_int, c_void_p]
		if idx==None:
			for i in range(nbu.size):
				self.nbu[i] = nbu[i]
				self.__hpipm.d_set_ocp_qp_dim_nbu(i, nbu[i], self.dim_struct)
		else:
			self.nbu[idx] = nbu
			self.__hpipm.d_set_ocp_qp_dim_nbu(idx, nbu, self.dim_struct)
		return


	def set_ng(self, ng, idx=None):
		self.__hpipm.d_set_ocp_qp_dim_ng.argtypes = [c_int, c_int, c_void_p]
		if idx==None:
			for i in range(ng.size):
				self.ng[i] = ng[i]
				self.__hpipm.d_set_ocp_qp_dim_ng(i, ng[i], self.dim_struct)
		else:
			self.ng[idx] = ng
			self.__hpipm.d_set_ocp_qp_dim_ng(idx, ng, self.dim_struct)
		return


	def set_ns(self, ns, idx=None):
		self.__hpipm.d_set_ocp_qp_dim_ns.argtypes = [c_int, c_int, c_void_p]
		if idx==None:
			for i in range(ns.size):
				self.ns[i] = ns[i]
				self.__hpipm.d_set_ocp_qp_dim_ns(i, ns[i], self.dim_struct)
		else:
			self.ns[idx] = ns
			self.__hpipm.d_set_ocp_qp_dim_ns(idx, ns, self.dim_struct)
		return




