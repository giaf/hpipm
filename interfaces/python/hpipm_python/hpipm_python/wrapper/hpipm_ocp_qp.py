from ctypes import *
import ctypes.util 
import numpy as np
#import faulthandler

#faulthandler.enable()



class hpipm_ocp_qp:
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

		# local qp memory
		self.A = []
		for i in range(N):
			self.A.append(np.zeros((nx[i+1], nx[i])))

		self.B = []
		for i in range(N):
			self.B.append(np.zeros((nx[i+1], nu[i])))

		self.b = []
		for i in range(N):
			self.b.append(np.zeros((nx[i+1], 1)))

		self.Q = []
		for i in range(N+1):
			self.Q.append(np.zeros((nx[i], nx[i])))

		self.R = []
		for i in range(N+1):
			self.R.append(np.zeros((nu[i], nu[i])))

		self.S = []
		for i in range(N+1):
			self.S.append(np.zeros((nu[i], nx[i])))

		self.q = []
		for i in range(N+1):
			self.q.append(np.zeros((nx[i], 1)))

		self.r = []
		for i in range(N+1):
			self.r.append(np.zeros((nu[i], 1)))

		self.Jx = []
		for i in range(N+1):
			self.Jx.append(np.zeros((nbx[i], nx[i])))

		self.lx = []
		for i in range(N+1):
			self.lx.append(np.zeros((nbx[i], 1)))

		self.ux = []
		for i in range(N+1):
			self.ux.append(np.zeros((nbx[i], 1)))

		self.Ju = []
		for i in range(N+1):
			self.Ju.append(np.zeros((nbu[i], nu[i])))

		self.lu = []
		for i in range(N+1):
			self.lu.append(np.zeros((nbu[i], 1)))

		self.uu = []
		for i in range(N+1):
			self.uu.append(np.zeros((nbu[i], 1)))

		self.idxb = []
		for i in range(N+1):
			self.idxb.append(np.zeros((nbx[i]+nbu[i], 1), dtype=np.int32))

		self.lb = []
		for i in range(N+1):
			self.lb.append(np.zeros((nbu[i]+nbx[i], 1)))

		self.ub = []
		for i in range(N+1):
			self.ub.append(np.zeros((nbu[i]+nbx[i], 1)))

		self.C = []
		for i in range(N+1):
			self.C.append(np.zeros((ng[i], nx[i])))

		self.D = []
		for i in range(N+1):
			self.D.append(np.zeros((ng[i], nu[i])))

		self.lg = []
		for i in range(N+1):
			self.lg.append(np.zeros((ng[i], 1)))

		self.ug = []
		for i in range(N+1):
			self.ug.append(np.zeros((ng[i], 1)))

		self.Zl = []
		for i in range(N+1):
			self.Zl.append(np.zeros((ns[i], 1)))

		self.Zu = []
		for i in range(N+1):
			self.Zu.append(np.zeros((ns[i], 1)))

		self.zl = []
		for i in range(N+1):
			self.zl.append(np.zeros((ns[i], 1)))

		self.zu = []
		for i in range(N+1):
			self.zu.append(np.zeros((ns[i], 1)))

		self.lls = []
		for i in range(N+1):
			self.lls.append(np.zeros((ns[i], 1)))

		self.lus = []
		for i in range(N+1):
			self.lus.append(np.zeros((ns[i], 1)))

		self.Jsu = []
		for i in range(N+1):
			self.Jsu.append(np.zeros((nbu[i], ns[i])))

		self.Jsx = []
		for i in range(N+1):
			self.Jsx.append(np.zeros((nbx[i], ns[i])))

		self.Jsg = []
		for i in range(N+1):
			self.Jsg.append(np.zeros((ng[i], ns[i])))

		self.idxs = []
		for i in range(N+1):
			self.idxs.append(np.zeros((ns[i], 1), dtype=np.int32))

		# load hpipm shared library
		__hpipm   = CDLL('libhpipm.so')
		self.__hpipm = __hpipm

		# C qp struct
		qp_struct_size = __hpipm.d_sizeof_ocp_qp()
		qp_struct = cast(create_string_buffer(qp_struct_size), c_void_p)
		self.qp_struct = qp_struct

		# C qp internal memory
		qp_mem_size = __hpipm.d_memsize_ocp_qp(dim.dim_struct)
		qp_mem = cast(create_string_buffer(qp_mem_size), c_void_p)
		self.qp_mem = qp_mem

		# create C qp
		__hpipm.d_create_ocp_qp(dim.dim_struct, qp_struct, qp_mem)



	def set_A(self, A, idx=None):
		self.__hpipm.d_cvt_colmaj_to_ocp_qp_A.argtypes = [c_int, POINTER(c_double), c_void_p]
		nx = self.dim.nx
		if idx==None:
			for i in range(len(A)):
				self.A[i] = A[i]
				self.A[i] = self.A[i].reshape((nx[i+1], nx[i]))
				self.A[i] = self.A[i].transpose()
				self.A[i] = np.ascontiguousarray(self.A[i], dtype=np.float64)
				tmp = cast(self.A[i].ctypes.data, POINTER(c_double))
				self.__hpipm.d_cvt_colmaj_to_ocp_qp_A(i, tmp, self.qp_struct)
		else:
			self.A[idx] = A
			self.A[idx] = self.A[idx].reshape((nx[idx+1], nx[idx]))
			self.A[idx] = self.A[idx].transpose()
			self.A[idx] = np.ascontiguousarray(self.A[idx], dtype=np.float64)
			tmp = cast(self.A[idx].ctypes.data, POINTER(c_double))
			self.__hpipm.d_cvt_colmaj_to_ocp_qp_A(idx, tmp, self.qp_struct)
		return


	def set_B(self, B, idx=None):
		self.__hpipm.d_cvt_colmaj_to_ocp_qp_B.argtypes = [c_int, POINTER(c_double), c_void_p]
		nx = self.dim.nx
		nu = self.dim.nu
		if idx==None:
			for i in range(len(B)):
				self.B[i] = B[i]
				self.B[i] = self.B[i].reshape((nx[i+1], nu[i]))
				self.B[i] = self.B[i].transpose()
				self.B[i] = np.ascontiguousarray(self.B[i], dtype=np.float64)
				tmp = cast(self.B[i].ctypes.data, POINTER(c_double))
				self.__hpipm.d_cvt_colmaj_to_ocp_qp_B(i, tmp, self.qp_struct)
		else:
			self.B[idx] = B
			self.B[idx] = self.B[idx].reshape((nx[idx+1], nu[idx]))
			self.B[idx] = self.B[idx].transpose()
			self.B[idx] = np.ascontiguousarray(self.B[idx], dtype=np.float64)
			tmp = cast(self.B[idx].ctypes.data, POINTER(c_double))
			self.__hpipm.d_cvt_colmaj_to_ocp_qp_B(idx, tmp, self.qp_struct)
		return


	def set_b(self, b, idx=None):
		self.__hpipm.d_cvt_colmaj_to_ocp_qp_b.argtypes = [c_int, POINTER(c_double), c_void_p]
		nx = self.dim.nx
		if idx==None:
			for i in range(len(b)):
				self.b[i] = b[i]
				self.b[i] = self.b[i].reshape((nx[i+1], 1))
				self.b[i] = np.ascontiguousarray(self.b[i], dtype=np.float64)
				tmp = cast(self.b[i].ctypes.data, POINTER(c_double))
				self.__hpipm.d_cvt_colmaj_to_ocp_qp_b(i, tmp, self.qp_struct)
		else:
			self.b[idx] = b
			self.b[idx] = self.b[idx].reshape((nx[idx+1], 1))
			self.b[idx] = np.ascontiguousarray(self.b[idx], dtype=np.float64)
			tmp = cast(self.b[idx].ctypes.data, POINTER(c_double))
			self.__hpipm.d_cvt_colmaj_to_ocp_qp_b(idx, tmp, self.qp_struct)
		return


	def set_Q(self, Q, idx=None):
		self.__hpipm.d_cvt_colmaj_to_ocp_qp_Q.argtypes = [c_int, POINTER(c_double), c_void_p]
		nx = self.dim.nx
		if idx==None:
			for i in range(len(Q)):
				self.Q[i] = Q[i]
				self.Q[i] = self.Q[i].reshape((nx[i], nx[i]))
				self.Q[i] = self.Q[i].transpose()
				self.Q[i] = np.ascontiguousarray(self.Q[i], dtype=np.float64)
				tmp = cast(self.Q[i].ctypes.data, POINTER(c_double))
				self.__hpipm.d_cvt_colmaj_to_ocp_qp_Q(i, tmp, self.qp_struct)
		else:
			self.Q[idx] = Q
			self.Q[idx] = self.Q[idx].reshape((nx[idx], nx[idx]))
			self.Q[idx] = self.Q[idx].transpose()
			self.Q[idx] = np.ascontiguousarray(self.Q[idx], dtype=np.float64)
			tmp = cast(self.Q[idx].ctypes.data, POINTER(c_double))
			self.__hpipm.d_cvt_colmaj_to_ocp_qp_Q(idx, tmp, self.qp_struct)
		return


	def set_R(self, R, idx=None):
		self.__hpipm.d_cvt_colmaj_to_ocp_qp_R.argtypes = [c_int, POINTER(c_double), c_void_p]
		nu = self.dim.nu
		if idx==None:
			for i in range(len(R)):
				self.R[i] = R[i]
				self.R[i] = self.R[i].reshape((nu[i], nu[i]))
				self.R[i] = self.R[i].transpose()
				self.R[i] = np.ascontiguousarray(self.R[i], dtype=np.float64)
				tmp = cast(self.R[i].ctypes.data, POINTER(c_double))
				self.__hpipm.d_cvt_colmaj_to_ocp_qp_R(i, tmp, self.qp_struct)
		else:
			self.R[idx] = R
			self.R[idx] = self.R[idx].reshape((nu[idx], nu[idx]))
			self.R[idx] = self.R[idx].transpose()
			self.R[idx] = np.ascontiguousarray(self.R[idx], dtype=np.float64)
			tmp = cast(self.R[idx].ctypes.data, POINTER(c_double))
			self.__hpipm.d_cvt_colmaj_to_ocp_qp_R(idx, tmp, self.qp_struct)
		return


	def set_S(self, S, idx=None):
		self.__hpipm.d_cvt_colmaj_to_ocp_qp_S.argtypes = [c_int, POINTER(c_double), c_void_p]
		nu = self.dim.nu
		nx = self.dim.nx
		if idx==None:
			for i in range(len(S)):
				self.S[i] = S[i]
				self.S[i] = self.S[i].reshape((nx[i], nu[i]))
				self.S[i] = self.S[i].transpose()
				self.S[i] = np.ascontiguousarray(self.S[i], dtype=np.float64)
				tmp = cast(self.S[i].ctypes.data, POINTER(c_double))
				self.__hpipm.d_cvt_colmaj_to_ocp_qp_S(i, tmp, self.qp_struct)
		else:
			self.S[idx] = S
			self.S[idx] = self.S[idx].reshape((nu[idx], nx[idx]))
			self.S[idx] = self.S[idx].transpose()
			self.S[idx] = np.ascontiguousarray(self.S[idx], dtype=np.float64)
			tmp = cast(self.S[idx].ctypes.data, POINTER(c_double))
			self.__hpipm.d_cvt_colmaj_to_ocp_qp_S(idx, tmp, self.qp_struct)
		return


	def set_q(self, q, idx=None):
		self.__hpipm.d_cvt_colmaj_to_ocp_qp_q.argtypes = [c_int, POINTER(c_double), c_void_p]
		nx = self.dim.nx
		if idx==None:
			for i in range(len(q)):
				self.q[i] = q[i]
				self.q[i] = self.q[i].reshape((nx[i], 1))
				self.q[i] = np.ascontiguousarray(self.q[i], dtype=np.float64)
				tmp = cast(self.q[i].ctypes.data, POINTER(c_double))
				self.__hpipm.d_cvt_colmaj_to_ocp_qp_q(i, tmp, self.qp_struct)
		else:
			self.q[idx] = q
			self.q[idx] = self.q[idx].reshape((nx[idx], 1))
			self.q[idx] = np.ascontiguousarray(self.q[idx], dtype=np.float64)
			tmp = cast(self.q[idx].ctypes.data, POINTER(c_double))
			self.__hpipm.d_cvt_colmaj_to_ocp_qp_q(idx, tmp, self.qp_struct)
		return


	def set_r(self, r, idx=None):
		self.__hpipm.d_cvt_colmaj_to_ocp_qp_r.argtypes = [c_int, POINTER(c_double), c_void_p]
		nu = self.dim.nu
		if idx==None:
			for i in range(len(r)):
				self.r[i] = r[i]
				self.r[i] = self.r[i].reshape((nx[i], 1))
				self.r[i] = np.ascontiguousarray(self.r[i], dtype=np.float64)
				tmp = cast(self.r[i].ctypes.data, POINTER(c_double))
				self.__hpipm.d_cvt_colmaj_to_ocp_qp_r(i, tmp, self.qp_struct)
		else:
			self.r[idx] = r
			self.r[idx] = self.r[idx].reshape((nx[idx], 1))
			self.r[idx] = np.ascontiguousarray(self.r[idx], dtype=np.float64)
			tmp = cast(self.r[idx].ctypes.data, POINTER(c_double))
			self.__hpipm.d_cvt_colmaj_to_ocp_qp_r(idx, tmp, self.qp_struct)
		return


	def set_Jx(self, Jx, idx=None):
		self.__hpipm.d_cvt_colmaj_to_ocp_qp_idxb.argtypes = [c_int, POINTER(c_int), c_void_p]
		nx = self.dim.nx
		nu = self.dim.nu
		nbx = self.dim.nbx
		nbu = self.dim.nbu
		if idx==None:
			for i in range(len(Jx)):
				self.Jx[i] = Jx[i]
				self.Jx[i] = self.Jx[i].reshape((nbx[i], nx[i]))
				self.Jx[i] = self.Jx[i].transpose()
				self.Jx[i] = np.ascontiguousarray(self.Jx[i], dtype=np.float64)
				for j in range(nbx[i]):
					k0 = -1
					for k in range(nx[i]):
						if (k0==-1) & (self.Jx[i][j][k]!=0):
							k0 = k
							self.idxb[i][nbu[i]+j] = nu[i]+k
				self.idxb[i] = np.ascontiguousarray(self.idxb[i], dtype=np.int32)
				tmp = cast(self.idxb[i].ctypes.data, POINTER(c_int))
				self.__hpipm.d_cvt_colmaj_to_ocp_qp_idxb(i, tmp, self.qp_struct)
		else:
			self.Jx[idx] = Jx
			self.Jx[idx] = self.Jx[idx].reshape((nbx[idx], nx[idx]))
			self.Jx[idx] = self.Jx[idx].transpose()
			self.Jx[idx] = np.ascontiguousarray(self.Jx[idx], dtype=np.float64)
			for j in range(nbx[idx]):
				k0 = -1
				for k in range(nx[idx]):
					if (k0==-1) & (self.Jx[idx][j][k]!=0):
						k0 = k
						self.idxb[idx][nbu[idx]+j] = nu[idx]+k
			self.idxb[idx] = np.ascontiguousarray(self.idxb[idx], dtype=np.int32)
			tmp = cast(self.idxb[idx].ctypes.data, POINTER(c_int))
			self.__hpipm.d_cvt_colmaj_to_ocp_qp_idxb(idx, tmp, self.qp_struct)
		return


	def set_lx(self, lx, idx=None):
		self.__hpipm.d_cvt_colmaj_to_ocp_qp_lb.argtypes = [c_int, POINTER(c_double), c_void_p]
		nbx = self.dim.nbx
		nbu = self.dim.nbu
		if idx==None:
			for i in range(len(lx)):
				self.lx[i] = lx[i]
				self.lx[i] = self.lx[i].reshape((nbx[i], 1))
				self.lx[i] = np.ascontiguousarray(self.lx[i], dtype=np.float64)
				for j in range(nbx[i]):
					self.lb[i][nbu[i]+j] = lx[i][j]
				self.lb[i] = np.ascontiguousarray(self.lb[i], dtype=np.float64)
				tmp = cast(self.lb[i].ctypes.data, POINTER(c_double))
				self.__hpipm.d_cvt_colmaj_to_ocp_qp_lb(i, tmp, self.qp_struct)
		else:
			self.lx[idx] = lx
			self.lx[idx] = self.lx[idx].reshape((nbx[idx], 1))
			self.lx[idx] = np.ascontiguousarray(self.lx[idx], dtype=np.float64)
			for j in range(nbx[idx]):
				self.lb[idx][nbu[idx]+j] = lx[j]
			self.lb[idx] = np.ascontiguousarray(self.lb[idx], dtype=np.float64)
			tmp = cast(self.lb[idx].ctypes.data, POINTER(c_double))
			self.__hpipm.d_cvt_colmaj_to_ocp_qp_lb(idx, tmp, self.qp_struct)
		return


	def set_ux(self, ux, idx=None):
		self.__hpipm.d_cvt_colmaj_to_ocp_qp_ub.argtypes = [c_int, POINTER(c_double), c_void_p]
		nbx = self.dim.nbx
		nbu = self.dim.nbu
		if idx==None:
			for i in range(len(ux)):
				self.ux[i] = ux[i]
				self.ux[i] = self.ux[i].reshape((nbx[i], 1))
				self.ux[i] = np.ascontiguousarray(self.ux[i], dtype=np.float64)
				for j in range(nbx[i]):
					self.ub[i][nbu[i]+j] = ux[i][j]
				self.ub[i] = np.ascontiguousarray(self.ub[i], dtype=np.float64)
				tmp = cast(self.ub[i].ctypes.data, POINTER(c_double))
				self.__hpipm.d_cvt_colmaj_to_ocp_qp_ub(i, tmp, self.qp_struct)
		else:
			self.ux[idx] = ux
			self.ux[idx] = self.ux[idx].reshape((nbx[idx], 1))
			self.ux[idx] = np.ascontiguousarray(self.ux[idx], dtype=np.float64)
			for j in range(nbx[idx]):
				self.ub[idx][nbu[idx]+j] = ux[j]
			self.ub[idx] = np.ascontiguousarray(self.ub[idx], dtype=np.float64)
			tmp = cast(self.ub[idx].ctypes.data, POINTER(c_double))
			self.__hpipm.d_cvt_colmaj_to_ocp_qp_ub(idx, tmp, self.qp_struct)
		return


	def set_Ju(self, Ju, idx=None):
		self.__hpipm.d_cvt_colmaj_to_ocp_qp_idxb.argtypes = [c_int, POINTER(c_int), c_void_p]
		nu = self.dim.nu
		nbu = self.dim.nbu
		if idx==None:
			for i in range(len(Ju)):
				self.Ju[i] = Ju[i]
				self.Ju[i] = self.Ju[i].reshape((nbu[i], nu[i]))
				self.Ju[i] = self.Ju[i].transpose()
				self.Ju[i] = np.ascontiguousarray(self.Ju[i], dtype=np.float64)
				for j in range(nbu[i]):
					k0 = -1
					for k in range(nu[i]):
						if (k0==-1) & (self.Ju[i][j][k]!=0):
							k0 = k
							self.idxb[i][j] = k
				self.idxb[i] = np.ascontiguousarray(self.idxb[i], dtype=np.int32)
				tmp = cast(self.idxb[i].ctypes.data, POINTER(c_int))
				self.__hpipm.d_cvt_colmaj_to_ocp_qp_idxb(i, tmp, self.qp_struct)
		else:
			self.Ju[idx] = Ju
			self.Ju[idx] = self.Ju[idx].reshape((nbu[idx], nu[idx]))
			self.Ju[idx] = self.Ju[idx].transpose()
			self.Ju[idx] = np.ascontiguousarray(self.Ju[idx], dtype=np.float64)
			for j in range(nbu[idx]):
				k0 = -1
				for k in range(nu[idx]):
					if (k0==-1) & (self.Ju[idx][j][k]!=0):
						k0 = k
						self.idxb[idx][j] = k
			self.idxb[idx] = np.ascontiguousarray(self.idxb[idx], dtype=np.int32)
			tmp = cast(self.idxb[idx].ctypes.data, POINTER(c_int))
			self.__hpipm.d_cvt_colmaj_to_ocp_qp_idxb(idx, tmp, self.qp_struct)
		return


	def set_lu(self, lu, idx=None):
		self.__hpipm.d_cvt_colmaj_to_ocp_qp_lb.argtypes = [c_int, POINTER(c_double), c_void_p]
		nbu = self.dim.nbu
		if idx==None:
			for i in range(len(lu)):
				self.lu[i] = lu[i]
				self.lu[i] = self.lu[i].reshape((nbu[i], 1))
				self.lu[i] = np.ascontiguousarray(self.lu[i], dtype=np.float64)
				for j in range(nbu[i]):
					self.lb[i][j] = lu[i][j]
				self.lb[i] = np.ascontiguousarray(self.lb[i], dtype=np.float64)
				tmp = cast(self.lb[i].ctypes.data, POINTER(c_double))
				self.__hpipm.d_cvt_colmaj_to_ocp_qp_lb(i, tmp, self.qp_struct)
		else:
			self.lu[idx] = lu
			self.lu[idx] = self.lu[idx].reshape((nbu[idx], 1))
			self.lu[idx] = np.ascontiguousarray(self.lu[idx], dtype=np.float64)
			for j in range(nbu[idx]):
				self.lb[idx][j] = lu[j]
			self.lb[idx] = np.ascontiguousarray(self.lb[idx], dtype=np.float64)
			tmp = cast(self.lb[idx].ctypes.data, POINTER(c_double))
			self.__hpipm.d_cvt_colmaj_to_ocp_qp_lb(idx, tmp, self.qp_struct)
		return


	def set_uu(self, uu, idx=None):
		self.__hpipm.d_cvt_colmaj_to_ocp_qp_ub.argtypes = [c_int, POINTER(c_double), c_void_p]
		nbu = self.dim.nbu
		if idx==None:
			for i in range(len(uu)):
				self.uu[i] = uu[i]
				self.uu[i] = self.uu[i].reshape((nbu[i], 1))
				self.uu[i] = np.ascontiguousarray(self.uu[i], dtype=np.float64)
				for j in range(nbu[i]):
					self.ub[i][j] = uu[i][j]
				self.ub[i] = np.ascontiguousarray(self.ub[i], dtype=np.float64)
				tmp = cast(self.ub[i].ctypes.data, POINTER(c_double))
				self.__hpipm.d_cvt_colmaj_to_ocp_qp_ub(i, tmp, self.qp_struct)
		else:
			self.uu[idx] = uu
			self.uu[idx] = self.uu[idx].reshape((nbu[idx], 1))
			self.uu[idx] = np.ascontiguousarray(self.uu[idx], dtype=np.float64)
			for j in range(nbu[idx]):
				self.ub[idx][j] = uu[j]
			self.ub[idx] = np.ascontiguousarray(self.ub[idx], dtype=np.float64)
			tmp = cast(self.ub[idx].ctypes.data, POINTER(c_double))
			self.__hpipm.d_cvt_colmaj_to_ocp_qp_ub(idx, tmp, self.qp_struct)
		return


	def set_C(self, C, idx=None):
		self.__hpipm.d_cvt_colmaj_to_ocp_qp_C.argtypes = [c_int, POINTER(c_double), c_void_p]
		nx = self.dim.nx
		ng = self.dim.ng
		if idx==None:
			for i in range(len(C)):
				self.C[i] = C[i]
				self.C[i] = self.C[i].reshape((ng[i], nx[i]))
				self.C[i] = self.C[i].transpose()
				self.C[i] = np.ascontiguousarray(self.C[i], dtype=np.float64)
				tmp = cast(self.C[i].ctypes.data, POINTER(c_double))
				self.__hpipm.d_cvt_colmaj_to_ocp_qp_C(i, tmp, self.qp_struct)
		else:
			self.C[idx] = C
			self.C[idx] = self.C[idx].reshape((ng[idx], nx[idx]))
			self.C[idx] = self.C[idx].transpose()
			self.C[idx] = np.ascontiguousarray(self.C[idx], dtype=np.float64)
			tmp = cast(self.C[idx].ctypes.data, POINTER(c_double))
			self.__hpipm.d_cvt_colmaj_to_ocp_qp_C(idx, tmp, self.qp_struct)
		return


	def set_D(self, D, idx=None):
		self.__hpipm.d_cvt_colmaj_to_ocp_qp_D.argtypes = [c_int, POINTER(c_double), c_void_p]
		nu = self.dim.nu
		ng = self.dim.ng
		if idx==None:
			for i in range(len(D)):
				self.D[i] = D[i]
				self.D[i] = self.D[i].reshape((ng[i], nu[i]))
				self.D[i] = self.D[i].transpose()
				self.D[i] = np.ascontiguousarray(self.D[i], dtype=np.float64)
				tmp = cast(self.D[i].ctypes.data, POINTER(c_double))
				self.__hpipm.d_cvt_colmaj_to_ocp_qp_D(i, tmp, self.qp_struct)
		else:
			self.D[idx] = D
			self.D[idx] = self.D[idx].reshape((ng[idx], nu[idx]))
			self.D[idx] = self.D[idx].transpose()
			self.D[idx] = np.ascontiguousarray(self.D[idx], dtype=np.float64)
			tmp = cast(self.D[idx].ctypes.data, POINTER(c_double))
			self.__hpipm.d_cvt_colmaj_to_ocp_qp_D(idx, tmp, self.qp_struct)
		return


	def set_lg(self, lg, idx=None):
		self.__hpipm.d_cvt_colmaj_to_ocp_qp_lg.argtypes = [c_int, POINTER(c_double), c_void_p]
		ng = self.dim.ng
		if idx==None:
			for i in range(len(lg)):
				self.lg[i] = lg[i]
				self.lg[i] = self.lg[i].reshape((ng[i], 1))
				self.lg[i] = np.ascontiguousarray(self.lg[i], dtype=np.float64)
				tmp = cast(self.lg[i].ctypes.data, POINTER(c_double))
				self.__hpipm.d_cvt_colmaj_to_ocp_qp_lg(i, tmp, self.qp_struct)
		else:
			self.lg[idx] = lg
			self.lg[idx] = self.lg[idx].reshape((ng[idx], 1))
			self.lg[idx] = np.ascontiguousarray(self.lg[idx], dtype=np.float64)
			tmp = cast(self.lg[idx].ctypes.data, POINTER(c_double))
			self.__hpipm.d_cvt_colmaj_to_ocp_qp_lg(idx, tmp, self.qp_struct)
		return


	def set_ug(self, ug, idx=None):
		self.__hpipm.d_cvt_colmaj_to_ocp_qp_ug.argtypes = [c_int, POINTER(c_double), c_void_p]
		ng = self.dim.ng
		if idx==None:
			for i in range(len(ug)):
				self.ug[i] = ug[i]
				self.ug[i] = self.ug[i].reshape((ng[i], 1))
				self.ug[i] = np.ascontiguousarray(self.ug[i], dtype=np.float64)
				tmp = cast(self.ug[i].ctypes.data, POINTER(c_double))
				self.__hpipm.d_cvt_colmaj_to_ocp_qp_ug(i, tmp, self.qp_struct)
		else:
			self.ug[idx] = ug
			self.ug[idx] = self.ug[idx].reshape((ng[idx], 1))
			self.ug[idx] = np.ascontiguousarray(self.ug[idx], dtype=np.float64)
			tmp = cast(self.ug[idx].ctypes.data, POINTER(c_double))
			self.__hpipm.d_cvt_colmaj_to_ocp_qp_ug(idx, tmp, self.qp_struct)
		return


	def set_Zl(self, Zl, idx=None):
		self.__hpipm.d_cvt_colmaj_to_ocp_qp_Zl.argtypes = [c_int, POINTER(c_double), c_void_p]
		ns = self.dim.ns
		if idx==None:
			for i in range(len(Zl)):
				tmp = Zl[i].copy()
				tmp = tmp.reshape((ns[i], ns[i]))
				for j in range(ns[i]):
					self.Zl[i][j] = tmp[j][j]
				self.Zl[i] = self.Zl[i].reshape((ns[i], 1))
				self.Zl[i] = np.ascontiguousarray(self.Zl[i], dtype=np.float64)
				tmp = cast(self.Zl[i].ctypes.data, POINTER(c_double))
				self.__hpipm.d_cvt_colmaj_to_ocp_qp_Zl(i, tmp, self.qp_struct)
		else:
			tmp = Zl.copy()
			tmp = tmp.reshape((ns[idx], ns[idx]))
			for j in range(ns[idx]):
				self.Zl[idx][j] = tmp[j][j]
			self.Zl[idx] = self.Zl[idx].reshape((ns[idx], 1))
			self.Zl[idx] = np.ascontiguousarray(self.Zl[idx], dtype=np.float64)
			tmp = cast(self.Zl[idx].ctypes.data, POINTER(c_double))
			self.__hpipm.d_cvt_colmaj_to_ocp_qp_Zl(idx, tmp, self.qp_struct)
		return


	def set_Zu(self, Zu, idx=None):
		self.__hpipm.d_cvt_colmaj_to_ocp_qp_Zu.argtypes = [c_int, POINTER(c_double), c_void_p]
		ns = self.dim.ns
		if idx==None:
			for i in range(len(Zu)):
				tmp = Zu[i].copy()
				tmp = tmp.reshape((ns[i], ns[i]))
				for j in range(ns[i]):
					self.Zu[i][j] = tmp[j][j]
				self.Zu[i] = self.Zu[i].reshape((ns[i], 1))
				self.Zu[i] = np.ascontiguousarray(self.Zu[i], dtype=np.float64)
				tmp = cast(self.Zu[i].ctypes.data, POINTER(c_double))
				self.__hpipm.d_cvt_colmaj_to_ocp_qp_Zu(i, tmp, self.qp_struct)
		else:
			tmp = Zu.copy()
			tmp = tmp.reshape((ns[idx], ns[idx]))
			for j in range(ns[idx]):
				self.Zu[idx][j] = tmp[j][j]
			self.Zu[idx] = self.Zu[idx].reshape((ns[idx], 1))
			self.Zu[idx] = np.ascontiguousarray(self.Zu[idx], dtype=np.float64)
			tmp = cast(self.Zu[idx].ctypes.data, POINTER(c_double))
			self.__hpipm.d_cvt_colmaj_to_ocp_qp_Zu(idx, tmp, self.qp_struct)
		return


	def set_zl(self, zl, idx=None):
		self.__hpipm.d_cvt_colmaj_to_ocp_qp_zl.argtypes = [c_int, POINTER(c_double), c_void_p]
		ns = self.dim.ns
		if idx==None:
			for i in range(len(zl)):
				self.zl[i] = zl[i]
				self.zl[i] = self.zl[i].reshape((ns[i], 1))
				self.zl[i] = np.ascontiguousarray(self.zl[i], dtype=np.float64)
				tmp = cast(self.zl[i].ctypes.data, POINTER(c_double))
				self.__hpipm.d_cvt_colmaj_to_ocp_qp_zl(i, tmp, self.qp_struct)
		else:
			self.zl[idx] = zl
			self.zl[idx] = self.zl[idx].reshape((ns[idx], 1))
			self.zl[idx] = np.ascontiguousarray(self.zl[idx], dtype=np.float64)
			tmp = cast(self.zl[idx].ctypes.data, POINTER(c_double))
			self.__hpipm.d_cvt_colmaj_to_ocp_qp_zl(idx, tmp, self.qp_struct)
		return


	def set_zu(self, zu, idx=None):
		self.__hpipm.d_cvt_colmaj_to_ocp_qp_zu.argtypes = [c_int, POINTER(c_double), c_void_p]
		ns = self.dim.ns
		if idx==None:
			for i in range(len(zu)):
				self.zu[i] = zu[i]
				self.zu[i] = self.zu[i].reshape((ns[i], 1))
				self.zu[i] = np.ascontiguousarray(self.zu[i], dtype=np.float64)
				tmp = cast(self.zu[i].ctypes.data, POINTER(c_double))
				self.__hpipm.d_cvt_colmaj_to_ocp_qp_zu(i, tmp, self.qp_struct)
		else:
			self.zu[idx] = zu
			self.zu[idx] = self.zu[idx].reshape((ns[idx], 1))
			self.zu[idx] = np.ascontiguousarray(self.zu[idx], dtype=np.float64)
			tmp = cast(self.zu[idx].ctypes.data, POINTER(c_double))
			self.__hpipm.d_cvt_colmaj_to_ocp_qp_zu(idx, tmp, self.qp_struct)
		return


	def set_lls(self, lls, idx=None):
		self.__hpipm.d_cvt_colmaj_to_ocp_qp_lls.argtypes = [c_int, POINTER(c_double), c_void_p]
		ns = self.dim.ns
		if idx==None:
			for i in range(len(lls)):
				self.lls[i] = lls[i]
				self.lls[i] = self.lls[i].reshape((ns[i], 1))
				self.lls[i] = np.ascontiguousarray(self.lls[i], dtype=np.float64)
				tmp = cast(self.lls[i].ctypes.data, POINTER(c_double))
				self.__hpipm.d_cvt_colmaj_to_ocp_qp_lls(i, tmp, self.qp_struct)
		else:
			self.lls[idx] = lls
			self.lls[idx] = self.lls[idx].reshape((ns[idx], 1))
			self.lls[idx] = np.ascontiguousarray(self.lls[idx], dtype=np.float64)
			tmp = cast(self.lls[idx].ctypes.data, POINTER(c_double))
			self.__hpipm.d_cvt_colmaj_to_ocp_qp_lls(idx, tmp, self.qp_struct)
		return


	def set_lus(self, lus, idx=None):
		self.__hpipm.d_cvt_colmaj_to_ocp_qp_lus.argtypes = [c_int, POINTER(c_double), c_void_p]
		ns = self.dim.ns
		if idx==None:
			for i in range(len(lus)):
				self.lus[i] = lus[i]
				self.lus[i] = self.lus[i].reshape((ns[i], 1))
				self.lus[i] = np.ascontiguousarray(self.lus[i], dtype=np.float64)
				tmp = cast(self.lus[i].ctypes.data, POINTER(c_double))
				self.__hpipm.d_cvt_colmaj_to_ocp_qp_lus(i, tmp, self.qp_struct)
		else:
			self.lus[idx] = lus
			self.lus[idx] = self.lus[idx].reshape((ns[idx], 1))
			self.lus[idx] = np.ascontiguousarray(self.lus[idx], dtype=np.float64)
			tmp = cast(self.lus[idx].ctypes.data, POINTER(c_double))
			self.__hpipm.d_cvt_colmaj_to_ocp_qp_lus(idx, tmp, self.qp_struct)
		return


	def set_Jsu(self, Jsu, idx=None):
		self.__hpipm.d_cvt_colmaj_to_ocp_qp_idxs.argtypes = [c_int, POINTER(c_int), c_void_p]
		nbx = self.dim.nbx
		nbu = self.dim.nbu
		ng = self.dim.ng
		ns = self.dim.ns
		if idx==None:
			for i in range(len(Jsu)):
				self.Jsu[i] = Jsu[i]
				self.Jsu[i] = self.Jsu[i].reshape((nbu[i], ns[i]))
				self.Jsu[i] = self.Jsu[i].transpose()
				self.Jsu[i] = np.ascontiguousarray(self.Jsu[i], dtype=np.float64)
				for j in range(ns[i]):
					k0 = -1
					for k in range(nbu[i]):
						if (k0==-1) & (self.Jsu[i][k][j]!=0):
							k0 = k
							self.idxs[i][j] = k
					for k in range(nbx[i]):
						if (k0==-1) & (self.Jsx[i][k][j]!=0):
							k0 = k
							self.idxs[i][j] = nbu[i]+k
					for k in range(ng[i]):
						if (k0==-1) & (self.Jsg[i][k][j]!=0):
							k0 = k
							self.idxs[i][j] = nbu[i]+ng[i]+k
				self.idxs[i] = np.ascontiguousarray(self.idxs[i], dtype=np.int32)
				tmp = cast(self.idxs[i].ctypes.data, POINTER(c_int))
				self.__hpipm.d_cvt_colmaj_to_ocp_qp_idxs(i, tmp, self.qp_struct)
		else:
			self.Jsu[idx] = Jsu
			self.Jsu[idx] = self.Jsu[idx].reshape((nbu[idx], ns[idx]))
			self.Jsu[idx] = self.Jsu[idx].transpose()
			self.Jsu[idx] = np.ascontiguousarray(self.Jsu[idx], dtype=np.float64)
			for j in range(ns[idx]):
				k0 = -1
				for k in range(nbu[idx]):
					if (k0==-1) & (self.Jsu[idx][k][j]!=0):
						k0 = k
						self.idxs[idx][j] = k
				for k in range(nbx[idx]):
					if (k0==-1) & (self.Jsx[idx][k][j]!=0):
						k0 = k
						self.idxs[idx][j] = nbu[idx]+k
				for k in range(ng[idx]):
					if (k0==-1) & (self.Jsg[idx][k][j]!=0):
						k0 = k
						self.idxs[idx][j] = nbu[idx]+ng[idx]+k
			self.idxs[idx] = np.ascontiguousarray(self.idxs[idx], dtype=np.int32)
			tmp = cast(self.idxs[idx].ctypes.data, POINTER(c_int))
			self.__hpipm.d_cvt_colmaj_to_ocp_qp_idxs(idx, tmp, self.qp_struct)
		return


	def set_Jsx(self, Jsx, idx=None):
		self.__hpipm.d_cvt_colmaj_to_ocp_qp_idxs.argtypes = [c_int, POINTER(c_int), c_void_p]
		nbx = self.dim.nbx
		nbu = self.dim.nbu
		ng = self.dim.ng
		ns = self.dim.ns
		if idx==None:
			for i in range(len(Jsx)):
				self.Jsx[i] = Jsx[i]
				self.Jsx[i] = self.Jsx[i].reshape((nbx[i], ns[i]))
				self.Jsx[i] = self.Jsx[i].transpose()
				self.Jsx[i] = np.ascontiguousarray(self.Jsx[i], dtype=np.float64)
				for j in range(ns[i]):
					k0 = -1
					for k in range(nbu[i]):
						if (k0==-1) & (self.Jsu[i][k][j]!=0):
							k0 = k
							self.idxs[i][j] = k
					for k in range(nbx[i]):
						if (k0==-1) & (self.Jsx[i][k][j]!=0):
							k0 = k
							self.idxs[i][j] = nbu[i]+k
					for k in range(ng[i]):
						if (k0==-1) & (self.Jsg[i][k][j]!=0):
							k0 = k
							self.idxs[i][j] = nbu[i]+ng[i]+k
				self.idxs[i] = np.ascontiguousarray(self.idxs[i], dtype=np.int32)
				tmp = cast(self.idxs[i].ctypes.data, POINTER(c_int))
				self.__hpipm.d_cvt_colmaj_to_ocp_qp_idxs(i, tmp, self.qp_struct)
		else:
			self.Jsx[idx] = Jsx
			self.Jsx[idx] = self.Jsx[idx].reshape((nbx[idx], ns[idx]))
			self.Jsx[idx] = self.Jsx[idx].transpose()
			self.Jsx[idx] = np.ascontiguousarray(self.Jsx[idx], dtype=np.float64)
			for j in range(ns[idx]):
				k0 = -1
				for k in range(nbu[idx]):
					if (k0==-1) & (self.Jsu[idx][k][j]!=0):
						k0 = k
						self.idxs[idx][j] = k
				for k in range(nbx[idx]):
					if (k0==-1) & (self.Jsx[idx][k][j]!=0):
						k0 = k
						self.idxs[idx][j] = nbu[idx]+k
				for k in range(ng[idx]):
					if (k0==-1) & (self.Jsg[idx][k][j]!=0):
						k0 = k
						self.idxs[idx][j] = nbu[idx]+ng[idx]+k
			self.idxs[idx] = np.ascontiguousarray(self.idxs[idx], dtype=np.int32)
			tmp = cast(self.idxs[idx].ctypes.data, POINTER(c_int))
			self.__hpipm.d_cvt_colmaj_to_ocp_qp_idxs(idx, tmp, self.qp_struct)
		return


	def set_Jsg(self, Jsg, idx=None):
		self.__hpipm.d_cvt_colmaj_to_ocp_qp_idxs.argtypes = [c_int, POINTER(c_int), c_void_p]
		nbx = self.dim.nbx
		nbu = self.dim.nbu
		ng = self.dim.ng
		ns = self.dim.ns
		if idx==None:
			for i in range(len(Jsg)):
				self.Jsg[i] = Jsg[i]
				self.Jsg[i] = self.Jsg[i].reshape((ng[i], ns[i]))
				self.Jsg[i] = self.Jsg[i].transpose()
				self.Jsg[i] = np.ascontiguousarray(self.Jsg[i], dtype=np.float64)
				for j in range(ns[i]):
					k0 = -1
					for k in range(nbu[i]):
						if (k0==-1) & (self.Jsu[i][k][j]!=0):
							k0 = k
							self.idxs[i][j] = k
					for k in range(nbx[i]):
						if (k0==-1) & (self.Jsx[i][k][j]!=0):
							k0 = k
							self.idxs[i][j] = nbu[i]+k
					for k in range(ng[i]):
						if (k0==-1) & (self.Jsg[i][k][j]!=0):
							k0 = k
							self.idxs[i][j] = nbu[i]+ng[i]+k
				self.idxs[i] = np.ascontiguousarray(self.idxs[i], dtype=np.int32)
				tmp = cast(self.idxs[i].ctypes.data, POINTER(c_int))
				self.__hpipm.d_cvt_colmaj_to_ocp_qp_idxs(i, tmp, self.qp_struct)
		else:
			self.Jsg[idx] = Jsg
			self.Jsg[idx] = self.Jsg[idx].reshape((ng[idx], ns[idx]))
			self.Jsg[idx] = self.Jsg[idx].transpose()
			self.Jsg[idx] = np.ascontiguousarray(self.Jsg[idx], dtype=np.float64)
			for j in range(ns[idx]):
				k0 = -1
				for k in range(nbu[idx]):
					if (k0==-1) & (self.Jsu[idx][k][j]!=0):
						k0 = k
						self.idxs[idx][j] = k
				for k in range(nbx[idx]):
					if (k0==-1) & (self.Jsx[idx][k][j]!=0):
						k0 = k
						self.idxs[idx][j] = nbu[idx]+k
				for k in range(ng[idx]):
					if (k0==-1) & (self.Jsg[idx][k][j]!=0):
						k0 = k
						self.idxs[idx][j] = nbu[idx]+ng[idx]+k
			self.idxs[idx] = np.ascontiguousarray(self.idxs[idx], dtype=np.int32)
			tmp = cast(self.idxs[idx].ctypes.data, POINTER(c_int))
			self.__hpipm.d_cvt_colmaj_to_ocp_qp_idxs(idx, tmp, self.qp_struct)
		return
	

	def print_C_struct(self):
		self.__hpipm.d_print_ocp_qp(self.qp_struct)




