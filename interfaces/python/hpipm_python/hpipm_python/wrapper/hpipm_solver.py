from ctypes import *
import ctypes.util 
import numpy as np
import faulthandler

faulthandler.enable()



class hpipm_solver:
	def __init__(self, qp_dims, qp_data):

		# load hpipm shared library
		__hpipm   = CDLL('libhpipm.so')
		self.__hpipm = __hpipm

		# extract dim
		N	= qp_dims.N

		# array of pointer
		A = (POINTER(c_double)*(N))()
		B = (POINTER(c_double)*(N))()
		b = (POINTER(c_double)*(N))()

		Q = (POINTER(c_double)*(N+1))()
		S = (POINTER(c_double)*(N+1))()
		R = (POINTER(c_double)*(N+1))()
		q = (POINTER(c_double)*(N+1))()
		r = (POINTER(c_double)*(N+1))()

		idxb = (POINTER(c_int)*(N+1))()
		lb = (POINTER(c_double)*(N+1))()
		ub = (POINTER(c_double)*(N+1))()

		C = (POINTER(c_double)*(N+1))()
		D = (POINTER(c_double)*(N+1))()
		lg = (POINTER(c_double)*(N+1))()
		ug = (POINTER(c_double)*(N+1))()

		Zl = (POINTER(c_double)*(N+1))()
		Zu = (POINTER(c_double)*(N+1))()
		zl = (POINTER(c_double)*(N+1))()
		zu = (POINTER(c_double)*(N+1))()

		idxs = (POINTER(c_int)*(N+1))()
		lls = (POINTER(c_double)*(N+1))()
		lus = (POINTER(c_double)*(N+1))()


		for i in range(N+1):

			if(i<N):
				# dynamics
				A[i] = cast(qp_data.A[i].ctypes.data, POINTER(c_double))
				B[i] = cast(qp_data.B[i].ctypes.data, POINTER(c_double))
				b[i] = cast(qp_data.b[i].ctypes.data, POINTER(c_double))

			# cost
			Q[i] = cast(qp_data.Q[i].ctypes.data, POINTER(c_double))
			S[i] = cast(qp_data.S[i].ctypes.data, POINTER(c_double))
			R[i] = cast(qp_data.R[i].ctypes.data, POINTER(c_double))
			q[i] = cast(qp_data.q[i].ctypes.data, POINTER(c_double))
			r[i] = cast(qp_data.r[i].ctypes.data, POINTER(c_double))

			# simple bounds
			if qp_dims.nbx[i]+qp_dims.nbu[i] > 0:
				lb[i] = cast(qp_data.lb[i].ctypes.data, POINTER(c_double))
				ub[i] = cast(qp_data.ub[i].ctypes.data, POINTER(c_double))
				idxb[i] = cast(qp_data.idxb[i].ctypes.data, POINTER(c_int))

			# polytopic constraints
			if qp_dims.ng[i] > 0:
				C[i] = cast(qp_data.C[i].ctypes.data, POINTER(c_double))
				D[i] = cast(qp_data.D[i].ctypes.data, POINTER(c_double))
				lg[i] = cast(qp_data.lg[i].ctypes.data, POINTER(c_double))
				ug[i] = cast(qp_data.ug[i].ctypes.data, POINTER(c_double))

			# slacks
			if qp_dims.ns[i] > 0:
				Zl[i] = cast(qp_data.Zl[i].ctypes.data, POINTER(c_double))
				Zu[i] = cast(qp_data.Zu[i].ctypes.data, POINTER(c_double))
				zl[i] = cast(qp_data.zl[i].ctypes.data, POINTER(c_double))
				zu[i] = cast(qp_data.zu[i].ctypes.data, POINTER(c_double))
				lls[i] = cast(qp_data.lls[i].ctypes.data, POINTER(c_double))
				lus[i] = cast(qp_data.lus[i].ctypes.data, POINTER(c_double))
				idxs[i] = cast(qp_data.idxs[i].ctypes.data, POINTER(c_int))

		# allocate memory for qp struct 
		qp_size = __hpipm.d_memsize_ocp_qp(qp_dims.dim)
		qp_mem = cast(create_string_buffer(qp_size), c_void_p)
		self.qp_mem = qp_mem

		# set up ocp_qp structure
		sizeof_d_ocp_qp = __hpipm.d_sizeof_ocp_qp()
		qp = cast(create_string_buffer(sizeof_d_ocp_qp), c_void_p)
		self.ocp_qp = qp

		__hpipm.d_create_ocp_qp(qp_dims.dim, qp, qp_mem)
		__hpipm.d_cvt_colmaj_to_ocp_qp(A, B, b, Q, S, R, q, r, idxb, lb, ub,
			C, D, lg, ug, Zl, Zu, zl, zu, idxs, lls, lus, qp)
		
		# allocate memory for ocp_qp_sol struct
		qp_sol_size = __hpipm.d_memsize_ocp_qp_sol(qp_dims.dim)
		qp_sol_mem = cast(create_string_buffer(qp_sol_size), c_void_p)
		self.qp_sol_mem = qp_sol_mem

		# set up ocp_qp_sol struct
		sizeof_d_ocp_qp_sol = __hpipm.d_sizeof_ocp_qp_sol()
		qp_sol = cast(create_string_buffer(sizeof_d_ocp_qp_sol), c_void_p)
		__hpipm.d_create_ocp_qp_sol(qp_dims.dim, qp_sol, qp_sol_mem)
		self.ocp_qp_sol = qp_sol

		# allocate memory for ipm_arg struct
		ipm_arg_size = __hpipm.d_memsize_ocp_qp_ipm_arg(qp_dims.dim)
		ipm_arg_mem = cast(create_string_buffer(ipm_arg_size), c_void_p)
		self.ipm_arg_mem = ipm_arg_mem

		# set up ipm_arg
		sizeof_d_ocp_qp_ipm_arg = __hpipm.d_sizeof_ocp_qp_ipm_arg()
		arg = cast(create_string_buffer(sizeof_d_ocp_qp_ipm_arg), c_void_p)
		self.ocp_qp_ipm_arg = arg

		__hpipm.d_create_ocp_qp_ipm_arg(qp_dims.dim, arg, ipm_arg_mem)
		__hpipm.d_set_default_ocp_qp_ipm_arg(1, arg)

		# allocate memory for ipm workspace 
		ipm_size = __hpipm.d_memsize_ocp_qp_ipm(qp_dims.dim, arg)
		ipm_mem = cast(create_string_buffer(ipm_size), c_void_p)
		self.ipm_mem = ipm_mem

		# set up ipm workspace
		sizeof_d_ocp_qp_ipm_workspace = __hpipm.d_sizeof_ocp_qp_ipm_workspace()
		workspace = cast(create_string_buffer(sizeof_d_ocp_qp_ipm_workspace), c_void_p)
		self.ocp_qp_ipm_workspace = workspace

		__hpipm.d_create_ocp_qp_ipm(qp_dims.dim, arg, workspace, ipm_mem)

		self.qp = qp
		self.qp_sol = qp_sol
		self.arg = arg
		self.workspace = workspace
		self.dim = qp_dims.dim
		
	def solve(self):
		return self.__hpipm.d_solve_ocp_qp_ipm(self.qp, self.qp_sol, 
			self.arg, self.workspace)

	def print_sol(self):
		self.__hpipm.d_print_ocp_qp_sol(self.ocp_qp_sol, self.dim)
		return 



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

		# allocate memory for C dimemsions struct
		dim_struct_size = __hpipm.d_sizeof_ocp_qp_dim()
		dim_struct = cast(create_string_buffer(dim_struct_size), c_void_p)
		self.dim = dim_struct

		dim_mem_size = __hpipm.d_memsize_ocp_qp_dim(N)
		dim_mem = cast(create_string_buffer(dim_mem_size), c_void_p)
		self.dim_mem = dim_mem

		# set up dimensions structure
		__hpipm.d_create_ocp_qp_dim(N, self.dim, self.dim_mem)

#		__hpipm.d_cvt_int_to_ocp_qp_dim(N,
#			cast(nx.ctypes.data, POINTER(c_double)),
#			cast(nu.ctypes.data, POINTER(c_double)),
#			cast(nbx.ctypes.data, POINTER(c_double)),
#			cast(nbu.ctypes.data, POINTER(c_double)),
#			cast(ng.ctypes.data, POINTER(c_double)),
#			cast(ns.ctypes.data, POINTER(c_double)),
#			self.dim)


	def set_nx(self, nx, idx=None):
#		self.__hpipm.d_set_ocp_qp_dim_nx.restypes = [None]
		self.__hpipm.d_set_ocp_qp_dim_nx.argtypes = [c_int, c_int, c_void_p]
		if idx==None:
			for i in range(nx.size):
				self.nx[i] = nx[i]
				self.__hpipm.d_set_ocp_qp_dim_nx(i, nx[i], self.dim)
		else:
			self.nx[idx] = nx
			self.__hpipm.d_set_ocp_qp_dim_nx(idx, nx, self.dim)
		return

	def set_nu(self, nu, idx=None):
		self.__hpipm.d_set_ocp_qp_dim_nu.argtypes = [c_int, c_int, c_void_p]
		if idx==None:
			for i in range(nu.size):
				self.nu[i] = nu[i]
				self.__hpipm.d_set_ocp_qp_dim_nu(i, nu[i], self.dim)
		else:
			self.nu[idx] = nu
			self.__hpipm.d_set_ocp_qp_dim_nu(idx, nu, self.dim)
		return

	def set_nbx(self, nbx, idx=None):
		self.__hpipm.d_set_ocp_qp_dim_nbx.argtypes = [c_int, c_int, c_void_p]
		if idx==None:
			for i in range(nbx.size):
				self.nbx[i] = nbx[i]
				self.__hpipm.d_set_ocp_qp_dim_nbx(i, nbx[i], self.dim)
		else:
			self.nbx[idx] = nbx
			self.__hpipm.d_set_ocp_qp_dim_nbx(idx, nbx, self.dim)
		return

	def set_nbu(self, nbu, idx=None):
		self.__hpipm.d_set_ocp_qp_dim_nbu.argtypes = [c_int, c_int, c_void_p]
		if idx==None:
			for i in range(nbu.size):
				self.nbu[i] = nbu[i]
				self.__hpipm.d_set_ocp_qp_dim_nbu(i, nbu[i], self.dim)
		else:
			self.nbu[idx] = nbu
			self.__hpipm.d_set_ocp_qp_dim_nbu(idx, nbu, self.dim)
		return

	def set_ng(self, ng, idx=None):
		self.__hpipm.d_set_ocp_qp_dim_ng.argtypes = [c_int, c_int, c_void_p]
		if idx==None:
			for i in range(ng.size):
				self.ng[i] = ng[i]
				self.__hpipm.d_set_ocp_qp_dim_ng(i, ng[i], self.dim)
		else:
			self.ng[idx] = ng
			self.__hpipm.d_set_ocp_qp_dim_ng(idx, ng, self.dim)
		return

	def set_ns(self, ns, idx=None):
		self.__hpipm.d_set_ocp_qp_dim_ns.argtypes = [c_int, c_int, c_void_p]
		if idx==None:
			for i in range(ns.size):
				self.ns[i] = ns[i]
				self.__hpipm.d_set_ocp_qp_dim_ns(i, ns[i], self.dim)
		else:
			self.ns[idx] = ns
			self.__hpipm.d_set_ocp_qp_dim_ns(idx, ns, self.dim)
		return



class hpipm_ocp_qp:
	def __init__(self, dims):
		N = dims.N
		nx = dims.nx
		nu = dims.nu
		nbx = dims.nbx
		nbu = dims.nbu
		ng = dims.ng
		ns = dims.ns

		self.dims = dims
#		self.dims = hpipm_ocp_qp_dim(N)
#		self.dims.set_nx(nx);
#		self.dims.set_nu(nu);
#		self.dims.set_nbx(nbx);
#		self.dims.set_nbu(nbu);
#		self.dims.set_ng(ng);
#		self.dims.set_ns(ns);

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



	def set_A(self, A, idx=None):
		nx = self.dims.nx
		if idx==None:
			for i in range(len(A)):
				self.A[i] = A[i]
				self.A[i] = np.ascontiguousarray(self.A[i], dtype=np.float64)
				self.A[i] = self.A[i].reshape((nx[i+1], nx[i]))
		else:
			self.A[idx] = A
			self.A[idx] = np.ascontiguousarray(self.A[idx], dtype=np.float64)
			self.A[idx] = self.A[i].reshape((nx[idx+1], nx[idx]))
		return

	def set_B(self, B, idx=None):
		nx = self.dims.nx
		nu = self.dims.nu
		if idx==None:
			for i in range(len(B)):
				self.B[i] = B[i]
				self.B[i] = np.ascontiguousarray(self.B[i], dtype=np.float64)
				self.B[i] = self.B[i].reshape((nx[i+1], nu[i]))
		else:
			self.B[idx] = B
			self.B[idx] = np.ascontiguousarray(self.B[idx], dtype=np.float64)
			self.B[idx] = self.B[idx].reshape((nx[idx+1], nu[idx]))
		return

	def set_b(self, b, idx=None):
		nx = self.dims.nx
		if idx==None:
			for i in range(len(b)):
				self.b[i] = b[i]
				self.b[i] = np.ascontiguousarray(self.b[i], dtype=np.float64)
				self.b[i] = self.b[i].reshape((nx[i+1], 1))
		else:
			self.b[idx] = b
			self.b[idx] = np.ascontiguousarray(self.b[idx], dtype=np.float64)
			self.b[idx] = self.b[idx].reshape((nx[idx+1], 1))
		return

	def set_Q(self, Q, idx=None):
		nx = self.dims.nx
		if idx==None:
			for i in range(len(Q)):
				self.Q[i] = Q[i]
				self.Q[i] = np.ascontiguousarray(self.Q[i], dtype=np.float64)
				self.Q[i] = self.Q[i].reshape((nx[i], nx[i]))
		else:
			self.Q[idx] = Q
			self.Q[idx] = np.ascontiguousarray(self.Q[idx], dtype=np.float64)
			self.Q[idx] = self.Q[idx].reshape((nx[idx], nx[idx]))
		return

	def set_R(self, R, idx=None):
		nu = self.dims.nu
		if idx==None:
			for i in range(len(R)):
				self.R[i] = R[i]
				self.R[i] = np.ascontiguousarray(self.R[i], dtype=np.float64)
				self.R[i] = self.R[i].reshape((nu[i], nu[i]))
		else:
			self.R[idx] = R
			self.R[idx] = np.ascontiguousarray(self.R[idx], dtype=np.float64)
			self.R[idx] = self.R[idx].reshape((nu[idx], nu[idx]))
		return

	def set_S(self, S, idx=None):
		nu = self.dims.nu
		nx = self.dims.nx
		if idx==None:
			for i in range(len(S)):
				self.S[i] = S[i]
				self.S[i] = np.ascontiguousarray(self.S[i], dtype=np.float64)
				self.S[i] = self.S[i].reshape((nx[i], nu[i]))
		else:
			self.S[idx] = S
			self.S[idx] = np.ascontiguousarray(self.S[idx], dtype=np.float64)
			self.S[idx] = self.S[idx].reshape((nu[idx], nx[idx]))
		return

	def set_q(self, q, idx=None):
		nx = self.dims.nx
		if idx==None:
			for i in range(len(q)):
				self.q[i] = q[i]
				self.q[i] = np.ascontiguousarray(self.q[i], dtype=np.float64)
				self.q[i] = self.q[i].reshape((nx[i], 1))
		else:
			self.q[idx] = q
			self.q[idx] = np.ascontiguousarray(self.q[idx], dtype=np.float64)
			self.q[idx] = self.q[idx].reshape((nx[idx], 1))
		return

	def set_r(self, r, idx=None):
		nu = self.dims.nu
		if idx==None:
			for i in range(len(r)):
				self.r[i] = r[i]
				self.r[i] = np.ascontiguousarray(self.r[i], dtype=np.float64)
				self.r[i] = self.r[i].reshape((nx[i], 1))
		else:
			self.r[idx] = r
			self.r[idx] = np.ascontiguousarray(self.r[idx], dtype=np.float64)
			self.r[idx] = self.r[idx].reshape((nx[idx], 1))
		return

	def set_Jx(self, Jx, idx=None):
		nx = self.dims.nx
		nu = self.dims.nu
		nbx = self.dims.nbx
		nbu = self.dims.nbu
		if idx==None:
			for i in range(len(Jx)):
				self.Jx[i] = Jx[i]
				self.Jx[i] = np.ascontiguousarray(self.Jx[i], dtype=np.float64)
				self.Jx[i] = self.Jx[i].reshape((nbx[i], nx[i]))
				for j in range(nbx[i]):
					k0 = -1
					for k in range(nx[i]):
						if (k0==-1) & (self.Jx[i][j][k]!=0):
							k0 = k
							self.idxb[i][nbu[i]+j] = nu[i]+k
				self.idxb[i] = np.ascontiguousarray(self.idxb[i], dtype=np.int32)
		else:
			self.Jx[idx] = Jx
			self.Jx[idx] = np.ascontiguousarray(self.Jx[idx], dtype=np.float64)
			self.Jx[idx] = self.Jx[idx].reshape((nbx[idx], nx[idx]))
			for j in range(nbx[idx]):
				k0 = -1
				for k in range(nx[idx]):
					if (k0==-1) & (self.Jx[idx][j][k]!=0):
						k0 = k
						self.idxb[idx][nbu[idx]+j] = nu[idx]+k
			self.idxb[idx] = np.ascontiguousarray(self.idxb[idx], dtype=np.int32)
		return

	def set_lx(self, lx, idx=None):
		nbx = self.dims.nbx
		nbu = self.dims.nbu
		if idx==None:
			for i in range(len(lx)):
				self.lx[i] = lx[i]
				self.lx[i] = np.ascontiguousarray(self.lx[i], dtype=np.float64)
				self.lx[i] = self.lx[i].reshape((nbx[i], 1))
				for j in range(nbx[i]):
					self.lb[i][nbu[i]+j] = lx[i][j]
				self.lb[i] = np.ascontiguousarray(self.lb[i], dtype=np.float64)
		else:
			self.lx[idx] = lx
			self.lx[idx] = np.ascontiguousarray(self.lx[idx], dtype=np.float64)
			self.lx[idx] = self.lx[idx].reshape((nbx[idx], 1))
			for j in range(nbx[idx]):
				self.lb[idx][nbu[idx]+j] = lx[j]
			self.lb[idx] = np.ascontiguousarray(self.lb[idx], dtype=np.float64)
		return

	def set_ux(self, ux, idx=None):
		nbx = self.dims.nbx
		nbu = self.dims.nbu
		if idx==None:
			for i in range(len(ux)):
				self.ux[i] = ux[i]
				self.ux[i] = np.ascontiguousarray(self.ux[i], dtype=np.float64)
				self.ux[i] = self.ux[i].reshape((nbx[i], 1))
				for j in range(nbx[i]):
					self.ub[i][nbu[i]+j] = ux[i][j]
				self.ub[i] = np.ascontiguousarray(self.ub[i], dtype=np.float64)
		else:
			self.ux[idx] = ux
			self.ux[idx] = np.ascontiguousarray(self.ux[idx], dtype=np.float64)
			self.ux[idx] = self.ux[idx].reshape((nbx[idx], 1))
			for j in range(nbx[idx]):
				self.ub[idx][nbu[idx]+j] = ux[j]
			self.ub[idx] = np.ascontiguousarray(self.ub[idx], dtype=np.float64)
		return

	def set_Ju(self, Ju, idx=None):
		nu = self.dims.nu
		nbu = self.dims.nbu
		if idx==None:
			for i in range(len(Ju)):
				self.Ju[i] = Ju[i]
				self.Ju[i] = np.ascontiguousarray(self.Ju[i], dtype=np.float64)
				self.Ju[i] = self.Ju[i].reshape((nbu[i], nu[i]))
				for j in range(nbu[i]):
					k0 = -1
					for k in range(nu[i]):
						if (k0==-1) & (self.Ju[i][j][k]!=0):
							k0 = k
							self.idxb[i][j] = k
				self.idxb[i] = np.ascontiguousarray(self.idxb[i], dtype=np.int32)
		else:
			self.Ju[idx] = Ju
			self.Ju[idx] = np.ascontiguousarray(self.Ju[idx], dtype=np.float64)
			self.Ju[idx] = self.Ju[idx].reshape((nbu[idx], nu[idx]))
			for j in range(nbu[idx]):
				k0 = -1
				for k in range(nu[idx]):
					if (k0==-1) & (self.Ju[idx][j][k]!=0):
						k0 = k
						self.idxb[idx][j] = k
			self.idxb[idx] = np.ascontiguousarray(self.idxb[idx], dtype=np.int32)
		return

	def set_lu(self, lu, idx=None):
		nbu = self.dims.nbu
		if idx==None:
			for i in range(len(lu)):
				self.lu[i] = lu[i]
				self.lu[i] = np.ascontiguousarray(self.lu[i], dtype=np.float64)
				self.lu[i] = self.lu[i].reshape((nbu[i], 1))
				for j in range(nbu[i]):
					self.lb[i][j] = lu[i][j]
				self.lb[i] = np.ascontiguousarray(self.lb[i], dtype=np.float64)
		else:
			self.lu[idx] = lu
			self.lu[idx] = np.ascontiguousarray(self.lu[idx], dtype=np.float64)
			self.lu[idx] = self.lu[idx].reshape((nbu[idx], 1))
			for j in range(nbu[idx]):
				self.lb[idx][j] = lu[j]
			self.lb[idx] = np.ascontiguousarray(self.lb[idx], dtype=np.float64)
		return

	def set_uu(self, uu, idx=None):
		nbu = self.dims.nbu
		if idx==None:
			for i in range(len(uu)):
				self.uu[i] = uu[i]
				self.uu[i] = np.ascontiguousarray(self.uu[i], dtype=np.float64)
				self.uu[i] = self.uu[i].reshape((nbu[i], 1))
				for j in range(nbu[i]):
					self.ub[i][j] = uu[i][j]
				self.ub[i] = np.ascontiguousarray(self.ub[i], dtype=np.float64)
		else:
			self.uu[idx] = uu
			self.uu[idx] = np.ascontiguousarray(self.uu[idx], dtype=np.float64)
			self.uu[idx] = self.uu[idx].reshape((nbu[idx], 1))
			for j in range(nbu[idx]):
				self.ub[idx][j] = uu[j]
			self.ub[idx] = np.ascontiguousarray(self.ub[idx], dtype=np.float64)
		return

	def set_C(self, C, idx=None):
		nx = self.dims.nx
		ng = self.dims.ng
		if idx==None:
			for i in range(len(C)):
				self.C[i] = C[i]
				self.C[i] = np.ascontiguousarray(self.C[i], dtype=np.float64)
				self.C[i] = self.C[i].reshape((ng[i], nx[i]))
		else:
			self.C[idx] = C
			self.C[idx] = np.ascontiguousarray(self.C[idx], dtype=np.float64)
			self.C[idx] = self.C[idx].reshape((ng[idx], nx[idx]))
		return

	def set_D(self, D, idx=None):
		nu = self.dims.nu
		ng = self.dims.ng
		if idx==None:
			for i in range(len(D)):
				self.D[i] = D[i]
				self.D[i] = np.ascontiguousarray(self.D[i], dtype=np.float64)
				self.D[i] = self.D[i].reshape((ng[i], nu[i]))
		else:
			self.D[idx] = D
			self.D[idx] = np.ascontiguousarray(self.D[idx], dtype=np.float64)
			self.D[idx] = self.D[idx].reshape((ng[idx], nu[idx]))
		return

	def set_lg(self, lg, idx=None):
		ng = self.dims.ng
		if idx==None:
			for i in range(len(lg)):
				self.lg[i] = lg[i]
				self.lg[i] = np.ascontiguousarray(self.lg[i], dtype=np.float64)
				self.lg[i] = self.lg[i].reshape((ng[i], 1))
		else:
			self.lg[idx] = lg
			self.lg[idx] = np.ascontiguousarray(self.lg[idx], dtype=np.float64)
			self.lg[idx] = self.lg[idx].reshape((ng[idx], 1))
		return

	def set_ug(self, ug, idx=None):
		ng = self.dims.ng
		if idx==None:
			for i in range(len(ug)):
				self.ug[i] = ug[i]
				self.ug[i] = np.ascontiguousarray(self.ug[i], dtype=np.float64)
				self.ug[i] = self.ug[i].reshape((ng[i], 1))
		else:
			self.ug[idx] = ug
			self.ug[idx] = np.ascontiguousarray(self.ug[idx], dtype=np.float64)
			self.ug[idx] = self.ug[idx].reshape((ng[idx], 1))
		return

	def set_Zl(self, Zl, idx=None):
		ns = self.dims.ns
		if idx==None:
			for i in range(len(Zl)):
				tmp = Zl[i].copy()
				tmp = tmp.reshape((ns[i], ns[i]))
				for j in range(ns[i]):
					self.Zl[i][j] = tmp[j][j]
				self.Zl[i] = np.ascontiguousarray(self.Zl[i], dtype=np.float64)
				self.Zl[i] = self.Zl[i].reshape((ns[i], 1))
		else:
			tmp = Zl.copy()
			tmp = tmp.reshape((ns[idx], ns[idx]))
			for j in range(ns[idx]):
				self.Zl[idx][j] = tmp[j][j]
			self.Zl[idx] = np.ascontiguousarray(self.Zl[idx], dtype=np.float64)
			self.Zl[idx] = self.Zl[idx].reshape((ns[idx], 1))
		return

	def set_Zu(self, Zu, idx=None):
		ns = self.dims.ns
		if idx==None:
			for i in range(len(Zu)):
				tmp = Zu[i].copy()
				tmp = tmp.reshape((ns[i], ns[i]))
				for j in range(ns[i]):
					self.Zu[i][j] = tmp[j][j]
				self.Zu[i] = np.ascontiguousarray(self.Zu[i], dtype=np.float64)
				self.Zu[i] = self.Zu[i].reshape((ns[i], 1))
		else:
			tmp = Zu.copy()
			tmp = tmp.reshape((ns[idx], ns[idx]))
			for j in range(ns[idx]):
				self.Zu[idx][j] = tmp[j][j]
			self.Zu[idx] = np.ascontiguousarray(self.Zu[idx], dtype=np.float64)
			self.Zu[idx] = self.Zu[idx].reshape((ns[idx], 1))
		return

	def set_zl(self, zl, idx=None):
		ns = self.dims.ns
		if idx==None:
			for i in range(len(zl)):
				self.zl[i] = zl[i]
				self.zl[i] = np.ascontiguousarray(self.zl[i], dtype=np.float64)
				self.zl[i] = self.zl[i].reshape((ns[i], 1))
		else:
			self.zl[idx] = zl
			self.zl[idx] = np.ascontiguousarray(self.zl[idx], dtype=np.float64)
			self.zl[idx] = self.zl[idx].reshape((ns[idx], 1))
		return

	def set_zu(self, zu, idx=None):
		ns = self.dims.ns
		if idx==None:
			for i in range(len(zu)):
				self.zu[i] = zu[i]
				self.zu[i] = np.ascontiguousarray(self.zu[i], dtype=np.float64)
				self.zu[i] = self.zu[i].reshape((ns[i], 1))
		else:
			self.zu[idx] = zu
			self.zu[idx] = np.ascontiguousarray(self.zu[idx], dtype=np.float64)
			self.zu[idx] = self.zu[idx].reshape((ns[idx], 1))
		return

	def set_lls(self, lls, idx=None):
		ns = self.dims.ns
		if idx==None:
			for i in range(len(lls)):
				self.lls[i] = lls[i]
				self.lls[i] = np.ascontiguousarray(self.lls[i], dtype=np.float64)
				self.lls[i] = self.lls[i].reshape((ns[i], 1))
		else:
			self.lls[idx] = lls
			self.lls[idx] = np.ascontiguousarray(self.lls[idx], dtype=np.float64)
			self.lls[idx] = self.lls[idx].reshape((ns[idx], 1))
		return

	def set_lus(self, lus, idx=None):
		ns = self.dims.ns
		if idx==None:
			for i in range(len(lus)):
				self.lus[i] = lus[i]
				self.lus[i] = np.ascontiguousarray(self.lus[i], dtype=np.float64)
				self.lus[i] = self.lus[i].reshape((ns[i], 1))
		else:
			self.lus[idx] = lus
			self.lus[idx] = np.ascontiguousarray(self.lus[idx], dtype=np.float64)
			self.lus[idx] = self.lus[idx].reshape((ns[idx], 1))
		return

	def set_Jsu(self, Jsu, idx=None):
		nbx = self.dims.nbx
		nbu = self.dims.nbu
		ng = self.dims.ng
		ns = self.dims.ns
		if idx==None:
			for i in range(len(Jsu)):
				self.Jsu[i] = Jsu[i]
				self.Jsu[i] = np.ascontiguousarray(self.Jsu[i], dtype=np.float64)
				self.Jsu[i] = self.Jsu[i].reshape((nbu[i], ns[i]))
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
		else:
			self.Jsu[idx] = Jsu
			self.Jsu[idx] = np.ascontiguousarray(self.Jsu[idx], dtype=np.float64)
			self.Jsu[idx] = self.Jsu[idx].reshape((nbu[idx], ns[idx]))
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
		return

	def set_Jsx(self, Jsx, idx=None):
		nbx = self.dims.nbx
		nbu = self.dims.nbu
		ng = self.dims.ng
		ns = self.dims.ns
		if idx==None:
			for i in range(len(Jsx)):
				self.Jsx[i] = Jsx[i]
				self.Jsx[i] = np.ascontiguousarray(self.Jsx[i], dtype=np.float64)
				self.Jsx[i] = self.Jsx[i].reshape((nbx[i], ns[i]))
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
		else:
			self.Jsx[idx] = Jsx
			self.Jsx[idx] = np.ascontiguousarray(self.Jsx[idx], dtype=np.float64)
			self.Jsx[idx] = self.Jsx[idx].reshape((nbx[idx], ns[idx]))
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
		return

	def set_Jsg(self, Jsg, idx=None):
		nbx = self.dims.nbx
		nbu = self.dims.nbu
		ng = self.dims.ng
		ns = self.dims.ns
		if idx==None:
			for i in range(len(Jsg)):
				self.Jsg[i] = Jsg[i]
				self.Jsg[i] = np.ascontiguousarray(self.Jsg[i], dtype=np.float64)
				self.Jsg[i] = self.Jsg[i].reshape((ng[i], ns[i]))
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
		else:
			self.Jsg[idx] = Jsg
			self.Jsg[idx] = np.ascontiguousarray(self.Jsg[idx], dtype=np.float64)
			self.Jsg[idx] = self.Jsg[idx].reshape((ng[idx], ns[idx]))
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
		return



class hpipm_ocp_qp_sol:
	def __init__(self):
		self.ux = None
		self.pi = None
		self.lam = None
		self.t = None

