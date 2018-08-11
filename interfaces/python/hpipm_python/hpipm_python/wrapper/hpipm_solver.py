from ctypes import *
import ctypes.util 
import numpy as np
import faulthandler

faulthandler.enable()



class hpipm_solver:
	def __init__(self, qp_dims, qp_data):

		# load blasfeo and hpipm shared libraries
		__blasfeo = CDLL('libblasfeo.so')
		__hpipm   = CDLL('libhpipm.so')

		# cast dimensions to contiguous int
		N	= qp_dims.N
		nx   = np.ascontiguousarray(qp_dims.nx,  dtype=np.int32)
		nu   = np.ascontiguousarray(qp_dims.nu,  dtype=np.int32)
		nbx  = np.ascontiguousarray(qp_dims.nbx, dtype=np.int32)
		nbu  = np.ascontiguousarray(qp_dims.nbu, dtype=np.int32)
		ng   = np.ascontiguousarray(qp_dims.ng,  dtype=np.int32)
		ns   = np.ascontiguousarray(qp_dims.ns,  dtype=np.int32)

		# allocate memory for dimemsions struct
		sizeof_d_ocp_qp_dim = __hpipm.d_sizeof_ocp_qp_dim()
		dim = cast(create_string_buffer(sizeof_d_ocp_qp_dim), c_void_p)
		self.ocp_qp_dim = dim

		dim_size = __hpipm.d_memsize_ocp_qp_dim(qp_dims.N)
		dim_mem = cast(create_string_buffer(dim_size), c_void_p)
		self.dim_mem = dim_mem

		# set up dimensions structure
		__hpipm.d_create_ocp_qp_dim(N, dim, dim_mem)
		__hpipm.d_cvt_int_to_ocp_qp_dim(N,
			cast(nx.ctypes.data, POINTER(c_double)),
			cast(nu.ctypes.data, POINTER(c_double)),
			cast(nbx.ctypes.data, POINTER(c_double)),
			cast(nbu.ctypes.data, POINTER(c_double)),
			cast(ng.ctypes.data, POINTER(c_double)),
			cast(ns.ctypes.data, POINTER(c_double)),
			dim)


		# new memory to handle change of API

		iidxb = []
		for i in range(N+1):
			iidxb.append(np.zeros((nbx[i]+nbu[i], 1), dtype=np.int32))
			for j in range(nbu[i]):
				k0 = -1
				for k in range(nu[i]):
					if (k0==-1) & (qp_data.Ju[i][j][k]!=0):
						k0 = k
						iidxb[i][j] = k
			for j in range(nbx[i]):
				k0 = -1
				for k in range(nx[i]):
					if (k0==-1) & (qp_data.Jx[i][j][k]!=0):
						k0 = k
						iidxb[i][nbu[i]+j] = nu[i]+k
			iidxb[i] = np.ascontiguousarray(iidxb[i], dtype=np.int32)

		llb = []
		for i in range(N+1):
			llb.append(np.zeros((nbx[i]+nbu[i], 1), dtype=np.float64))
			for j in range(nbu[i]):
				llb[i][j] = qp_data.lu[i][j]
			for j in range(nbx[i]):
				llb[i][nbu[i]+j] = qp_data.lx[i][j]
			llb[i] = np.ascontiguousarray(llb[i], dtype=np.float64)

		uub = []
		for i in range(N+1):
			uub.append(np.zeros((nbx[i]+nbu[i], 1), dtype=np.float64))
			for j in range(nbu[i]):
				uub[i][j] = qp_data.uu[i][j]
			for j in range(nbx[i]):
				uub[i][nbu[i]+j] = qp_data.ux[i][j]
			uub[i] = np.ascontiguousarray(uub[i], dtype=np.float64)


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
		d_ls = (POINTER(c_double)*(N+1))()
		d_us = (POINTER(c_double)*(N+1))()


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
				lb[i] = cast(llb[i].ctypes.data, POINTER(c_double))
				ub[i] = cast(uub[i].ctypes.data, POINTER(c_double))
				idxb[i] = cast(iidxb[i].ctypes.data, POINTER(c_int))

			# polytopic constraints
			if qp_dims.ng[i] > 0:
				C[i] = cast(qp_data.C[i].ctypes.data, POINTER(c_double))
				D[i] = cast(qp_data.D[i].ctypes.data, POINTER(c_double))
				lg[i] = cast(qp_data.lg[i].ctypes.data, POINTER(c_double))
				ug[i] = cast(qp_data.ug[i].ctypes.data, POINTER(c_double))

			# slacks
			if qp_dims.ns[i] > 0:
				qp_data.Zl[i] = np.ascontiguousarray(qp_data.Zl[i], dtype=np.float64)
				Zl[i] = cast(qp_data.Zl[i].ctypes.data, POINTER(c_double))
				qp_data.Zu[i] = np.ascontiguousarray(qp_data.Zu[i], dtype=np.float64)
				Zu[i] = cast(qp_data.Zu[i].ctypes.data, POINTER(c_double))

				qp_data.zl[i] = np.ascontiguousarray(qp_data.zl[i], dtype=np.float64)
				zl[i] = cast(qp_data.zl[i].ctypes.data, POINTER(c_double))
				qp_data.zu[i] = np.ascontiguousarray(qp_data.zu[i], dtype=np.float64)
				zu[i] = cast(qp_data.zu[i].ctypes.data, POINTER(c_double))

				qp_data.d_ls[i] = np.ascontiguousarray(qp_data.d_ls[i], dtype=np.float64)
				d_ls[i] = cast(qp_data.d_ls[i].ctypes.data, POINTER(c_double))
				qp_data.d_us[i] = np.ascontiguousarray(qp_data.d_us[i], dtype=np.float64)
				d_us[i] = cast(qp_data.d_us[i].ctypes.data, POINTER(c_double))

				# slack indeces
				qp_data.idxs[i] = np.ascontiguousarray(qp_data.idxs[i], dtype=np.int32)
				idxs[i] = cast(qp_data.idxs[i].ctypes.data, POINTER(c_int))

		# allocate memory for qp struct 
		qp_size = __hpipm.d_memsize_ocp_qp(dim)
		qp_mem = cast(create_string_buffer(qp_size), c_void_p)
		self.qp_mem = qp_mem

		# set up ocp_qp structure
		sizeof_d_ocp_qp = __hpipm.d_sizeof_ocp_qp()
		qp = cast(create_string_buffer(sizeof_d_ocp_qp), c_void_p)
		self.ocp_qp = qp

		__hpipm.d_create_ocp_qp(dim, qp, qp_mem)
		__hpipm.d_cvt_colmaj_to_ocp_qp(A, B, b, Q, S, R, q, r, idxb, lb, ub,
			C, D, lg, ug, Zl, Zu, zl, zu, idxs, d_ls, d_us, qp)
		
		# allocate memory for ocp_qp_sol struct
		qp_sol_size = __hpipm.d_memsize_ocp_qp_sol(dim)
		qp_sol_mem = cast(create_string_buffer(qp_sol_size), c_void_p)
		self.qp_sol_mem = qp_sol_mem

		# set up ocp_qp_sol struct
		sizeof_d_ocp_qp_sol = __hpipm.d_sizeof_ocp_qp_sol()
		qp_sol = cast(create_string_buffer(sizeof_d_ocp_qp_sol), c_void_p)
		__hpipm.d_create_ocp_qp_sol(dim, qp_sol, qp_sol_mem)
		self.ocp_qp_sol = qp_sol

		# allocate memory for ipm_arg struct
		ipm_arg_size = __hpipm.d_memsize_ocp_qp_ipm_arg(dim)
		ipm_arg_mem = cast(create_string_buffer(ipm_arg_size), c_void_p)
		self.ipm_arg_mem = ipm_arg_mem

		# set up ipm_arg
		sizeof_d_ocp_qp_ipm_arg = __hpipm.d_sizeof_ocp_qp_ipm_arg()
		arg = cast(create_string_buffer(sizeof_d_ocp_qp_ipm_arg), c_void_p)
		self.ocp_qp_ipm_arg = arg

		__hpipm.d_create_ocp_qp_ipm_arg(dim, arg, ipm_arg_mem)
		__hpipm.d_set_default_ocp_qp_ipm_arg(1, arg)

		# allocate memory for ipm workspace 
		ipm_size = __hpipm.d_memsize_ocp_qp_ipm(dim, arg)
		ipm_mem = cast(create_string_buffer(ipm_size), c_void_p)
		self.ipm_mem = ipm_mem

		# set up ipm workspace
		sizeof_d_ocp_qp_ipm_workspace = __hpipm.d_sizeof_ocp_qp_ipm_workspace()
		workspace = cast(create_string_buffer(sizeof_d_ocp_qp_ipm_workspace), c_void_p)
		self.ocp_qp_ipm_workspace = workspace

		__hpipm.d_create_ocp_qp_ipm(dim, arg, workspace, ipm_mem)

		self.qp = qp
		self.qp_sol = qp_sol
		self.arg = arg
		self.workspace = workspace
		self.dim = dim
		
		self.__hpipm = __hpipm
		self.__blasfeo = __blasfeo

	def solve(self):
		return self.__hpipm.d_solve_ocp_qp_ipm(self.qp, self.qp_sol, 
			self.arg, self.workspace)

	def print_sol(self):
		self.__hpipm.d_print_ocp_qp_sol(self.ocp_qp_sol, self.ocp_qp_dim)
		return 



def hpipm_solve(qp_dims, qp_data):
	# load blasfeo and hpipm shared libraries
	__blasfeo = CDLL('libblasfeo.so')
	__hpipm   = CDLL('libhpipm.so')
	
	# cast dimensions to contiguous int
	# TODO(andrea): int32 might not be portable to Windows
	nx   = np.ascontiguousarray(qp_dims.nx,  dtype=np.int32)
	nu   = np.ascontiguousarray(qp_dims.nu,  dtype=np.int32)
	nbx  = np.ascontiguousarray(qp_dims.nbx, dtype=np.int32)
	nbu  = np.ascontiguousarray(qp_dims.nbu, dtype=np.int32)
	ng   = np.ascontiguousarray(qp_dims.ng,  dtype=np.int32)
	ns   = np.ascontiguousarray(qp_dims.ns,  dtype=np.int32)
	N	= qp_dims.N

	# allocate memory for dimemsions struct
	sizeof_d_ocp_qp_dim = __hpipm.d_sizeof_ocp_qp_dim()
	dim = cast(create_string_buffer(sizeof_d_ocp_qp_dim), c_void_p)
	ocp_qp_dim = dim

	dim_size = __hpipm.d_memsize_ocp_qp_dim(qp_dims.N)
	dim_mem = cast(create_string_buffer(dim_size), c_void_p)
	dim_mem = dim_mem

	# set up dimensions structure
	__hpipm.d_create_ocp_qp_dim(N, dim, dim_mem)
	__hpipm.d_cvt_int_to_ocp_qp_dim(N, 
		cast(nx.ctypes.data, POINTER(c_int)), 
		cast(nu.ctypes.data, POINTER(c_int)), 
		cast(nbx.ctypes.data, POINTER(c_int)), 
		cast(nbu.ctypes.data, POINTER(c_int)), 
		cast(ng.ctypes.data, POINTER(c_int)), 
		cast(ns.ctypes.data, POINTER(c_int)), 
		dim)

	A = (POINTER(c_double)*(N))()
	B = (POINTER(c_double)*(N))()
	b = (POINTER(c_double)*(N))()  

	Q = (POINTER(c_double)*(N+1))()	
	S = (POINTER(c_double)*(N+1))()   
	R = (POINTER(c_double)*(N+1))()   
	q = (POINTER(c_double)*(N+1))()   
	r = (POINTER(c_double)*(N+1))()   
   
	d_lb = (POINTER(c_double)*(N+1))()
	d_ub = (POINTER(c_double)*(N+1))()

	C = (POINTER(c_double)*(N+1))()   
	D = (POINTER(c_double)*(N+1))()   

	d_lg = (POINTER(c_double)*(N+1))()
	d_ug = (POINTER(c_double)*(N+1))()

	Zl = (POINTER(c_double)*(N+1))()  
	Zu = (POINTER(c_double)*(N+1))()  

	zl = (POINTER(c_double)*(N+1))()  
	zu = (POINTER(c_double)*(N+1))()  

	d_ls = (POINTER(c_double)*(N+1))()
	d_us = (POINTER(c_double)*(N+1))()

	idxb = (POINTER(c_int)*(N+1))()
	idxs = (POINTER(c_int)*(N+1))()

	x0 = (POINTER(c_double)*1)()  

	for i in range(N):
		# dynamics
		qp_data.A[i] = np.ascontiguousarray(qp_data.A[i], dtype= np.float64)
		A[i] = cast(qp_data.A[i].ctypes.data, POINTER(c_double))
		qp_data.B[i] = np.ascontiguousarray(qp_data.B[i], dtype=np.float64)
		B[i] = cast(qp_data.B[i].ctypes.data, POINTER(c_double))
		qp_data.b[i] = np.ascontiguousarray(qp_data.b[i], dtype=np.float64)
		b[i] = cast(qp_data.b[i].ctypes.data, POINTER(c_double))
		 
		# cost
		qp_data.Q[i] = np.ascontiguousarray(qp_data.Q[i], dtype=np.float64)
		Q[i] = cast(qp_data.Q[i].ctypes.data, POINTER(c_double))
		qp_data.S[i] = np.ascontiguousarray(qp_data.S[i], dtype=np.float64)
		S[i] = cast(qp_data.S[i].ctypes.data, POINTER(c_double))
		qp_data.R[i] = np.ascontiguousarray(qp_data.R[i], dtype=np.float64)
		R[i] = cast(qp_data.R[i].ctypes.data, POINTER(c_double))

		qp_data.q[i] = np.ascontiguousarray(qp_data.q[i], dtype=np.float64)
		q[i] = cast(qp_data.q[i].ctypes.data, POINTER(c_double))
		qp_data.r[i] = np.ascontiguousarray(qp_data.r[i], dtype=np.float64)
		r[i] = cast(qp_data.r[i].ctypes.data, POINTER(c_double))

		# simple bounds
		if qp_dims.nbx[i]+qp_dims.nbu[i] > 0:
			qp_data.d_lb[i] = np.ascontiguousarray(qp_data.d_lb[i], dtype=np.float64)
			d_lb[i] = cast(qp_data.d_lb[i].ctypes.data, POINTER(c_double))
			qp_data.d_ub[i] = np.ascontiguousarray(qp_data.d_ub[i], dtype=np.float64)
			d_ub[i] = cast(qp_data.d_ub[i].ctypes.data, POINTER(c_double))
		 
			qp_data.idxb[i] = np.ascontiguousarray(qp_data.idxb[i], dtype=np.int32)
			idxb[i] = cast(qp_data.idxb[i].ctypes.data, POINTER(c_int))

		# polytopic constraints
		if qp_dims.ng[i] > 0:
			qp_data.C[i] = np.ascontiguousarray(qp_data.C[i], dtype=np.float64)
			C[i] = cast(qp_data.C[i].ctypes.data, POINTER(c_double))
			qp_data.D[i] = np.ascontiguousarray(qp_data.D[i], dtype=np.float64)
			D[i] = cast(qp_data.D[i].ctypes.data, POINTER(c_double))

			qp_data.d_lg[i] = np.ascontiguousarray(qp_data.d_lg[i], dtype=np.float64)
			d_lg[i] = cast(qp_data.d_lg[i].ctypes.data, POINTER(c_double))
			qp_data.d_ug[i] = np.ascontiguousarray(qp_data.d_ug[i], dtype=np.float64)
			d_ug[i] = cast(qp_data.d_ug[i].ctypes.data, POINTER(c_double))
		 

		# slacks
		if qp_dims.ns[i] > 0:
			qp_data.Zl[i] = np.ascontiguousarray(qp_data.Zl[i], dtype=np.float64)
			Zl[i] = cast(qp_data.Zl[i].ctypes.data, POINTER(c_double))
			qp_data.Zu[i] = np.ascontiguousarray(qp_data.Zu[i], dtype=np.float64)
			Zu[i] = cast(qp_data.Zu[i].ctypes.data, POINTER(c_double))
		 
			qp_data.zl[i] = np.ascontiguousarray(qp_data.zl[i], dtype=np.float64)
			zl[i] = cast(qp_data.zl[i].ctypes.data, POINTER(c_double))
			qp_data.zu[i] = np.ascontiguousarray(qp_data.zu[i], dtype=np.float64)
			zu[i] = cast(qp_data.zu[i].ctypes.data, POINTER(c_double))
		 
			qp_data.d_ls[i] = np.ascontiguousarray(qp_data.d_ls[i], dtype=np.float64)
			d_ls[i] = cast(qp_data.d_ls[i].ctypes.data, POINTER(c_double))
			qp_data.d_us[i] = np.ascontiguousarray(qp_data.d_us[i], dtype=np.float64)
			d_us[i] = cast(qp_data.d_us[i].ctypes.data, POINTER(c_double))
		 
			# slack indeces
			qp_data.idxs[i] = np.ascontiguousarray(qp_data.idxs[i], dtype=np.int32)
			idxs[i] = cast(qp_data.idxs[i].ctypes.data, POINTER(c_int))
	
	i = N
	
	# cost
	qp_data.Q[i] = np.ascontiguousarray(qp_data.Q[i], dtype=np.float64)
	Q[i] = cast(qp_data.Q[i].ctypes.data, POINTER(c_double))
	qp_data.S[i] = np.ascontiguousarray(qp_data.S[i], dtype=np.float64)
	S[i] = cast(qp_data.S[i].ctypes.data, POINTER(c_double))
	qp_data.R[i] = np.ascontiguousarray(qp_data.R[i], dtype=np.float64)
	R[i] = cast(qp_data.R[i].ctypes.data, POINTER(c_double))

	qp_data.q[i] = np.ascontiguousarray(qp_data.q[i], dtype=np.float64)
	q[i] = cast(qp_data.q[i].ctypes.data, POINTER(c_double))
	qp_data.r[i] = np.ascontiguousarray(qp_data.r[i], dtype=np.float64)
	r[i] = cast(qp_data.r[i].ctypes.data, POINTER(c_double))

	# simple bounds
	if qp_dims.nbx[i]+qp_dims.nbu[i] > 0:
		qp_data.d_lb[i] = np.ascontiguousarray(qp_data.d_lb[i], dtype=np.float64)
		d_lb[i] = cast(qp_data.d_lb[i].ctypes.data, POINTER(c_double))
		qp_data.d_ub[i] = np.ascontiguousarray(qp_data.d_ub[i], dtype=np.float64)
		d_ub[i] = cast(qp_data.d_ub[i].ctypes.data, POINTER(c_double))
	 
		qp_data.idxb[i] = np.ascontiguousarray(qp_data.idxb[i], dtype=np.int32)
		idxb[i] = cast(qp_data.idxb[i].ctypes.data, POINTER(c_int))

	# polytopic constraints
	if qp_dims.ng[i] > 0:
		qp_data.C[i] = np.ascontiguousarray(qp_data.C[i], dtype=np.float64)
		C[i] = cast(qp_data.C[i].ctypes.data, POINTER(c_double))
		qp_data.D[i] = np.ascontiguousarray(qp_data.D[i], dtype=np.float64)
		D[i] = cast(qp_data.D[i].ctypes.data, POINTER(c_double))

		qp_data.d_lg[i] = np.ascontiguousarray(qp_data.d_lg[i], dtype=np.float64)
		d_lg[i] = cast(qp_data.d_lg[i].ctypes.data, POINTER(c_double))
		qp_data.d_ug[i] = np.ascontiguousarray(qp_data.d_ug[i], dtype=np.float64)
		d_ug[i] = cast(qp_data.d_ug[i].ctypes.data, POINTER(c_double))
	 

	# slacks
	if qp_dims.ns[i] > 0:
		qp_data.Zl[i] = np.ascontiguousarray(qp_data.Zl[i], dtype=np.float64)
		Zl[i] = cast(qp_data.Zl[i].ctypes.data, POINTER(c_double))
		qp_data.Zu[i] = np.ascontiguousarray(qp_data.Zu[i], dtype=np.float64)
		Zu[i] = cast(qp_data.Zu[i].ctypes.data, POINTER(c_double))
	 
		qp_data.zl[i] = np.ascontiguousarray(qp_data.zl[i], dtype=np.float64)
		zl[i] = cast(qp_data.zl[i].ctypes.data, POINTER(c_double))
		qp_data.zu[i] = np.ascontiguousarray(qp_data.zu[i], dtype=np.float64)
		zu[i] = cast(qp_data.zu[i].ctypes.data, POINTER(c_double))
	 
		qp_data.d_ls[i] = np.ascontiguousarray(qp_data.d_ls[i], dtype=np.float64)
		d_ls[i] = cast(qp_data.d_ls[i].ctypes.data, POINTER(c_double))
		qp_data.d_us[i] = np.ascontiguousarray(qp_data.d_us[i], dtype=np.float64)
		d_us[i] = cast(qp_data.d_us[i].ctypes.data, POINTER(c_double))
	 
		# slack indeces
		qp_data.idxs[i] = np.ascontiguousarray(qp_data.idxs[i], dtype=np.int32)
		idxs[i] = cast(qp_data.idxs[i].ctypes.data, POINTER(c_int))

	# allocate memory for qp struct 
	qp_size = __hpipm.d_memsize_ocp_qp(dim)
	qp_mem = cast(create_string_buffer(qp_size), c_void_p)
	qp_mem = qp_mem

	# set up ocp_qp structure
	sizeof_d_ocp_qp = __hpipm.d_sizeof_ocp_qp()
	qp = cast(create_string_buffer(sizeof_d_ocp_qp), c_void_p)
	ocp_qp = qp

	__hpipm.d_create_ocp_qp(dim, qp, qp_mem)
	__hpipm.d_cvt_colmaj_to_ocp_qp(A, B, b, Q, S, R, q, r, idxb, d_lb, 
		d_ub, C, D, d_lg, d_ug, Zl, Zu, zl, zu, idxs, d_ls, d_us, qp)
	
	# allocate memory for ocp_qp_sol struct
	qp_sol_size = __hpipm.d_memsize_ocp_qp_sol(dim)
	qp_sol_mem = cast(create_string_buffer(qp_sol_size), c_void_p)
	qp_sol_mem = qp_sol_mem

	# set up ocp_qp_sol struct
	sizeof_d_ocp_qp_sol = __hpipm.d_sizeof_ocp_qp_sol()
	qp_sol = cast(create_string_buffer(sizeof_d_ocp_qp_sol), c_void_p)
	__hpipm.d_create_ocp_qp_sol(dim, qp_sol, qp_sol_mem)
	ocp_qp_sol = qp_sol

	# allocate memory for ipm_arg struct
	ipm_arg_size = __hpipm.d_memsize_ocp_qp_ipm_arg(dim)
	ipm_arg_mem = cast(create_string_buffer(ipm_arg_size), c_void_p)
	ipm_arg_mem = ipm_arg_mem

	# set up ipm_arg
	sizeof_d_ocp_qp_ipm_arg = __hpipm.d_sizeof_ocp_qp_ipm_arg()
	arg = cast(create_string_buffer(sizeof_d_ocp_qp_ipm_arg), c_void_p)
	ocp_qp_ipm_arg = arg

	__hpipm.d_create_ocp_qp_ipm_arg(dim, arg, ipm_arg_mem)
	__hpipm.d_set_default_ocp_qp_ipm_arg(1, arg)

	# allocate memory for ipm workspace 
	ipm_size = __hpipm.d_memsize_ocp_qp_ipm(dim, arg)
	ipm_mem = cast(create_string_buffer(ipm_size), c_void_p)
	ipm_mem = ipm_mem

	# set up ipm workspace
	sizeof_d_ocp_qp_ipm_workspace = __hpipm.d_sizeof_ocp_qp_ipm_workspace()
	workspace = cast(create_string_buffer(sizeof_d_ocp_qp_ipm_workspace), c_void_p)
	ocp_qp_ipm_workspace = workspace

	__hpipm.d_create_ocp_qp_ipm(dim, arg, workspace, ipm_mem)

	return __hpipm.d_solve_ocp_qp_ipm(qp, qp_sol, 
		arg, workspace)



class hpipm_ocp_qp_dims:
	def __init__(self, N):
		self.N	= N
		self.nx   = np.zeros(N+1, dtype=int)
		self.nu   = np.zeros(N+1, dtype=int)
		self.nbx  = np.zeros(N+1, dtype=int)
		self.nbu  = np.zeros(N+1, dtype=int)
		self.ng   = np.zeros(N+1, dtype=int)
		self.ns   = np.zeros(N+1, dtype=int)

	def set_nx(self, nx, idx=None):
		if idx==None:
			for i in range(nx.size):
				self.nx[i] = nx[i]
		else:
			self.nx[idx] = nx
		return

	def set_nu(self, nu, idx=None):
		if idx==None:
			for i in range(nu.size):
				self.nu[i] = nu[i]
		else:
			self.nu[idx] = nu
		return

	def set_nbx(self, nbx, idx=None):
		if idx==None:
			for i in range(nbx.size):
				self.nbx[i] = nbx[i]
		else:
			self.nbx[idx] = nbx
		return

	def set_nbu(self, nbu, idx=None):
		if idx==None:
			for i in range(nbu.size):
				self.nbu[i] = nbu[i]
		else:
			self.nbu[idx] = nbu
		return

	def set_ng(self, ng, idx=None):
		if idx==None:
			for i in range(ng.size):
				self.ng[i] = ng[i]
		else:
			self.ng[idx] = ng
		return

	def set_ns(self, ns, idx=None):
		if idx==None:
			for i in range(ns.size):
				self.ns[i] = ns[i]
		else:
			self.ns[idx] = ns
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

		self.dims = hpipm_ocp_qp_dims(N)
		self.dims.set_nx(nx);
		self.dims.set_nu(nu);
		self.dims.set_nbx(nbx);
		self.dims.set_nbu(nbu);
		self.dims.set_ng(ng);
		self.dims.set_ns(ns);

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


		# old interface

		self.Zl   = None
		self.Zu   = None

		self.zl   = None
		self.zu   = None

		self.d_ls = None
		self.d_us = None

		self.idxs = None

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
		nbx = self.dims.nbx
		if idx==None:
			for i in range(len(Jx)):
				self.Jx[i] = Jx[i]
				self.Jx[i] = np.ascontiguousarray(self.Jx[i], dtype=np.float64)
				self.Jx[i] = self.Jx[i].reshape((nbx[i], nx[i]))
		else:
			self.Jx[idx] = Jx
			self.Jx[idx] = np.ascontiguousarray(self.Jx[idx], dtype=np.float64)
			self.Jx[idx] = self.Jx[idx].reshape((nbx[idx], nx[idx]))
		return

	def set_lx(self, lx, idx=None):
		nbx = self.dims.nbx
		if idx==None:
			for i in range(len(lx)):
				self.lx[i] = lx[i]
				self.lx[i] = np.ascontiguousarray(self.lx[i], dtype=np.float64)
				self.lx[i] = self.lx[i].reshape((nbx[i], 1))
		else:
			self.lx[idx] = lx
			self.lx[idx] = np.ascontiguousarray(self.lx[idx], dtype=np.float64)
			self.lx[idx] = self.lx[idx].reshape((nbx[idx], 1))
		return

	def set_ux(self, ux, idx=None):
		nbx = self.dims.nbx
		if idx==None:
			for i in range(len(ux)):
				self.ux[i] = ux[i]
				self.ux[i] = np.ascontiguousarray(self.ux[i], dtype=np.float64)
				self.ux[i] = self.ux[i].reshape((nbx[i], 1))
		else:
			self.ux[idx] = ux
			self.ux[idx] = np.ascontiguousarray(self.ux[idx], dtype=np.float64)
			self.ux[idx] = self.ux[idx].reshape((nbx[idx], 1))
		return

	def set_Ju(self, Ju, idx=None):
		nu = self.dims.nu
		nbu = self.dims.nbu
		if idx==None:
			for i in range(len(Ju)):
				self.Ju[i] = Ju[i]
				self.Ju[i] = np.ascontiguousarray(self.Ju[i], dtype=np.float64)
				self.Ju[i] = self.Ju[i].reshape((nbu[i], nu[i]))
		else:
			self.Ju[idx] = Ju
			self.Ju[idx] = np.ascontiguousarray(self.Ju[idx], dtype=np.float64)
			self.Ju[idx] = self.Ju[idx].reshape((nbu[idx], nu[idx]))
		return

	def set_lu(self, lu, idx=None):
		nbu = self.dims.nbu
		if idx==None:
			for i in range(len(lu)):
				self.lu[i] = lu[i]
				self.lu[i] = np.ascontiguousarray(self.lu[i], dtype=np.float64)
				self.lu[i] = self.lu[i].reshape((nbu[i], 1))
		else:
			self.lu[idx] = lu
			self.lu[idx] = np.ascontiguousarray(self.lu[idx], dtype=np.float64)
			self.lu[idx] = self.lu[idx].reshape((nbu[idx], 1))
		return

	def set_uu(self, uu, idx=None):
		nbu = self.dims.nbu
		if idx==None:
			for i in range(len(uu)):
				self.uu[i] = uu[i]
				self.uu[i] = np.ascontiguousarray(self.uu[i], dtype=np.float64)
				self.uu[i] = self.uu[i].reshape((nbu[i], 1))
		else:
			self.uu[idx] = uu
			self.uu[idx] = np.ascontiguousarray(self.uu[idx], dtype=np.float64)
			self.uu[idx] = self.uu[idx].reshape((nbu[idx], 1))
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



class hpipm_ocp_qp_sol:
	def __init__(self):
		self.ux = None
		self.pi = None
		self.lam = None
		self.t = None

