from hpipm_python.external.blasfeo.wrapper import *
from hpipm_python.hpipm.core.wrapper import *
from hpipm_python.hpipm.d_ocp_qp.wrapper import *

from ctypes import *
import numpy as np

class hpipm_solver:
    def __init__(self, qp_dims, qp_data):
      
        # load blasfeo and hpipm shared libraries
        __blasfeo = CDLL('libblasfeo.so')
        __hpipm   = CDLL('libhpipm.so')

        # cast dimensions to contiguous int
        nx  = np.ascontiguousarray(qp_dims.nx, dtype=np.int64)
        nu  = np.ascontiguousarray(qp_dims.nu, dtype=np.int64)
        nbx = np.ascontiguousarray(qp_dims.nbx, dtype=np.int64)
        nbu = np.ascontiguousarray(qp_dims.nbu, dtype=np.int64)
        ng  = np.ascontiguousarray(qp_dims.ng, dtype=np.int64)
        ns  = np.ascontiguousarray(qp_dims.ns, dtype=np.int64)

        # allocate memory for dimemsions struct
        dim = d_ocp_qp_dim()
        dim_size = __hpipm.d_memsize_ocp_qp_dim(qp_dims.N)
        dim_mem = c_void_p()
        __blasfeo.v_zeros(byref(dim_mem), dim_size)
        self.dim_mem = dim_mem

        # set up dimensions structure
        __hpipm.d_create_ocp_qp_dim(qp_dims.N, byref(dim), dim_mem)
        __hpipm.d_cvt_int_to_ocp_qp_dim(qp_dims.N, c_void_p(nx.ctypes.data), c_void_p(nu.ctypes.data), 
                c_void_p(nbx.ctypes.data), c_void_p(nbu.ctypes.data), c_void_p(ng.ctypes.data), c_void_p(ns.ctypes.data), byref(dim))
        
        N = qp_dims.N

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
            A[i] = POINTER(c_double)(c_void_p(qp_data.A[i].ctypes.data))
            qp_data.B[i] = np.ascontiguousarray(qp_data.B[i], dtype=np.float64)
            B[i] = POINTER(c_double)(c_void_p(qp_data.B[i].ctypes.data)) 
            qp_data.b[i] = np.ascontiguousarray(qp_data.b[i], dtype=np.float64)
            b[i] = POINTER(c_double)(c_void_p(qp_data.b[i].ctypes.data)) 
             
            # cost
            qp_data.Q[i] = np.ascontiguousarray(qp_data.Q[i], dtype=np.float64)
            Q[i] = POINTER(c_double)(c_void_p(qp_data.Q[i].ctypes.data)) 
            qp_data.S[i] = np.ascontiguousarray(qp_data.S[i], dtype=np.float64)
            S[i] = POINTER(c_double)(c_void_p(qp_data.S[i].ctypes.data)) 
            qp_data.R[i] = np.ascontiguousarray(qp_data.R[i], dtype=np.float64)
            R[i] = POINTER(c_double)(c_void_p(qp_data.R[i].ctypes.data)) 

            qp_data.q[i] = np.ascontiguousarray(qp_data.q[i], dtype=np.float64)
            q[i] = POINTER(c_double)(c_void_p(qp_data.q[i].ctypes.data)) 
            qp_data.r[i] = np.ascontiguousarray(qp_data.r[i], dtype=np.float64)
            r[i] = POINTER(c_double)(c_void_p(qp_data.r[i].ctypes.data)) 

            # simple bounds
            if qp_dims.nb[i] > 0:
                qp_data.d_lb[i] = np.ascontiguousarray(qp_data.d_lb[i], dtype=np.float64)
                d_lb[i] = POINTER(c_double)(c_void_p(qp_data.d_lb[i].ctypes.data)) 
                qp_data.d_ub[i] = np.ascontiguousarray(qp_data.d_ub[i], dtype=np.float64)
                d_ub[i] = POINTER(c_double)(c_void_p(qp_data.d_ub[i].ctypes.data)) 
             
                qp_data.idxb[i] = np.ascontiguousarray(qp_data.idxb[i], dtype=np.int64)
                idxb[i] = POINTER(c_int)(c_void_p(qp_data.idxb[i].ctypes.data)) 

            # polytopic constraints
            if qp_dims.ng[i] > 0:
                qp_data.C[i] = np.ascontiguousarray(qp_data.C[i], dtype=np.float64)
                C[i] = POINTER(c_double)(c_void_p(qp_data.C[i].ctypes.data)) 
                qp_data.D[i] = np.ascontiguousarray(qp_data.D[i], dtype=np.float64)
                D[i] = POINTER(c_double)(c_void_p(qp_data.D[i].ctypes.data)) 

                qp_data.d_lg[i] = np.ascontiguousarray(qp_data.d_lg[i], dtype=np.float64)
                d_lg[i] = POINTER(c_double)(c_void_p(qp_data.d_lg[i].ctypes.data)) 
                qp_data.d_ug[i] = np.ascontiguousarray(qp_data.d_ug[i], dtype=np.float64)
                d_ug[i] = POINTER(c_double)(c_void_p(qp_data.d_ug[i].ctypes.data)) 
             

            # slacks
            if qp_dims.ns[i] > 0:
                qp_data.Zl[i] = np.ascontiguousarray(qp_data.Zl[i], dtype=np.float64)
                Zl[i] = POINTER(c_double)(c_void_p(qp_data.Zl[i].ctypes.data)) 
                qp_data.Zu[i] = np.ascontiguousarray(qp_data.Zu[i], dtype=np.float64)
                Zu[i] = POINTER(c_double)(c_void_p(qp_data.Zu[i].ctypes.data)) 
             
                qp_data.zl[i] = np.ascontiguousarray(qp_data.zl[i], dtype=np.float64)
                zl[i] = POINTER(c_double)(c_void_p(qp_data.zl[i].ctypes.data)) 
                qp_data.zu[i] = np.ascontiguousarray(qp_data.zu[i], dtype=np.float64)
                zu[i] = POINTER(c_double)(c_void_p(qp_data.zu[i].ctypes.data)) 
             
                qp_data.d_ls[i] = np.ascontiguousarray(qp_data.d_ls[i], dtype=np.float64)
                d_ls[i] = POINTER(c_double)(c_void_p(qp_data.d_ls[i].ctypes.data)) 
             
                qp_data.d_us[i] = np.ascontiguousarray(qp_data.d_us[i], dtype=np.float64)
                d_us[i] = POINTER(c_double)(c_void_p(qp_data.d_us[i].ctypes.data)) 
             
                # slack indeces
                qp_data.idxs[i] = np.ascontiguousarray(qp_data.idxs[i], dtype=np.int64)
                idxs[i] = POINTER(c_int)(c_void_p(qp_data.idxs[i].ctypes.data)) 
        
        i = N
        
        # cost
        qp_data.Q[i] = np.ascontiguousarray(qp_data.Q[i], dtype=np.float64)
        Q[i] = POINTER(c_double)(c_void_p(qp_data.Q[i].ctypes.data)) 
        qp_data.S[i] = np.ascontiguousarray(qp_data.S[i], dtype=np.float64)
        S[i] = POINTER(c_double)(c_void_p(qp_data.S[i].ctypes.data)) 
        qp_data.R[i] = np.ascontiguousarray(qp_data.R[i], dtype=np.float64)
        R[i] = POINTER(c_double)(c_void_p(qp_data.R[i].ctypes.data)) 

        qp_data.q[i] = np.ascontiguousarray(qp_data.q[i], dtype=np.float64)
        q[i] = POINTER(c_double)(c_void_p(qp_data.q[i].ctypes.data)) 
        qp_data.r[i] = np.ascontiguousarray(qp_data.r[i], dtype=np.float64)
        r[i] = POINTER(c_double)(c_void_p(qp_data.r[i].ctypes.data)) 

        # simple bounds
        if qp_dims.nb[i] > 0:
            qp_data.d_lb[i] = np.ascontiguousarray(qp_data.d_lb[i], dtype=np.float64)
            d_lb[i] = POINTER(c_double)(c_void_p(qp_data.d_lb[i].ctypes.data)) 
            qp_data.d_ub[i] = np.ascontiguousarray(qp_data.d_ub[i], dtype=np.float64)
            d_ub[i] = POINTER(c_double)(c_void_p(qp_data.d_ub[i].ctypes.data)) 
         
            qp_data.idxb[i] = np.ascontiguousarray(qp_data.idxb[i], dtype=np.int64)
            idxb[i] = POINTER(c_int)(c_void_p(qp_data.idxb[i].ctypes.data)) 

        # polytopic constraints
        if qp_dims.ng[i] > 0:
            qp_data.C[i] = np.ascontiguousarray(qp_data.C[i], dtype=np.float64)
            C[i] = POINTER(c_double)(c_void_p(qp_data.C[i].ctypes.data)) 
            qp_data.D[i] = np.ascontiguousarray(qp_data.D[i], dtype=np.float64)
            D[i] = POINTER(c_double)(c_void_p(qp_data.D[i].ctypes.data)) 

            qp_data.d_lg[i] = np.ascontiguousarray(qp_data.d_lg[i], dtype=np.float64)
            d_lg[i] = POINTER(c_double)(c_void_p(qp_data.d_lg[i].ctypes.data)) 
            qp_data.d_ug[i] = np.ascontiguousarray(qp_data.d_ug[i], dtype=np.float64)
            d_ug[i] = POINTER(c_double)(c_void_p(qp_data.d_ug[i].ctypes.data)) 
         

        # slacks
        if qp_dims.ns[i] > 0:
            qp_data.Zl[i] = np.ascontiguousarray(qp_data.Zl[i], dtype=np.float64)
            Zl[i] = POINTER(c_double)(c_void_p(qp_data.Zl[i].ctypes.data)) 
            qp_data.Zu[i] = np.ascontiguousarray(qp_data.Zu[i], dtype=np.float64)
            Zu[i] = POINTER(c_double)(c_void_p(qp_data.Zu[i].ctypes.data)) 
         
            qp_data.zl[i] = np.ascontiguousarray(qp_data.zl[i], dtype=np.float64)
            zl[i] = POINTER(c_double)(c_void_p(qp_data.zl[i].ctypes.data)) 
            qp_data.zu[i] = np.ascontiguousarray(qp_data.zu[i], dtype=np.float64)
            zu[i] = POINTER(c_double)(c_void_p(qp_data.zu[i].ctypes.data)) 
         
            qp_data.d_ls[i] = np.ascontiguousarray(qp_data.d_ls[i], dtype=np.float64)
            d_ls[i] = POINTER(c_double)(c_void_p(qp_data.d_ls[i].ctypes.data)) 
         
            qp_data.d_us[i] = np.ascontiguousarray(qp_data.d_us[i], dtype=np.float64)
            d_us[i] = POINTER(c_double)(c_void_p(qp_data.d_us[i].ctypes.data)) 
         
            # slack indeces
            qp_data.idxs[i] = np.ascontiguousarray(qp_data.idxs[i], dtype=np.int64)
            idxs[i] = POINTER(c_int)(c_void_p(qp_data.idxs[i].ctypes.data)) 

        # allocate memory for qp struct 
        qp_size = __hpipm.d_memsize_ocp_qp(byref(dim))
        qp_mem = c_void_p()
        __blasfeo.v_zeros(byref(qp_mem), qp_size)
        self.qp_mem = qp_mem

        # set up ocp_qp structure
        qp = d_ocp_qp()
        __hpipm.d_create_ocp_qp(byref(dim), byref(qp), qp_mem)
        __hpipm.d_cvt_rowmaj_to_ocp_qp(A, B, b, Q, S, R, q, r, idxb, d_lb, d_ub, C, D, d_lg, d_ug, Zl, Zu, zl, zu, idxs, d_ls, d_us, byref(qp))
        
        # allocate memory for ocp_qp_sol struct
        qp_sol_size = __hpipm.d_memsize_ocp_qp_sol(byref(dim))
        qp_sol_mem = c_void_p()
        __blasfeo.v_zeros(byref(qp_sol_mem), qp_sol_size)
        self.qp_sol_mem = qp_sol_mem

        # set up ocp_qp_sol struct
        qp_sol = d_ocp_qp_sol() 
        __hpipm.d_create_ocp_qp_sol(byref(dim), byref(qp_sol), qp_sol_mem)

        # allocate memory for ipm_arg struct
        ipm_arg_size = __hpipm.d_memsize_ocp_qp_ipm_arg(byref(dim))
        ipm_arg_mem = c_void_p()
        __blasfeo.v_zeros(byref(ipm_arg_mem), ipm_arg_size)
        self.ipm_arg_mem = ipm_arg_mem
    
        # set up ipm_arg
        arg = d_ocp_qp_ipm_arg()
        __hpipm.d_create_ocp_qp_ipm_arg(byref(dim), byref(arg), ipm_arg_mem)
        __hpipm.d_set_default_ocp_qp_ipm_arg(byref(arg))

        # allocate memory for ipm workspace 
        ipm_size = __hpipm.d_memsize_ocp_qp_ipm(byref(dim), byref(arg))
        ipm_mem = c_void_p()
        __blasfeo.v_zeros(byref(ipm_mem), ipm_size)
        self.ipm_mem = ipm_mem

        # set up ipm workspace
        workspace = d_ocp_qp_ipm_workspace()
        __hpipm.d_create_ocp_qp_ipm(byref(dim), byref(arg), byref(workspace), ipm_mem)

        self.qp = qp
        self.qp_sol = qp_sol
        self.arg = arg
        self.workspace = workspace

        self.__hpipm = __hpipm
        self.__blasfeo = __blasfeo

    def solve(self):
        return self.__hpipm.d_solve_ocp_qp_ipm(byref(self.qp), byref(self.qp_sol), byref(self.arg), byref(self.workspace))

class hpipm_dims:
    def __init__(self):
        self.nx   = None
        self.nu   = None
        self.nb   = None
        self.ng   = None
        self.ns   = None
        self.nbx  = None
        self.nbu  = None
        
        self.N    = None


class hpipm_data:
    def __init__(self):
        self.A    = None
        self.B    = None
        self.b    = None

        self.Q    = None
        self.S    = None
        self.R    = None
        self.q    = None
        self.r    = None
       
        self.d_lb = None
        self.d_ub = None

        self.C    = None
        self.D    = None

        self.d_lg = None
        self.d_ug = None

        self.idxb = None

        self.Zl   = None
        self.Zu   = None

        self.zl   = None
        self.zu   = None

        self.d_ls = None
        self.d_us = None

        self.idxs = None

        self.x0   = None
