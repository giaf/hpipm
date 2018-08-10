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
        N    = qp_dims.N
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
        self.qp_mem = qp_mem

        # set up ocp_qp structure
        sizeof_d_ocp_qp = __hpipm.d_sizeof_ocp_qp()
        qp = cast(create_string_buffer(sizeof_d_ocp_qp), c_void_p)
        self.ocp_qp = qp

        __hpipm.d_create_ocp_qp(dim, qp, qp_mem)
        __hpipm.d_cvt_colmaj_to_ocp_qp(A, B, b, Q, S, R, q, r, idxb, d_lb, 
            d_ub, C, D, d_lg, d_ug, Zl, Zu, zl, zu, idxs, d_ls, d_us, qp)
        
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
    N    = qp_dims.N

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
		self.N    = N
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


		# old interface
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

	def set_A(self, A, idx=None):
		if idx==None:
			for i in range(len(A)):
				self.A[i] = A[i]
		else:
			self.A[idx] = A
		return

	def set_B(self, B, idx=None):
		if idx==None:
			for i in range(len(B)):
				self.B[i] = B[i]
		else:
			self.B[idx] = B
		return

	def set_b(self, b, idx=None):
		if idx==None:
			for i in range(len(b)):
				self.b[i] = b[i]
		else:
			self.b[idx] = b
		return

	def set_Q(self, Q, idx=None):
		if idx==None:
			for i in range(len(Q)):
				self.Q[i] = Q[i]
		else:
			self.Q[idx] = Q
		return

	def set_R(self, R, idx=None):
		if idx==None:
			for i in range(len(R)):
				self.R[i] = R[i]
		else:
			self.R[idx] = R
		return

	def set_S(self, S, idx=None):
		if idx==None:
			for i in range(len(S)):
				self.S[i] = S[i]
		else:
			self.S[idx] = S
		return

	def set_q(self, q, idx=None):
		if idx==None:
			for i in range(len(q)):
				self.q[i] = q[i]
		else:
			self.q[idx] = q
		return

	def set_r(self, r, idx=None):
		if idx==None:
			for i in range(len(r)):
				self.r[i] = r[i]
		else:
			self.r[idx] = r
		return

	def set_Jx(self, Jx, idx=None):
		if idx==None:
			for i in range(len(Jx)):
				self.Jx[i] = Jx[i]
		else:
			self.Jx[idx] = Jx
		return

	def set_lx(self, lx, idx=None):
		if idx==None:
			for i in range(len(lx)):
				self.lx[i] = lx[i]
		else:
			self.lx[idx] = lx
		return

	def set_ux(self, ux, idx=None):
		if idx==None:
			for i in range(len(ux)):
				self.ux[i] = ux[i]
		else:
			self.ux[idx] = ux
		return

	def set_Ju(self, Ju, idx=None):
		if idx==None:
			for i in range(len(Ju)):
				self.Ju[i] = Ju[i]
		else:
			self.Ju[idx] = Ju
		return

	def set_lu(self, lu, idx=None):
		if idx==None:
			for i in range(len(lu)):
				self.lu[i] = lu[i]
		else:
			self.lu[idx] = lu
		return

	def set_uu(self, uu, idx=None):
		if idx==None:
			for i in range(len(uu)):
				self.uu[i] = uu[i]
		else:
			self.uu[idx] = uu
		return



class hpipm_ocp_qp_sol:
    def __init__(self):
        self.ux = None
        self.pi = None
        self.lam = None
        self.t = None

