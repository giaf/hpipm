from ctypes import *
import ctypes.util 
import numpy as np

class hpipm_solver:
    def __init__(self, qp_dims, qp_data):
      
        # load blasfeo and hpipm shared libraries
        __blasfeo = CDLL('libblasfeo.so')
        __hpipm   = CDLL('libhpipm.so')
        
        # cast dimensions to contiguous int
        nx   = np.ascontiguousarray(qp_dims.nx,  dtype=np.int32)
        nu   = np.ascontiguousarray(qp_dims.nu,  dtype=np.int32)
        nb   = np.ascontiguousarray(qp_dims.nb,  dtype=np.int32)
        nbx  = np.ascontiguousarray(qp_dims.nbx, dtype=np.int32)
        nbu  = np.ascontiguousarray(qp_dims.nbu, dtype=np.int32)
        ng   = np.ascontiguousarray(qp_dims.ng,  dtype=np.int32)
        ns   = np.ascontiguousarray(qp_dims.ns,  dtype=np.int32)
        nsbx = np.ascontiguousarray(qp_dims.nsbx,  dtype=np.int32)
        nsbu = np.ascontiguousarray(qp_dims.nsbu,  dtype=np.int32)
        nsg  = np.ascontiguousarray(qp_dims.nsg,  dtype=np.int32)
        N    = qp_dims.N

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
            cast(nsbx.ctypes.data, POINTER(c_double)), 
            cast(nsbu.ctypes.data, POINTER(c_double)), 
            cast(nsg.ctypes.data, POINTER(c_double)), 
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
            if qp_dims.nb[i] > 0:
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
        if qp_dims.nb[i] > 0:
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
    nb   = np.ascontiguousarray(qp_dims.nb,  dtype=np.int32)
    nbx  = np.ascontiguousarray(qp_dims.nbx, dtype=np.int32)
    nbu  = np.ascontiguousarray(qp_dims.nbu, dtype=np.int32)
    ng   = np.ascontiguousarray(qp_dims.ng,  dtype=np.int32)
    ns   = np.ascontiguousarray(qp_dims.ns,  dtype=np.int32)
    nsbx = np.ascontiguousarray(qp_dims.nsbx,  dtype=np.int32)
    nsbu = np.ascontiguousarray(qp_dims.nsbu,  dtype=np.int32)
    nsg  = np.ascontiguousarray(qp_dims.nsg,  dtype=np.int32)
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
        cast(nsbx.ctypes.data, POINTER(c_int)), 
        cast(nsbu.ctypes.data, POINTER(c_int)), 
        cast(nsg.ctypes.data, POINTER(c_int)), 
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
        if qp_dims.nb[i] > 0:
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
    if qp_dims.nb[i] > 0:
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

class hpipm_dims:
    def __init__(self):
        self.nx   = None
        self.nu   = None
        self.nb   = None
        self.nbx  = None
        self.nbu  = None
        self.ng   = None
        self.ns   = None
        self.nsbx = None
        self.nsbu = None
        self.nsg  = None
        
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

class hpipm_sol:
    def __init__(self):
        self.ux = None
        self.pi = None
        self.lam = None
        self.t = None

