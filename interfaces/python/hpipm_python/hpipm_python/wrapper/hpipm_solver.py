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

        # cast dimensions to int
        nx  = qp_dims.nx.astype(int)
        nu  = qp_dims.nu.astype(int)
        nbx = qp_dims.nbx.astype(int)
        nbu = qp_dims.nbu.astype(int)
        ng  = qp_dims.ng.astype(int)
        ns  = qp_dims.ns.astype(int)

        # allocate memory for dimemsions struct
        dim = d_ocp_qp_dim()
        dim_size = __hpipm.d_memsize_ocp_qp_dim(qp_dims.N)
        dim_mem = c_void_p()
        __blasfeo.v_zeros(byref(dim_mem), dim_size)
        self.dim_mem = dim_mem

        # set up dimension structure
        __hpipm.d_create_ocp_qp_dim(qp_dims.N, byref(dim), dim_mem)
        __hpipm.d_cvt_int_to_ocp_qp_dim(qp_dims.N, c_void_p(nx.ctypes.data), c_void_p(nu.ctypes.data), 
                c_void_p(nbx.ctypes.data), c_void_p(nbu.ctypes.data), c_void_p(ng.ctypes.data), c_void_p(ns.ctypes.data), byref(dim))
        
        # cast data to double
        nx  = qp_data.nx.astype(float64)
        nu  = qp_data.nu.astype(float64)
        nbx = qp_data.nbx.astype(float64)
        nbu = qp_data.nbu.astype(float64)
        ng  = qp_data.ng.astype(float64)
        ns  = qp_data.ns.astype(float64)

        # allocate memory for qp struct 
	qp_size = d_memsize_ocp_qp(byref(dim));
	qp_mem = c_void_p()
        __blasfeo.v_zeros(byref(qp_mem), qp_size);
        self.qp_mem = qp_mem

        # setup ocp_qp structure
        qp = __hpipm.d_ocp_qp()
        d_create_ocp_qp(byref(dim), byref(qp), qp_mem);
        d_cvt_colmaj_to_ocp_qp(A, B, b, Q, S, R, q, r, idxb, d_lb, d_ub, C, D, d_lg, d_ug, Zl, Zu, zl, zu, idxs, d_ls, d_us, &qp)

class hpipm_dims:
    def __init__(self):
        self.nx   = None
        self.nu   = None
        self.nb   = None
        self.ng   = None
        self.ns   = None
        self.nbx  = None
        self.nbu  = None
        self.
        self.N    = None

        self.idbx = None
        self.idbs = None

class hpipm_data:
    def __init__(self):
        self.A  = None
        self.B  = None
        self.b  = None

        self.Q  = None
        self.S  = None
        self.R  = None
        self.q  = None
        self.r  = None
       
        self.d_lb = None
        self.d_ub = None

        self.C = None
        self.D = None

        self.d_lg = None
        self.d_ug = None

        self.Zl = None
        self.Zu = None

        self.zl = None
        self.zu = None

        self.d_ls = None
        self.d_us = None

        self.x0 = None
