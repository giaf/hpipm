from hpipm_python.external.blasfeo.wrapper import *
from hpipm_python.hpipm.core.wrapper import *
from hpipm_python.hpipm.d_ocp_qp.wrapper import *

from ctypes import *
import numpy as np

class hpipm_solver:
    def __init__(self, N, _nx, _nu, _nbx, _nbu, _ng, _ns):
       
        
        nx  = _nx.astype(int)
        nu  = _nu.astype(int)
        nbx = _nbx.astype(int)
        nbu = _nbu.astype(int)
        ng  = _ng.astype(int)
        ns  = _ns.astype(int)

        __blasfeo = CDLL('libblasfeo.so')
        __hpipm   = CDLL('libhpipm.so')

        dim = d_ocp_qp_dim()
        dim_size = __hpipm.d_memsize_ocp_qp_dim(N)
        dim_mem = c_void_p()
        
        __blasfeo.v_zeros(byref(dim_mem), dim_size)
        __hpipm.d_create_ocp_qp_dim(N, byref(dim), dim_mem)
        __hpipm.d_cvt_int_to_ocp_qp_dim(N, c_void_p(nx.ctypes.data), c_void_p(nu.ctypes.data), 
                c_void_p(nbx.ctypes.data), c_void_p(nbu.ctypes.data), c_void_p(ng.ctypes.data), c_void_p(ns.ctypes.data), byref(dim))


