from hpipm_python.external.blasfeo.wrapper import *
from hpipm_python.hpipm.core.wrapper import *
from hpipm_python.hpipm.d_ocp_qp.wrapper import *

from ctypes import *

import numpy as np

class hpipm_solver:
    def __init__(self, dim, N):
        
        _blasfeo = CDLL('libblasfeo.so')
        _hpipm = CDLL('libhpipm.so')
        dim_size = _hpipm.d_memsize_ocp_qp_dim(N)
        dim_mem = c_void_p()
        _blasfeo.v_zeros(byref(dim_mem), dim_size)



