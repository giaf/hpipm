from hpipm_python.external.blasfeo.wrapper import *
from hpipm_python.hpipm.core.wrapper import *
from hpipm_python.hpipm.d_ocp_qp.wrapper import *

from ctypes import *

import numpy as np

class hpipm_solver:
    def __init__(self, dim, N):
        
        _h = CDLL('libhpipm.so')
        dim_size = _h.d_memsize_ocp_qp_dim(N)
        dim_mem = c_void_p()
        v_zeros(byref(dim_size), dim_mem)



