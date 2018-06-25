
from ctypes import *
from os import *
from hpipm_python.external.blasfeo.wrapper import *

# d_ocp_qp - include/hpipm_d_ocp_qp_ipm_dim.h
class d_core_qp_ipm_workspace(Structure):
    _fields_ = [
	("v",           POINTER(c_double)), # primal variables
	("pi",          POINTER(c_double)), # equality constraints multipliers
	("lam",         POINTER(c_double)), # inequality constraints multipliers
	("t",           POINTER(c_double)), # inequality constraints slacks
	("t_inv",       POINTER(c_double)), # inverse of t
	("dv",          POINTER(c_double)), # step in v
	("dpi",         POINTER(c_double)), # step in pi
	("dlam",        POINTER(c_double)), # step in lam
	("dt",          POINTER(c_double)), # step in t
	("res_g",       POINTER(c_double)), # q-residuals
	("res_b",       POINTER(c_double)), # b-residuals
        ("res_d",       POINTER(c_double)), # d-residuals
	("res_m",       POINTER(c_double)), # m-residuals
	("res_m_bkp",   POINTER(c_double)), # m-residuals
	("Gamma",       POINTER(c_double)), # Hessian update
	("gamma",       POINTER(c_double)), # gradient update
	("alpha",       c_double), # step length
	("alpha_prim",  c_double), # step length
	("alpha_dual",  c_double), # step length
	("sigma",       c_double), # centering XXX
	("mu",          c_double), # duality measuere
	("mu_aff",      c_double), # affine duality measuere
	("nc_inv",      c_double), # 1.0/nt, where nt is the total number of constraints
	("nv",          c_int),    # number of primal variables
	("ne",          c_int),    # number of equality constraints
	("nc",          c_int),    # number of (two-sided) inequality constraints
	("memsize",     c_int),    # memory size (in bytes) of workspace
        ]

