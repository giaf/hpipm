
from ctypes import *
from os import *
from hpipm_python.external.blasfeo.wrapper import *

# d_ocp_qp - include/hpipm_d_ocp_qp_ipm_dim.h
class d_ocp_qp_dim(Structure):
    _fields_ = [
        ("nx",      POINTER(c_int)),            # number of states
        ("nu",      POINTER(c_int)),            # number of inputs
        ("nb",      POINTER(c_int)),            # number of box constraints
        ("nbx",     POINTER(c_int)),            # number of state box constraints
        ("nbu",     POINTER(c_int)),            # number of input box constraints
        ("ng",      POINTER(c_int)),            # number of general constraints
        ("ns",      POINTER(c_int)),            # number of soft constraints
        ("nsbx",    POINTER(c_int)),            # number of soft state box constraints
        ("nsbu",    POINTER(c_int)),            # number of soft input box constraints
        ("nsg",     POINTER(c_int)),            # number of soft general constraints
        ("N",       c_int),                     # horizon length
        ("memsize", c_int)                      # memory size
        ]

class d_core_qp_ipm_workspace(Structure):
    _fields_ = [
	double *v; # primal variables
	double *pi; # equality constraints multipliers
	double *lam; # inequality constraints multipliers
	double *t; # inequality constraints slacks
	double *t_inv; # inverse of t
	double *dv; # step in v
	double *dpi; # step in pi
	double *dlam; # step in lam
	double *dt; # step in t
	double *res_g; # q-residuals
	double *res_b; # b-residuals
	double *res_d; # d-residuals
	double *res_m; # m-residuals
	double *res_m_bkp; # m-residuals
	double *Gamma; # Hessian update
	double *gamma; # gradient update
	double alpha; # step length
	double alpha_prim; # step length
	double alpha_dual; # step length
	double sigma; # centering XXX
	double mu; # duality measuere
	double mu_aff; # affine duality measuere
	double nc_inv; # 1.0/nt, where nt is the total number of constraints
	int nv; # number of primal variables
	int ne; # number of equality constraints
	int nc; # number of (two-sided) inequality constraints
	int memsize; # memory size (in bytes) of workspace
        ]

