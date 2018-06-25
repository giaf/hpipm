from ctypes import *
from os import *

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




# d_ocp_qp - include/hpipm_d_ocp.h
class d_ocp_qp(Structure):
    _fields_ = [
        ("dim",     POINTER(d_ocp_qp_dim)),     # ocp_qp dimensios struct
        ("BAbt",    POINTER(blasfeo_dmat)),     # dynamics matrix & vector work space
        ("RSQrq",   POINTER(blasfeo_dmat)),     # hessian of cost & vector work space
        ("DCt",     POINTER(blasfeo_dmat)),     # inequality constraints matrix
        ("b",       POINTER(blasfeo_dvec)),     # dynamics vector
        ("rqz",     POINTER(blasfeo_dvec)),     # gradient of cost & gradient of slacks
        ("d",       POINTER(blasfeo_dvec)),     # inequality constraints vector
        ("m",       POINTER(blasfeo_dvec)),     # rhs of complementarity condition
        ("Z",       POINTER(blasfeo_dvec)),     # (diagonal) hessian of slacks
        ("idxb",    POINTER(POINTER(c_int)),    # index of box constraints
        ("idxs",    POINTER(POINTER(c_int)),    # index of soft constraints
        ("memsize", c_int)                      # memory size in bytes
    



# d_ocp_qp - include/hpipm_d_ocp_qp_sol.h
class d_ocp_qp_sol(Structure):
    _fields_ = [
        ("dim",     POINTER(d_ocp_qp_dim),      # dimensions
        ("ux",      POINTER(blasfeo_dvec),      # input-states
        ("pi",      POINTER(blasfeo_dvec),      # eq. multipliers
        ("lam",     POINTER(blasfeo_dvec),      # ineq. multipliers
        ("t",       POINTER(blasfeo_dvec),      # slack variables 
        ("misc",    c_void_p)                   # miscellaneous 
        ("memsize", c_int)                      # memory size in bytes
        ]
        



# d_ocp_qp - hpipm_d_ocp_qp_ipm.h
class d_ocp_qp_ipm_arg(Structure):
    _fields_ = [
        ("mu0",             c_double), # initial value for duality measure
        ("alpha_min",       c_double), # exit cond on step length
        ("res_g_max",       c_double), # exit cond on inf norm of residuals
        ("res_b_max",       c_double), # exit cond on inf norm of residuals
        ("res_d_max",       c_double), # exit cond on inf norm of residuals
        ("res_m_max",       c_double), # exit cond on inf norm of residuals
        ("reg_prim",        c_double), # reg of primal hessian
        ("lam_min",         c_double), # min value in lam vector
        ("t_min",           c_double), # min value in t vector
        ("iter_max",        c_int),    # exit cond in iter number
        ("stat_max",        c_int),    # iterations saved in stat
        ("pred_corr",       c_int),    # use Mehrotra's predictor-corrector IPM algirthm
        ("cond_pred_corr",  c_int),    # conditional Mehrotra's predictor-corrector
        ("itref_pred_max",  c_int),    # max number of iterative refinement steps for predictor step
        ("itref_corr_max",  c_int),    # max number of iterative refinement steps for corrector step
        ("warm_start",      c_int),    # 0 no warm start, 1 warm start primal sol
        ("lq_fact",         c_int),    # 0 syrk+potrf, 1 mix, 2 lq
        ("abs_form",        c_int),    # absolute IPM formulation
        ("comp_dual_sol",   c_int),    # dual solution (only for abs_form==1)
        ("comp_res_exit",   c_int),    # compute res. on exit (only for abs_form==1 and comp_dual_sol==1)
        ("memsize",         c_int)     # memory size
        ]       




class d_ocp_qp_ipm_workspace(Structure):
    _fields_ [
        ("core_workspace",  POINTER(d_core_qp_ipm_workspace)),
        ("res_workspace",   POINTER(d_ocp_qp_res_workspace)),
        ("sol_step",        POINTER(d_ocp_qp_sol)),
        ("sol_itref",       POINTER(d_ocp_qp_sol),
        ("qp_step",         POINTER(d_ocp_qp),
        ("qp_itref",        POINTER(d_ocp_qp),
        ("res_itref",       POINTER(d_ocp_qp_res),
        ("res",             POINTER(d_ocp_qp_res),
        ("Gamma",           POINTER(blasfeo_dvec), # hessian update
        ("gamma",           POINTER(blasfeo_dvec), # hessian update
        ("tmp_nsM",         POINTER(blasfeo_dvec), # work space of size nxM
        ("tmp_nbgM",        POINTER(blasfeo_dvec), # work space of size nbM+ngM
        ("tmp_nsM",         POINTER(blasfeo_dvec), # work space of size nsM
        ("Pb",              POINTER(blasfeo_dvec), # Pb
        ("Zs_inv",          POINTER(blasfeo_dvec),
        ("L",               POINTER(blasfeo_dmat),
        ("Lh",              POINTER(blasfeo_dmat),
        ("AL",              POINTER(blasfeo_dmat),
        ("lq0",             POINTER(blasfeo_dmat),
        ("tmp_m"            POINTER(blasfeo_dvec),

        ("stat",            POINTER(c_double)), # convergence statistics
        ("use_hess_fact",   POINTER(c_int)),
        ("lq_work0",        c_void_p),
        ("qp_res",          c_int*4), # infinity norm of residuals
        ("iter",            c_int), # iteration number
        ("stat_max",        c_int), # iterations saved in stat
        ("use_Pb",          c_int),
        ("memsize",         c_int),
        ]
