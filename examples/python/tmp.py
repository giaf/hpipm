from hpipm_python import *
from hpipm_python.common import *
import numpy as np
import time

nv = 8 #number of variables
ne = 8 # number of equality constraints
ng = 0
dim = hpipm_dense_qp_dim()
dim.set('nv', nv)
dim.set('ne', ne)
dim.set('ng', ng)
H = np.array([[ 2, -2,  0,  0,  0,  0,  0,  0],
[-2,  4, -2,  0,  0,  0,  0,  0],
[ 0, -2,  4, -2,  0,  0,  0,  0],
[ 0,  0, -2,  2,  0,  0,  0,  0],
[ 0,  0,  0,  0,  2, -2,  0,  0],
[ 0,  0,  0,  0, -2,  4, -2,  0],
[ 0,  0,  0,  0,  0, -2,  4, -2],
[ 0,  0,  0,  0,  0,  0, -2,  2]])
g = np.array([[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.]])
A = np.array([[ 1,  0,  0,  0,  0,  0,  0,  0],
[ 0,  0,  0,  0,  1,  0,  0,  0],
[ 0,  0,  0,  1,  0,  0,  0,  0],
[ 0,  0,  0,  0,  0,  0,  0,  1],
[-3,  3,  0,  0,  0,  0,  0,  0],
[ 0,  0,  0,  0, -3,  3,  0,  0],
[ 0,  0, -3,  3,  0,  0,  0,  0],
[ 0,  0,  0,  0,  0,  0, -3,  3]])
b = np.array([[ 1. ],[ 2. ],[10. ],[12. ],[ 1.5],[ 3. ],[ 0. ],[ 1.5]])
qp = hpipm_dense_qp(dim)
qp.set('H', H)
qp.set('g', g)
qp.set('A', A)
qp.set('b', b)
qp.print_C_struct()

qp_sol = hpipm_dense_qp_sol(dim)
mode = 'speed'
arg = hpipm_dense_qp_solver_arg(dim, mode)
arg.set('mu0', 1e4)
arg.set('iter_max', 30)
arg.set('tol_stat', 1e-4)
arg.set('tol_eq', 1e-5)
arg.set('tol_ineq', 1e-5)
arg.set('tol_comp', 1e-5)
arg.set('reg_prim', 1e-12)
solver = hpipm_dense_qp_solver(dim, arg)
optimization_start_time = time.time()
solver.solve(qp, qp_sol)
optimization_end_time = time.time()
v = qp_sol.get('v')

print(v)

print("A @ v[:,None]: " , np.dot(A, v))
