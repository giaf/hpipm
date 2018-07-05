from hpipm_python import *
import numpy as np

qp_data = hpipm_data()

A = np.array([1, 0, 1, 1])
B = np.array([0, 1])
b = np.array([0, 0])

Q = np.array([1, 0, 0, 1])
S = np.array([0, 0])
R = np.array([1])
q = np.array([1, 1])
r = np.array([0])

qp_data.A = [A, A, A, A, A]
qp_data.B = [B, B, B, B, B]
qp_data.b = [b, b, b, b, b]
qp_data.Q = [Q, Q, Q, Q, Q, Q]
qp_data.S = [S, S, S, S, S, S]
qp_data.R = [R, R, R, R, R, R]
qp_data.q = [q, q, q, q, q, q]
qp_data.r = [r, r, r, r, r, r]

x0 = np.array([1, 1])

qp_data.d_lb = [x0]
qp_data.d_ub = [x0]

qp_data.idxb = [np.array([1, 2])]

qp_dims = hpipm_dims()

qp_dims.nx = np.array([2, 2, 2, 2, 2, 2])
qp_dims.nu = np.array([1, 1, 1, 1, 1, 0])
qp_dims.nb = np.array([2, 0, 0, 0, 0, 0])
qp_dims.ng = np.array([0, 0, 0, 0, 0, 0])
qp_dims.ns = np.array([0, 0, 0, 0, 0, 0])
qp_dims.nbx = np.array([2, 0, 0, 0, 0, 0])
qp_dims.nbu = np.array([0, 0, 0, 0, 0, 0])
qp_dims.N = 5

# set up solver
solver = hpipm_solver(qp_dims, qp_data)

# solve qp
return_flag = solver.solve()

print('hpipm returned with flag {0:2d}'.format(return_flag))
solver.print_sol()
