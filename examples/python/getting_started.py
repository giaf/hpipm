from hpipm_python import *
from hpipm_python.common import *
import numpy as np



# dims
qp_dims = hpipm_dims(5)

qp_dims.set_nx(np.array([2, 2, 2, 2, 2, 2]))
qp_dims.set_nu(np.array([1, 1, 1, 1, 1]))
qp_dims.set_nbx(2, 0)
#qp_dims.set_nbu(np.array([0, 0, 0, 0, 0, 0]))
#qp_dims.set_ng(np.array([0, 0, 0, 0, 0, 0]))
#qp_dims.set_ns(np.array([0, 0, 0, 0, 0, 0]))

print(qp_dims.nx)
print(qp_dims.nu)
print(qp_dims.nbx)
print(qp_dims.nbu)
print(qp_dims.ng)
print(qp_dims.ns)



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

# set up solver
solver = hpipm_solver(qp_dims, qp_data)

# solve qp
return_flag = solver.solve()

print('HPIPM returned with flag {0:1d}.'.format(return_flag))

if return_flag == 0:
    print('-> QP solved! Solution:\n')
    solver.print_sol()
else:
    print('-> Solver failed!')

