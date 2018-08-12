from hpipm_python import *
from hpipm_python.common import *
import numpy as np
import time



# dims
dims = hpipm_ocp_qp_dims(5)

dims.set_nx(np.array([2, 2, 2, 2, 2, 2]))
dims.set_nu(np.array([1, 1, 1, 1, 1]))
dims.set_nbx(2, 0)
#dims.set_ng(2, 0)
dims.set_nbx(2, 5)
dims.set_ns(2, 5)


# data
A = np.array([1, 0, 1, 1])
B = np.array([0, 1])
#b = np.array([0, 0])

Q = np.array([1, 0, 0, 1])
S = np.array([0, 0])
R = np.array([1])
q = np.array([1, 1])
#r = np.array([0])

Jx = np.array([1, 0, 0, 1])
x0 = np.array([1, 1])
Jsx = np.array([1, 0, 0, 1])
Zl = np.array([0e5, 0, 0, 0e5])
Zu = np.array([0e5, 0, 0, 0e5])



# qp
qp = hpipm_ocp_qp(dims)

qp.set_A([A, A, A, A, A])
qp.set_B([B, B, B, B, B])
#qp.set_b([b, b, b, b, b])
qp.set_Q([Q, Q, Q, Q, Q, Q])
qp.set_S([S, S, S, S, S])
qp.set_R([R, R, R, R, R])
qp.set_q([q, q, q, q, q, q])
#qp.set_r([r, r, r, r, r])
qp.set_Jx(Jx, 0)
qp.set_lx(x0, 0)
qp.set_ux(x0, 0)
#qp.set_C(Jx, 0)
#qp.set_lg(x0, 0)
#qp.set_ug(x0, 0)
qp.set_Jx(Jx, 5)
qp.set_Jsx(Jsx, 5)
qp.set_Zl(Zl, 5)
qp.set_Zu(Zu, 5)


# set up solver
start_time = time.time()
solver = hpipm_solver(dims, qp)
end_time = time.time()

print('create solver time {:e}'.format(end_time-start_time))

# solve qp
start_time = time.time()
return_flag = solver.solve()
end_time = time.time()

print('solve time {:e}'.format(end_time-start_time))

print('HPIPM returned with flag {0:1d}.'.format(return_flag))

if return_flag == 0:
    print('-> QP solved! Solution:\n')
    solver.print_sol()
else:
    print('-> Solver failed!')

