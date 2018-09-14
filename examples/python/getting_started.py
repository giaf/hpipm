from hpipm_python import *
from hpipm_python.common import *
import numpy as np
import time



# dims
N = 5

start_time = time.time()
dims = hpipm_ocp_qp_dim(N)
end_time = time.time()
print('create dim time {:e}'.format(end_time-start_time))

start_time = time.time()
dims.set_nx(np.array([2, 2, 2, 2, 2, 2]))
end_time = time.time()
print('set nx time {:e}'.format(end_time-start_time))
dims.set_nu(np.array([1, 1, 1, 1, 1]))
dims.set_nbx(2, 0)
#dims.set_ng(2, 0)
dims.set_nbx(2, 5)
#dims.set_ns(2, 5)


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
Zl = np.array([1e5, 0, 0, 1e5])
Zu = np.array([1e5, 0, 0, 1e5])
zl = np.array([1e5, 1e5])
zu = np.array([1e5, 1e5])



# qp
start_time = time.time()
qp = hpipm_ocp_qp(dims)
end_time = time.time()
print('create qp time {:e}'.format(end_time-start_time))

start_time = time.time()
qp.set_A([A, A, A, A, A])
end_time = time.time()
print('set A time {:e}'.format(end_time-start_time))
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
#qp.set_Jsx(Jsx, 5)
#qp.set_Zl(Zl, 5)
#qp.set_Zu(Zu, 5)
#qp.set_zl(zl, 5)
#qp.set_zu(zu, 5)

#qp.print_C_struct()


# qp sol
start_time = time.time()
qp_sol = hpipm_ocp_qp_sol(dims)
end_time = time.time()
print('create qp_sol time {:e}'.format(end_time-start_time))


# set up solver
start_time = time.time()
solver = hpipm_ocp_qp_solver(dims)
end_time = time.time()
print('create solver time {:e}'.format(end_time-start_time))


# solve qp
start_time = time.time()
return_flag = solver.solve(qp, qp_sol)
end_time = time.time()
print('solve time {:e}'.format(end_time-start_time))

print('HPIPM returned with flag {0:1d}.'.format(return_flag))

if return_flag == 0:
    print('-> QP solved! Solution:\n')
    qp_sol.print()
else:
    print('-> Solver failed!')

# extract and print sol
print('u =')
u = qp_sol.get_u()
for i in range(N+1):
	print(u[i])

print('x =')
for i in range(N+1):
	tmp = qp_sol.get_x(i)
	print(tmp)

