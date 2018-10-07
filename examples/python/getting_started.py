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
dims.set('nx', np.array([2, 2, 2, 2, 2, 2]))
end_time = time.time()
print('set nx time {:e}'.format(end_time-start_time))
dims.set('nu', np.array([1, 1, 1, 1, 1]))
dims.set('nbx', 2, 0)
#dims.set('ng', 2, 0)
dims.set('nbx', 2, 5)
#dims.set('ns', 2, 5)


# data
if 0:
	# data as a contiguous array (interpreted as row-major)
	A = np.array([1, 1, 0, 1])
else:
	# data as a matrix
	A = np.zeros((2,2))
	A[0][0] = 1.0
	A[0][1] = 1.0
	A[1][1] = 1.0
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
qp.set('A', [A, A, A, A, A])
end_time = time.time()
print('set A time {:e}'.format(end_time-start_time))
qp.set('B', [B, B, B, B, B])
#qp.set('b', [b, b, b, b, b])
qp.set('Q', [Q, Q, Q, Q, Q, Q])
qp.set('S', [S, S, S, S, S])
qp.set('R', [R, R, R, R, R])
qp.set('q', [q, q, q, q, q, q])
#qp.set('r', [r, r, r, r, r])
qp.set('Jx', Jx, 0)
qp.set('lx', x0, 0)
qp.set('ux', x0, 0)
#qp.set('C', Jx, 0)
#qp.set('lg', x0, 0)
#qp.set('ug', x0, 0)
qp.set('Jx', Jx, 5)
#qp.set('Jsx', Jsx, 5)
#qp.set('Zl', Zl, 5)
#qp.set('Zu', Zu, 5)
#qp.set('zl', zl, 5)
#qp.set('zu', zu, 5)

qp.print_C_struct()


# qp sol
start_time = time.time()
qp_sol = hpipm_ocp_qp_sol(dims)
end_time = time.time()
print('create qp_sol time {:e}'.format(end_time-start_time))


# set up solver arg
start_time = time.time()
arg = hpipm_ocp_qp_solver_arg(dims)
end_time = time.time()
print('create solger arguments time {:e}'.format(end_time-start_time))

arg.set_mu0(1e4)
arg.set_iter_max(30)
arg.set_tol_stat(1e-4)
arg.set_tol_eq(1e-5)
arg.set_tol_ineq(1e-5)
arg.set_tol_comp(1e-5)
arg.set_reg_prim(1e-12)


# set up solver
start_time = time.time()
solver = hpipm_ocp_qp_solver(dims, arg)
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
    qp_sol.print_C_struct()
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


# print solver statistics
print('\nsolver statistics:\n')
print('ipm return = {0:1d}\n'.format(return_flag))
res_stat = solver.get_res_stat()
print('ipm max res stat = {:e}\n'.format(res_stat))
res_eq = solver.get_res_eq()
print('ipm max res eq   = {:e}\n'.format(res_eq))
res_ineq = solver.get_res_ineq()
print('ipm max res ineq = {:e}\n'.format(res_ineq))
res_comp = solver.get_res_comp()
print('ipm max res comp = {:e}\n'.format(res_comp))
iters = solver.get_iter()
print('ipm iter = {0:1d}\n'.format(iters))
stat = solver.get_stat()
print('stat =')
print('\talpha_aff\tmu_aff\t\tsigma\t\talpha\t\tmu')
for ii in range(iters):
	print('\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}'.format(stat[ii][0], stat[ii][1], stat[ii][2], stat[ii][3], stat[ii][4]))
print('')
