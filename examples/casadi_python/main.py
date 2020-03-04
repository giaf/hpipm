import casadi as ca
import numpy as np
import scipy as sp
#import matplotlib.pyplot as plt
from PMSM_model import *
from reference_solver import *
import datetime

N  = 12
nx = 2
nu = 2
nv = N*nu
Ts = 0.00005
x0 = np.array([[-10.0], [-5.0]])

Qc  = 10*np.eye(nx)
QNc = 0.0001*np.eye(nx)
Rc  = 0.001*np.eye(nu)

torque_ref = 10.0

# rho parameter
rho = 0.0000002
# alpha parameter
alpha = 1.0
print(alpha)

par = PMSM_par()
u_max = par.u_max
w_ref = par.w_ref
Ad, Bd, cd, xplus = get_PMSM_dynamics(w_ref, Ts, par)
print(Ad)

# build A and B matrices
A = np.zeros((nx*N, nx)) 
B = np.zeros((nx*N, nu*N)) 
A_hat = np.zeros((nx*N, nx))

for i in range(N):
    A[i*nx:(i+1)*nx,:] = np.linalg.matrix_power(Ad, i+1)
    A_hat[i*nx:(i+1)*nx,:] = np.eye(nx)
    for j in range(i):
        A_hat[i*nx:(i+1)*nx,:] =  A_hat[i*nx:(i+1)*nx,:] + \
            np.linalg.matrix_power(Ad, j+1)

for i in range(N):
    for j in range(N-i):
        B[(i+j)*nx:(i+j+1)*nx, (j)*nu:(j+1)*nu]  = \
            np.dot(np.linalg.matrix_power(Ad, i), Bd)

Q = np.kron(Ts*Qc, np.eye(N))
R = np.kron(Ts*Rc, np.eye(N))
Q[(N-1)*nx:N*nx, (N-1)*nx:N*nx] = QNc

# compute steady-state
ss_calc = ss_calculator(w_ref, Ts, par)
ss = ss_calc.compute(torque_ref)
ssf = ss.full()

xss = ssf[0:2]
uss = ssf[2:4]
uss_lin = np.linalg.solve(Bd, xss - np.dot(Ad, xss))
x_ref = np.kron(np.ones((N,1)), xss)
u_ref = np.kron(np.ones((N,1)), uss)

# check steady state
err = np.linalg.norm(xss - np.dot(Ad, xss) - np.dot(Bd, uss) - cd)
if err > 1e-8:
    raise Exception('Computed steady state is not a steady state for the linear system')

# optimization variables
x = ca.MX.sym('x', nv)
z = ca.MX.sym('z', nv)
u = ca.MX.sym('u', nv)

# define functions g an f: min_x f(x) s.t g(x) \leq 0 (different notation from Nishihara2015)

H = np.dot(np.dot(B.transpose(), Q), B) + R

a_tilde = ( np.dot(x0.transpose(), A.transpose()) + \
        np.dot(cd.transpose(), A_hat.transpose()) ).transpose()

h = (+ np.dot(a_tilde.transpose(), np.dot(Q, B)) \
    - np.dot(x_ref.transpose(), np.dot(Q, B)) \
    - np.dot(u_ref.transpose(), R)).transpose()

f = 1.0/2.0*ca.mtimes(ca.mtimes(ca.transpose(x), H), x) + ca.mtimes(ca.transpose(h), x)

g = []
for i in range(N):
    g.append(x[i*2]*x[i*2] + x[i*2+1]*x[i*2+1] - u_max**2)

g = ca.vertcat(*g)

# create reference solver (Ipopt)
ref_solver = reference_solver(nv, x, g, h, H)

ref_solver.solve()
ref_sol = ref_solver.x.full()





### HPIPM OCP QCQP ###

from hpipm_python import *
from hpipm_python.common import *
import time
import sys
import os

# check that env.sh has been run
env_run = os.getenv('ENV_RUN')
if env_run!='true':
	print('ERROR: env.sh has not been sourced! Before executing this example, run:')
	print('source env.sh')
	sys.exit(1)

travis_run = os.getenv('TRAVIS_RUN')
#travis_run = 'true'

# dim
nq = 1

dim = hpipm_ocp_qcqp_dim(N)

dim.set('nx', nx, 0, N) # number of states
dim.set('nu', nu, 0, N-1) # number of inputs
dim.set('nbx', nx, 0) # number of state bounds
dim.set('nq', nq, 0, N-1)

# print to shell
dim.print_C_struct()


# qp
Qd = Ts*Qc
qd = -np.dot(Qd,xss)
QdN = QNc
qdN = -np.dot(QdN,xss)
Rd = Ts*Rc
rd = -np.dot(Rd,uss)

Jx0 = np.array([1, 0, 0, 1]).reshape(nx,nx)

Rq = 2*np.eye(nx)
uq = np.array([u_max*u_max])

qp = hpipm_ocp_qcqp(dim)

qp.set('A', Ad, 0, N-1)
qp.set('B', Bd, 0, N-1)
qp.set('b', cd, 0, N-1)
qp.set('Q', Qd, 0, N-1)
qp.set('q', qd, 0, N-1)
qp.set('Q', QdN, N)
qp.set('q', qdN, N)
qp.set('R', Rd, 0, N-1)
qp.set('r', rd, 0, N-1)
qp.set('Jx', Jx0, 0)
qp.set('lx', x0, 0)
qp.set('ux', x0, 0)
qp.set('Rq', Rq, 0, N-1)
qp.set('uq', uq, 0, N-1)

# print to shell
qp.print_C_struct()



# qp sol
qp_sol = hpipm_ocp_qcqp_sol(dim)


# set up solver arg
#mode = 'speed_abs'
mode = 'speed'
#mode = 'balance'
#mode = 'robust'
# create and set default arg based on mode
arg = hpipm_ocp_qcqp_solver_arg(dim, mode)

# create and set default arg based on mode
arg.set('mu0', 1e1)
arg.set('iter_max', 40)
arg.set('tol_stat', 1e-8)
arg.set('tol_eq', 1e-8)
arg.set('tol_ineq', 1e-8)
arg.set('tol_comp', 1e-8)
arg.set('reg_prim', 1e-12)

# set up solver
solver = hpipm_ocp_qcqp_solver(dim, arg)


# solve qp
start_time = time.time()
solver.solve(qp, qp_sol)
end_time = time.time()
if(travis_run!='true'):
	print('solve time {:e}'.format(end_time-start_time))


if(travis_run!='true'):
	qp_sol.print_C_struct()

if(travis_run!='true'):
	print('quadratic constr')
for i in range(N):
	ui = qp_sol.get('u', i)
	if(travis_run!='true'):
		print(0.5*np.dot(ui.transpose(),np.dot(Rq,ui)))

# get solver statistics
status = solver.get('status')
res_stat = solver.get('max_res_stat')
res_eq = solver.get('max_res_eq')
res_ineq = solver.get('max_res_ineq')
res_comp = solver.get('max_res_comp')
iters = solver.get('iter')
stat = solver.get('stat')
if(travis_run!='true'):
	print('\nsolver statistics:\n')
	print('ipm return = {0:1d}\n'.format(status))
	print('ipm max res stat = {:e}\n'.format(res_stat))
	print('ipm max res eq   = {:e}\n'.format(res_eq))
	print('ipm max res ineq = {:e}\n'.format(res_ineq))
	print('ipm max res comp = {:e}\n'.format(res_comp))
	print('ipm iter = {0:1d}\n'.format(iters))
	print('stat =')
	print('\titer\talpha_aff\tmu_aff\t\tsigma\t\talpha_prim\talpha_dual\tmu\t\tres_stat\tres_eq\t\tres_ineq\tres_comp')
	for ii in range(iters+1):
		print('\t{:d}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}'.format(ii, stat[ii][0], stat[ii][1], stat[ii][2], stat[ii][3], stat[ii][4], stat[ii][5], stat[ii][6], stat[ii][7], stat[ii][8], stat[ii][9]))
	print('')


