###################################################################################################
#                                                                                                 #
# This file is part of HPIPM.                                                                     #
#                                                                                                 #
# HPIPM -- High-Performance Interior Point Method.                                                #
# Copyright (C) 2019 by Gianluca Frison.                                                          #
# Developed at IMTEK (University of Freiburg) under the supervision of Moritz Diehl.              #
# All rights reserved.                                                                            #
#                                                                                                 #
# The 2-Clause BSD License                                                                        #
#                                                                                                 #
# Redistribution and use in source and binary forms, with or without                              #
# modification, are permitted provided that the following conditions are met:                     #
#                                                                                                 #
# 1. Redistributions of source code must retain the above copyright notice, this                  #
#    list of conditions and the following disclaimer.                                             #
# 2. Redistributions in binary form must reproduce the above copyright notice,                    #
#    this list of conditions and the following disclaimer in the documentation                    #
#    and/or other materials provided with the distribution.                                       #
#                                                                                                 #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND                 #
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED                   #
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE                          #
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR                 #
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES                  #
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;                    #
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND                     #
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT                      #
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS                   #
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                                    #
#                                                                                                 #
# Author: Andrea Zanelli, Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de             #
#                                                                                                 #
###################################################################################################

import casadi as ca
import numpy as np
#import scipy as sp
#import matplotlib.pyplot as plt
from PMSM_model import *
from reference_solver import *
import datetime
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


### problem setup ###

# PMSM parameters
par = PMSM_par()
w_ref = par.w_ref # speed reference
Ts = 5e-5

# dim
N  = 3
nx = 2
nu = 2
nq = 1

# dynamics
Ad, Bd, cd, xplus = PMSM_dynamics_disc(w_ref, Ts, par)
x0 = np.array([[-10.0], [-5.0]])

# cost: quadratic penalty w.r.t. x and u references
Qc  = 10*np.eye(nx)
QcN = 0.0001*np.eye(nx)
Rc  = 0.001*np.eye(nu)

Qd = Ts*Qc
QdN = QcN
Rd = Ts*Rc

# constr
u_max = par.u_max

# compute steady-state for x and u reference
torque_ref = 10.0
ss_calc = PMSM_ss_calculator(w_ref, Ts, par)
ss = ss_calc.compute(torque_ref)
ssf = ss.full()

xss = ssf[0:2]
uss = ssf[2:4]
#uss_lin = np.linalg.solve(Bd, xss - np.dot(Ad, xss))

# check steady state
err = np.linalg.norm(xss - np.dot(Ad, xss) - np.dot(Bd, uss) - cd)
if err > 1e-8:
	raise Exception('Computed steady state is not a steady state for the linear system')


### condensing ###

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

Q = np.kron(Qd, np.eye(N))
R = np.kron(Rd, np.eye(N))
Q[(N-1)*nx:N*nx, (N-1)*nx:N*nx] = QdN

x_ref = np.kron(np.ones((N,1)), xss)
u_ref = np.kron(np.ones((N,1)), uss)

H = np.dot(np.dot(B.transpose(), Q), B) + R

a_tilde = ( np.dot(x0.transpose(), A.transpose()) + \
		np.dot(cd.transpose(), A_hat.transpose()) ).transpose()

h = (+ np.dot(a_tilde.transpose(), np.dot(Q, B)) \
	- np.dot(x_ref.transpose(), np.dot(Q, B)) \
	- np.dot(u_ref.transpose(), R)).transpose()


## ipopt solver for condensed problem ##

if(travis_run!='true'):
	print('\n\nipopt solver\n')

# number of ipopt optmization variables
nv = N*nu

# optimization variables
x = ca.MX.sym('x', nv)

# define functions g an f: min_x f(x) s.t g(x) \leq 0 (different notation from Nishihara2015)

# consider h as parametric
hp = ca.MX.sym('hp', nv, 1)

f = 1.0/2.0*ca.mtimes(ca.mtimes(ca.transpose(x), H), x) + ca.mtimes(ca.transpose(hp), x)

g = []
for i in range(N):
	g.append(x[i*2]*x[i*2] + x[i*2+1]*x[i*2+1] - u_max**2)
g = ca.vertcat(*g)

# create reference solver (Ipopt)
ref_solver = reference_solver(x, f, g, hp)
ref_solver.update_gradient(h)

start_time = time.time()
ref_solver.solve()
end_time = time.time()
ipopt_time = end_time - start_time
if(travis_run!='true'):
	print('\nipopt solve time {:e}'.format(ipopt_time))

ref_u_sol = ref_solver.x_sol.full()
ref_x_sol = a_tilde + np.dot(B, ref_u_sol)
ref_u_sol = ref_u_sol.reshape(N,nu)
ref_x_sol = ref_x_sol.reshape(N,nx)
if(travis_run!='true'):
	print('\nu')
	print(ref_u_sol)
	print('\nx')
	print(ref_x_sol)


### HPIPM OCP QCQP solver for full-space problem ###

from hpipm_python import *
from hpipm_python.common import *

if(travis_run!='true'):
	print('\n\nHPIPM solver')

codegen_data = 1; # export qp data in the file ocp_qcqp_data.c for use from C examples

# dim
dim = hpipm_ocp_qcqp_dim(N)

dim.set('nx', nx, 0, N) # number of states
dim.set('nu', nu, 0, N-1) # number of inputs
dim.set('nbx', nx, 0) # number of state bounds
dim.set('nq', nq, 0, N-1)

# print to shell
#dim.print_C_struct()
# codegen
if codegen_data:
	dim.codegen('ocp_qcqp_data.c', 'w')


# qp
qd = -np.dot(Qd,xss)
qdN = -np.dot(QdN,xss)
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
#qp.print_C_struct()
# codegen
if codegen_data:
	qp.codegen('ocp_qcqp_data.c', 'a')



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

# codegen
if codegen_data:
	arg.codegen('ocp_qcqp_data.c', 'a')

# set up solver
solver = hpipm_ocp_qcqp_solver(dim, arg)


# solve qp
start_time = time.time()
solver.solve(qp, qp_sol)
end_time = time.time()
hpipm_time = end_time - start_time
if(travis_run!='true'):
	print('\nhpipm solve time {:e}'.format(hpipm_time))


if(travis_run!='true'):
	print('')
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



### return HPIPM status ###

if status==0:
	print('\nsuccess!\n')
else:
	print('\nSolution failed, solver returned status {0:1d}\n'.format(status))



sys.exit(int(status))


