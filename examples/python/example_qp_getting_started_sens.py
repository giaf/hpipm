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
# Author: Gianluca Frison, gianluca.frison (at) imtek.uni-freiburg.de                             #
#                                                                                                 #
###################################################################################################
from hpipm_python import *
from hpipm_python.common import *
import numpy as np
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



# define flags
codegen_data = 0 # export qp data in the file ocp_qp_data.c for use from C examples
reduce_eq_dof = 1 # eliminates e.g. equality constraint on x0
warm_start = 0 # set to 1 to warm-start the primal variable

if reduce_eq_dof:
	codegen_data = 0 # not supported yet


# dim
N = 5
nx = 2
nu = 1

nbx = nx
#ng = nx
#ns = nx

dim = hpipm_ocp_qp_dim(N)

dim.set('nx', nx, 0, N) # number of states
dim.set('nu', nu, 0, N-1) # number of inputs
dim.set('nbx', nbx, 0) # number of state bounds
#dim.set('ng', nx, 0)
dim.set('nbx', nbx, N)
#dim.set('ns', nx, N)
dim.set('nbxe', nx, 0)

# print to shell
#dim.print_C_struct()
# codegen
if codegen_data:
	dim.codegen('ocp_qp_data.c', 'w')



# data
if 0:
	# data as a contiguous array (interpreted as row-major)
	A = np.array([1, 1, 0, 1]).reshape(nx,nx)
else:
	# data as a matrix
	A = np.zeros((2,2))
	A[0][0] = 1.0
	A[0][1] = 1.0
	A[1][1] = 1.0
B = np.array([0, 1]).reshape(nx,nu)
#b = np.array([0, 0]).reshape(nx,1)

Q = np.array([1, 0, 0, 1]).reshape(nx,nx)
S = np.array([0, 0]).reshape(nu,nx)
R = np.array([1]).reshape(nu,nu)
q = np.array([1, 1]).reshape(nx,1)
#r = np.array([0]).reshape(nu,1)

Jx = np.array([1, 0, 0, 1]).reshape(nbx,nx)
x0 = np.array([1, 1]).reshape(nx,1)
#Jsx = np.array([1, 0, 0, 1]).reshape(nbx,ns)
#Zl = np.array([1e5, 1e5]).reshape(ns,1)
#Zu = np.array([1e5, 1e5]).reshape(ns,1)
#zl = np.array([1e5, 1e5]).reshape(ns,1)
#zu = np.array([1e5, 1e5]).reshape(ns,1)

if reduce_eq_dof:
	Jbxe = np.array([1, 1]).reshape(nbx,1)


# qp
qp = hpipm_ocp_qp(dim)

qp.set('A', A, 0, N-1)
qp.set('B', B, 0, N-1)
#qp.set('b', [b, b, b, b, b])
qp.set('Q', Q, 0, N)
qp.set('S', S, 0, N-1)
qp.set('R', R, 0, N-1)
qp.set('q', q, 0, N)
#qp.set('r', r, 0, N-1)
qp.set('Jx', Jx, 0)
qp.set('lx', x0, 0)
qp.set('ux', x0, 0)
#qp.set('C', Jx, 0)
#qp.set('lg', x0, 0)
#qp.set('ug', x0, 0)
qp.set('Jx', Jx, N)
#qp.set('Jsx', Jsx, N)
#qp.set('Zl', Zl, N)
#qp.set('Zu', Zu, N)
#qp.set('zl', zl, N)
#qp.set('zu', zu, N)
if reduce_eq_dof:
	qp.set('Jbxe', Jbxe, 0)

# print to shell
#qp.print_C_struct()
# codegen
if codegen_data:
	qp.codegen('ocp_qp_data.c', 'a')


# qp sol
qp_sol = hpipm_ocp_qp_sol(dim)


# set up solver arg
#mode = 'speed_abs'
mode = 'speed'
#mode = 'balance'
#mode = 'robust'
# create and set default arg based on mode
if reduce_eq_dof:
	arg = hpipm_ocp_qp_solver_arg2(dim, mode)
else:
	arg = hpipm_ocp_qp_solver_arg(dim, mode)

# create and set default arg based on mode
arg.set('mu0', 1e4)
arg.set('iter_max', 30)
arg.set('tol_stat', 1e-4)
arg.set('tol_eq', 1e-5)
arg.set('tol_ineq', 1e-5)
arg.set('tol_comp', 1e-5)
arg.set('reg_prim', 1e-12)

if reduce_eq_dof:
	arg.set('reduce_eq_dof', 1)

# if warm_start=1, then the primal variable is initialized from qp_sol
arg.set('warm_start', warm_start)
if warm_start:
	x_guess  = np.array([[1.0, 1.0],[2., -1.025], [ 0.975, -0.85 ], [ 0.125, -0.275], [-0.15, 0.15], [0.0, 0.0]])
	u_guess  = np.array([[-2.025], [0.175], [0.575], [0.425], [-0.15]])
	for i in range(N+1):
		qp_sol.set('x', i, x_guess[i])
	for i in range(N):
		qp_sol.set('u', i, u_guess[i])

# codegen
if codegen_data:
	arg.codegen('ocp_qp_data.c', 'a')

# set up solver
if reduce_eq_dof:
	solver = hpipm_ocp_qp_solver2(dim, arg)
else:
	solver = hpipm_ocp_qp_solver(dim, arg)


# solve qp
start_time = time.time()
solver.solve(qp, qp_sol)
end_time = time.time()
if(travis_run!='true'):
	print('solve time {:e}'.format(end_time-start_time))


if(travis_run!='true'):
	qp_sol.print_C_struct()

# extract and print sol
if(travis_run!='true'):
	print('u =')
u = qp_sol.get('u', 0, N)
for i in range(N+1):
	if(travis_run!='true'):
		print(u[i])


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



if status==0:
	print('\nsuccess!\n')
else:
	print('\nSolution failed, solver returned status {0:1d}\n'.format(status))



# solution sensitivities
# qp seeds
qp_seed = hpipm_ocp_qp_seed(dim)
seed_x0 = np.zeros((nx,1))
seed_x0[0] = 1.0
qp_seed.set('seed_lbx', seed_x0, 0);
qp_seed.set('seed_ubx', seed_x0, 0);
qp_seed.print_C_struct()
# qp sens
qp_sens = hpipm_ocp_qp_sol(dim)
solver.sens_frw(qp, qp_seed, qp_sens)
qp_sens.print_C_struct()

# finite differences
h = 0.001
idx = 0
dx = np.zeros([nx,1])
dx[idx] = 1.0
x0h = x0 + h*dx
#print(x0h)
qp.set('lx', x0h, 0)
qp.set('ux', x0h, 0)
#qp.print_C_struct()
qp_solh = hpipm_ocp_qp_sol(dim)
solver.solve(qp, qp_solh)
# build finite differences sol
qp_sens_fd = hpipm_ocp_qp_sol(dim)
# check u
u = qp_sol.get('u', 0, N)
uh = qp_solh.get('u', 0, N)
for i in range(N+1):
	fd = (uh[i]-u[i])/h
	qp_sens_fd.set('u', fd, i)
# check x
x = qp_sol.get('x', 0, N)
xh = qp_solh.get('x', 0, N)
for i in range(N+1):
	fd = (xh[i]-x[i])/h
	qp_sens_fd.set('x', fd, i)
# check pi
pi = qp_sol.get('pi', 0, N-1)
pih = qp_solh.get('pi', 0, N-1)
for i in range(N):
	fd = (pih[i]-pi[i])/h
	qp_sens_fd.set('pi', fd, i)
# check lam_lbx
lam_lbx = qp_sol.get('lam_lbx', 0, N)
lam_lbxh = qp_solh.get('lam_lbx', 0, N)
for i in range(N+1):
	fd = (lam_lbxh[i]-lam_lbx[i])/h
	qp_sens_fd.set('lam_lbx', fd, i)
# check lam_ubx
lam_ubx = qp_sol.get('lam_ubx', 0, N)
lam_ubxh = qp_solh.get('lam_ubx', 0, N)
for i in range(N+1):
	fd = (lam_ubxh[i]-lam_ubx[i])/h
	qp_sens_fd.set('lam_ubx', fd, i)
qp_sens_fd.print_C_struct()




#exit
sys.exit(int(status))
