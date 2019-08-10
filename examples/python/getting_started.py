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



# define flags
codegen_data = 1; # export qp data in the file qp_data.c for use from C examples



# dims
N = 5

start_time = time.time()
dims = hpipm_ocp_qp_dim(N)
end_time = time.time()
print('create dim time {:e}'.format(end_time-start_time))

start_time = time.time()
dims.set('nx', np.array([2, 2, 2, 2, 2, 2])) # number of states
end_time = time.time()
print('set nx time {:e}'.format(end_time-start_time))
dims.set('nu', np.array([1, 1, 1, 1, 1])) # number of inputs
dims.set('nbx', 2, 0) # number of state bounds
#dims.set('ng', 2, 0)
dims.set('nbx', 2, 5)
#dims.set('ns', 2, 5)

# print to shell
#dims.print_C_struct()
# codegen
if codegen_data:
	dims.codegen('qp_data.c', 'w')



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

# print to shell
#qp.print_C_struct()
# codegen
if codegen_data:
	qp.codegen('qp_data.c', 'a')


# qp sol
start_time = time.time()
qp_sol = hpipm_ocp_qp_sol(dims)
end_time = time.time()
print('create qp_sol time {:e}'.format(end_time-start_time))


# set up solver arg
start_time = time.time()
arg = hpipm_ocp_qp_solver_arg(dims)
end_time = time.time()
print('create solver arguments time {:e}'.format(end_time-start_time))

arg.set_mu0(1e4)
arg.set_iter_max(30)
arg.set_tol_stat(1e-4)
arg.set_tol_eq(1e-5)
arg.set_tol_ineq(1e-5)
arg.set_tol_comp(1e-5)
arg.set_reg_prim(1e-12)

# codegen
if codegen_data:
	arg.codegen('qp_data.c', 'a')

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
res_stat = solver.get_max_res_stat()
print('ipm max res stat = {:e}\n'.format(res_stat))
res_eq = solver.get_max_res_eq()
print('ipm max res eq   = {:e}\n'.format(res_eq))
res_ineq = solver.get_max_res_ineq()
print('ipm max res ineq = {:e}\n'.format(res_ineq))
res_comp = solver.get_max_res_comp()
print('ipm max res comp = {:e}\n'.format(res_comp))
iters = solver.get_iter()
print('ipm iter = {0:1d}\n'.format(iters))
stat = solver.get_stat()
print('stat =')
print('\titer\talpha_aff\tmu_aff\t\tsigma\t\talpha\t\tmu\t\tres_stat\tres_eq\t\tres_ineq\tres_comp')
for ii in range(iters+1):
	print('\t{:d}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}'.format(ii, stat[ii][0], stat[ii][1], stat[ii][2], stat[ii][3], stat[ii][4], stat[ii][5], stat[ii][6], stat[ii][7], stat[ii][8]))
print('')
