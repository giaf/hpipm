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
from utils import check_env, is_travis_run


def main(travis_run):

	# define flags
	codegen_data = 1 # export qp data in the file ocp_qp_data.c for use from C examples
	print_structs = 1

	reduce_eq_dof = 0 # eliminates e.g. equality constraint on x0

	if reduce_eq_dof:
		codegen_data = 0 # not supported yet

	# dim
	N = 5
	nx = 2
	nu = 1
	nbx = nx

	dim = hpipm_ocp_qp_dim(N)

	dim.set('nx', nx, 0, N) # number of states
	dim.set('nu', nu, 0, N-1) # number of inputs
	dim.set('nbxe', nx, 0)
	dim.set('nbx', nbx, N)

	# data
	if 0:
		# data as a contiguous array (interpreted as row-major)
		A = np.array([1, 1, 0, 1]).reshape(nx,nx)
	else:
		# data as a matrix
		A = np.zeros((2,2))
		A[0, 0] = 1.0
		A[0, 1] = 1.0
		A[1, 1] = 1.0

	B = np.array([0, 1]).reshape(nx,nu)
	# b = np.array([0, 0]).reshape(nx,1)

	Q = np.array([1, 0, 0, 1]).reshape(nx,nx)
	S = np.array([0, 0]).reshape(nu,nx)
	R = np.array([1]).reshape(nu,nu)
	q = np.array([1, 1]).reshape(nx,1)
	# r = np.array([0]).reshape(nu,1)

	Jx = np.eye(nx)
	x0 = np.array([1, 1]).reshape(nx,1)

	# qp
	qp = hpipm_ocp_qp(dim)

	qp.set('A', A, 0, N-1)
	qp.set('B', B, 0, N-1)
	# qp.set('b', [b, b, b, b, b])
	qp.set('Q', Q, 0, N)
	qp.set('S', S, 0, N-1)
	qp.set('R', R, 0, N-1)
	qp.set('q', q, 0, N)
	# qp.set('r', r, 0, N-1)

	qp.set('Jx', Jx, 0)
	qp.set('lx', x0, 0)
	qp.set('ux', x0, 0)

	qp.set('Jx', Jx, N)

	if reduce_eq_dof:
		Jbxe = np.array([1, 1]).reshape(nbx, 1)
		qp.set('Jbxe', Jbxe, 0)

	# qp sol
	qp_sol = hpipm_ocp_qp_sol(dim)

	# set up solver arg
	# mode = 'speed_abs'
	mode = 'speed'
	# mode = 'balance'
	# mode = 'robust'

	# create and set default arg based on mode
	if reduce_eq_dof:
		arg = hpipm_ocp_qp_solver_arg2(dim, mode)
		arg.set('reduce_eq_dof', 1)
	else:
		arg = hpipm_ocp_qp_solver_arg(dim, mode)

	arg.set('mu0', 1e4)
	arg.set('iter_max', 30)
	arg.set('tol_stat', 1e-4)
	arg.set('tol_eq', 1e-5)
	arg.set('tol_ineq', 1e-5)
	arg.set('tol_comp', 1e-5)
	arg.set('reg_prim', 1e-12)

	# set up solver
	if reduce_eq_dof:
		solver = hpipm_ocp_qp_solver2(dim, arg)
	else:
		solver = hpipm_ocp_qp_solver(dim, arg)

	# solve qp
	start_time = time.time()
	solver.solve(qp, qp_sol)
	end_time = time.time()

	# extract and print sol
	u_traj = qp_sol.get('u', 0, N-1)
	x_traj = qp_sol.get('x', 0, N)

	K_traj = solver.get_feedback(qp, 'ric_K', 0, N-1)

	if not travis_run:
		print('solve time {:e}'.format(end_time-start_time))

		for n in range(N):
			print(f'\nu_{n} =')
			print(u_traj[n].flatten())

		for n in range(N):
			print(f'\nK_{n} =')
			print(K_traj[n].flatten())

		for n in range(N+1):
			print(f'\nx_{n} =')
			print(x_traj[n].flatten())

	# get solver statistics
	status = solver.get('status')
	# res_stat = solver.get('max_res_stat')
	# res_eq = solver.get('max_res_eq')
	# res_ineq = solver.get('max_res_ineq')
	# res_comp = solver.get('max_res_comp')
	# iters = solver.get('iter')

	if status == 0:
		print('\nsuccess!\n')
	else:
		print('\nSolution failed, solver returned status {0:1d}\n'.format(status))

	if not travis_run:
		solver.print_stats()

	if not travis_run and print_structs:
		dim.print_C_struct()
		qp.print_C_struct()
		qp_sol.print_C_struct()

	# codegen
	if codegen_data:
		dim.codegen('ocp_qp_data.c', 'w')
		qp.codegen('ocp_qp_data.c', 'a')
		arg.codegen('ocp_qp_data.c', 'a')

	return status


if __name__ == '__main__':

	check_env()
	travis_run = is_travis_run()

	status = main(travis_run)

	sys.exit(int(status))