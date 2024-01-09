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


def main(travis_run: bool):

	# define flags
	codegen_data = 1; # export qp data in the file ocp_qcqp_data.c for use from C examples

	# dim
	N = 5
	nx = 2
	nu = 1

	nbx = nx
	nq = 1

	dim = hpipm_ocp_qcqp_dim(N)

	dim.set('nx', nx, 0, N) # number of states
	dim.set('nu', nu, 0, N-1) # number of inputs
	dim.set('nbx', nbx, 0) # number of state bounds
	dim.set('nq', nq, N)

	# print to shell
	# dim.print_C_struct()

	# data
	A = np.array([1, 1, 0, 1]).reshape(nx,nx)
	B = np.array([0, 1]).reshape(nx,nu)

	Q = np.array([1, 0, 0, 1]).reshape(nx,nx)
	R = np.array([1]).reshape(nu,nu)
	q = np.array([1, 1]).reshape(nx,1)

	Jx = np.array([1, 0, 0, 1]).reshape(nbx,nx)
	x0 = np.array([1, 1]).reshape(nx,1)

	Qq = 2*np.eye(nx)
	uq = np.array([0.1])

	# qp
	qp = hpipm_ocp_qcqp(dim)

	qp.set('A', A, 0, N-1)
	qp.set('B', B, 0, N-1)
	qp.set('Q', Q, 0, N)
	qp.set('R', R, 0, N-1)
	qp.set('q', q, 0, N)
	qp.set('Jx', Jx, 0)
	qp.set('lx', x0, 0)
	qp.set('ux', x0, 0)
	qp.set('Qq', Qq, N)
	qp.set('uq', uq, N)

	# print to shell
	# qp.print_C_struct()

	# qp sol
	qp_sol = hpipm_ocp_qcqp_sol(dim)

	# set up solver arg
	# mode = 'speed_abs'
	mode = 'speed'
	# mode = 'balance'
	# mode = 'robust'

	# create and set default arg based on mode
	arg = hpipm_ocp_qcqp_solver_arg(dim, mode)

	arg.set('mu0', 1e1)
	arg.set('iter_max', 30)
	arg.set('tol_stat', 1e-4)
	arg.set('tol_eq', 1e-5)
	arg.set('tol_ineq', 1e-5)
	arg.set('tol_comp', 1e-5)
	arg.set('reg_prim', 1e-12)

	# codegen
	if codegen_data:
		dim.codegen('ocp_qcqp_data.c', 'w')
		arg.codegen('ocp_qcqp_data.c', 'a')
		qp.codegen('ocp_qcqp_data.c', 'a')

	# set up solver
	solver = hpipm_ocp_qcqp_solver(dim, arg)

	# solve qp
	start_time = time.time()
	solver.solve(qp, qp_sol)
	end_time = time.time()

	if not travis_run:
		print('solve time {:e}'.format(end_time-start_time))
		qp_sol.print_C_struct()

	# extract and print sol
	u_traj = qp_sol.get('u', 0, N)
	x_traj = qp_sol.get('x', 0, N)

	xN = x_traj[-1]

	if not travis_run:
		print('u =')
		for u in u_traj:
				print(u)
		print('x =')
		for x in x_traj:
				print(x)

	# get solver statistics
	status = solver.get('status')
	res_stat = solver.get('max_res_stat')
	res_eq = solver.get('max_res_eq')
	res_ineq = solver.get('max_res_ineq')
	res_comp = solver.get('max_res_comp')
	iters = solver.get('iter')

	if not travis_run:

		print('quadratic constr')
		print(0.5*np.dot(xN.transpose(),np.dot(Qq,xN)))

		print('ipm return = {0:1d}\n'.format(status))
		print('ipm max res stat = {:e}\n'.format(res_stat))
		print('ipm max res eq   = {:e}\n'.format(res_eq))
		print('ipm max res ineq = {:e}\n'.format(res_ineq))
		print('ipm max res comp = {:e}\n'.format(res_comp))
		print('ipm iter = {0:1d}\n'.format(iters))

		solver.print_stats()

	if status == 0:
		print('\nsuccess!\n')
	else:
		print('\nSolution failed, solver returned status {0:1d}\n'.format(status))

	return status


if __name__ == '__main__':

	check_env()
	travis_run = is_travis_run()

	status = main(travis_run)

	sys.exit(int(status))