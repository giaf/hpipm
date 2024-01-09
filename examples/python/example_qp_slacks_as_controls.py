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

	x0_elim = True
	# x0_elim = False

	# define flags
	codegen_data = 1; # export qp data in the file ocp_qp_data.c for use from C examples

	bigM = 1e2

	# dim
	N = 1
	nx = 2
	nu = 1

	nbx = nx
	nbu = nu
	#ng = nx
	ns = nx

	dim = hpipm_ocp_qp_dim(N)

	if not x0_elim:
		dim.set('nx', nx, 0)    # number of states
		dim.set('nbx', nbx, 0)  # number of state bounds

	dim.set('nx', nx, 1, N) # number of states
	dim.set('nu', nu, 0)    # number of inputs
	dim.set('nbu', nbu, 0)  # number of input bounds

	dim.set('nu', 2*ns, N)
	dim.set('nbu', 2*ns, N)
	# dim.set('nbx', nbx, N)
	dim.set('ng', 2*ns, N)  # number of inequality constraints


	# print to shell
	print("HPIPM print C dims")
	dim.print_C_struct()

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
	b = np.array([0, 0]).reshape(nx,1)

	Q = np.array([1, 0, 0, 1]).reshape(nx,nx)
	# S = np.array([0, 0]).reshape(nu,nx)
	R = np.array([1]).reshape(nu,nu)
	# q = np.array([1, 1]).reshape(nx,1)
	# r = np.array([0]).reshape(nu,1)

	C = np.vstack([np.eye(nx), -np.eye(nx)])
	D = np.block([[np.eye(nx), 0*np.eye(nx)],
				[0*np.eye(nx), np.eye(nx)]])

	Jx = np.array([1, 0, 0, 1]).reshape(nbx,nx)
	x0 = np.array([1, 1]).reshape(nx,1)
	# Jsx = np.array([1, 0, 0, 1]).reshape(nbx,ns)

	Zl = 0 * np.array([1e2, 1e2]).reshape(2,1)
	Zu = 0 * np.array([1e2, 1e2]).reshape(2,1)
	zl = 1 * np.array([1e2, 1e2]).reshape(2,1)
	zu = 1 * np.array([1e2, 1e2]).reshape(2,1)

	Ju = np.array([1]).reshape(nbu,nu)
	lbu = np.array([-1.0]).reshape(nbu,1)
	ubu = np.array([1.0]).reshape(nbu,1)

	# lbs = np.array([0, 0]).reshape(ns,1)
	# ubs = np.array([1.0]).reshape(nbu,1)

	# qp
	qp = hpipm_ocp_qp(dim)

	qp.set('B', B, 0, N-1)

	if x0_elim:
		b0 = A @ x0 + b
		qp.set('b', b0, 0)
	else:
		qp.set('A', A, 0)
		qp.set('b', b, 0)

	qp.set('A', A, 1, N-1)
	qp.set('b', b, 1, N-1)

	if not x0_elim:
		qp.set('Q', Q, 0)

	qp.set('Q', Q, 1, N)
	qp.set('R', R, 0)

	if not x0_elim:
		qp.set('Jx', Jx, 0)
		qp.set('lx', x0, 0)
		qp.set('ux', x0, 0)

	qp.set('Ju', Ju, 0)
	qp.set('lbu', lbu, 0)
	qp.set('ubu', ubu, 0)

	# Time instant 1
	qp.set('R', np.diag(np.concatenate([Zl.squeeze(), Zu.squeeze()])), N)
	qp.set('r', np.concatenate([zl, zu]), N)

	# qp.set('Jx', Jx, N)
	# qp.set('lx', -bigM, N)
	# qp.set('ux', bigM, N)

	qp.set('Ju', np.eye(2*ns), N)
	qp.set('lbu', np.zeros((2*ns,1)), N)
	qp.set('ubu', bigM*np.ones((2*ns,1)), N)
	qp.set('ubu_mask', 0*np.ones((2*ns,1)), N)

	qp.set('C', C, N)
	qp.set('D', D, N)
	# qp.set('lg', np.array([0, 0, -bigM, -bigM]).reshape(2*nx, 1), N)
	# qp.set('ug', np.array([bigM, bigM, 0, 0]).reshape(2*nx, 1), N)
	qp.set('lg', np.array([0, 0, 0, 0]).reshape(2*ns, 1), N)
	qp.set('ug', np.array([bigM, bigM, bigM, bigM]).reshape(2*ns, 1), N)
	qp.set('ug_mask', 0*np.ones((2*ns, 1)), N)

	# print to shell
	print("QP formulation HPIPM C print")
	qp.print_C_struct()

	# qp sol
	qp_sol = hpipm_ocp_qp_sol(dim)

	# set up solver arg
	# mode = 'speed_abs'
	mode = 'speed'
	# mode = 'balance'
	# mode = 'robust'

	# create and set default arg based on mode
	arg = hpipm_ocp_qp_solver_arg(dim, mode)

	# create and set default arg based on mode
	arg.set('mu0', 1e4)
	arg.set('iter_max', 30)
	arg.set('tol_stat', 1e-10)
	arg.set('tol_eq', 1e-10)
	arg.set('tol_ineq', 1e-10)
	arg.set('tol_comp', 1e-10)
	arg.set('reg_prim', 1e-12)

	# codegen
	if codegen_data:
		dim.codegen('ocp_qp_data.c', 'w')
		qp.codegen('ocp_qp_data.c', 'a')
		arg.codegen('ocp_qp_data.c', 'a')

	# set up solver
	solver = hpipm_ocp_qp_solver(dim, arg)

	# solve qp
	start_time = time.time()
	solver.solve(qp, qp_sol)
	end_time = time.time()

	if not travis_run:
		print('solve time {:e}'.format(end_time-start_time))
		print("HPIPM solution C print:")
		qp_sol.print_C_struct()

		# inputs
		print('u =')
		u_traj = qp_sol.get('u', 0, N)
		for u in u_traj:
			print(u)

		# states
		print('x =')
		x_traj = qp_sol.get('x', 0, N)
		for x in x_traj:
			print(x)

		# slack of lower constraints
		print('sl =')
		sl_traj = qp_sol.get('sl', 0, N)
		for sl in sl_traj:
			print(sl)

		# slack of upper constraints
		print('su =')
		su_traj = qp_sol.get('su', 0, N)
		for su in su_traj:
			print(su)


	# get solver statistics
	status = solver.get('status')
	# res_stat = solver.get('max_res_stat')
	# res_eq = solver.get('max_res_eq')
	# res_ineq = solver.get('max_res_ineq')
	# res_comp = solver.get('max_res_comp')
	# iters = solver.get('iter')

	if not travis_run:
		solver.print_stats()

	if status==0:
		print('\nsuccess!\n')
	else:
		print('\nSolution failed, solver returned status {0:1d}\n'.format(status))

	return status


if __name__ == '__main__':

	check_env()
	travis_run = is_travis_run()

	status = main(travis_run)
	sys.exit(int(status))