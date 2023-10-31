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


def main():

	# check that env.sh has been run
	env_run = os.getenv('ENV_RUN')
	if env_run != 'true':
		print('ERROR: env.sh has not been sourced! Before executing this example, run:')
		print('source env.sh')
		sys.exit(1)

	travis_run = os.getenv('TRAVIS_RUN')

	# flags
	codegen_data = 1 # export qp data in the file dense_qp_data.c for use from C examples
	print_structs = 0
	warm_start = 0 # set to 1 to warm-start the primal variable

	# dim
	nv = 3  # number of variables
	ne = 1  # number of equality constraints
	nb = 0  # number of box constraints
	ng = 3  # number of general (inequality) constraints

	dim = hpipm_dense_qp_dim()

	dim.set('nv', nv)
	dim.set('nb', nb)
	dim.set('ne', ne)
	dim.set('ng', ng)

	# this data is taken the from qpsolvers example
	# <https://github.com/qpsolvers/qpsolvers/blob/6f33bbc4ed797d92be1b88aa493a4398a8d8f73f/examples/quadratic_program.py>
	M = np.array([[1.0, 2.0, 0.0], [-8.0, 3.0, 2.0], [0.0, 1.0, 1.0]])
	H = np.dot(M.T, M)
	g = np.dot(np.array([3.0, 2.0, 3.0]), M)

	C = np.array([[1.0, 2.0, 1.0], [2.0, 0.0, 1.0], [-1.0, 2.0, -1.0]])
	d = np.array([3.0, 2.0, -2.0])
	A = np.array([[1.0, 1.0, 1.0]])
	b = np.array([1.0])

	# QP
	qp = hpipm_dense_qp(dim)

	qp.set('H', H)
	qp.set('g', g)
	qp.set('C', C)
	qp.set('ug', d)
	qp.set('lg', d)  # arbitrary
	qp.set('lg_mask', np.zeros(3))  # disable lower bound
	qp.set('A', A)
	qp.set('b', b)

	# qp sol
	qp_sol = hpipm_dense_qp_sol(dim)

	# set up solver arg
	# mode = 'speed_abs'
	mode = 'speed'
	# mode = 'balance'
	# mode = 'robust'

	# create and set default arg based on mode
	arg = hpipm_dense_qp_solver_arg(dim, mode)

	arg.set('mu0', 1e4)
	arg.set('iter_max', 30)
	arg.set('tol_stat', 1e-4)
	arg.set('tol_eq', 1e-5)
	arg.set('tol_ineq', 1e-5)
	arg.set('tol_comp', 1e-5)
	arg.set('reg_prim', 1e-12)

	# if warm_start=1, then the primal variable is initialized from qp_sol
	arg.set('warm_start', warm_start)
	qp_sol.set('v', np.array([0.30769244, -0.6923054, 1.38461296]))

	# set up solver
	solver = hpipm_dense_qp_solver(dim, arg)

	start_time = time.time()
	solver.solve(qp, qp_sol)
	end_time = time.time()

	if travis_run != 'true':
		print('solve time {:e}'.format(end_time-start_time))

	v = qp_sol.get('v')
	pi = qp_sol.get('pi')
	lam_lb = qp_sol.get('lam_lb')
	lam_ub = qp_sol.get('lam_ub')
	lam_lg = qp_sol.get('lam_lg')
	lam_ug = qp_sol.get('lam_ug')

	print('v      = {}'.format(v.flatten()))
	print('pi     = {}'.format(pi.flatten()))
	print('lam_lb = {}'.format(lam_lb.flatten()))
	print('lam_ub = {}'.format(lam_ub.flatten()))
	print('lam_lg = {}'.format(lam_lg.flatten()))
	print('lam_ug = {}'.format(lam_ug.flatten()))

	# get solver statistics
	status = solver.get('status')
	# res_stat = solver.get('max_res_stat')
	# res_eq = solver.get('max_res_eq')
	# res_ineq = solver.get('max_res_ineq')
	# res_comp = solver.get('max_res_comp')
	# iters = solver.get('iter')

	if travis_run != 'true':
		solver.print_stats()

	if status==0:
		print('\nsuccess!\n')
	else:
		print('\nSolution failed, solver returned status {0:1d}\n'.format(status))

	if print_structs:
		dim.print_C_struct()
		qp.print_C_struct()

	if codegen_data:
		dim.codegen('dense_qp_data.c', 'w')
		qp.codegen('dense_qp_data.c', 'a')
		arg.codegen('dense_qp_data.c', 'a')

	return status


if __name__ == '__main__':
	status = main()
	sys.exit(int(status))
