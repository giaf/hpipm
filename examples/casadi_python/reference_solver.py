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
import os


class reference_solver:

	def __init__(self, x, f_exp, g_exp, p):

		nx = x.size(1)
		ng = g_exp.size(1)

		self.x_sol = np.zeros((nx,1))

		# create reference solver 
#		x0 = np.zeros(nx)
#		lbx = []
#		ubx = []
		J = f_exp 
		g = g_exp
		self.lbg = -1e12*np.ones((ng,1))
		self.ubg = np.zeros((ng,1))
		
		travis_run = os.getenv('TRAVIS_RUN')
		if(travis_run!='true'):
			print_level = 3
		else:
			print_level = 0

		prob = {'f': f_exp, 'x': x, 'g': g_exp, 'p': p}
		opts = {'ipopt': {'print_level': print_level, 'dual_inf_tol':1e-10, 'constr_viol_tol': 1e-10, 'compl_inf_tol': 1e-10}}
		self.ipopt_solver = ca.nlpsol('solver', 'ipopt', prob, opts)

	def solve(self):
		# call Ipopt
#		print(self.h)
		sol = self.ipopt_solver(p=self.h, lbg=self.lbg, ubg=self.ubg)

		stats = self.ipopt_solver.stats()
		if stats['return_status']  != 'Solve_Succeeded':
			raise Exception('reference solver failed!')
		self.x_sol = sol['x']
#		print(self.x_sol)
		
	def update_gradient(self, h):
		self.h = h
