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

TOL = 1e-12

def main(travis_run):

    N = 5
    nx = 3
    nu = 2

    dim = hpipm_ocp_qp_dim(N)

    dim.set('nx', nx, 0, N) # number of states
    dim.set('nu', nu, 0, N-1) # number of inputs
    dim.set('nbx', nx, 0)

    A = np.diag([1., -0.9, 2])
    B = np.arange(nx*nu).reshape((nx, nu))

    Q = np.eye(nx)
    S = np.zeros((nu, nx))
    R = 0.1*np.eye(nu)
    q = -0.001*np.ones((nx,1))

    x0 = 5*np.ones((nx, 1))

    # QP data
    qp = hpipm_ocp_qp(dim)

    qp.set('A', A, 0, N-1)
    qp.set('B', B, 0, N-1)
    qp.set('Q', Q, 0, N)
    qp.set('S', S, 0, N-1)
    qp.set('R', R, 0, N-1)
    qp.set('q', q, 0, N)

    qp.set('idxbx', np.arange(nx), 0)
    qp.set('lbx', x0, 0)
    qp.set('ubx', x0, 0)

    # qp solution
    qp_sol = hpipm_ocp_qp_sol(dim)

    # set up solver arg
    mode = 'robust'
    arg = hpipm_ocp_qp_solver_arg(dim, mode)

    solver = hpipm_ocp_qp_solver(dim, arg)

    # solve qp
    start_time = time.time()
    solver.solve(qp, qp_sol)
    end_time = time.time()

    status = solver.get('status')

    # extract and print sol
    u_traj = qp_sol.get('u', 0, N-1)
    x_traj = qp_sol.get('x', 0, N)
    pi_traj = qp_sol.get('pi', 0, N-1)

    lam_lbx = qp_sol.get('lam_lbx', 0)
    lam_ubx = qp_sol.get('lam_ubx', 0)
    pi0 = lam_lbx - lam_ubx

    # extract linear feedback law
    K_traj = solver.get_feedback(qp, 'ric_K', 0, N-1)
    k_traj = solver.get_feedback(qp, 'ric_k', 0, N-1)

    Ls_traj = solver.get_feedback(qp, 'ric_Ls', 0, N-1)

    Lr_traj = solver.get_feedback(qp, 'ric_Lr', 0, N-1)
    lr_traj = solver.get_feedback(qp, 'ric_lr', 0, N-1)

    P_traj = solver.get_feedback(qp, 'ric_P', 0, N)
    p_traj = solver.get_feedback(qp, 'ric_p', 0, N)

    print("\ncheck u")
    for n in range(N):
        diff_u = k_traj[n] + K_traj[n] @ x_traj[n] - u_traj[n]
        max_diff_u = np.linalg.norm(diff_u, ord=np.inf)
        print(max_diff_u)
        assert max_diff_u <= TOL

    print("\ncheck multiplier")
    # NOTE the first multiplier is a bit off,
    # this is expected (?) due to the initial state constraints being treated as upper and lower bounds
    diff_multiplier = p_traj[0] + P_traj[0] @ x_traj[0] - pi0

    max_diff_multiplier = np.linalg.norm(diff_multiplier, ord=np.inf)
    print(max_diff_multiplier)

    for n in range(1, N):
        diff_multiplier = p_traj[n] + P_traj[n] @ x_traj[n] - pi_traj[n-1]
        max_diff_multiplier = np.linalg.norm(diff_multiplier, ord=np.inf)
        print(max_diff_multiplier)
        assert max_diff_multiplier <= TOL


    if not travis_run:

        solver.print_stats()
        print('solve time {:e}'.format(end_time-start_time))

        for n in range(N):
            print(f'\nu_{n} =')
            print(u_traj[n].flatten())

        for n in range(N+1):
            print(f'\nx_{n} =')
            print(x_traj[n].flatten())

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