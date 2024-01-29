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
codegen_data = 1; # export qp data in the file dense_qp_data.c for use from C examples
warm_start = 0; # set to 1 to warm-start the primal variable


# dim
nv = 3  # number of variables
ne = 1  # number of equality constraints
nb = 3  # number of box constraints
ng = 0  # number of general (inequality) constraints

dim = hpipm_dense_qp_dim()
dim.set('nv', nv)
dim.set('nb', nb)
dim.set('ne', ne)
dim.set('ng', ng)

# codegen
if codegen_data:
	dim.codegen('dense_qp_data.c', 'w')

# data
H = np.array([[ 65., -22., -16.], [-22., 14., 7.], [-16., 7., 5.]])
g = np.array([-13.,  15.,   7.])
A = np.array([[1., 0, 0]])
b = np.array([-0.5])
lb = np.array([-0.5, -2, -0.8])

# qp
qp = hpipm_dense_qp(dim)
qp.set('H', H)
qp.set('g', g)
qp.set('A', A)
qp.set('b', b)
qp.set('idxb', np.array([0, 1, 2]))
qp.set('lb', lb)
qp.set('ub_mask', np.zeros(3))

# codegen
if codegen_data:
	qp.codegen('dense_qp_data.c', 'a')

arg = hpipm_dense_qp_solver_arg(dim, 'robust')
qp_sol = hpipm_dense_qp_sol(dim)
solver = hpipm_dense_qp_solver(dim, arg)
solver.solve(qp, qp_sol)

# codegen
if codegen_data:
	arg.codegen('dense_qp_data.c', 'a')

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



sys.exit(int(status))
