from hpipm_python import *
import numpy as np

N = 10

nx = np.array([2, 2, 2, 2, 2, 2])
nu = np.array([1, 1, 1, 1, 1, 0])
nb = np.array([2, 0, 0, 0, 0, 0])
ng = np.array([0, 0, 0, 0, 0, 0])
ns = np.array([0, 0, 0, 0, 0, 0])
nbx = np.array([2, 0, 0, 0, 0, 0])
nbu = np.array([0, 0, 0, 0, 0, 0])

solver = hpipm_solver(N, nx, nu, nbx, nbu, ng, ns)

