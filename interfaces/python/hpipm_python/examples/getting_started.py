from hpipm_python import *
import numpy as np

qp_data = hpipm_data()

qp_data.A = np.array([1, 0, 1, 1])
qp_data.B = np.array([0, 1])
qp_data.b = np.array([0, 0])

qp_data.Q = np.array([1, 0, 0, 1])
qp_data.S = np.array([0, 0])
qp_data.R = np.array([1])
qp_data.q = np.array([1, 1])
qp_data.r = np.array([0])

qp_data.x0 = np.array([1, 1])

qp_dims = hpipm_dims()

qp_dims.idxb0 = np.array([1, 2])

qp_dims.nx = np.array([2, 2, 2, 2, 2, 2])
qp_dims.nu = np.array([1, 1, 1, 1, 1, 0])
qp_dims.nb = np.array([2, 0, 0, 0, 0, 0])
qp_dims.ng = np.array([0, 0, 0, 0, 0, 0])
qp_dims.ns = np.array([0, 0, 0, 0, 0, 0])
qp_dims.nbx = np.array([2, 0, 0, 0, 0, 0])
qp_dims.nbu = np.array([0, 0, 0, 0, 0, 0])
qp_dims.N = 10

solver = hpipm_solver(qp_dims, qp_data)

