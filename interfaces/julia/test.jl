# ENV["PYTHON"] = /usr/bin/python3
# Pkg.build("PyCall")

using PyCall
@pyimport hpipm_python as hp
@pyimport numpy as np

include("py2ju_obj.jl")
include("ju2py_obj.jl")
include("hpipm_solver.jl")

qp_data = hpipm_data()

A = [1.0 0; 1 1]
B = [0, 1]
b = [0, 0]

Q = [1 0; 0 1]
S = [0, 0]
R = [1]
q = [1, 1]
r = [0]

qp_data["A"] = [A, A, A, A, A]
qp_data["B"] = [B, B, B, B, B]
qp_data["b"] = [b, b, b, b, b]
qp_data["Q"] = [Q, Q, Q, Q, Q, Q]
qp_data["S"] = [S, S, S, S, S, S]
qp_data["R"] = [R, R, R, R, R, R]
qp_data["q"] = [q, q, q, q, q, q]
qp_data["r"] = [r, r, r, r, r, r]

x0 = [1, 1]

qp_data["d_lb"] = [x0]
qp_data["d_ub"] = [x0]

qp_data["idxb"] = [[1, 2]]
py_qp_data = hp.hpipm_data()
py_qp_data = ju2py_obj(qp_data, py_qp_data) 
#qp_dims = hp.hpipm_dims()
#
#qp_dims[:nx]   = np.array([2, 2, 2, 2, 2, 2])
#qp_dims[:nu]   = np.array([1, 1, 1, 1, 1, 0])
#qp_dims[:nb]   = np.array([2, 0, 0, 0, 0, 0])
#qp_dims[:nbx]  = np.array([2, 0, 0, 0, 0, 0])
#qp_dims[:nbu]  = np.array([0, 0, 0, 0, 0, 0])
#qp_dims[:ng]   = np.array([0, 0, 0, 0, 0, 0])
#qp_dims[:ns]   = np.array([0, 0, 0, 0, 0, 0])
#qp_dims[:nsbx] = np.array([0, 0, 0, 0, 0, 0])
#qp_dims[:nsbu] = np.array([0, 0, 0, 0, 0, 0])
#qp_dims[:nsg]  = np.array([0, 0, 0, 0, 0, 0])
#qp_dims[:N]    = 5
#
## set up solver
#solver = hp.hpipm_solver(qp_dims, qp_data)
#
## solve qp
#return_flag = solver[:solve]()
#
#println("HPIPM returned with flag $(return_flag)\n")
#
#if return_flag == 0
#    println("-> QP solved Solution:\n")
#    solver[:print_sol]()
#else
#    println("-> Solver failed.\n")
#end
