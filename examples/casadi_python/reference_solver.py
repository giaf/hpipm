import casadi as ca
import numpy as np

class reference_solver:
    def __init__(self, nv, x, g_exp, h, H):
        self.x = np.zeros((nv,1))

        hp = ca.MX.sym('hp', nv, 1)
        f_exp = 1.0/2.0*ca.mtimes(ca.mtimes(ca.transpose(x), H), x) \
            + ca.mtimes(ca.transpose(hp), x)
        self.h = h 

        self.ipopt_solver = None

        # create reference solver 
        w = x
        w0 = np.zeros(nv)
        lbw = []
        ubw = []
        J = f_exp 
        g = g_exp
        self.lbg = -1e12*np.ones((int(nv/2),1))
        self.ubg = np.zeros((int(nv/2),1))
        
        p = ca.vertcat(hp)

        prob = {'f': J, 'x': w, 'g': g, 'p': p}
        opts = {'ipopt': {'print_level': 3, 'dual_inf_tol':1e-10, 'constr_viol_tol': 1e-10, 'compl_inf_tol': 1e-10}}
        self.ipopt_solver = ca.nlpsol('solver', 'ipopt', prob, opts)

    def solve(self):
        # call Ipopt
        print(self.h)
        sol = self.ipopt_solver(p=self.h, lbg=self.lbg, ubg=self.ubg)

        stats = self.ipopt_solver.stats()
        if stats['return_status']  != 'Solve_Succeeded':
            raise Exception('reference solver failed!')
        self.x = sol['x']
        print(self.x)
        
    def update_gradient(self, h):
        self.h = h
