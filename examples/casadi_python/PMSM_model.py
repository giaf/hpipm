import casadi as ca
import numpy as np
import scipy as sp
from scipy import signal

def PMSM_dynamics(i_d, i_q, u_d, u_q, par, w):
    # parameters
    p = par.p 
    theta = par.theta 
    Rs = par.Rs 
    Ld = par.Ld 
    Lq = par.Lq 
    psi_pm = par.psi_pm 
    m_load = par.m_load
    u_max = par.u_max

    xdot = ca.vertcat(-Rs/Ld*i_d + Lq/Ld*w*i_q + u_d/Ld, \
        -Rs/Lq*i_q - Ld/Lq*w*i_d - psi_pm/Lq*w + u_q/Lq)

    return xdot

class PMSM_par:
    def __init__(self):
        self.p = 3
        self.theta = 342e-4
        self.Rs = 0.11
        self.Ld = 3.35e-3
        self.Lq = 1.2*self.Ld
        self.psi_pm = 0.3765
        self.w_ref = 3*209.4
        self.m_load = 25
        self.u_max = 290
        self.J = np.array([[0, -1], [1, 0]])

def PMSM_ode_ss(x, w, par):
    # parameters
    p = par.p 
    theta = par.theta 
    Rs = par.Rs 
    Ld = par.Ld 
    Lq = par.Lq 
    psi_pm = par.psi_pm 
    m_load = par.m_load
    u_max = par.u_max

    # variables
    i_d = x[0]
    i_q = x[1]
    u_d = x[2]
    u_q = x[3]
    
    xdot = PMSM_dynamics(i_d, i_q, u_d, u_q, par, w)
    return xdot

class ss_calculator():
    def __init__(self, w, Ts, par):
        # parameters
        p = par.p 
        theta = par.theta 
        Rs = par.Rs 
        Ld = par.Ld 
        Lq = par.Lq 
        psi_pm = par.psi_pm 
        m_load = par.m_load
        u_max = par.u_max
        J_par = par.J

        i_d = ca.MX.sym('i_d', 1, 1)
        i_q = ca.MX.sym('i_q', 1, 1)
        u_d = ca.MX.sym('u_d', 1, 1)
        u_q = ca.MX.sym('u_q', 1, 1)

        x = ca.vertcat(i_d, i_q)
        u = ca.vertcat(u_d, u_q)

        xdot = PMSM_dynamics(i_d, i_q, u_d, u_q, par, w)

        # create solver steady-state computation 
        w = ca.vertcat(x, u)
        w0 = np.ones(4)
        lbw = []
        ubw = []
        J = 1/2*(i_d**2 + i_q**2)
        g = ca.vertcat(xdot, 3/2*p*ca.mtimes(ca.mtimes(ca.horzcat(i_d, i_q), \
            J_par), ca.vertcat(psi_pm, psi_pm)))
        
        prob = {'f': J, 'x': w, 'g': g}
        opts = {'ipopt': {'print_level': 0}}
        self.ss_solver = ca.nlpsol('solver', 'ipopt', prob, opts);
        xdot_fun = ca.Function('xdot_fun', [ca.vertcat(i_d, i_q, u_d, u_q)], [xdot])
        self.xdot_fun = xdot_fun
    
    def compute(self, torque_ref): 
        lbg = np.zeros((3,1))
        ubg = np.zeros((3,1))
        lbg[2] = torque_ref
        ubg[2] = torque_ref
        sol = self.ss_solver(lbg = lbg, ubg = ubg)
        stats = self.ss_solver.stats()
        if stats['return_status']  != 'Solve_Succeeded':
            raise Exception('steady-state calculation failed!')

        # check steady-state
        err = np.linalg.norm(self.xdot_fun(sol['x']))
        if err > 1e-6:
            raise Exception('Solution found does not seem to be a steady-state solution')
        return sol['x']

def get_PMSM_dynamics(w, Ts, par, Qc=None, Rc=None):
    # parameters
    p = par.p 
    # import pdb; pdb.set_trace()
    theta = par.theta 
    Rs = par.Rs 
    Ld = par.Ld 
    Lq = par.Lq 
    psi_pm = par.psi_pm 
    m_load = par.m_load
    u_max = par.u_max
    J_par = par.J
    
    i_d = ca.MX.sym('i_d', 1, 1)
    i_q = ca.MX.sym('i_q', 1, 1)
    u_d = ca.MX.sym('u_d', 1, 1)
    u_q = ca.MX.sym('u_q', 1, 1)

    x = ca.vertcat(i_d, i_q)
    u = ca.vertcat(u_d, u_q)

    xdot = PMSM_dynamics(i_d, i_q, u_d, u_q, par, w)

    # fixed step Runge-Kutta 4 integrator
    if Qc is None or Rc is None:
        M  = 100 # RK4 steps per interval
        DT = Ts/M
        f  = ca.Function('f', [x, u], [xdot])
        X0 = ca.MX.sym('X0', 2, 1)
        U  = ca.MX.sym('U', 2, 1)
        X  = X0

        for j in range(M):
            k1 = f(X, U)
            k2 = f(X + DT/2 * k1, U)
            k3 = f(X + DT/2 * k2, U)
            k4 = f(X + DT * k3, U)
            X=X+DT/6*(k1 +2*k2 +2*k3 +k4)

        # x_{k+1} = Ax_k + Bu_k
        A_exp = ca.jacobian(X, X0)
        B_exp = ca.jacobian(X, U)
        A_fun = ca.Function('A_fun', [X0, U], [A_exp])
        B_fun = ca.Function('B_fun', [X0, U], [B_exp])

        x_plus = ca.Function('x_plus', [X0, U], [X])

        A_exp_c = ca.jacobian(xdot, x)
        B_exp_c = ca.jacobian(xdot, u)
        A_fun_c = ca.Function('A_fun', [x, u], [A_exp_c])
        B_fun_c = ca.Function('B_fun', [x, u], [B_exp_c])
        A_c = A_fun_c(np.zeros((2,1)), np.zeros((2,1)))
        B_c = B_fun_c(np.zeros((2,1)), np.zeros((2,1)))
        Ad, Bd, Cd, Dd, dt = sp.signal.cont2discrete((A_c.full(),B_c.full(), \
            np.zeros(2), np.zeros(2)), dt = Ts)

        A = A_fun(np.zeros((2,1)), np.zeros((2,1)))
        B = B_fun(np.zeros((2,1)), np.zeros((2,1)))
        c = x_plus(np.zeros((2,1)), np.zeros((2,1))).full() 
        return A, B, c, x_plus
    else:
        M  = 100 # RK4 steps per interval
        DT = Ts/M
        X0    = ca.MX.sym('X0', 2, 1)
        X_ref = ca.MX.sym('X_ref', 2, 1)
        U  = ca.MX.sym('U', 2, 1)
        U_ref  = ca.MX.sym('U_ref', 2, 1)
        X  = X0
        stage_cost = 1/2*ca.mtimes(ca.mtimes((X_ref - x).T, Qc), \
            (X_ref - x)) + 1/2*ca.mtimes(ca.mtimes((U_ref - u).T, Rc), (U_ref - u))

        f = ca.Function('f', [x, u], [xdot, stage_cost])
        J = 0

        for j in range(M):
            k1, k1_j = f(X, U)
            k2, k2_j  = f(X + DT/2 * k1, U)
            k3, k3_j = f(X + DT/2 * k2, U)
            k4, k4_j = f(X + DT * k3, U)
            X=X+DT/6*(k1 +2*k2 +2*k3 +k4)
            J = J + DT/6*(k1_j + 2*k2_j + 2*k3_j + k4_j)

        x_plus = ca.Function('x_plus', [X0, U, X_ref, U_ref], [X, J])

        return x_plus
