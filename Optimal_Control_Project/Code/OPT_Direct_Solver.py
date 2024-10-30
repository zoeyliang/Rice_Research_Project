from scipy.integrate import solve_ivp

from pyomo.dae import *
from pyomo.environ import *

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from math import pi
from math import nan
from pyomo.environ import exp




# def add(a,c):
#     b = a+c
#     return b

# z = add(x,y)










def hypersonic_problem(RandNo):
    # test control for simulation: 11 degrees AOA
    def alpha_test(t):
        return rad(11)

    # generating simulation using scipy
    t0     = 0
    tfinal = 3000

    x_init = 0
    y_init = 80
    v_init = 5
    gamma_init = rad(-5)

    y_final = 0

    # used for initial simulation
    def hit_ground(t, state):
        return state[1] - y_final
    hit_ground.terminal = True

    # 2D flight dynamics with flap rate control
    def flight2D(t, state, control):
        x = state[0]
        y = state[1]
        v = state[2]
        gamma = state[3]

        alpha = control
        
        result = np.zeros(4)
        result[0] = v * cos(gamma)
        result[1] = v * sin(gamma)
        result[2] = -1 / mass * (D(y, v, alpha, RandNo) + mass * g(y) * sin(gamma))
        result[3] = 1 / (mass * v) * (L(y, v, alpha, RandNo) 
                                    - mass * g(y) * cos(gamma) 
                                    + mass * (v ** 2 * cos(gamma)) / (R_E + y) )
        
        return result

    sol  = solve_ivp(lambda t, y: flight2D(t, y, alpha_test(t)), [t0, tfinal], [x_init, y_init, v_init, gamma_init],
                                events=hit_ground, max_step=1, dense_output=True)
    
    # time
    t_sim_pts = sol.t

    # duration & normalized time
    T_sim = t_sim_pts[-1]
    tau_sim_pts = t_sim_pts / T_sim

    # discrete states
    x_sim_pts = sol.y[0, :]
    y_sim_pts = sol.y[1, :]
    v_sim_pts = new_func(sol)
    gamma_sim_pts = sol.y[3, :]

    # discrete controls
    alpha_sim_pts = np.array([alpha_test(t) for t in sol.t])

    #--------------------------------------------------------------------------------------------------
    # continuous states in normalized time (linear interpolation)
    def x_sim(tau):
        return np.interp(tau, tau_sim_pts, x_sim_pts)

    def y_sim(tau):
        return np.interp(tau, tau_sim_pts, y_sim_pts)

    def v_sim(tau):
        return np.interp(tau, tau_sim_pts, v_sim_pts)

    def gamma_sim(tau):
        return np.interp(tau, tau_sim_pts, gamma_sim_pts)

    # continuous controls in normalized time (linear interpolation)
    def alpha_sim(tau):
        return np.interp(tau, tau_sim_pts, alpha_sim_pts)
    #--------------------------------------------------------------------------------------------------

    #--------------------------------------------------------------------------------------------------
    # state derivatives
    def dx_sim(tau):
        return T_sim * v_sim(tau) * cos(gamma_sim(tau))

    def dy_sim(tau):
        return T_sim * v_sim(tau) * sin(gamma_sim(tau))

    def dv_sim(tau):
        return -T_sim / mass * (D(y_sim(tau), v_sim(tau), alpha_sim(tau), RandNo) + mass * g(y_sim(tau)) * sin(gamma_sim(tau)))

    def dgamma_sim(tau):
        return T_sim / (mass * v_sim(tau)) * (L(y_sim(tau), v_sim(tau), alpha_sim(tau), RandNo) 
                                        - mass * g(y_sim(tau)) * cos(gamma_sim(tau)) 
                                        + mass * v_sim(tau) ** 2 * cos(gamma_sim(tau)) / (R_E + y_sim(tau)))
    #--------------------------------------------------------------------------------------------------


    #--------------------------------------------------------------------------------------------------
    # creating model
    m = ConcreteModel()

    v_min = 1e-3

    # noise
    m.RandNo = RandNo

    # normalized time horizon
    m.tau = ContinuousSet(bounds = (0, 1))

    # duration
    m.T = Var(within=NonNegativeReals, initialize = T_sim)

    # state variables + derivatives
    m.x = Var(m.tau, within=NonNegativeReals, initialize = lambda m, tau: x_sim(tau))
    m.dx = DerivativeVar(m.x, wrt = m.tau, initialize = lambda m, tau: dx_sim(tau))

    m.y = Var(m.tau, initialize = lambda m, tau: y_sim(tau))
    m.dy = DerivativeVar(m.y, wrt = m.tau, initialize = lambda m, tau: dy_sim(tau))

    m.v = Var(m.tau, bounds = (v_min, None), initialize = lambda m, tau: v_sim(tau))
    m.dv = DerivativeVar(m.v, wrt = m.tau, initialize = lambda m, tau: dv_sim(tau))

    m.gamma = Var(m.tau, bounds = (gamma_min, gamma_max), initialize = lambda m, tau: gamma_sim(tau))
    m.dgamma = DerivativeVar(m.gamma, wrt = m.tau, initialize = lambda m, tau: dgamma_sim(tau))

    # control
    m.alpha = Var(m.tau, bounds = (alpha_min, alpha_max), initialize = lambda m, tau: alpha_sim(tau))

    # dynamics constraints
    m.ode_x = Constraint(m.tau, rule = lambda m, tau: 
        m.dx[tau] == m.T * m.v[tau] * cos(m.gamma[tau]))
    m.ode_x[m.tau.first()].deactivate()
                        
    m.ode_y = Constraint(m.tau, rule = lambda m, tau: 
        m.dy[tau] == m.T * m.v[tau] * sin(m.gamma[tau]))
    m.ode_y[m.tau.first()].deactivate()
                        
    m.ode_v = Constraint(m.tau, rule = lambda m, tau: 
        m.dv[tau] == - m.T / mass * (D(m.y[tau], m.v[tau], m.alpha[tau], RandNo) + mass * g(m.y[tau]) * sin(m.gamma[tau])))
    m.ode_v[m.tau.first()].deactivate()

    m.ode_gamma = Constraint(m.tau, rule = lambda m, tau: 
        m.dgamma[tau] == m.T / (mass * m.v[tau]) * (L(m.y[tau], m.v[tau], m.alpha[tau], RandNo) - mass * g(m.y[tau]) * cos(m.gamma[tau])
                                + mass * m.v[tau] ** 2 * cos(m.gamma[tau]) / (R_E + m.y[tau])))
    m.ode_gamma[m.tau.first()].deactivate()
    

    # initial conditions
    x_init = 0
    y_init = 80
    v_init = 5
    gamma_init = rad(-5)

    # final conditions
    y_final = 0

    # initial conditions
    m.x[0].fix(x_init)
    m.y[0].fix(y_init)
    m.v[0].fix(v_init)
    m.gamma[0].fix(gamma_init)

    # final conditions
    m.y[1].fix(y_final)

    # dynamic pressure constraint
    m.dynamic_pressure = Constraint(m.tau, rule = lambda m, tau: 
        q_bar(m.y[tau], m.v[tau]) <= q_bar_max)

    # stagnation heat rate constraint
    m.stag_heat_rate = Constraint(m.tau, rule = lambda m, tau: 
        dQ_stag(m.y[tau], m.v[tau]) <= dQ_stag_max)

    # objective function
    def _obj(m):
        return -m.x[1]
    m.obj = Objective(rule=_obj)

    return m
    #--------------------------------------------------------------------------------------------------

    
def new_func(sol):
    v_sim_pts = sol.y[2, :]
    return v_sim_pts

def solve(m, nfe=100, max_iter=None, tee=True):
    # discretize the problem using collocation on 100 subintervals and 3 collocation points per subinterval
    discretizer = TransformationFactory('dae.collocation')
    discretizer.apply_to(m, nfe=nfe, ncp=3, scheme='LAGRANGE-RADAU')

    # solve model
    # solver = SolverFactory('ipopt')
    solver = SolverFactory('ipopt', executable="/Applications/anaconda3/bin/ipopt")
    solver.set_executable("/Applications/anaconda3/bin/ipopt")

    
    if(max_iter is not None):
        solver.options['max_iter'] = max_iter
    results = solver.solve(m, tee=tee)
    # results = 1
    return results

# wrapper class for optimal solution of control problem
class problem_results:
    def __init__(self, m):
        self.t_opt = np.array([m.T() * tau for tau in m.tau])
        self.tau_opt = np.array([tau for tau in m.tau])
        self.T_opt = m.T()

        # states
        self.x_opt = np.array([m.x[tau]() for tau in m.tau])
        self.y_opt = np.array([m.y[tau]() for tau in m.tau])
        self.v_opt = np.array([m.v[tau]() for tau in m.tau])
        self.gamma_opt = np.array([m.gamma[tau]() for tau in m.tau])

        # controls
        self.alpha_opt = np.array([m.alpha[tau]() for tau in m.tau if tau != 0])

        # path constraints
        self.q_bar_opt = np.array([q_bar(m.y[tau](), m.v[tau]()) for tau in m.tau])
        self.dQ_stag_opt = np.array([dQ_stag(m.y[tau](), m.v[tau]()) for tau in m.tau])

        # lift, drag, L/D
        self.L_opt = np.array([L(m.y[tau](), m.v[tau](), m.alpha[tau](), m.RandNo) for tau in m.tau if tau != 0])
        self.D_opt = np.array([D(m.y[tau](), m.v[tau](), m.alpha[tau](), m.RandNo) for tau in m.tau if tau != 0])
        self.LD_opt = self.L_opt / self.D_opt



# radians to degrees
def deg(r):
    return 180 / pi * r

# degrees to radians
def rad(d):
    return pi / 180 * d

# physical constants
R_E         = 6371                    # km                          radius of earth
mu          = 3.986e5                 # km^3/s^2                    gravitational parameter

# vehicle-specific constants
mass        = 1200                    # kg                          mass of vehicle
A_w         = 6e-6                    # km^2                        reference area
C           = 37.356                  #                             heat rate constant
q_bar_max   = 40                      # kPa                         maximum dynamic pressure
dQ_stag_max = 6                       # MW/m^2                      maximum stagnation heat rate

# bounds for states/controls
v_min = 1e-3 # km/s

# all angles in radians during optimization, but plotted in degrees

gamma_min = rad(-30) # -30 deg
gamma_max = rad(30) # 30 deg

alpha_min = rad(-10) # -10 deg
alpha_max = rad(30) # 30 deg

# atmospheric density
# y in km, rho in kg/km^3
def rho(y):
    return 1.225e9 * exp(-0.14 * y)

# lift coefficient
# y in km, v in km/s, alpha in radians
# def c_L(y, v, alpha, RandNo):
#     return -0.04*RandNo[0] + 0.8*RandNo[1]*alpha
def c_L(y, v, alpha, RandNo):
    return RandNo[0] + RandNo[1]*alpha

# drag coefficient
# y in km, v in km/s, alpha in radians
# def c_D(y, v, alpha, RandNo):
#     return 0.012*RandNo[2] - 0.01*RandNo[3]*alpha + 0.6*RandNo[4]*alpha**2
def c_D(y, v, alpha, RandNo):
    return RandNo[2] + RandNo[3]*alpha + RandNo[4]*alpha**2

# stagnation heat rate
# y in km, v in km/s, dQ_stag in MW/m^2
def dQ_stag(y, v):
    return 1e-6 * C * rho(y) ** 0.5 * v ** 3.05

# gravitational model
# y in km, g in km/s^2
def g(y):
    return mu / (R_E + y) ** 2

# dynamic pressure
# y in km, v in km/s, q_bar in kPa
def q_bar(y, v):
    return 1e-6 * 0.5 * rho(y) * v ** 2

# lift force
# y in km, v in km/s, alpha in radians, L in kN
def L(y, v, alpha, RandNo):
    return 1e6 * q_bar(y, v) * c_L(y, v, alpha, RandNo) * A_w

# drag force
# y in km, v in km/s, alpha in radians, D in kN
def D(y, v, alpha, RandNo):
    return 1e6 * q_bar(y, v) * c_D(y, v, alpha, RandNo) * A_w



def OPT_Direct_Solver(sample):
    m = hypersonic_problem(sample)
    solve(m)
    results = problem_results(m)
    return results

z = OPT_Direct_Solver(x) # this line is for Matlab usage