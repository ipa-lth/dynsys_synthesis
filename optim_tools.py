#!/usr/bin/env python
import cvxpy
import numpy as np
from numpy import linalg as LA
import collections
import math

def accurate_solve(problem, solver, **kwargs_solver):
    try:
        problem.solve(solver=solver, **kwargs_solver)
    except Exception, e:
        print "Caught Exception: {}".format(str(e))
        pass
    
    while problem.status not in ["infeasible", "unbounded", "optimal"]:
        if solver is not cvxpy.MOSEK:
            if 'max_iters' in kwargs_solver:
                kwargs_solver['max_iters'] += kwargs_solver['max_iters']
            else:
                kwargs_solver['max_iters'] = 5000
            print "Status with {} iterations: {}\nExtending to {} iterations".format(kwargs_solver['max_iters']/2, problem.status, kwargs_solver['max_iters'])

            problem.solve(solver=solver, **kwargs_solver)
        else:
            print "'{}' Solutions with MOSEK not implemented".format(problem.status)
            break
    return problem

def checked_solve(problem, solvers, **kwargs_solver):
    repeat_solve = True
    #legacy check
    if type(solvers) is list:
        for solver, kwargs in solvers:
            if repeat_solve == True:
                #print "Solving with {} with {}".format(solver, *list(kwargs))
                try:
                    problem.solve(solver=solver, **kwargs)
                    repeat_solve = False
                except Exception, e:
                    print " -----> Exception by {}: {}".format(solver, str(e))
                    repeat_solve = True
                else:
                    # Check constraints
                    if 'inaccurate' in problem.status:
                        if 'optimal' in problem.status:
                            print "Violation for optimal_inaccurate"
                            print "Max violation:", max([c.violation for c in problem.constraints])

                        elif 'unbounded' in problem.status:
                            print "Status is unbounded! Possible ERROR in objective?"
                            print "Max violation:", max([c.violation for c in problem.constraints])

                        else: # 'infeasible' in problem.status:
                            print "Violation for infeasible_inaccurate"
                            print "Max violation:", max([c.violation for c in problem.constraints])

                        repeat_solve = False
        return problem, not repeat_solve
    else:
        #legacy solve
        try:
            problem.solve(solver=solvers, **kwargs_solver)
        except Exception, e:
            print " -----> Exception: {}".format(str(e))
            return problem, False

        repeat_solve = False
        # Check constraints
        if 'inaccurate' in problem.status:

            if 'optimal' in problem.status:
                print "Violation for optimal_inaccurate"
                print "Max violation:", max([c.violation for c in problem.constraints])
                repeat_solve = True

            elif 'unbounded' in problem.status:
                print "Status is unbounded! Possible ERROR in objective?"
                print "Max violation:", max([c.violation for c in problem.constraints])

            else: # 'infeasible' in problem.status:
                print "Violation for infeasible_inaccurate"
                print "Max violation:", max([c.violation for c in problem.constraints])
                repeat_solve = True

        return problem, True

'''
## Example bisection code (MATLAB)
## https://see.stanford.edu/materials/lsocoee364a/hw6sol.pdf (p.3)

## Example bisection code
## https://www.coursehero.com/file/p3eovvo/To-solve-the-problem-we-can-use-a-bisection-method-solving-an-LP-feasibility/

import numpy as np
import matplotlib.pyplot as plt
import cvxpy as cvx

k = 201
t = -3 + 6 * np.arange(k) / (k - 1)
y = np.exp(t)

Tpowers = np.vstack((np.ones(k), t, t**2)).T

a = cvx.Variable(3)
b = cvx.Variable(2)
gamma = cvx.Parameter(sign='positive')
lhs = cvx.abs(Tpowers * a - (y[:, np.newaxis] * Tpowers) * cvx.vstack(1, b))
rhs = gamma * Tpowers * cvx.vstack(1, b)
problem = cvx.Problem(cvx.Minimize(0), [lhs <= rhs])

l, u = 0, np.exp(3) # initial upper and lower bounds
bisection_tol = 1e-3 # bisection tolerance
while u - l >= bisection_tol:
    gamma.value = (l + u) / 2
    # solve the feasibility problem
    problem.solve()
    if problem.status == 'optimal':
        u = gamma.value
        a_opt = a.value
        b_opt = b.value
        objval_opt = gamma.value
    else:
        l = gamma.value

y_fit = (Tpowers * a_opt / (Tpowers * np.vstack((1, b_opt)))).A1

plt.figure()
plt.plot(t, y, 'b', t, y_fit, 'r+')
plt.xlabel('t')
plt.ylabel('y')
plt.show()

plt.figure()
plt.plot(t, y_fit - y)
plt.xlabel('t')
plt.ylabel('err')
plt.show()
'''

''' Bisection function
Tries to find the max value of parameter of the given feasibility problem by bisection


l: lower bound with infeasible solution
u: upper bound with optimal solution (if None, see note2!)
bisection_tol: bisection tolerance
sol: solver for problem
problem: cvxpy Problem
variables: array of variables in problems
parameter:

Note: make sure solution with l is infeasible, with u is optimal
Example: [[o_Q, o_z0, o_z1], o_g] = bisect_max(
                                        0, 2,
                                        prob, g,
                                        [Q, z0, z1],
                                        solver=cvxpy.SCS,
                                        verbose=False)


Note2: If upper bound is None, the optimization tries to find it by doubling the value until it is not feasible any more. lower bound is shifted if problem is found feasible, however. Initial u is l+1.0
'''
def bisect_max(l, u, problem, parameter, variables,
           bisection_tol=1e-3, solver=cvxpy.CVXOPT, bisect_verbose=False, **kwargs_solver):

    if solver is cvxpy.SCS:
        print "Setting 'Warmstart' option to True"
        kwargs_solver['warm_start'] = True
            
    if (u is not None):
        # cross check bound
        if (u < l):
            #print "upperBound < lowerBound"
            raise ValueError("upperBound({}) < lowerBound({})".format(u, l))

    elif (u is None) and (l is not None):
        # First iteration
        u = l + 1.0
        if bisect_verbose:
            print "processing upper bound: {}".format(u)
        parameter.value = u
        problem, ok = checked_solve(problem, solver, **kwargs_solver)
        if ok:
            uStatus = problem.status
        else:
            uStatus = 'unknown'
        
        while 'optimal' in uStatus:
            #if u >= l: # shift upper bound if found feasible, this condition is always true
                #l = u
            l = u
            u = 2.0*u
            if bisect_verbose:
                print "processing upper bound: {}".format(u)
            parameter.value = u
            problem, ok = checked_solve(problem, solver, **kwargs_solver)
            if ok:
                uStatus = problem.status
            else:
                uStatus = 'unknown'
        if bisect_verbose:
            print 'found bounds: [{}-{}]'.format(l, u)
    else:
        raise ValueError("Not implemented")
    
    # check validity solution of l is optimal, solution of u is infeasible
    parameter.value = l
    problem, ok = checked_solve(problem, solver, **kwargs_solver)
    if ok:
        lStatus = problem.status
    else:
        lStatus = 'unknown'

    parameter.value = u
    problem, ok = checked_solve(problem, solver, **kwargs_solver)
    if ok:
        uStatus = problem.status
    else:
        uStatus = 'unknown'

    if not ('optimal' in lStatus and uStatus in ['infeasible', 'unknown']):
        #print "UpperBound({})={}, LowerBound({})={}".format(u, uStatus, l, lStatus)
        raise ValueError("UpperBound({})={}, LowerBound({})={}".format(u, uStatus, l, lStatus))

    variables_opt = [None] * len(variables)
    
    #temp_iters = kwargs_solver['max_iters']
    while u - l >= bisection_tol:
        parameter.value = (l + u) / 2.0
        ## solve the feasibility problem
        problem, ok = checked_solve(problem, solver, **kwargs_solver)
        if ok:
            status = problem.status
        else:
            status = 'unknown'

        if bisect_verbose:
            print "Range: {}-{}; parameter {} -> {}".format(l, u, parameter.value, status)
            
        if 'infeasible' in status or 'unknown' in status:
            u = parameter.value
            #kwargs_solver['max_iters'] = temp_iters
        elif 'optimal' in status:
            l = parameter.value
            # update Variables
            for i in range(len(variables)):
                variables_opt[i] = variables[i].value
            # update Parameters
            objval_opt = parameter.value
        else:
            raise ValueError("Problem may not bisectionable. Found Intermediate solution of '{}'".format(status))

    # Solve problem again for last feasible value (To ensure solved problem in prob instance at the end)
    parameter.value = objval_opt
    problem, ok = checked_solve(problem, solver, **kwargs_solver)

    return [variables_opt, objval_opt]

def bisect_max2(l, u, problem, parameter, variables,
           bisection_tol=1e-3, solvers=[cvxpy.CVXOPT, {}], bisect_verbose=False):
          
    if (u is not None):
        # cross check bound
        if (u < l):
            #print "upperBound < lowerBound"
            raise ValueError("upperBound({}) < lowerBound({})".format(u, l))

    elif (u is None) and (l is not None):
        # First iteration
        u = l + 1.0
        if bisect_verbose:
            print "processing upper bound: {}".format(u)
        parameter.value = u
        problem, ok = checked_solve(problem, solvers)
        if ok:
            uStatus = problem.status
        else:
            uStatus = 'unknown'
        
        while 'optimal' in uStatus:
            #if u >= l: # shift upper bound if found feasible, this condition is always true
                #l = u
            l = u
            u = 2.0*u
            if bisect_verbose:
                print "processing upper bound: {}".format(u)
            parameter.value = u
            problem, ok = checked_solve(problem, solvers)
            if ok:
                uStatus = problem.status
            else:
                uStatus = 'unknown'
        if bisect_verbose:
            print 'found bounds: [{}-{}]'.format(l, u)
    else:
        raise ValueError("Not implemented")
    
    # check validity solution of l is optimal, solution of u is infeasible
    parameter.value = l
    problem, ok = checked_solve(problem, solvers)
    if ok:
        lStatus = problem.status
    else:
        lStatus = 'unknown'

    parameter.value = u
    problem, ok = checked_solve(problem, solvers)
    if ok:
        uStatus = problem.status
    else:
        uStatus = 'unknown'

    if not ('optimal' in lStatus and uStatus in ['infeasible', 'unknown']):
        #print "UpperBound({})={}, LowerBound({})={}".format(u, uStatus, l, lStatus)
        raise ValueError("UpperBound({})={}, LowerBound({})={}".format(u, uStatus, l, lStatus))

    variables_opt = [None] * len(variables)
    
    #temp_iters = kwargs_solver['max_iters']
    while u - l >= bisection_tol:
        parameter.value = (l + u) / 2.0
        ## solve the feasibility problem
        problem, ok = checked_solve(problem, solvers)
        if ok:
            status = problem.status
        else:
            status = 'unknown'

        if bisect_verbose:
            print "Range: {}-{}; parameter {} -> {}".format(l, u, parameter.value, status)
            
        if status in ['infeasible', 'unknown'] :
            u = parameter.value
            #kwargs_solver['max_iters'] = temp_iters
        elif 'optimal' in status:
            l = parameter.value
            # update Variables
            for i in range(len(variables)):
                variables_opt[i] = variables[i].value
            # update Parameters
            objval_opt = parameter.value
        else:
            raise ValueError("Problem may not bisectionable. Found Intermediate solution of '{}'".format(status))

    # Solve problem again for last feasible value (To ensure solved problem in prob instance at the end)
    parameter.value = objval_opt
    problem, ok = checked_solve(problem, solvers)

    return [variables_opt, objval_opt]
''' Bisection function
Tries to find the min value of parameter of the given feasibility problem by bisection


l: lower bound with optimal solution
u: upper bound with infeasable solution
bisection_tol: bisection tolerance
sol: solver for problem
problem: cvxpy Problem
variables: array of variables in problems
parameter:

Note: make sure solution with l is infeasible, with u is optimal
Example: [[o_Q, o_z0, o_z1], o_g] = bisect_min(
                                        0, 2,
                                        prob, g,
                                        [Q, z0, z1],
                                        solver=cvxpy.SCS,
                                        verbose=False)
'''
def bisect_min(l, u, problem, parameter, variables,
           bisection_tol=1e-3, solver=cvxpy.CVXOPT, verbose=False):

    # cross check bound
    if (u < l):
        #print "upperBound < lowerBound"
        raise ValueError("upperBound({}) < lowerBound({})".format(u, l))

    # check validity solution of l is optimal, solution of u is infeasible
    parameter.value = l
    problem.solve(solver=solver)
    lStatus = problem.status

    parameter.value = u
    problem.solve(solver=solver)
    uStatus = problem.status

    if not ('optimal' in uStatus and 'optimal' not in lStatus):
        #print "UpperBound({})={}, LowerBound({})={}".format(u, uStatus, l, lStatus)
        raise ValueError("UpperBound({})={}, LowerBound({})={}".format(u, uStatus, l, lStatus))

    variables_opt = [None] * len(variables)
    while u - l >= bisection_tol:
        parameter.value = (l + u) / 2.0
        ## solve the feasibility problem
        problem.solve(solver=solver)
        if verbose:
            print "Range: {}-{}; parameter {} -> {}".format(l, u, parameter.value, problem.status)

        if 'optimal' in problem.status:
            u = parameter.value
            # update Variables
            for i in range(len(variables)):
                variables_opt[i] = variables[i].value
            #a_opt = a.value
            #b_opt = b.value
            # update Parameters
            objval_opt = parameter.value
        else:
            l = parameter.value
    return [variables_opt, objval_opt]

''' Helper functions'''
### Seite 37 (im Text zwischen (4.12) und (4.13))
def _N(n):
    return np.diag([p for p in xrange(-n, 0, 1)])

### Seite 55; (4.64)
def _M(n):
    return np.diag([p for p in xrange(0, n, 1)])

### Seite 35; (4.6)
# D = diag(v^n, ... v^2, v)
def _D(v, n):
    return np.diag([v**x for x in range(n, 0, -1)])

# D^-1 = (diag(v^n, ... v^2, v))^-1 = diag(v^-n, ..., v^-2, v^-1)
def _D_inv(v, n):
    return np.diag([v**-x for x in range(n, 0, -1)])

assert(np.allclose(_D_inv(0.123, 15), LA.inv(_D(0.123, 15))))

def _H(k, n):
    H = np.zeros((n, n))
    if k>n or k<=0:
        print "k = {}, smaller 0 and bigger {} is not really defined!".format(k, n-1)
        #raise ArithmeticError("k = {}, k smaller 0 and bigger n = {} is not defined".format(k, n))
    else:
        H[k-1, k-1] = 1
    return H
    
### Seite 55; (4.64)
def _P(l, k, n):
    I = np.eye(n)
    Mn = _M(n)
    P = I
    if k == 0:
        pass # P=I
    else:
        for q in xrange(0, k, 1):
            P = P * ((l-q)*I + Mn)
    return np.matrix(P)

def sat(v, u_max):
    return np.clip(v, -u_max, u_max)

''' State Space Simulator '''
def simulate(A, B, C, D, regulator_func, s, T, delay=None, umax=None, x0=0):
    #intitialize y, u
    y = np.matrix(np.zeros((C.shape[0],len(T))))
    u = np.zeros((len(T),np.size(x0,1)))
    u_sat = np.zeros((len(T),np.size(x0,1)))
    if type(x0) is int:
        xt = np.matrix([x0]*len(A)).T
        print "x0 = \n{}".format(xt)
    else:
        xt = x0

    if delay:
        s_queue = collections.deque(maxlen=int(math.ceil(delay/(T[1]-T[0]))))
        
    for i, t in enumerate(T):
        if delay:
            s_queue.append(s[i])
            if len(s_queue) == int(math.ceil(delay/(T[1]-T[0]))):
                s_delay = s_queue[0]
            else:
                s_delay = 0
        else:
            s_delay = s[i]
        
        u[[i],:] = regulator_func(y[:,i-1], s_delay, xt)

        if umax is not None:
            u_sat[[i],:] = sat(u[[i],:], umax)
        else:
            u_sat[[i],:] = u[[i],:]

        x_dot = A.dot(xt) + B.dot(u_sat[[i],:])
        
        y[:,i] = C.dot(xt) + D.dot(u_sat[[i],:])

        if i < len(T)-1:
            xt = xt + x_dot*(T[i+1]-T[i])
    return y, u, u_sat

# example Regulator function
def exampleRegulator(y, s, x):
    # fill-in K matrix euation. Below is just a controller matrix for
    # the Inverted Pendulum pendulum problem
    K = np.array([-70.7107  ,-37.8345  ,105.5298   ,20.9238])
    return s-K.dot(x)

# no controller just forwarding setpoint
def openLoop(y, s, x):
    #print s
    return s

''' Steuerungsnormalform aus python control '''
# Flipping matrixes to fit Adamy definition
def reverse_x_order(T):
    return np.flipud(np.fliplr(T))

'''
Returns transformed A, b, c, d and Transformation-Matrix T (x_trans = T*x) and Steuerbarkeitsmatrix Q
'''
def get_Steuerungsnormalform(A, b, c, d):
    #https://www.eit.hs-karlsruhe.de/mesysto/teil-a-zeitkontinuierliche-signale-und-systeme/darstellung-von-systemen-im-zustandsraum/transformation-auf-eine-bestimmte-darstellungsform/transformation-einer-zustandsgleichung-in-regelungsnormalform.html

    # Image(url = "https://www.eit.hs-karlsruhe.de/mesysto/fileadmin/images/Skript_SYS_V_10_0_2/Kapitel_10_3/Grafik_10_3_7_HQ.png")


    # Berechnung der inversen Steuerbarkeitsmatrix
    n = A.shape[0]
    Q = b #
    for i in range(1, n):
        Q = np.hstack([Q, LA.matrix_power(A,i).dot(b)])
    Q_inv = LA.inv(Q)

    #Zeilenvektor t_1.T entspricht der letzten Zeile der inversen Steuerbarkeitsmatrix
    #Image(url="https://www.eit.hs-karlsruhe.de/mesysto/fileadmin/images/Skript_SYS_V_10_0_2/Kapitel_10_3/Formel_10_3_51_HQ.png")

    t1 = Q_inv[-1,:]

    # Berechnung der Transformationsmatrix 
    #Image(url="https://www.eit.hs-karlsruhe.de/mesysto/fileadmin/images/Skript_SYS_V_10_0_2/Kapitel_10_3/Grafik_10_3_8_HQ.png")
    T = t1
    for i in range(1, n):
        T = np.vstack([T, t1.dot(LA.matrix_power(A,i))])

    #Bestimmung der Zustandsraumdarstellung in Regelungsnormalform 
    #Image(url="https://www.eit.hs-karlsruhe.de/mesysto/fileadmin/images/Skript_SYS_V_10_0_2/Kapitel_10_3/Grafik_10_3_9_HQ.png")
    #Image(url="https://www.eit.hs-karlsruhe.de/mesysto/fileadmin/images/Skript_SYS_V_10_0_2/Kapitel_10_3/Grafik_10_3_10_HQ.png")

    A0 = T.dot(A).dot(LA.inv(T))
    b0 = T.dot(b)
    c0 = (c.T.dot(LA.inv(T))).T

    return (A0, b0, c0, d), T, Q

# Helper function to take the calculation time
import time
class Timer(object):
    def __init__(self, name=None):
        self.name = name

    def __enter__(self):
        self.tstart = time.time()

    def __exit__(self, type, value, traceback):
        if self.name:
            print '[%s]' % self.name,
        print 'Elapsed: %s' % (time.time() - self.tstart) 
