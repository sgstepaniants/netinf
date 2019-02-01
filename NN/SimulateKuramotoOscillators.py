import numpy as np
from InitFunctions import *
from MakeNetworks import tridiag
from scipy.integrate import ode
from matplotlib import pyplot as plt


def simulateKuramoto(A, tspan, K, p0, w0):
    n = A.shape[1]

    def kur(t, p):
        r = np.array([p,]*n)
        dist = r - r.transpose()
    
        v = w0 + np.sum(K/n * np.multiply(A, np.sin(dist)), axis=1)
        return v
    
    
    # Create an `ode` instance to solve the system of differential
    # equations defined by `kur`, and set the solver method to 'dop853'.
    solver = ode(kur)
    solver.set_integrator('dop853')
    
    # Set the initial value p(0) = p0.
    t0 = tspan[0]
    solver.set_initial_value(p0, t0)
    
    # Create an array to hold the solution.
    # Put the initial value in the solution array.
    t1 = tspan[-1]
    sol = np.empty((tspan.shape[0], n))
    sol[0] = p0
    
    # Repeatedly call the `integrate` method to advance the
    # solution to time tspan[k], and save the solution in sol[k].
    k = 1
    while solver.successful() and solver.t < t1:
        solver.integrate(tspan[k])
        sol[k] = solver.y[0:n]
        k += 1
    
    return sol
