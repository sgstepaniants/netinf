import numpy as np
from InitFunctions import *
from MakeNetworks import tridiag
from scipy.integrate import ode
from matplotlib import pyplot as plt


def simulateSpringMass(n, tspan, K, p0, v0, mass, damp, bc):
    def springMass(t, z):
        p = z[0:n]
        p = p - unifpfn(n, bc)
        v = z[n:]
    
        if bc == 'circ':
            if K.shape[0] == n and K.shape[1] == n:
                r = np.array([p,]*n)
                dist = r - r.transpose()
                force = np.sum((K + K.transpose()) * dist, axis=1)
            else:
                raise Exception('Spring constants for circular boundary conditions must have dimension n x n')
        elif bc == 'fixed':
            if K.shape[0] == n+2 and K.shape[1] == n+2:
                p = np.insert(p, 0, 0)
                p = np.insert(p, n+1, 0)
                r = np.array([p,]*(n+2))
                dist = r - r.transpose()
                force = np.sum(K * dist, axis=1)
                force = force[1:(n + 1)]
            else:
                raise Exception('Spring constants for fixed boundary conditions must have dimension (n + 2) x (n + 2)')
        elif bc == 'free':
            if K.shape[0] == n and K.shape[1] == n:
                r = np.array([p,]*n)
                dist = r - r.transpose()
                force = np.sum(K * dist, axis=1)
            else:
                raise Exception('Spring constants for free boundary conditions must have dimension n x n')
    
        a = np.divide(np.subtract(force, np.multiply(damp, v)), mass)
        f = np.concatenate((v, a))
        return f
    
    
    # Create an `ode` instance to solve the system of differential
    # equations defined by `springMass`, and set the solver method to 'dop853'.
    solver = ode(springMass)
    solver.set_integrator('dop853')
    
    # Set the initial value z(0) = z0.
    t0 = tspan[0]
    
    z0 = np.concatenate((p0, v0))
    solver.set_initial_value(z0, t0)
    
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
    
    # Plot the solution.
    #for i in range(n):
    #    plt.plot(tspan, sol[:, i], label='node '+str(i+1))
    #plt.xlabel('t')
    #plt.grid(True)
    #plt.legend()
    #plt.show()
    
    return sol
