import numpy as np
import os
import errno
import shelve
from SimulateSpringMassOscillators import simulateSpringMass
from MakeNetworks import tridiag, erdos_reyni
from InitFunctions import *

#******************************************************************************
# Create Dataset of Spring Mass Simulations and their Connectivity Matrices
#******************************************************************************
# arguments to ode solver

ns = [5]
ps = [0.5]
tspan = np.linspace(0, 10, 200)
bc = 'fixed'

for n in ns:
    for p in ps:
        mass = constmfn(n, 1)
        spring = 1
        damp = zerocfn(n)
        
        trials = 10000
        di = n * tspan.shape[0]
        df = n * n
        X = np.zeros((trials, di))
        Y = np.zeros((trials, df))
        for i in range(trials):
            # Use for fixed boundary conditions
            K = tridiag(n+2, corners=False)
            adj = erdos_reyni(n, p)
            K[1:n+1, 1:n+1] = adj
            K = spring * K
            
            # Use for free boundary conditions
            #adj = symmetric_erdos_reyni(n, p)
            #K = spring * adj
            
            # add true connectivity matrix to dataset
            Y[i] = adj.flatten()
            
            # simulate data and add it to the dataset
            p0 = randpfn(n)
            v0 = randvfn(n)
            X[i] = simulateSpringMass(n, tspan, K, p0, v0, mass, damp, bc).flatten()
            print i
        
        filename = 'data/' + str(bc) + '/n=' + str(n) + '/p=' + str(p) + '/exp'
        if not os.path.exists(os.path.dirname(filename)):
            try:
                os.makedirs(os.path.dirname(filename))
            except OSError as exc:  # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise

        data = shelve.open(filename)
        data['sims'] = X               # simulations
        data['adjs'] = Y               # true connectivity matrix
        data['n'] = n                  # number of nodes
        data['p'] = p                  # probability connection is made
        data['tspan'] = tspan          # simulation time
        data['pfn'] = 'randpfn'        # initial positions function
        data['vfn'] = 'randvfn'        # initial velocities function
        data['mfn'] = 'constmfn(1)'    # masses function
        data['cfn'] = 'zerocfn'        # damping coefficients function
        data['k'] = spring             # spring constant (same for all springs)
        data['bc'] = bc                # simulation boundary conditions
        data.close()
