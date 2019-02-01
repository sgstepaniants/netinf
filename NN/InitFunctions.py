import numpy as np
import math
from numpy.random import random_sample, randn

# functions that initialize the node positions
def unifpfn(n, bc):
    if bc == 'circ':
        p = np.arange(0.0, n) / n - 0.5
    elif bc == 'fixed':
        p = np.arange(1.0, n+1) / (n+1) - 0.5
    elif bc == 'free':
        p = np.arange(0.0, n) / (n-1) - 0.5
    else:
        raise Exception('invalid boundary conditions')
    return p


def randpfn(n):
    return random_sample(n) - 0.5


def randcircpfn(n):
    return 2 * math.pi * random_sample(n)


def randwfn(n):
    return 2 * random_sample(n) - 1


def smallrandpfn(n, bc):
    p = unifpfn(n, bc)
    if bc == 'circ':
        p += randn(n) / (2*n)
    elif bc == 'fixed':
        p += randn(n) / (2*(n+1))
    elif bc == 'free':
        p += randn(n) / (2*(n-1))
    else:
        raise Exception('invalid boundary conditions')
    return p


# functions that initialize the node velocities
def randvfn(n):
    return randn(n) 


def zerovfn(n):
    return np.zeros(n)


# functions that initialize the node masses
def constmfn(n, const):
    return const * np.ones(n)


def randmfn(n):
    return random_sample(n)


# functions that initialize the damping coefficients
def constcfn(n, const):
    return const * np.ones(n)


def zerocfn(n):
    return np.zeros(n)
