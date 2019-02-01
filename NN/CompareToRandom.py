def fitVelocities(pos, vel, modelMat):
    m, n = pos.shape
    adj = np.zeros([n, n])
    vel_hat = np.zeros([m, n])
    interactionMats = np.zeros([n, n, m])
    
    for i in range(n):
        for j in range(n):
            interactionMats[i, j] = modelMat[i, j]
        interactionMats[i, i] = 0
        w = np.linalg.lstsq(interactionMats[i].transpose(), vel[:, i], rcond=None)[0]
        vel_hat[:, i] = w.dot(interactionMats[i])
        adj[i] = w

    return vel_hat, adj


#******************************************************************************
# Create Dataset of Kuramoto Oscillator Positions
#******************************************************************************
n = 3
p = 0.5
K = 10
A = erdos_reyni(n, p)
#A = np.array([[0, 1, 1, 0, 0], [0, 0, 1, 0, 0], [0, 1, 0, 1, 0], [0, 0, 1, 0, 0], [0, 1, 1, 1, 0]])
print A

deltat = 0.01
endtime = 5
reps = 2
nobs = int(endtime/deltat)
tspan = np.linspace(0, endtime, nobs)
trials = reps * (nobs - 2)

noise = 0.01

X = noise * np.random.randn(trials, n)
Y = noise * np.random.randn(trials, n)
for r in range(reps):
    p0 = randcircpfn(n)
    w0 = 0 * randwfn(n)
    
    data = simulateKuramoto(A, tspan, K, p0, w0)
    # record node positions as predictor variables
    X[r*(nobs-2):(r+1)*(nobs-2)] += data[1:nobs-1]
    # record node velocities as output variables
    #Y[r*(nobs-2):(r+1)*(nobs-2)] = (data[2:nobs] - data[0:nobs-2]) / (2*deltat)
    Y[r*(nobs-2):(r+1)*(nobs-2)] += np.zeros([nobs-2, n])
    for i in range(nobs-2):
        r = np.array([X[i],]*n)
        dist = r - r.transpose()
        Y[i] = w0 + np.sum(K/n * np.multiply(A, np.sin(dist)), axis=1)

# Plot the solution.
plt.figure(1)
plt.title('Node Positions')
for i in range(n):
    plt.plot(range(trials), X[:, i], label='node '+str(i+1))
plt.xlabel('t')
plt.grid(True)
plt.legend()
plt.show()

plt.figure(2)
plt.title('Node Velocities')
for i in range(n):
    plt.plot(range(trials), Y[:, i], label='node '+str(i+1))
plt.xlabel('t')
plt.grid(True)
plt.legend()
plt.show()


#******************************************************************************
# Compare True Model Velocity Predictions to Random Model Predictions
#******************************************************************************
trueModelMat = 
pred_vel, adj = fitVelocities(X, Y, trueModelMat)
