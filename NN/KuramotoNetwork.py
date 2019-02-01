import numpy as np
import random
import time
import torch
import torch.nn as nn
from torch.autograd import Variable
import torch.nn.functional as F
import torch.utils.data
import torch.optim as optim
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split
from matplotlib import pyplot as plt

from SimulateKuramotoOscillators import simulateKuramoto
from MakeNetworks import tridiag, erdos_reyni
from InitFunctions import *


def fitVelocities(pos, vel, model):
    m, n = pos.shape
    adj = np.zeros([n, n])
    vel_hat = np.zeros([m, n])
    interactionMats = np.zeros([n, n, m])
    
    for i in range(n):
        for j in range(n):
            interactionMats[i, j] = model(pos[:, [i, j]])
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
# Build Feed-Forward Neural Network
#******************************************************************************

# Hyper Parameters
batch_size = 100
momentum = 0.9

# Positions
X = torch.Tensor(X)
# Velocities
Y = torch.Tensor(Y)

# Dataset Loader (Input Pipline)
trainDataset = torch.utils.data.TensorDataset(X, Y)
train_loader = torch.utils.data.DataLoader(dataset=trainDataset, batch_size=batch_size, shuffle=True)


class Net(torch.nn.Module):

    def __init__(self):
        super(Net, self).__init__()
        self.fc1 = torch.nn.Linear(2, 100)
        self.fc2 = torch.nn.Linear(100, 100)
        self.fc3 = torch.nn.Linear(100, 100)
        self.fc4 = torch.nn.Linear(100, 100)
        self.fc5 = torch.nn.Linear(100, 1)
        #torch.nn.init.xavier_uniform(self.fc1.weight)

    def forward(self, x):
        x = self.fc1(x)
        x = F.leaky_relu(x, negative_slope=0.03)
        x = self.fc2(x)
        x = F.leaky_relu(x, negative_slope=0.03)
        x = self.fc3(x)
        x = F.leaky_relu(x, negative_slope=0.03)
        x = self.fc4(x)
        x = F.leaky_relu(x, negative_slope=0.03)
        x = self.fc5(x)
        return x

model = Net()

#******************************************************************************
# Train Feed-Forward Neural Network
#******************************************************************************

# Loss and Optimizer
# Softmax is internally computed.
# Set parameters to be updated.
loss_fn = nn.MSELoss()

# number of epochs/iterations for training the neural net
num_epochs = 5
freq = 1

# store train squared error for each epoch
train_err = np.zeros(num_epochs * trials / batch_size / freq)

# store predicted velocities and inferred adjacency matrices for each epoch
pred_mats = np.zeros([n, n, num_epochs * trials / batch_size / freq])

learning_rate = [1000 * 0.5**i for i in range(num_epochs * trials / batch_size)]
weight_decay = 1

# Training the Model
for epoch in range(num_epochs):
    for i, (pos, vel) in enumerate(train_loader):
        pos = Variable(pos)
        vel = Variable(vel)

        # Forward + Backward + Optimize
        ind = epoch * trials / batch_size + i
        #optimizer = torch.optim.SGD(model.parameters(), lr=learning_rate[ind],
        #                            weight_decay=weight_decay[ind], momentum=momentum)
        optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate[ind],
                                    weight_decay=weight_decay)
        optimizer.zero_grad()
        
        vel_hat, adj = fitVelocities(pos, vel.data, model)
        l1_reg = torch.norm(torch.Tensor(adj), 1)
        loss = loss_fn(Variable(torch.Tensor(vel_hat), requires_grad=True), vel) + l1_reg
        loss.backward()
        optimizer.step()

        if (i + 1) % freq == 0:
            err = loss.data[0]
            ind = epoch * trials / batch_size / freq + (i + 1) / freq - 1
            train_err[ind] = err
            pred_mats[:, :, ind] = adj
            print ('Epoch: [%d/%d], Step: [%d/%d], Loss: %.4f'
                   % (epoch + 1, num_epochs, i + 1, math.ceil(1.0*trials/batch_size), err))

pred_vel, adj = fitVelocities(Variable(X), Y, model)


#******************************************************************************
# Plot Train and Test Errors
#******************************************************************************
xs = np.arange(0, len(train_err))
plt.figure(3)
plt.title('Mean Squared Error over Epochs')
plt.plot(xs, train_err, 'b-', ms=2, label='train')
plt.legend(loc='upper right', title='Legend')
plt.xlabel('Effective Epoch')
plt.ylabel('Mean Squared Error')
plt.show()

print "Min Train Error: " + str(min(train_err))


plt.figure(4)
plt.title('True and Predicted Node Velocities')
for i in range(n):
    plt.plot(range(trials), Y[:, i].numpy(), label='true node '+str(i+1))
    plt.plot(range(trials), pred_vel[:, i], label='pred node '+str(i+1))
plt.xlabel('t')
plt.grid(True)
plt.legend()
plt.show()


#******************************************************************************
# Plot Neural Network Model Output
#******************************************************************************
npts = 100
x = np.linspace(0, 2*math.pi, npts)
y = np.linspace(0, 2*math.pi, npts)
xv, yv = np.meshgrid(x, y)
interactionFunc = np.sin(xv - yv)
plt.matshow(interactionFunc)


modelInteractionFunc = np.zeros([npts, npts])
for i in range(npts):
    for j in range(npts):
        modelInteractionFunc[i, j] = model(Variable(torch.Tensor([x[i], y[j]])))

plt.matshow(modelInteractionFunc)


print(adj)

#row_sums = np.abs(adj).sum(axis=1)
#A_pred = adj / row_sums[:, np.newaxis]
#print(A_pred)
