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
from sklearn import preprocessing
from matplotlib import pyplot as plt

from SimulateKuramotoOscillators import simulateKuramoto
from MakeNetworks import tridiag, erdos_reyni
from InitFunctions import *

def calculateCoupling(pos, model):
    m, n = pos.shape
    coupling = np.zeros([m, n**2])
    for i in range(n):
        for j in range(n):
            coupling[:, n*i+j] = model(pos[:, [i, j]])
    return coupling


#******************************************************************************
# Create Dataset of Kuramoto Oscillator Positions
#******************************************************************************
n = 3
p = 0.5
K = 10
A = erdos_reyni(n, p)
#A = np.array([[0, 1, 1, 0, 0], [0, 0, 1, 0, 0], [0, 1, 0, 1, 0], [0, 0, 1, 0, 0], [0, 1, 1, 1, 0]])
print A

deltat = 0.1
endtime = 5
reps = 10
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
    Y[r*(nobs-2):(r+1)*(nobs-2)] = (data[2:nobs] - data[0:nobs-2]) / (2*deltat)
    #Y[r*(nobs-2):(r+1)*(nobs-2)] += np.zeros([nobs-2, n])
    #for i in range(nobs-2):
    #    r = np.array([X[i],]*n)
    #    dist = r - r.transpose()
    #    Y[i] = w0 + np.sum(K/n * np.multiply(A, np.sin(dist)), axis=1)

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

# Normalize the data
X = preprocessing.normalize(X)
Y = preprocessing.normalize(Y)

# Positions
X = torch.Tensor(X)
# Velocities
Y = torch.Tensor(Y)

# Dataset Loader (Input Pipline)
trainDataset = torch.utils.data.TensorDataset(X, Y)
train_loader = torch.utils.data.DataLoader(dataset=trainDataset, batch_size=batch_size, shuffle=True)


class Net1(torch.nn.Module):
    
    def __init__(self):
        super(Net1, self).__init__()
        self.apply(init_weights)
        self.fc1 = torch.nn.Linear(2, 10)
        self.fc2 = torch.nn.Linear(10, 10)
        self.fc3 = torch.nn.Linear(10, 10)
        self.fc4 = torch.nn.Linear(10, 10)
        self.fc5 = torch.nn.Linear(10, 1)

    def forward(self, x):
        x = self.fc1(x)
        x = F.leaky_relu(x, negative_slope=0.03)
        x = self.fc2(x)
        x = F.sigmoid(x)
        x = self.fc3(x)
        x = F.leaky_relu(x, negative_slope=0.03)
        x = self.fc4(x)
        x = F.sigmoid(x)
        x = self.fc5(x)
        return x


class Net2(torch.nn.Module):

    def __init__(self):
        super(Net2, self).__init__()
        self.fc1 = torch.nn.Linear(n**2, n, bias=False)

    def forward(self, x):
        x = self.fc1(x)
        return x
    
    def zeroWeights(self):
        #print(self.fc1.weight.data)
        for i in range(n):
            if i > 0:
                self.fc1.weight[i, 0:n*i].data.clamp_(0, 0)
            if i < n-1:
                self.fc1.weight[i, n*i+n:n**2].data.clamp_(0, 0)
        #print(self.fc1.weight.data)
    
    def adjWeights(self):
        adj = np.zeros([n, n])
        for i in range(n):
            adj[i, :] = self.fc1.weight[i, n*i:n*i+n].data
        return adj
    
def init_weights(m):
        if type(m) == nn.Linear:
            torch.nn.init.xavier_uniform(m.weight)

model1 = Net1()
#model1.apply(init_weights)

model2 = Net2()
#model2.apply(init_weights)

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

learning_rate = [10 * 0.1**i for i in range(num_epochs * trials / batch_size)]
weight_decay = [10 * 0.1**i for i in range(num_epochs * trials / batch_size)]

# Training the Model
for epoch in range(num_epochs):
    for i, (pos, vel) in enumerate(train_loader):
        pos = Variable(pos)
        vel = Variable(vel)
        
        model2.zeroWeights()

        # Forward + Backward + Optimize
        ind = epoch * trials / batch_size + i
        optimizer1 = torch.optim.SGD(model1.parameters(), lr=learning_rate[ind],
                                    weight_decay=weight_decay[ind], momentum=momentum)
        optimizer2 = torch.optim.SGD(model2.parameters(), lr=learning_rate[ind],
                                    weight_decay=weight_decay[ind], momentum=momentum)
        #optimizer1 = torch.optim.Adam(model1.parameters(), lr=learning_rate[ind],
        #                            weight_decay=weight_decay)
        #optimizer2 = torch.optim.Adam(model2.parameters(), lr=learning_rate[ind],
        #                            weight_decay=weight_decay)
        
        optimizer1.zero_grad()
        optimizer2.zero_grad()
        
        coupling = calculateCoupling(pos, model1)
        vel_hat = model2(Variable(torch.Tensor(coupling)))
        
        loss = loss_fn(vel_hat, vel)
        loss.backward()
        
        optimizer1.step()
        optimizer2.step()

        if (i + 1) % freq == 0:
            err = loss.data[0]
            ind = epoch * trials / batch_size / freq + (i + 1) / freq - 1
            train_err[ind] = err
            pred_mats[:, :, ind] = adj
            print ('Epoch: [%d/%d], Step: [%d/%d], Loss: %.4f'
                   % (epoch + 1, num_epochs, i + 1, math.ceil(1.0*trials/batch_size), err))


coupling = calculateCoupling(Variable(X), model1)
pred_vel = model2(Variable(torch.Tensor(coupling))).data
adj = model2.adjWeights()


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
    plt.plot(range(trials), pred_vel[:, i].numpy(), label='pred node '+str(i+1))
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
couplingFunc = np.sin(xv - yv)
plt.matshow(couplingFunc)


modelCouplingFunc = np.zeros([npts, npts])
for i in range(npts):
    for j in range(npts):
        modelCouplingFunc[i, j] = model1(Variable(torch.Tensor([x[i], y[j]])))

plt.matshow(modelCouplingFunc)


print(adj.transpose)

#row_sums = np.abs(adj).sum(axis=1)
#A_pred = adj / row_sums[:, np.newaxis]
#print(A_pred)
