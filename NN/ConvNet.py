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


def inputImage(pos):
    m, n = pos.shape
    image = np.zeros([m, 1, n, 2*n])
    for i in range(n):
        for j in range(n):
            image[:, 0, i, 2*j] = pos[:, i];
            image[:, 0, i, 2*j+1] = pos[:, j];
    return image


#******************************************************************************
# Create Dataset of Kuramoto Oscillator Positions
#******************************************************************************
n = 2
p = 0.5
K = 4
A = erdos_reyni(n, p)
A = np.array([[0, 1], [1, 0]])
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

# Positions
X = torch.Tensor(X)
# Velocities
Y = torch.Tensor(Y)

# Dataset Loader (Input Pipline)
trainDataset = torch.utils.data.TensorDataset(X, Y)
train_loader = torch.utils.data.DataLoader(dataset=trainDataset, batch_size=batch_size, shuffle=True)


num_filters = 2
class Net(torch.nn.Module):

    def __init__(self):
        super(Net, self).__init__()
        self.pairconv = torch.nn.Conv2d(1, num_filters, kernel_size=(1,2), stride=(1,2))
        self.depthconv1 = torch.nn.Conv3d(1, num_filters, kernel_size=(num_filters,1,1))
        self.depthconv2 = torch.nn.Conv3d(1, num_filters, kernel_size=(num_filters,1,1))
        self.flattenconv = torch.nn.Conv3d(1, 1, kernel_size=(num_filters,1,1))
        self.adjweights = torch.nn.Parameter(torch.randn([n, n]))
    
    def computeCoupling(self, x):
        x = self.pairconv(x)
        x = x.unsqueeze(1)
        x = F.relu(x)
        
        # First depth convolution
        x = self.depthconv1(x)
        x = x.transpose(1, 2)
        x = F.relu(x)
        
        # Second depth convolution
        x = self.depthconv2(x)
        x = x.transpose(1, 2)
        x = F.relu(x)
        
        x = self.flattenconv(x)
        x = x.squeeze()
        return x
    
    def forward(self, x):
        x = self.computeCoupling(x)
        x = x * self.adjweights
        x = torch.sum(x, dim=1)
        return x


def init_weights(m):
    if (type(m) == torch.nn.Conv2d) | (type(m) == torch.nn.Conv3d):
        torch.nn.init.xavier_uniform(m.weight)


model = Net()
model.apply(init_weights)

#******************************************************************************
# Train Feed-Forward Neural Network
#******************************************************************************

# Loss and Optimizer
# Softmax is internally computed.
# Set parameters to be updated.
loss_fn = nn.MSELoss()

# number of epochs/iterations for training the neural net
num_epochs = 2
freq = 1

# store train squared error for each epoch
train_err = np.zeros(num_epochs * trials / batch_size / freq)

# store predicted velocities and inferred adjacency matrices for each epoch
pred_mats = np.zeros([n, n, num_epochs * trials / batch_size / freq])

##LBFGS

learning_rate = [1 * 0.5**i for i in range(num_epochs * trials / batch_size)]
weight_decay = [0 * 0.1**i for i in range(num_epochs * trials / batch_size)]

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
                                    weight_decay=weight_decay[ind])
        optimizer.zero_grad()
        
        image = inputImage(pos)
        vel_hat = model(Variable(torch.Tensor(image)))
        
        #l1_reg = torch.norm(torch.Tensor(adj), 1)
        loss = loss_fn(vel_hat, vel)
        loss.backward()
        optimizer.step()

        if (i + 1) % freq == 0:
            err = loss.data[0]
            ind = epoch * trials / batch_size / freq + (i + 1) / freq - 1
            train_err[ind] = err
            print ('Epoch: [%d/%d], Step: [%d/%d], Loss: %.4f'
                   % (epoch + 1, num_epochs, i + 1, math.ceil(1.0*trials/batch_size), err))

image = inputImage(X)
pred_vel = model(Variable(torch.Tensor(image)))
adj = model.adjweights


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
    plt.plot(range(trials), pred_vel[:, i].data.numpy(), label='pred node '+str(i+1))
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


image = np.zeros([npts**2, 1, n, 2*n])
image[:, 0, 0, 0] = np.repeat(x, npts)
image[:, 0, 0, 1] = np.tile(x, npts)
coupling = model.computeCoupling(Variable(torch.Tensor(image))).data
modelInteractionFunc = np.reshape(coupling[:, 0, 0], (npts, npts))

plt.matshow(modelInteractionFunc)


print(adj)

#row_sums = np.abs(adj).sum(axis=1)
#A_pred = adj / row_sums[:, np.newaxis]
#print(A_pred)
