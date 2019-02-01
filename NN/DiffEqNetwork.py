import numpy as np
import random
import torch
from torch.autograd import Variable
import torch.nn.functional as F
import torch.optim as optim
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split
from matplotlib import pyplot as plt

from SimulateSpringMassOscillators import simulateSpringMass
from MakeNetworks import tridiag, erdos_reyni
from InitFunctions import *

#******************************************************************************
# Create Dataset of Spring Mass Positions and Velocities
#******************************************************************************
n = 2

trials = 100000
tspan = np.linspace(0, 10, trials+2)
deltat = tspan[1] - tspan[0]
bc = 'fixed'

mass = constmfn(n, 1)
spring = 1
damp = zerocfn(n)

K = tridiag(n+2, corners=False)
prob = 0.5
adj = erdos_reyni(n, prob)
K[1:n+1, 1:n+1] = adj
K = spring * K

p0 = randpfn(n)
v0 = randvfn(n)

# record node positions as predictor variables
data = simulateSpringMass(n, tspan, K, p0, v0, mass, damp, bc)
# record node velocities as output variables
Y = 100*(data[2:trials+2] - data[0:trials]) / (2*deltat)
X = data[1:trials+1]

#******************************************************************************
# Preprocess Data
#******************************************************************************
X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.33)
X_train = torch.Tensor(X_train)
X_test = torch.Tensor(X_test)
Y_train = torch.Tensor(Y_train)
Y_test = torch.Tensor(Y_test)

#******************************************************************************
# Build Feed-Forward Neural Network
#******************************************************************************
# 1 hidden layer
di = X_train.shape[1]
d1 = 2
df = Y_train.shape[1]

class Net(torch.nn.Module):

    def __init__(self):
        super(Net, self).__init__()
        self.fc1 = torch.nn.Linear(di, d1)
        torch.nn.init.xavier_uniform(self.fc1.weight)
        
        self.fc2 = torch.nn.Linear(d1, df)
        torch.nn.init.xavier_uniform(self.fc2.weight)

    def forward(self, x):
        x = self.fc1(x)
        x = F.sigmoid(x)
        x = self.fc2(x)
        return x

network_model = Net()

#******************************************************************************
# Train Feed-Forward Neural Network
#******************************************************************************
# store train and test squared errors for each iteration
train_err = []
test_err = []

# number of epochs/iterations for training the neural net
num_epochs = 1000
batch_size = 500
freq = 10

learning_rate = [1 * 0.5**i for i in range(num_epochs / freq)]
weight_decay = [1 * 0**i for i in range(num_epochs / freq)]
momentum = 0.9

# for each epoch, do a step of SGD and recompute the weights of the neural net
i = 0
for idx in range(num_epochs):
    print idx
    
    if idx % freq == 0:
        # update the optimizer
        opt = optim.SGD(network_model.parameters(), lr=learning_rate[i],
                 weight_decay=weight_decay[i], momentum = 0.9)
        
        # save the train and test errors every time we have sampled freq batches
        Y_train_hat = network_model(Variable(X_train)).data
        train_err.append(mean_squared_error(Y_train, Y_train_hat))
        Y_test_hat = network_model(Variable(X_test)).data
        test_err.append(mean_squared_error(Y_test, Y_test_hat))
        i += 1

    # train the neural net using mini-batch SGD
    batch_idx = random.sample(range(0, X_train.shape[0]), batch_size)
    network_model.zero_grad()

    # get a batch sample of the output
    Y_batch_train = Variable(Y_train[batch_idx, :])
    # predict the batch sample of the output
    Y_batch_train_hat = network_model(Variable(X_train[batch_idx, :]))
    # compute the mean squared error as our loss
    loss = F.mse_loss(Y_batch_train, Y_batch_train_hat)

    # this computes the gradient for us!
    loss.backward()
    # this does the parameter for us!
    opt.step()


# get network weights
W = list(network_model.parameters())
#A = np.matrix(W[4].data) * np.matrix(W[2].data) * np.matrix(W[0].data)

#******************************************************************************
# Plot Train and Test Errors
#******************************************************************************
xs = np.arange(0, num_epochs / freq)
plt.figure(1)
plt.title('Mean Squared Error over Effective Epochs')
plt.plot(xs, train_err, 'bo', ms=2, label='train')
plt.plot(xs, test_err, 'ro', ms=2, label='test')
plt.legend(loc='upper right', title='Legend')
plt.xlabel('Effective Epoch')
plt.ylabel('Mean Squared Error')
plt.show()

print "Min Train Error: " + str(min(train_err))
print "Min Test Error: " + str(min(test_err))


#******************************************************************************
# Plot Prediction of Node Velocities
#******************************************************************************
plt.figure(2)
plt.title('Prediction of Node Velocities')
plt.plot(tspan[1:trials+1], X, ms=1, label='true')
plt.plot(tspan[1:trials+1], network_model(Variable(torch.Tensor(X))).data.numpy(), ms=1, label='pred')
plt.legend(loc='upper right', title='Legend')
plt.xlabel('Time')
plt.ylabel('Node Velocities')
plt.show()
