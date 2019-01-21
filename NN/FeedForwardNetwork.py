import os
import errno
import shelve
import numpy as np
import random
import torch
from torch.autograd import Variable
import torch.nn.functional as F
import torch.optim as optim
from sklearn.metrics import mean_squared_error
from matplotlib import pyplot as plt

#******************************************************************************
# Read Dataset of Spring Mass Simulations and their Connectivity Matrices
#******************************************************************************
filepath = 'data/fixed/n=5/p=0.5/'
data = shelve.open(filepath + 'exp')
X = data['sims']
Y = data['adjs']
n = data['n']
data.close()


#******************************************************************************
# Preprocess Data
#******************************************************************************
trials = X.shape[0]
train_size = int(0.6 * trials)
X_train = torch.Tensor(X[0:train_size])
Y_train = torch.Tensor(Y[0:train_size])
X_test = torch.Tensor(X[train_size:])
Y_test = torch.Tensor(Y[train_size:])


#******************************************************************************
# Build Feed-Forward Neural Network
#******************************************************************************
# 2 hidden layers
di = X_train.shape[1]
d1 = 100
d2 = 100
df = Y_train.shape[1]
# build the computational graph
network_model = torch.nn.Sequential(
                torch.nn.Linear(di, d1),
                torch.nn.ReLU(),
                torch.nn.Linear(d1, d2),
                torch.nn.ReLU(),
                torch.nn.Linear(d2, df),
                torch.nn.Sigmoid()
            )


#******************************************************************************
# Train Feed-Forward Neural Network
#******************************************************************************
# store train and test squared errors for each iteration
train_err = []
test_err = []

# number of epochs/iterations for training the neural net
num_epochs = 1000
batch_size = 600
freq = 100
learning_rate = [0.1 * 1**i for i in range(num_epochs / freq)]
weight_decay = [0 * 0.8**i for i in range(num_epochs / freq)]

# for each epoch, do a step of SGD and recompute the weights of the neural net
i = 0
for idx in range(num_epochs):
    print idx
    
    if idx % freq == 0:
        # update the optimizer
        opt = optim.Adam(network_model.parameters(), lr=learning_rate[i],
                 weight_decay=weight_decay[i])
        
        # save the train and test errors every time we have sampled freq batches
        Y_train_hat = network_model(Variable(X_train)).data
        train_err.append(mean_squared_error(Y_train, Y_train_hat))
        Y_test_hat = network_model(Variable(X_test)).data
        test_err.append(mean_squared_error(Y_test, Y_test_hat))
        i += 1

    # train the neural net using mini-batch SGD
    batch_idx = random.sample(range(0, train_size), batch_size)
    network_model.zero_grad()

    # get a batch sample of the output
    Y_batch_train = Variable(Y_train[batch_idx, :], requires_grad=True)
    # predict the batch sample of the output
    Y_batch_train_hat = Variable(network_model(Variable(X_train[batch_idx, :])).data, requires_grad=False)
    # compute the mean squared error as our loss
    loss = F.binary_cross_entropy(Y_batch_train, Y_batch_train_hat)

    # this computes the gradient for us!
    loss.backward()
    # this does the parameter for us!
    opt.step()


#******************************************************************************
# Plot Train and Test Errors
#******************************************************************************
xs = np.arange(0, num_epochs / freq)
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
# Find Optimal Threshold for Connectivity Matrix
#******************************************************************************
thresh = 0.5
aveTrainFalsePos = 0.0
aveTrainFalseNeg = 0.0
for i in range(train_size):
    y = np.array(Y_train[i])
    y_hat = np.array(network_model(Variable(X_train[i])).data) > thresh
    aveTrainFalsePos += np.count_nonzero(y_hat - y == 1)
    aveTrainFalseNeg += np.count_nonzero(y - y_hat == 1)
aveTrainFalsePos /= train_size
aveTrainFalseNeg /= train_size

print "Average Train False Positives: " + str(aveTrainFalsePos)  # predicted per simulation
print "Average Train False Negatives: " + str(aveTrainFalseNeg)  # predicted per simulation


aveTestFalsePos = 0.0
aveTestFalseNeg = 0.0
for i in range(trials - train_size):
    y = np.array(Y_test[i])
    y_hat = np.array(network_model(Variable(X_test[i])).data) > thresh
    aveTestFalsePos += np.count_nonzero(y_hat - y == 1)
    aveTestFalseNeg += np.count_nonzero(y - y_hat == 1)
aveTestFalsePos /= (trials - train_size)
aveTestFalseNeg /= (trials - train_size)

print "Average Test False Positives: " + str(aveTestFalsePos)  # predicted per simulation
print "Average Test False Negatives: " + str(aveTestFalseNeg)  # predicted per simulation


#******************************************************************************
# Plot Examples of True and Predicted Connectivity Matrices
#******************************************************************************
i = 1
true = Y[i].reshape(n, n)
pred = np.array(network_model(Variable(torch.Tensor(X)[i])).data).reshape(n, n)

plt.matshow(true)
for (i, j), z in np.ndenumerate(true):
    plt.text(j, i, '{:0.2f}'.format(z), ha='center', va='center')
plt.show()

plt.matshow(pred)
for (i, j), z in np.ndenumerate(pred):
    plt.text(j, i, '{:0.2f}'.format(z), ha='center', va='center')
plt.show()


#******************************************************************************
# Save Network Parameters and Results
#******************************************************************************
filename = filepath + 'results'
if not os.path.exists(os.path.dirname(filename)):
    try:
        os.makedirs(os.path.dirname(filename))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise

results = shelve.open(filename)
results['trials'] = trials
results['train_size'] = train_size
results['test_size'] = trials - train_size
results['network_model'] = network_model
results['num_epochs'] = num_epochs
results['batch_size'] = batch_size
results['learning_rate'] = learning_rate
results['weight_decay'] = weight_decay
results['freq'] = freq
results['loss'] = 'MSE'
results['train_err'] = train_err
results['min_train_err'] = min(train_err)
results['test_err'] = test_err
results['min_test_err'] = min(test_err)
results['ave_train_false_pos'] = aveTrainFalsePos
results['ave_train_false_neg'] = aveTrainFalseNeg
results['ave_test_false_pos'] = aveTestFalsePos
results['ave_test_false_neg'] = aveTestFalseNeg
results.close()
