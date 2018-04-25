% This script saves a .mat file with the parameters we usually want to use in
% our experiments. We can then load that .mat file and change the parts that we
% want to experiment with.

clear all; close all; clc;

% history of node positions
global prev_pos;
prev_pos = [];

% what value to threshold the amplitudes of the nodes at
global thresh;
thresh = 0.1;

% how many iterations to look back to calculate amplitude of nodes
global past;
past = 30;


% nvars is the number of variables/nodes (oscillators for us)
nvars = 12;

% set distributions for random initial conditions, and noise
measParam = 0.2;
stacks = 5;
randpfn = @(n) rand(n, 1) - 0.5; % random [-0.5, 0.5]
randvfn = @(n) zeros(n, 1); % masses start at rest (0 velocity)
randmfn = @(n) ones(n, 1);
randkfn = @(n) ones(n, 1);
randcfn = @(n) zeros(n, 1); % undamped oscillations

% set preprocessing step for data
noisefn  = @(data) AddGaussianNoise(data, measParam);
haskelfn = @(data) MakeHaskel(data, stacks);
preprocfn = @(data) noisefn(haskelfn(data));

% set the sampling in time
deltat = 0.1; % space between time points
endtime = 25 + (stacks - 1) * deltat; % solve NNCoupled model from 0 to T (endtime)
nobs = round(endtime / deltat); % number of time points (observations)

% Decide how many times to solve NNCoupled model, randomizing the
% frequencies, initial conditions, and noise each time.
% We generate ntrials * reps random instances.
ntrials = 100; % number of instances to pass to GC algorithm at once
reps = 1; % number of times to run GC algorithm (each with a different ...
% set of ntrials instances)

% if splitting the data in time, where do we split it?
% standard experiment: no split, so use time index 1 to nobs
tsplits = nobs - stacks + 1;

% How often save workspace (after freq runs of GC)
freq = 1;

% A, j, p were temporary, so no need to save them
clearvars A Asmall j p
save('UsualParams.mat')
