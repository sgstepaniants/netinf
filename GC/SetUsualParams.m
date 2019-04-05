% This script saves a .mat file with the parameters we usually want to use in
% our experiments. We can then load that .mat file and change the parts that we
% want to experiment with.

clear all; close all; clc;
addpath('./InitFunctions/')

% nvars is the number of variables/nodes (oscillators for us)
nvars = 12;

bc = 'fixed';

% set distributions for random initial conditions, and noise
measParam = 0.01;
stacks = 1;

% set preprocessing step for data
noisefn  = @(data) AddGaussianNoise(data, measParam);

rank = 5;
hankelfn = @(data) MakeHankel(data, stacks, rank);

forcingFunc = zeros([nvars, length(tSpan)]);

%preprocfn = @(data) noisefn(data);
preprocfn = @(data) NoiseThenDetrend(data, noisefn);
%preprocfn = @(data) NoiseThenSine(data, noisefn);
%preprocfn = @(data) hankelfn(mydetrend(noisefn(data)));

% set the sampling in time
deltat = 0.05; % space between time points
endtime = 10;
%endtime = 25 + (stacks - 1) * deltat; % solve NNCoupled model from 0 to T (endtime)
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

save('UsualParams.mat')
