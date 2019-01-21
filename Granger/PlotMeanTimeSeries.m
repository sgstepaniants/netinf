% Take mean of time series data over several trials for n masses
% connected to their nearest neighbors by springs.

clear all; close all; clc;

nvars = 4;

deltat = 0.1;
endtime = 10;
nobs = endtime / deltat;
tSpan = linspace(0, endtime, nobs);

ntrials = 500;

%randpfn = @(n) rand(n, 1) - 0.5; % random [-0.5, 0.5]
randpfn = @(n) [-0.5; -0.25; 0; 0.25] + 0.25 * rand(4, 1);
randvfn = @(n) zeros(n, 1); % masses start at rest (0 velocity)
randmfn = @(n) ones(n, 1); % all masses are equal
randkfn = @(n) ones(n, 1); % all spring constants are equal
randcfn = @(n) zeros(n, 1); % undamped oscillations
bc = 'circ';

X = GenerateNNCoupledData(nvars, tSpan, ntrials, randpfn, randvfn, randmfn, randkfn, randcfn, bc, 0);

% plot the average time series simulation
plot(mean(X, 3)')
