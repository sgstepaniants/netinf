% Generate stacked Hankel matrix of time series data for n masses
% connected to their nearest neighbors by springs.

clear all; close all; clc;

nvars = 3;
deltat = 0.1;
endtime = 20;
nobs = endtime / deltat;
tSpan = linspace(0, endtime, nobs);
ntrials = 1000;

randpfn = @(n) rand(n, 1) - 0.5; % random [-0.5, 0.5]
%randpfn = @(n) [rand(1) / 2 - 0.25, 0.25, rand(1) / 2 + 0.25, 0.75]';
randvfn = @(n) zeros(n, 1); % masses start at rest (0 velocity)
randmfn = @(n) ones(n, 1); % all masses are equal
randkfn = @(n) ones(n, 1); % all spring constants are equal
randcfn = @(n) zeros(n, 1); % undamped oscillations
bc = 'fixed';

stacks = 2;
X = GenerateNNCoupledData(nvars, tSpan, ntrials, randpfn, randvfn, randmfn, randkfn, randcfn, bc, 0);
%H = MakeHaskel(X, stacks);



% plot the final time series simulation
%plot(H(:, :, 1)')
plot(mean(X, 3)')
