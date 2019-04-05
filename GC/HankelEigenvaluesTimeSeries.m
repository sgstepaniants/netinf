% Stack the nncoupled system into a Hankel
% matrix s times and compute its V eigenvalue time series.

clear all; close all; clc;

nvars = 4;
ntrials = 1;
deltat = 0.1;

randpfn = @(n) rand(n, 1) - 0.5; % random [-0.5, 0.5]
randvfn = @(n) zeros(n, 1); % masses start at rest (0 velocity)
randmfn = @(n) ones(n, 1); % all masses are equal
randkfn = @(n) ones(n, 1); % all spring constants are equal
randcfn = @(n) zeros(n, 1); % undamped oscillations
bc = 'circ';

stacks = 1 : 10 : 100;

for s = stacks
    endtime = 25 + (s - 1) * deltat;
    nobs = round(endtime / deltat);
    tSpan = linspace(0, endtime, nobs);

    X = GenerateNNCoupledData(nvars, tSpan, ntrials, randpfn, randvfn, randmfn, randkfn, randcfn, bc, 0);

    [U, S, V] = MakeHankel(X, s);
    plot(V)
    pause(0.1)
end
