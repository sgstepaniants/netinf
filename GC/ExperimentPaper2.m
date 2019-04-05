% Simulates three connected spring mass oscillators with fixed boundary
% conditions. They have random initial positions and are released to
% oscillate.

clear all; close all; clc;

expnum = 'Paper2';

% Run MVGC toolbox 10 times and vote for final connectivity matrix.
nvars = 3;
mats = [0, 0, 0; 0, 0, 0; 0, 0, 0];
mats(:, :, 2) = [0, 0, 1; 1, 0, 1; 1, 0, 0];
mats(:, :, 3) = [0, 1, 1; 0, 0, 1; 1, 1, 0];

Ks = repmat(MakeNetworkTriDiag(nvars+2, false), [1, 1, 3]);
Ks(2:nvars+1, 2:nvars+1, :) = mats;

% Define time sampling.
endtime = 10;
nobs = 500;
tSpan = linspace(0, endtime, nobs);

% Number of times to split time interval.
tsplits = nobs;
% How often to save worspace (after runs of MVGC).
freq = 1;

% Specify noise and prepocessing for data.
noisefn  = @(data) WhiteGaussianNoise(data, 0.1);
preprocfn = @(data) NoiseThenDetrend(data, noisefn);

% Specify focring function for oscillators.
forcingFunc = zeros([nvars, nobs]);

% Specify boundary conditions.
bc = 'fixed';

% Number of simulation trials.
ntrials = 100;
% Number of times to repeat MVGC toolbox runs.
reps = 10;

BaseNNExperiment(expnum, mats, Ks, @(n)randpfn(n), @(n)zerovfn(n), @(n)constmfn(n,1), ...
    @(n)zerocfn(n), bc, forcingFunc, preprocfn, tSpan, ntrials, reps, ...
    tsplits, freq)

save(strcat('exp', expnum, '/', 'exp', expnum, '-params.mat'))
