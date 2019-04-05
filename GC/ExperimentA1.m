clear all; close all; clc;

expnum = 'A1';

% Run MVGC toolbox 10 times and vote for final connectivity matrix.
nvars = 4;
mats = MakeNetworkTriDiag(nvars, 1);
Ks = MakeNetworkTriDiag(nvars+2, 1);

% Define sampling time.
deltat = 0.05; % space between time points
endtime = 10;
nobs = round(endtime / deltat); % number of time points (observations)

% Specify noise and prepocessing for data.
noisefn  = @(data) WhiteGaussianNoise(data, 0.1);
preprocfn = @(data) NoiseThenDetrend(data, noisefn);

% Specify focring function for oscillators.
forcingFunc = zeros([nvars, nobs]);

% Specify boundary conditions.
bc = 'fixed';

% Number of times to repeat MVGC toolbox runs.
reps = 10;

ntrials = 100;
BaseNNExperiment(expnum, mats, Ks, @(n)randpfn(n), @(n)zerovfn(n), @(n)constmfn(n,1), ...
    @(n)zerocfn(n), bc, forcingFunc, preprocfn, deltat, endtime, ntrials, reps, ...
    tsplits, freq)

save(strcat('exp', expnum, '/', 'exp', expnum, '-params.mat'))
