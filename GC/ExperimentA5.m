load('UsualParams.mat')
expnum = 'A5';

% Run two nodes 10 times with a slight shift in one of the masses
% and get back estimated MVGC connectivity matrix.
nvars = 2;
mats = MakeNetworkTriDiag(nvars, 1);
reps = 1;

noisefn  = @(data) WhiteGaussianNoise(data, 0.005);
preprocfn = @(data) noisefn(data);
%randpfn = @(n) [0.1, 1 : n - 1]' / n; % evenly spaced
%randvfn = @(n) [0.5, zeros(1, n - 1)]'; % one mass perturbed

BaseExperiment(expnum, mats, randpfn, randvfn, randmfn, randkfn, ...
    preprocfn, deltat, endtime, ntrials, reps, tsplits, freq)
