load('UsualParams.mat')
expnum = 'A6';

% Run three nodes 1 time with ran
% and get back estimated MVGC connectivity matrix.
nvars = 3;
mats = MakeNetworkTriDiag(nvars, 1, true);
reps = 1;

noisefn  = @(data) WhiteGaussianNoise(data, 0.01);
%preprocfn = @(data) NoiseThenDetrend(data, noisefn);
preprocfn = @(data) noisefn(data);

BaseExperiment(expnum, mats, randpfn, randvfn, randmfn, randkfn, ...
    preprocfn, deltat, endtime, ntrials, reps, tsplits, freq)
