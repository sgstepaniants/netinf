load('UsualParams.mat')
expnum = 'A2';

% Run small matrices 10 times with evenly spaced masses and get back
% a fully disconnected connectivity matrix.
nvars = 6;
mats = MakeNetworkTriDiag(nvars, 1);
reps = 10;

randpfn = @(n) [0 : 1 / n : (n - 1) / n]'; % evenly spaced

BaseExperiment(expnum, smallMats, randpfn, randvfn, randmfn, randkfn, ...
    preprocfn, deltat, endtime, ntrials, reps, tsplits, freq)
