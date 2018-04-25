load('UsualParams.mat')
expnum = 'A4';

% Run small matrices 10 times with a slight shift in one of the masses
% and get back estimated MVGC connectivity matrix.
nvars = 6;
mats = MakeNetworkTriDiag(nvars, 1);
reps = 10;

randpfn = @(n) [0 : 1 / n : (n - 1) / n]'; % evenly spaced
randvfn = @(n) [0.1, zeros(1, n - 1)]'; % one mass perturbed

BaseExperiment(expnum, smallMats, randpfn, randvfn, randmfn, randkfn, ...
    preprocfn, deltat, endtime, ntrials, reps, tsplits, freq)
