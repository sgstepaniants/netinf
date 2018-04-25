load('UsualParams.mat')
expnum = 'A1';

% Run MVGC toolbox 10 times and vote for final connectivity matrix.
nvars = 4;
mats = MakeNetworkTriDiag(nvars, 1);
reps = 10;

BaseExperiment(expnum, mats, randpfn, randvfn, randmfn, randkfn, ...
    randcfn, preprocfn, deltat, endtime, ntrials, reps, tsplits, freq, stacks)
