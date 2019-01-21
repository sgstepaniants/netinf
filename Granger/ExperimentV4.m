load('UsualParams.mat')
expnum = 'V5_100s_rank=5_fixed';

nvars = 5;
mats = MakeNetworkTriDiag(rank);
reps = 10;
ntrials = 100;

% Save all initial conditions in experiment folder.
save('UsualParams.mat')

BaseHankelExperiment(expnum, nvars, mats, randpfn, randvfn, randmfn, randkfn, ...
    randcfn, preprocfn, deltat, endtime, ntrials, reps, tsplits, freq)

copyfile('UsualParams.mat', strcat('exp', expnum, '/', 'exp', expnum, '-params.mat'))
