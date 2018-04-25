load('UsualParams.mat')
expnum = 'H4_5s';

nvars = 4;
mats = MakeNetworkTriDiag(nvars, stacks);
reps = 6;
ntrials = 100;

% Save all initial conditions in experiment folder.
save('UsualParams.mat')

BaseExperiment(expnum, mats, randpfn, randvfn, randmfn, randkfn, ...
    randcfn, preprocfn, deltat, endtime, ntrials, reps, tsplits, freq, stacks)

copyfile('UsualParams.mat', strcat('exp', expnum, '/', 'exp', expnum, '-params.mat'))
