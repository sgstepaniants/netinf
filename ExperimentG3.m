load('UsualParams.mat')
expnum = 'G3';

% Run MVGC toolbox several times and vote for final connectivity matrix.
nvars = 3;
mats = MakeNetworkTriDiag(nvars, stacks);
reps = 30;
ntrials = 1;

% Save all initial conditions in experiment folder.
save('UsualParams.mat')

BaseExperiment(expnum, mats, randpfn, randvfn, randmfn, randkfn, ...
    randcfn, preprocfn, deltat, endtime, ntrials, reps, tsplits, freq, stacks)

copyfile('UsualParams.mat', strcat('exp', expnum, '/', 'exp', expnum, '-params.mat'))
