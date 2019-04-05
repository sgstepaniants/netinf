load('UsualParams.mat')
expnum = 'GD4';

% Run MVGC toolbox several times and vote for final connectivity matrix.
nvars = 4;
mats = MakeNetworkTriDiag(nvars, stacks);
reps = 10;
ntrials = 10;

% Save all initial conditions in experiment folder.
save('UsualParams.mat')

BaseExperiment(expnum, mats, randpfn, randvfn, randmfn, randkfn, ...
    randcfn, preprocfn, deltat, endtime, ntrials, reps, tsplits, freq, stacks)

copyfile('UsualParams.mat', strcat('exp', expnum, '/', 'exp', expnum, '-params.mat'))
