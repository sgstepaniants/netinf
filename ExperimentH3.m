load('UsualParams.mat')
expnum = 'H3_2s';

% Run MVGC toolbox 10 times and vote for final connectivity matrix.
nvars = 3;
mats = MakeNetworkTriDiag(nvars, stacks);
reps = 1;
ntrials = 100;

% Save all initial conditions in experiment folder.
save('UsualParams.mat')

BaseExperiment(expnum, mats, randpfn, randvfn, randmfn, randkfn, ...
    randcfn, preprocfn, deltat, endtime, ntrials, reps, tsplits, freq, stacks)

copyfile('UsualParams.mat', strcat('exp', expnum, '/', 'exp', expnum, '-params.mat'))
