load('UsualParams.mat')
expnum = 'G2';

% Run MVGC toolbox and vote for final connectivity matrix.
nvars = 2;
mats = MakeNetworkTriDiag(nvars);
reps = 20;
ntrials = 1;

% Save all initial conditions in experiment folder.
save('UsualParams.mat')

BaseExperiment(expnum, mats, randpfn, randvfn, randmfn, randkfn, ...
    randcfn, preprocfn, deltat, endtime, ntrials, reps, tsplits, freq)

copyfile('UsualParams.mat', strcat('exp', expnum, '/', 'exp', expnum, '-params.mat'))
