load('UsualParams.mat')
expnum = 'S4_smallMass=5';

% Run MVGC toolbox several times and vote for final connectivity matrix.
nvars = 4;
mats = MakeNetworkTriDiag(nvars);
reps = 10;
ntrials = 100;

% kick the first mass but keep the rest at their equilibrium positions
pert = 0.1;
randpfn = @(n) [pert; [1 : n - 1]' / n]; % evenly spaced

% make one mass small and the rest large
smallMass = 5;
randmfn = @(n) [smallMass; 100 * ones(n - 1, 1)];

% Save all initial conditions in experiment folder.
save('UsualParams.mat')

BaseExperiment(expnum, mats, randpfn, randvfn, randmfn, randkfn, ...
    randcfn, preprocfn, deltat, endtime, ntrials, reps, tsplits, freq, stacks)

copyfile('UsualParams.mat', strcat('exp', expnum, '/', 'exp', expnum, '-params.mat'))
