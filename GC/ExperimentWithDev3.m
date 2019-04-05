load('UsualParams.mat')
expnum = 'WithDev3';

% Run MVGC toolbox several times and vote for final connectivity matrix.
nvars = 6;
p = 0.5;

reps = 100;
ntrials = 20;

% Specify parameters of the spring mass system.
pert = 1;
mass = 1;

spring = 1;
mat = MakeNetworkSymmER(nvars/2, p, true);

K = MakeNetworkTriDiag(nvars/2 + 2, false);
K(2:nvars/2+1, 2:nvars/2+1) = mat;
K = spring * K;

damp = 0;
bc = 'fixed';

% Save all initial conditions in experiment folder.
save('UsualParams.mat')

BaseExperiment(expnum, mat, K, @(n)randpfn(n), @(n)zerovfn(n), @(n)constmfn(n, mass), ...
    @(n)constcfn(n, damp), preprocfn, deltat, endtime, ntrials, ...
    reps, tsplits, freq)

copyfile('UsualParams.mat', strcat('exp', expnum, '/', 'exp', expnum, '-params.mat'))
