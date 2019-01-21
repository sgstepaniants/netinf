load('UsualParams.mat')
expnum = 'P2';

% Run MVGC toolbox several times and vote for final connectivity matrix.
nvars = 2;
p = 0.5;

mat = MakeNetworkSymmER(nvars, p, true);

reps = 100;
ntrials = 20;

% Specify parameters of the spring mass system.
pert = 1;
mass = 1;

spring = 0.01;
K = MakeNetworkTriDiag(nvars + 2, false);
K(2:nvars+1, 2:nvars+1) = mat;
K = spring * K;

damp = 0.5;
bc = 'fixed';

% Save all initial conditions in experiment folder.
save('UsualParams.mat')

BaseExperiment(expnum, mat, K, @(n)onepfn(n, pert, bc), @(n)zerovfn(n), @(n)constmfn(n, mass), ...
    @(n)constcfn(n, damp), preprocfn, deltat, endtime, ntrials, ...
    reps, tsplits, freq, stacks)

copyfile('UsualParams.mat', strcat('exp', expnum, '/', 'exp', expnum, '-params.mat'))
