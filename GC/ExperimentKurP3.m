load('UsualParams.mat')
expnum = 'KurP3';

% Run MVGC toolbox several times and vote for final connectivity matrix.
nvars = 3;

reps = 100;
ntrials = 20;

pfn = @(n) 2*pi*rand([n, 1]); % uniform [0, 2pi]
wfn = @(n) 2*rand([n, 1]) - ones(n,1); % uniform [-1, 1]
cfn = @(n) constcfn(n, 10);

A = MakeNetworkER(nvars, 0.5, true) % adjacency matrix
K = 10; % connection strength

kickTime = 50;
pert = 50;

% Save all initial conditions in experiment folder.
save('UsualParams.mat')

BaseKuramotoExperiment(expnum, A, K, pfn, wfn, cfn, kickTime, pert, ...
    preprocfn, deltat, endtime, ntrials, reps, tsplits, freq)

copyfile('UsualParams.mat', strcat('exp', expnum, '/', 'exp', expnum, '-params.mat'))
