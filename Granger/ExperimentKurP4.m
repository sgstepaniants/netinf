load('UsualParams.mat')
expnum = 'KurP4';

% Run MVGC toolbox several times and vote for final connectivity matrix.
nvars = 3;

reps = 1;
ntrials = 100;

pfn = @(n) 2*pi*rand([n, 1]); % uniform [0, 2pi]
wfn = @(n) 2*rand([n, 1]) - ones(n,1); % uniform [-1, 1]
cfn = @(n) constcfn(n, 0.4);

%A = MakeNetworkER(nvars, 0.5, true) % adjacency matrix
A = [[0, 1, 1]; [1, 0, 1]; [0, 1, 0]]
K = 1; % connection strength

kickTime = nobs/4;
pert = 10;

% Save all initial conditions in experiment folder.
save('UsualParams.mat')

BaseKuramotoExperiment(expnum, A, K, pfn, wfn, cfn, kickTime, pert, ...
    preprocfn, deltat, endtime, ntrials, reps, tsplits, freq)

copyfile('UsualParams.mat', strcat('exp', expnum, '/', 'exp', expnum, '-params.mat'))
