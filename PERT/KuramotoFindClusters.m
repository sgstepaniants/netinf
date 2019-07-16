% Generate time series data for n Kuramoto oscillators with certain
% initial conditions and boundary conditions.

clear all; close all; clc;
addpath('../DataScripts/SimulateData/')
addpath('../DataScripts/SimulateData/InitFunctions/')
addpath('../BAgraph_dir/')

nvars = 10;

endtime = 10;
deltat = 0.1;
nobs = round(endtime / deltat);
tSpan = linspace(0, endtime, nobs);

noisefn = @(data) WhiteGaussianNoise(data, 0.1);

prob = 0.1;
mat = MakeNetworkER(nvars, prob, true);
mat(1:floor(nvars/2), 1:floor(nvars/2)) = BAgraph_dir(floor(nvars/2), floor(3*nvars/8), floor(3*nvars/8)).';
mat(floor(nvars/2)+1:end, floor(nvars/2)+1:end) = BAgraph_dir(floor(nvars/2), floor(3*nvars/8), floor(3*nvars/8)).';
K = 20; % connection strength

pfn = @(n) randfn(n, 0, 2*pi);
wfn = @(n) randfn(n, -1, 1);

forcingFunc = zeros([nvars, nobs]);

numTrials = 100;
data = GenerateKuramotoData(mat, tSpan, numTrials, K, pfn, wfn, forcingFunc);

repData = permute(repmat(data, [1, 1, 1, nvars]), [1, 4, 2, 3]);
aveCorr = mean(cos(repData - permute(repData, [2, 1, 3, 4])), 4);

G = digraph(mat.');
figure(1)
plot(G)

figure(2)
aves = aveCorr(:, :, 100);
aves(1 : nvars + 1 : nvars^2) = [];
%h = histogram(aves, 20);
imagesc(normalize(aveCorr(:, :, 100)) > 0.85)

figure(3)
imagesc(aveCorr(:, :, 100))
